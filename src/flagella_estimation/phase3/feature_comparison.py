"""Grouped 3D/2D feature comparison for Phase 3 clip datasets.

The 3D features are inherited from the Phase 2 source summary.  The 2D
features are measured from the actual body-only clip pixels, rather than from
the body-centred source trajectory.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


THREED_FEATURES = (
    "cell_displacement",
    "cell_path_length",
    "cell_mean_speed",
    "cell_speed_std",
    "cell_straightness",
    "cell_axis_angle_change",
    "cell_axis_angle_std",
    "cell_angular_velocity_rms",
    "cell_wobble",
    "flagella_axis_alignment",
    "flagella_axis_spread",
    "cell_flagella_axis_angle",
    "cell_flagella_axis_stability",
)

PIXEL_FEATURES = (
    "2d_area_px_mean",
    "2d_area_px_std",
    "2d_major_axis_px_mean",
    "2d_minor_axis_px_mean",
    "2d_aspect_ratio_mean",
    "2d_aspect_ratio_std",
    "2d_axis_angle_change_rad",
    "2d_axis_angular_velocity_rms_rad_s",
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _shape(frame: np.ndarray) -> tuple[float, float, float, float]:
    ys, xs = np.nonzero(frame < 255)
    if len(xs) < 3:
        return (float("nan"),) * 4
    centered = np.column_stack((xs - xs.mean(), ys - ys.mean()))
    values, vectors = np.linalg.eigh(np.cov(centered, rowvar=False))
    values = np.maximum(values, 0.0)
    major, minor = np.sqrt(values[::-1]) * 4.0
    direction = vectors[:, -1]
    angle = math.atan2(float(direction[1]), float(direction[0]))
    return float(len(xs)), float(major), float(minor), angle


def body_axis_angles_rad(clip: np.ndarray) -> np.ndarray:
    """Return unwrapped, pi-periodic long-axis angles from body silhouettes."""
    values = np.asarray([_shape(frame) for frame in clip], dtype=float)
    # Capsule orientation is axial (pi-periodic), hence unwrap doubled angles.
    return np.unwrap(2.0 * values[:, 3]) / 2.0


def pixel_features(clip: np.ndarray, fps: float) -> dict[str, float]:
    values = np.asarray([_shape(frame) for frame in clip], dtype=float)
    area, major, minor, _ = values.T
    aspect = np.divide(major, minor, out=np.full_like(major, np.nan), where=minor > 0)
    angle = body_axis_angles_rad(clip)
    angular_velocity = np.diff(angle) * fps
    return {
        "2d_area_px_mean": float(np.nanmean(area)),
        "2d_area_px_std": float(np.nanstd(area)),
        "2d_major_axis_px_mean": float(np.nanmean(major)),
        "2d_minor_axis_px_mean": float(np.nanmean(minor)),
        "2d_aspect_ratio_mean": float(np.nanmean(aspect)),
        "2d_aspect_ratio_std": float(np.nanstd(aspect)),
        "2d_axis_angle_change_rad": float(angle[-1] - angle[0]),
        "2d_axis_angular_velocity_rms_rad_s": float(
            np.sqrt(np.nanmean(angular_velocity**2))
        )
        if len(angular_velocity)
        else 0.0,
    }


def grouped_nearest_centroid_pixel_baseline(
    rows: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Evaluate pixel features with source-run-held-out nearest centroids.

    ``group_key`` is held out as a whole, so sibling windows from a source run
    never contribute to the centroid used for that source's prediction.
    """
    candidate_rows = [
        row
        for row in rows
        if bool(row["training_candidate"]) and int(row["n_flagella"]) in {1, 2, 3}
    ]
    groups = sorted({str(row["group_key"]) for row in candidate_rows})
    predictions: list[dict[str, Any]] = []
    for held_out_group in groups:
        train = [row for row in candidate_rows if row["group_key"] != held_out_group]
        test = [row for row in candidate_rows if row["group_key"] == held_out_group]
        if not train or not test:
            continue
        train_x = np.asarray(
            [[float(row[key]) for key in PIXEL_FEATURES] for row in train], dtype=float
        )
        center = np.nanmean(train_x, axis=0)
        scale = np.nanstd(train_x, axis=0)
        usable = np.isfinite(center) & np.isfinite(scale) & (scale > 1.0e-12)
        if not np.any(usable):
            continue
        centroids = {
            n: np.nanmean(
                np.asarray(
                    [
                        [float(row[key]) for key in PIXEL_FEATURES]
                        for row in train
                        if int(row["n_flagella"]) == n
                    ],
                    dtype=float,
                )[:, usable],
                axis=0,
            )
            for n in (1, 2, 3)
        }
        for row in test:
            x = np.asarray([float(row[key]) for key in PIXEL_FEATURES], dtype=float)
            x = (x[usable] - center[usable]) / scale[usable]
            distances = {
                n: float(
                    np.linalg.norm(x - (centroid - center[usable]) / scale[usable])
                )
                for n, centroid in centroids.items()
            }
            predicted = min(distances, key=distances.get)
            predictions.append(
                {
                    "clip_id": row["clip_id"],
                    "group_key": held_out_group,
                    "n_flagella": int(row["n_flagella"]),
                    "predicted_n_flagella": predicted,
                    "correct": predicted == int(row["n_flagella"]),
                    "nearest_centroid_distance": distances[predicted],
                }
            )
    accuracy = (
        float(np.mean([bool(row["correct"]) for row in predictions]))
        if predictions
        else float("nan")
    )
    return predictions, {
        "feature_set": "2d_pixel",
        "evaluation_unit": "clip",
        "group_holdout": True,
        "eligible_group_count": len(groups),
        "prediction_count": len(predictions),
        "accuracy": accuracy,
    }


def evaluate(dataset_dir: Path, output_dir: Path) -> Path:
    """Write per-clip 3D+2D features and grouped distribution summaries.

    n=4 diagnostic rows remain in distributions, while nearest-centroid
    accuracy is restricted to Phase 3 training candidates (n=1..3).
    """
    manifest = json.loads((dataset_dir / "manifest.json").read_text(encoding="utf-8"))
    source = Path(manifest["input_dataset"])
    source_rows = {row["sample_id"]: row for row in _read_csv(source / "summary.csv")}
    records = [
        json.loads(line)
        for line in (dataset_dir / "clip_metadata.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
        if line
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for record in records:
        clip = record["clip"]
        source_row = source_rows[record["provenance"]["run_id"]]
        row: dict[str, Any] = {
            "clip_id": clip["clip_id"],
            "sample_id": record["provenance"]["run_id"],
            "group_key": record["track"]["group_key"],
            "n_flagella": record["labels"]["n_flagella"],
            "training_candidate": record["qc"]["training_candidate"],
            "diagnostic_only": record["qc"]["diagnostic_only"],
        }
        for key in THREED_FEATURES:
            row[f"3d_{key}"] = source_row.get(key, "")
        array = np.load(Path(clip["output_path"]))
        row.update(pixel_features(array, float(clip["frame_rate_hz"])))
        rows.append(row)
    fields = list(rows[0]) if rows else ["clip_id"]
    _write_csv(output_dir / "features_3d_2d.csv", rows, fields)
    numeric = [key for key in fields if key.startswith(("2d_", "3d_"))]
    summary: list[dict[str, Any]] = []
    for n in sorted({int(row["n_flagella"]) for row in rows}):
        subset = [row for row in rows if int(row["n_flagella"]) == n]
        for key in numeric:
            data = np.asarray(
                [float(row[key]) for row in subset if str(row[key]).strip() != ""],
                dtype=float,
            )
            data = data[np.isfinite(data)]
            summary.append(
                {
                    "n_flagella": n,
                    "feature": key,
                    "count": len(data),
                    "mean": float(np.mean(data)) if len(data) else float("nan"),
                    "std": float(np.std(data)) if len(data) else float("nan"),
                }
            )
    _write_csv(
        output_dir / "feature_summary_by_n.csv",
        summary,
        list(summary[0]) if summary else ["n_flagella"],
    )
    predictions, baseline = grouped_nearest_centroid_pixel_baseline(rows)
    _write_csv(
        output_dir / "grouped_nearest_centroid_pixel.csv",
        predictions,
        list(predictions[0]) if predictions else ["clip_id"],
    )
    (output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "dataset_dir": str(dataset_dir),
                "feature_sets": {
                    "3d": list(THREED_FEATURES),
                    "2d": [key for key in numeric if key.startswith("2d_")],
                },
                "evaluation_policy": "n=4 distribution-only; grouped downstream evaluation must hold out group_key",
                "grouped_nearest_centroid_pixel": baseline,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return output_dir
