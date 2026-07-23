#!/usr/bin/env python3
"""Analyze whether 2D projected motion preserves flagella-count differences."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
import math
from pathlib import Path
import sys
from typing import Any
from zoneinfo import ZoneInfo

import numpy as np


DEFAULT_DATASET_DIR = Path(
    "outputs/phase2_analysis/flagella_count_behavior/datasets/v1"
)
DEFAULT_N_FLAGELLA = (1, 2, 3)
FEATURE_COLUMNS = [
    "xy_displacement_um",
    "xy_path_length_um",
    "xy_mean_speed_um_s",
    "xy_speed_std_um_s",
    "xy_speed_cv",
    "xy_straightness",
    "xy_range_x_um",
    "xy_range_y_um",
    "xy_extent_um",
    "xy_body_axis_angle_change_deg",
    "xy_body_axis_angle_std_deg",
    "xy_body_axis_angular_velocity_rms_rad_s",
    "xy_velocity_heading_change_deg",
    "xy_velocity_heading_std_deg",
]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _to_float(value: Any, default: float = float("nan")) -> float:
    try:
        if value in (None, ""):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _to_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _finite(values: np.ndarray) -> np.ndarray:
    return values[np.isfinite(values)]


def _mean(values: np.ndarray) -> float:
    values = _finite(values.astype(float))
    return float(np.mean(values)) if values.size else float("nan")


def _std(values: np.ndarray) -> float:
    values = _finite(values.astype(float))
    if values.size == 0:
        return float("nan")
    if values.size == 1:
        return 0.0
    return float(np.std(values, ddof=1))


def _unwrap_degrees(values_deg: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values_deg)
    out = np.full_like(values_deg, np.nan, dtype=float)
    if not np.any(finite):
        return out
    out[finite] = np.rad2deg(np.unwrap(np.deg2rad(values_deg[finite])))
    return out


def _quat_to_xy_axis_angle_deg(qx: float, qy: float, qz: float, qw: float) -> float:
    q = np.asarray([qx, qy, qz, qw], dtype=float)
    norm = float(np.linalg.norm(q))
    if not math.isfinite(norm) or norm <= 0.0:
        return float("nan")
    x, y, z, w = q / norm
    # Rotate the model body x-axis by q and keep the observed XY heading.
    axis_x = 1.0 - 2.0 * (y * y + z * z)
    axis_y = 2.0 * (x * y + z * w)
    if abs(axis_x) < 1.0e-15 and abs(axis_y) < 1.0e-15:
        return float("nan")
    return float(np.rad2deg(np.arctan2(axis_y, axis_x)))


def _trajectory_arrays(rows: list[dict[str, str]]) -> dict[str, np.ndarray]:
    return {
        "t": np.asarray([_to_float(row.get("t")) for row in rows], dtype=float),
        "x": np.asarray([_to_float(row.get("x")) for row in rows], dtype=float),
        "y": np.asarray([_to_float(row.get("y")) for row in rows], dtype=float),
        "qx": np.asarray([_to_float(row.get("qx")) for row in rows], dtype=float),
        "qy": np.asarray([_to_float(row.get("qy")) for row in rows], dtype=float),
        "qz": np.asarray([_to_float(row.get("qz")) for row in rows], dtype=float),
        "qw": np.asarray([_to_float(row.get("qw")) for row in rows], dtype=float),
    }


def summarize_trajectory_2d(trajectory_csv: Path) -> dict[str, float | int]:
    rows = _read_csv(trajectory_csv)
    if len(rows) < 2:
        return {key: float("nan") for key in FEATURE_COLUMNS} | {
            "trajectory_step_count": len(rows)
        }

    arr = _trajectory_arrays(rows)
    t = arr["t"]
    x = arr["x"]
    y = arr["y"]
    valid = np.isfinite(t) & np.isfinite(x) & np.isfinite(y)
    if int(np.count_nonzero(valid)) < 2:
        return {key: float("nan") for key in FEATURE_COLUMNS} | {
            "trajectory_step_count": len(rows)
        }

    t = t[valid]
    # The v1 2D renderer centers each frame on the body.  Body-center
    # translation is therefore absent from projection.mp4, so evaluate body
    # motion features in that camera frame instead of raw world XY.
    x = np.zeros_like(x[valid])
    y = np.zeros_like(y[valid])
    dx = np.diff(x)
    dy = np.diff(y)
    dt = np.diff(t)
    segment = np.hypot(dx, dy)
    valid_dt = np.isfinite(dt) & (dt > 0.0)
    speed = np.divide(
        segment,
        dt,
        out=np.full_like(segment, np.nan, dtype=float),
        where=valid_dt,
    )
    path_length = float(np.nansum(segment))
    displacement = float(math.hypot(float(x[-1] - x[0]), float(y[-1] - y[0])))

    q_valid = valid
    angle = np.asarray(
        [
            _quat_to_xy_axis_angle_deg(qx, qy, qz, qw)
            for qx, qy, qz, qw in zip(
                arr["qx"][q_valid],
                arr["qy"][q_valid],
                arr["qz"][q_valid],
                arr["qw"][q_valid],
                strict=True,
            )
        ],
        dtype=float,
    )
    angle_unwrapped = _unwrap_degrees(angle)
    d_angle = np.diff(angle_unwrapped)
    angular_velocity = np.divide(
        np.deg2rad(d_angle),
        dt,
        out=np.full_like(d_angle, np.nan, dtype=float),
        where=valid_dt & np.isfinite(d_angle),
    )

    heading = np.rad2deg(np.arctan2(dy, dx))
    heading[segment <= 1.0e-12] = np.nan
    heading_unwrapped = _unwrap_degrees(heading)
    heading_change = (
        float(heading_unwrapped[-1] - heading_unwrapped[0])
        if np.isfinite(heading_unwrapped[0]) and np.isfinite(heading_unwrapped[-1])
        else float("nan")
    )

    speed_mean = _mean(speed)
    speed_std = _std(speed)
    return {
        "trajectory_step_count": len(rows),
        "xy_displacement_um": displacement,
        "xy_path_length_um": path_length,
        "xy_mean_speed_um_s": speed_mean,
        "xy_speed_std_um_s": speed_std,
        "xy_speed_cv": speed_std / speed_mean
        if math.isfinite(speed_std) and abs(speed_mean) > 1.0e-30
        else float("nan"),
        "xy_straightness": displacement / path_length
        if path_length > 0.0
        else float("nan"),
        "xy_range_x_um": float(np.nanmax(x) - np.nanmin(x)),
        "xy_range_y_um": float(np.nanmax(y) - np.nanmin(y)),
        "xy_extent_um": float(
            math.hypot(
                float(np.nanmax(x) - np.nanmin(x)), float(np.nanmax(y) - np.nanmin(y))
            )
        ),
        "xy_body_axis_angle_change_deg": float(angle_unwrapped[-1] - angle_unwrapped[0])
        if np.isfinite(angle_unwrapped[0]) and np.isfinite(angle_unwrapped[-1])
        else float("nan"),
        "xy_body_axis_angle_std_deg": _std(angle_unwrapped),
        "xy_body_axis_angular_velocity_rms_rad_s": float(
            math.sqrt(np.nanmean(angular_velocity**2))
        )
        if np.isfinite(angular_velocity).any()
        else float("nan"),
        "xy_velocity_heading_change_deg": heading_change,
        "xy_velocity_heading_std_deg": _std(heading_unwrapped),
    }


def _group_id(row: dict[str, Any]) -> str:
    return (
        f"as{int(float(row['attach_seed'])):03d}_ps{int(float(row['phase_seed'])):03d}"
    )


def _standardize(train: np.ndarray, test: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = np.zeros(train.shape[1], dtype=float)
    std = np.ones(train.shape[1], dtype=float)
    for col_idx in range(train.shape[1]):
        finite = _finite(train[:, col_idx])
        if finite.size:
            mean[col_idx] = float(np.mean(finite))
            raw_std = float(np.std(finite, ddof=0))
            if math.isfinite(raw_std) and raw_std > 0.0:
                std[col_idx] = raw_std
    train_out = (train - mean) / std
    test_out = (test - mean) / std
    train_out[~np.isfinite(train_out)] = 0.0
    test_out[~np.isfinite(test_out)] = 0.0
    return train_out, test_out


def grouped_nearest_centroid(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups = sorted({_group_id(row) for row in rows})
    out: list[dict[str, Any]] = []
    for group in groups:
        train_rows = [row for row in rows if _group_id(row) != group]
        test_rows = [row for row in rows if _group_id(row) == group]
        train_x = np.asarray(
            [
                [_to_float(row.get(feature)) for feature in FEATURE_COLUMNS]
                for row in train_rows
            ],
            dtype=float,
        )
        test_x = np.asarray(
            [
                [_to_float(row.get(feature)) for feature in FEATURE_COLUMNS]
                for row in test_rows
            ],
            dtype=float,
        )
        train_x, test_x = _standardize(train_x, test_x)
        labels = sorted({int(float(row["n_flagella"])) for row in train_rows})
        centroids = {
            label: np.mean(
                train_x[
                    np.asarray(
                        [int(float(row["n_flagella"])) == label for row in train_rows]
                    )
                ],
                axis=0,
            )
            for label in labels
        }
        for row, features in zip(test_rows, test_x, strict=True):
            distances = {
                label: float(np.nansum((features - centroid) ** 2))
                for label, centroid in centroids.items()
            }
            predicted = min(distances, key=distances.get)
            actual = int(float(row["n_flagella"]))
            out.append(
                {
                    "sample_id": row["sample_id"],
                    "group_id": group,
                    "n_flagella": actual,
                    "predicted_n_flagella": predicted,
                    "correct": actual == predicted,
                    "nearest_centroid_distance": distances[predicted],
                }
            )
    return out


def _feature_stats(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    labels = sorted({int(float(row["n_flagella"])) for row in rows})
    for feature in FEATURE_COLUMNS:
        for label in labels:
            values = np.asarray(
                [
                    _to_float(row.get(feature))
                    for row in rows
                    if int(float(row["n_flagella"])) == label
                ],
                dtype=float,
            )
            finite = _finite(values)
            out.append(
                {
                    "feature": feature,
                    "n_flagella": label,
                    "count": int(finite.size),
                    "mean": float(np.mean(finite)) if finite.size else float("nan"),
                    "std": float(np.std(finite, ddof=1))
                    if finite.size > 1
                    else (0.0 if finite.size == 1 else float("nan")),
                    "min": float(np.min(finite)) if finite.size else float("nan"),
                    "max": float(np.max(finite)) if finite.size else float("nan"),
                }
            )
    return out


def analyze_2d_separability(
    *,
    dataset_dir: Path,
    output_dir: Path | None = None,
    n_flagella: tuple[int, ...] = DEFAULT_N_FLAGELLA,
    ml_candidates_only: bool = True,
    overwrite: bool = False,
) -> Path:
    summary_csv = dataset_dir / "summary.csv"
    if not summary_csv.is_file():
        raise FileNotFoundError(summary_csv)
    output_dir = output_dir or dataset_dir / "analysis" / "projection_2d"
    if output_dir.exists() and not overwrite:
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    missing_trajectories: list[tuple[str, Path]] = []
    for sample in _read_csv(summary_csv):
        n = int(_to_float(sample.get("n_flagella"), default=-1.0))
        if n not in n_flagella:
            continue
        if ml_candidates_only and not _to_bool(sample.get("use_for_ml_candidate")):
            continue
        raw_dir = Path(str(sample.get("raw_dir", "")))
        trajectory_csv = raw_dir / "trajectory.csv"
        if not trajectory_csv.is_file():
            missing_trajectories.append(
                (str(sample.get("sample_id", "")), trajectory_csv)
            )
            continue
        rows.append(
            {
                "sample_id": sample.get("sample_id", ""),
                "dataset_id": sample.get("dataset_id", ""),
                "n_flagella": n,
                "attach_seed": sample.get("attach_seed", ""),
                "phase_seed": sample.get("phase_seed", ""),
                "quality_class": sample.get("quality_class", ""),
                "use_for_ml_candidate": sample.get("use_for_ml_candidate", ""),
                "raw_dir": str(raw_dir),
                "trajectory_csv": str(trajectory_csv),
                **summarize_trajectory_2d(trajectory_csv),
            }
        )
    if missing_trajectories:
        missing = ", ".join(
            f"{sample_id or '<missing sample_id>'}: {path}"
            for sample_id, path in missing_trajectories
        )
        raise FileNotFoundError(
            "Missing trajectory.csv for eligible samples: " + missing
        )
    if not rows:
        raise ValueError(f"No matching trajectory rows found under {dataset_dir}")

    features_path = output_dir / "features_2d.csv"
    stats_path = output_dir / "feature_summary_by_n_flagella.csv"
    baseline_path = output_dir / "grouped_nearest_centroid_baseline.csv"
    _write_csv(
        features_path,
        rows,
        [
            "sample_id",
            "dataset_id",
            "n_flagella",
            "attach_seed",
            "phase_seed",
            "quality_class",
            "use_for_ml_candidate",
            "raw_dir",
            "trajectory_csv",
            "trajectory_step_count",
            *FEATURE_COLUMNS,
        ],
    )
    _write_csv(
        stats_path,
        _feature_stats(rows),
        ["feature", "n_flagella", "count", "mean", "std", "min", "max"],
    )
    baseline_rows = grouped_nearest_centroid(rows)
    _write_csv(
        baseline_path,
        baseline_rows,
        [
            "sample_id",
            "group_id",
            "n_flagella",
            "predicted_n_flagella",
            "correct",
            "nearest_centroid_distance",
        ],
    )
    accuracy = sum(1 for row in baseline_rows if row["correct"]) / max(
        len(baseline_rows), 1
    )
    manifest = {
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "dataset_dir": str(dataset_dir),
        "summary_csv": str(summary_csv),
        "n_flagella": list(n_flagella),
        "ml_candidates_only": ml_candidates_only,
        "sample_count": len(rows),
        "grouped_nearest_centroid_accuracy": accuracy,
        "outputs": {
            "features_2d_csv": str(features_path),
            "feature_summary_by_n_flagella_csv": str(stats_path),
            "grouped_nearest_centroid_baseline_csv": str(baseline_path),
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return output_dir


def _parse_n_flagella(value: str) -> tuple[int, ...]:
    return tuple(int(part.strip()) for part in value.split(",") if part.strip())


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", type=Path, default=DEFAULT_DATASET_DIR)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--n-flagella", default="1,2,3")
    parser.add_argument("--include-non-ml-candidates", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)
    output_dir = analyze_2d_separability(
        dataset_dir=args.dataset_dir,
        output_dir=args.output_dir,
        n_flagella=_parse_n_flagella(args.n_flagella),
        ml_candidates_only=not args.include_non_ml_candidates,
        overwrite=bool(args.overwrite),
    )
    print(output_dir)


if __name__ == "__main__":
    main()
