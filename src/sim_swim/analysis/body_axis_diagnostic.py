"""Read completed Phase 2 archives for body-axis direction diagnostics.

This is descriptive analysis: it adds no QC gate or stationarity threshold.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from copy import deepcopy
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from sim_swim.analysis.model_development_evaluation import collect_rows
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.helix_axis import estimate_body_axis
from sim_swim.sim.params import SimulationConfig


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _angle_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.rad2deg(np.arccos(np.clip(np.sum(a * b, axis=-1), -1.0, 1.0)))


def axis_rows(
    *,
    condition_id: str,
    times: np.ndarray,
    axes: np.ndarray,
    wobble_final_deg: float,
    window_s: float = 0.1,
) -> tuple[list[dict], list[dict], dict]:
    """Return signed 3D directions and nonoverlapping physical-time windows."""
    if times.ndim != 1 or axes.shape != (len(times), 3) or len(times) < 2:
        raise ValueError("Body-axis archive must have at least two aligned states")
    if not np.isfinite(times).all() or not np.isfinite(axes).all():
        raise ValueError("Body-axis archive contains non-finite data")
    if not np.all(np.diff(times) > 0) or not np.allclose(
        np.linalg.norm(axes, axis=1), 1.0, atol=1e-8
    ):
        raise ValueError("Body-axis times or unit vectors are invalid")
    initial = np.broadcast_to(axes[0], axes.shape)
    initial_angle = _angle_deg(axes, initial)
    step_angle = np.r_[0.0, _angle_deg(axes[1:], axes[:-1])]
    series = [
        {
            "condition_id": condition_id,
            "t_s": float(t),
            "axis_x": float(axis[0]),
            "axis_y": float(axis[1]),
            "axis_z": float(axis[2]),
            "angle_from_initial_deg": float(angle),
            "step_angle_deg": float(step),
        }
        for t, axis, angle, step in zip(times, axes, initial_angle, step_angle)
    ]
    windows = []
    # Integrator accumulation can put a nominal 2.0 s endpoint a few ps past
    # the boundary. Keep that sample in the final full window.
    end_s = round(float(times[-1]), 9)
    for index in range(int(np.ceil(end_s / window_s))):
        start, end = index * window_s, min((index + 1) * window_s, end_s)
        mask = (times >= start - 1e-10) & (times < end - 1e-10)
        if np.isclose(end, end_s):
            mask |= np.isclose(times, end_s, atol=1e-9)
        segment = axes[mask]
        if len(segment) < 2:
            continue
        mean_axis = segment.mean(axis=0)
        mean_norm = float(np.linalg.norm(mean_axis))
        if mean_norm < 1e-12:
            mean_axis = segment[0]
        else:
            mean_axis /= mean_norm
        deviations = _angle_deg(segment, np.broadcast_to(mean_axis, segment.shape))
        windows.append(
            {
                "condition_id": condition_id,
                "window_start_s": float(start),
                "window_end_s": float(end),
                "sample_count": len(segment),
                "direction_rms_deg": float(np.sqrt(np.mean(deviations**2))),
                "direction_max_deg": float(np.max(deviations)),
                "endpoint_drift_deg": float(_angle_deg(segment[:1], segment[-1:])[0]),
                "max_angle_from_initial_deg": float(np.max(initial_angle[mask])),
                "mean_direction_degenerate": mean_norm < 1e-12,
            }
        )
    summary = {
        "condition_id": condition_id,
        "archive_states": len(times),
        "final_t_s": end_s,
        "max_angle_from_initial_deg": float(np.max(initial_angle)),
        "body_axis_wobble_rms_deg_final": wobble_final_deg,
    }
    return series, windows, summary


def _verified_hashes(
    run_dir: Path, verification: Path | None, source_commit: str
) -> int | None:
    if verification is None:
        return None
    record = json.loads(verification.read_text(encoding="utf-8"))
    if record.get("sha256_match") is not True:
        raise ValueError("Source synchronization did not pass SHA-256 verification")
    if record.get("source_commit") != source_commit:
        raise ValueError("Source synchronization commit mismatch")
    files = dict(record.get("files", {}))
    if not files or len(files) != int(record.get("file_count", -1)):
        raise ValueError("Source synchronization file count mismatch")
    for relative, expected in files.items():
        path = run_dir / relative
        if not path.is_file() or _sha256(path) != expected:
            raise ValueError(f"Source SHA-256 mismatch: {relative}")
    return len(files)


def analyze(
    *,
    config_path: Path,
    run_dir: Path,
    output_dir: Path,
    sync_verification: Path | None = None,
) -> dict:
    config = load_yaml(config_path)
    rows, provenance = collect_rows(config=config, run_dirs=[run_dir])
    verified_count = _verified_hashes(
        run_dir, sync_verification, str(provenance[0]["git"]["commit"])
    )
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    records = {str(item["condition_id"]): item for item in manifest["conditions"]}
    base = SimulationConfig.from_dict(load_yaml(Path(config["base_config"])))
    all_series, all_windows, all_summaries = [], [], []
    output_dir.mkdir(parents=True, exist_ok=True)
    for row in rows:
        cid = str(row["condition_id"])
        record = records[cid]
        cfg = base.with_overrides(record["config_overrides"])
        model = ModelBuilder(cfg).build()
        condition_dir = run_dir / "conditions" / cid
        summary = json.loads(
            (condition_dir / "run_summary.json").read_text(encoding="utf-8")
        )
        expected_steps = (
            config["development_evaluation"]
            .get("long_duration", {})
            .get("expected_steps")
        )
        if expected_steps is not None and (
            summary["execution"].get("row_count") != expected_steps
            or summary["execution"].get("expected_total_steps") != expected_steps
            or summary["execution"].get("step_indices_contiguous_from_zero") is not True
        ):
            raise ValueError(f"Incomplete step coverage: {cid}")
        with np.load(
            condition_dir / "state_archive.npz", allow_pickle=False
        ) as archive:
            times = np.asarray(archive["t"], dtype=float)
            positions = np.asarray(archive["bead_positions_um"], dtype=float) * 1e-6
        if len(times) != len(positions):
            raise ValueError(f"Archive time/position mismatch: {cid}")
        axes = []
        first_layer, last_layer = (
            model.body_layer_indices[0],
            model.body_layer_indices[-1],
        )
        for p in positions:
            if (
                np.linalg.norm(p[last_layer].mean(axis=0) - p[first_layer].mean(axis=0))
                <= 1e-18
            ):
                raise ValueError(f"Degenerate body axis: {cid}")
            axes.append(
                estimate_body_axis(p, model.body_layer_indices, model.body_indices).axis
            )
        axes = np.asarray(axes)
        wobble = float(summary["all_step_metrics"]["body_axis_wobble_rms_deg"]["final"])
        series, windows, compact = axis_rows(
            condition_id=cid, times=times, axes=axes, wobble_final_deg=wobble
        )
        compact["strict_qc"] = row["screen_status"]
        all_series.extend(series)
        all_windows.extend(windows)
        all_summaries.append(compact)
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.plot(times, [item["angle_from_initial_deg"] for item in series])
        ax.set(xlabel="time [s]", ylabel="body axis from initial [deg]", title=cid)
        fig.tight_layout()
        fig.savefig(output_dir / f"{cid}_axis.png", dpi=140)
        plt.close(fig)
    _write_csv(output_dir / "axis_timeseries.csv", all_series)
    _write_csv(output_dir / "axis_windows.csv", all_windows)
    _write_csv(output_dir / "axis_summary.csv", all_summaries)
    result = {
        "source_run": str(run_dir.resolve()),
        "source_config": str(config_path),
        "source_git": provenance[0]["git"],
        "condition_count": len(rows),
        "source_sha256_file_count": verified_count,
        "window_s": 0.1,
        "interpretation": "diagnostic only; no stationarity threshold",
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return result


def compare_axis_outputs(
    *,
    reference_config: Path,
    off_config: Path,
    reference_axis_dir: Path,
    off_axis_dir: Path,
    output_dir: Path,
) -> Path:
    """Pair windows only after proving the campaign condition differs by OFF."""
    before = normalize_campaign_config(load_yaml(reference_config))
    after = normalize_campaign_config(load_yaml(off_config))
    if before["base_config"] != after["base_config"]:
        raise ValueError("Reference and OFF base configs differ")
    expected_base = deepcopy(before["base_overrides"])
    expected_base.setdefault("potentials", {}).setdefault(
        "spring_spring_repulsion", {}
    )["flagella_flagella_enabled"] = False
    if after["base_overrides"] != expected_base:
        raise ValueError("Reference and OFF base overrides differ beyond repulsion")
    old_conditions = {
        str(item["condition_id"]): item for item in build_campaign_conditions(before)
    }
    new_conditions = {
        str(item["condition_id"]): item for item in build_campaign_conditions(after)
    }
    if set(old_conditions) != set(new_conditions) or len(old_conditions) != 13:
        raise ValueError("Reference and OFF condition IDs are not the same 13")
    for cid, old in old_conditions.items():
        expected = deepcopy(old["config_overrides"])
        expected.setdefault("potentials", {}).setdefault("spring_spring_repulsion", {})[
            "flagella_flagella_enabled"
        ] = False
        if new_conditions[cid]["config_overrides"] != expected:
            raise ValueError(f"Condition differs beyond repulsion: {cid}")
    old_manifest = json.loads((reference_axis_dir / "manifest.json").read_text())
    new_manifest = json.loads((off_axis_dir / "manifest.json").read_text())
    if old_manifest["condition_count"] != 13 or new_manifest["condition_count"] != 13:
        raise ValueError("Axis analysis does not cover 13 conditions")
    if str(old_manifest["source_config"]) != str(reference_config) or str(
        new_manifest["source_config"]
    ) != str(off_config):
        raise ValueError("Axis analysis source configs differ from comparison inputs")
    old_rows = {
        (r["condition_id"], r["window_start_s"]): r
        for r in _read_csv(reference_axis_dir / "axis_windows.csv")
    }
    new_rows = {
        (r["condition_id"], r["window_start_s"]): r
        for r in _read_csv(off_axis_dir / "axis_windows.csv")
    }
    if set(old_rows) != set(new_rows):
        raise ValueError("Reference and OFF time windows do not match")
    metrics = ("direction_rms_deg", "direction_max_deg", "endpoint_drift_deg")
    pairs = []
    for cid, start in sorted(old_rows):
        old, new = old_rows[(cid, start)], new_rows[(cid, start)]
        if old["window_end_s"] != new["window_end_s"]:
            raise ValueError(f"Window end differs: {cid} {start}")
        pair = {
            "condition_id": cid,
            "window_start_s": start,
            "window_end_s": old["window_end_s"],
        }
        for metric in metrics:
            on_value, off_value = float(old[metric]), float(new[metric])
            pair[f"on_{metric}"] = on_value
            pair[f"off_{metric}"] = off_value
            pair[f"off_minus_on_{metric}"] = off_value - on_value
        pairs.append(pair)
    output_dir.mkdir(parents=True, exist_ok=True)
    destination = output_dir / "paired_axis_windows.csv"
    _write_csv(destination, pairs)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sync-verification", type=Path)
    parser.add_argument("--reference-config", type=Path)
    parser.add_argument("--reference-axis-dir", type=Path)
    args = parser.parse_args()
    result = analyze(
        config_path=args.config,
        run_dir=args.run_dir,
        output_dir=args.output_dir,
        sync_verification=args.sync_verification,
    )
    if args.reference_config or args.reference_axis_dir:
        if not args.reference_config or not args.reference_axis_dir:
            parser.error("Both reference options are required for paired comparison")
        result["paired_axis_windows"] = str(
            compare_axis_outputs(
                reference_config=args.reference_config,
                off_config=args.config,
                reference_axis_dir=args.reference_axis_dir,
                off_axis_dir=args.output_dir,
                output_dir=args.output_dir,
            )
        )
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
