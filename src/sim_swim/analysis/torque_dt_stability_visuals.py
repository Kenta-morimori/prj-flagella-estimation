"""Offline qualitative-review artifacts for Issue #61 torque-dt campaigns."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")


FEATURES: tuple[tuple[str, str, str], ...] = (
    ("qc_pass", "QC pass (1=pass)", "status"),
    ("max_flag_bond_rel_err", "max flag bond relative error", "max"),
    ("max_hook_len_rel_err", "max hook length relative error", "max"),
    ("max_axis_mean_deviation_deg", "max axis mean deviation [deg]", "max"),
    ("min_axis_rearward_projection", "min axis rearward projection", "min"),
    ("final_bundle_participation_ratio", "final bundle participation ratio", "final"),
    ("max_bundle_radius_um", "max bundle radius [um]", "max"),
)

METRIC_KEYS = {
    "max_flag_bond_rel_err": "flag_bond_rel_err_max",
    "max_hook_len_rel_err": "hook_len_rel_err_max",
    "max_axis_mean_deviation_deg": "flag_helix_axis_mean_deviation_deg_max",
    "min_axis_rearward_projection": "flag_helix_axis_rearward_projection_min",
    "final_bundle_participation_ratio": "bundle_participation_ratio",
    "max_bundle_radius_um": "flag_helix_bundle_radius_max_um",
}


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _metric(summary: dict[str, Any], name: str, statistic: str) -> float:
    value = (
        summary.get("all_step_metrics", {}).get(METRIC_KEYS[name], {}).get(statistic)
    )
    return float(value) if value is not None else float("nan")


def feature_rows(run_dir: Path) -> list[dict[str, Any]]:
    """Extract compact shape and bundling diagnostics without restarting a run."""

    manifest = _read_json(run_dir / "run_manifest.json")
    rows: list[dict[str, Any]] = []
    for record in manifest["conditions"]:
        output_dir = Path(str(record["output_dir"]))
        summary = _read_json(output_dir / "run_summary.json")
        gates = dict(summary.get("gates", {}) or {})
        passes = summary.get("execution", {}).get("status") == "completed" and all(
            not bool(gates.get(name, {}).get("any_fail", True))
            for name in ("finite", "shape_nonbody", "shape_body")
        )
        row: dict[str, Any] = {
            "condition_id": record["condition_id"],
            "torque_Nm_per_flagellum": float(record["torque_Nm_per_flagellum"]),
            "dt_star": float(record["dt_star"]),
            "tau_s": float(record["time"]["tau_s"]),
            "duration_s": float(record["time"]["duration_s"]),
            "qc_pass": int(passes),
        }
        for key, _, statistic in FEATURES[1:]:
            row[key] = _metric(summary, key, statistic)
        rows.append(row)
    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _heatmap(rows: list[dict[str, Any]], *, output_path: Path, title: str) -> None:
    import matplotlib.pyplot as plt

    torques = sorted({float(row["torque_Nm_per_flagellum"]) for row in rows})
    dt_stars = sorted({float(row["dt_star"]) for row in rows}, reverse=True)
    by_cell = {
        (float(row["torque_Nm_per_flagellum"]), float(row["dt_star"])): row
        for row in rows
    }
    figure, axes = plt.subplots(3, 3, figsize=(14, 12), constrained_layout=True)
    for axis, (key, label, _) in zip(axes.flat, FEATURES, strict=False):
        matrix = np.full((len(torques), len(dt_stars)), np.nan)
        for row_index, torque in enumerate(torques):
            for column_index, dt_star in enumerate(dt_stars):
                row = by_cell.get((torque, dt_star))
                if row is not None:
                    matrix[row_index, column_index] = float(row[key])
        image = axis.imshow(matrix, origin="lower", aspect="auto", cmap="viridis")
        axis.set_title(label, fontsize=10)
        axis.set_xticks(range(len(dt_stars)), [f"{value:.0e}" for value in dt_stars])
        axis.set_yticks(range(len(torques)), [f"{value:.2g}" for value in torques])
        axis.set_xlabel("dt_star")
        axis.set_ylabel("torque [N m / flagellum]")
        for (row_index, column_index), value in np.ndenumerate(matrix):
            axis.text(
                column_index,
                row_index,
                "—" if not np.isfinite(value) else f"{value:.3g}",
                ha="center",
                va="center",
                color="white"
                if np.isfinite(value) and value > np.nanmean(matrix)
                else "black",
                fontsize=9,
            )
        figure.colorbar(image, ax=axis, shrink=0.8)
    axes.flat[-1].axis("off")
    figure.suptitle(title)
    figure.savefig(output_path, dpi=220)
    plt.close(figure)


def build_visuals(run_dir: Path, output_dir: Path) -> dict[str, Path]:
    """Write a feature table and heatmaps from a completed torque-dt campaign."""

    rows = feature_rows(run_dir)
    if not rows:
        raise ValueError(f"No conditions found in {run_dir / 'run_manifest.json'}")
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "torque_dt_feature_summary.csv"
    heatmap_path = output_dir / "torque_dt_feature_heatmaps.png"
    _write_csv(csv_path, rows)
    _heatmap(
        rows,
        output_path=heatmap_path,
        title="Issue #61: 2010 project torque-linked 1 tau diagnostics",
    )
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "kind": "phase2_2010_torque_dt_stability_visuals",
                "input_run_dir": str(run_dir),
                "outputs": {
                    "feature_summary_csv": str(csv_path),
                    "feature_heatmaps_png": str(heatmap_path),
                },
                "interpretation": "Qualitative diagnostics only; 1 tau does not establish long-time bundling or dt adoption.",
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return {"csv": csv_path, "heatmap": heatmap_path, "manifest": manifest_path}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args(argv)
    output_dir = args.output_dir or args.run_dir / "analysis" / "torque_dt_visuals"
    outputs = build_visuals(args.run_dir, output_dir)
    print(outputs["heatmap"])
