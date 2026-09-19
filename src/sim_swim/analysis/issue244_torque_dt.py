"""Combine and visualize the Issue #244 local torque--dt short screens.

The initial ``hook`` failure is deliberately interpreted here, rather than in
the simulator's strict generic gate.  This keeps the canonical QC contract
unchanged while making the approved short-screen exception explicit and
reproducible.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")

HOOK_LENGTH_LIMIT = 1.0
METRICS: tuple[tuple[str, str], ...] = (
    ("hook_angle_err_max_deg", "max hook angle error [deg]"),
    ("flag_bond_rel_err_max", "max flag bond relative error"),
    ("body_spring_max_stretch_ratio", "max body spring stretch ratio"),
    ("motor_force_balance_residual_ratio", "max motor force residual ratio"),
    ("motor_torque_balance_residual_ratio", "max motor torque residual ratio"),
    ("wall_time_s", "wall time [s]"),
    ("steps_per_s", "steps/s"),
)


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _max_metric(summary: dict[str, Any], name: str) -> float:
    value = dict(summary.get("all_step_metrics", {}) or {}).get(name, {})
    try:
        return float(dict(value or {}).get("max"))
    except (TypeError, ValueError):
        return float("nan")


def _bool(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"1", "true", "yes"}


def _initial_hook_warning(summary: dict[str, Any], *, dt_internal_s: float) -> bool:
    """Return true only for the documented one-step hook-angle transient."""

    gates = dict(summary.get("gates", {}) or {})
    finite = dict(gates.get("finite", {}) or {})
    nonbody = dict(gates.get("shape_nonbody", {}) or {})
    body = dict(gates.get("shape_body", {}) or {})
    snapshot = dict(nonbody.get("first_failure_snapshot", {}) or {})
    metrics = dict(snapshot.get("metrics", {}) or {})
    try:
        first_t_s = float(nonbody.get("first_observed_fail_t_s"))
        hook_angle = float(metrics.get("hook_angle_err_max_deg"))
        local_hook = float(metrics.get("local_attach_first_rel_err"))
        hook_length = float(metrics.get("hook_len_rel_err_max"))
    except (TypeError, ValueError):
        return False
    motor_metrics = (
        _max_metric(summary, "motor_force_balance_residual_ratio"),
        _max_metric(summary, "motor_torque_balance_residual_ratio"),
    )
    return (
        not bool(finite.get("any_fail", True))
        and not bool(body.get("any_fail", True))
        and _bool(nonbody.get("final_pass"))
        and bool(nonbody.get("any_fail", False))
        and str(nonbody.get("first_failure_category") or "") == "hook"
        and int(nonbody.get("observed_fail_sample_count") or 0) == 1
        and int(snapshot.get("step", -1)) == 0
        and math.isclose(first_t_s, dt_internal_s, rel_tol=1e-8, abs_tol=1e-15)
        and hook_angle > 30.0
        and local_hook <= HOOK_LENGTH_LIMIT
        and hook_length <= HOOK_LENGTH_LIMIT
        and all(math.isfinite(value) for value in motor_metrics)
    )


def _screen_status(summary: dict[str, Any], *, dt_internal_s: float) -> str:
    gates = dict(summary.get("gates", {}) or {})
    finite = dict(gates.get("finite", {}) or {})
    nonbody = dict(gates.get("shape_nonbody", {}) or {})
    body = dict(gates.get("shape_body", {}) or {})
    if _initial_hook_warning(summary, dt_internal_s=dt_internal_s):
        return "warning"
    if (
        not bool(finite.get("any_fail", True))
        and not bool(nonbody.get("any_fail", True))
        and not bool(body.get("any_fail", True))
    ):
        return "pass"
    return "fail"


def _row(record: dict[str, Any]) -> dict[str, Any]:
    output_dir = Path(str(record["output_dir"]))
    summary = _read_json(output_dir / "run_summary.json")
    performance = _read_json(output_dir / "performance.json")
    time = dict(record.get("time", {}) or {})
    axes = dict(record.get("axis_values", {}) or {})
    dt_star = float(time["dt_star"])
    values: dict[str, Any] = {
        "condition_id": record["condition_id"],
        "n_flagella": int(axes["n_flagella"]),
        "torque_Nm": float(axes["motor_torque"]),
        "dt_star": dt_star,
        "dt_internal_s": float(time["dt_internal_s"]),
        "screen_status": _screen_status(
            summary, dt_internal_s=float(time["dt_internal_s"])
        ),
        "raw_nonbody_any_fail": bool(
            dict(summary.get("gates", {}).get("shape_nonbody", {}) or {}).get(
                "any_fail", False
            )
        ),
        "first_failure_category": str(
            dict(summary.get("gates", {}).get("shape_nonbody", {}) or {}).get(
                "first_failure_category"
            )
            or ""
        ),
        "wall_time_s": float(performance["wall_time_s"]),
        "steps_per_s": float(performance["steps_per_s"]),
    }
    for metric, _ in METRICS[:5]:
        values[metric] = _max_metric(summary, metric)
    return values


def feature_rows(
    *, baseline_run_dir: Path, coarse_run_dir: Path
) -> list[dict[str, Any]]:
    """Load the two 10-condition runs and require one complete 2×5 grid per n."""

    rows: list[dict[str, Any]] = []
    for run_dir in (baseline_run_dir, coarse_run_dir):
        manifest = _read_json(run_dir / "run_manifest.json")
        for record in manifest.get("conditions", []):
            rows.append(_row(dict(record)))
    expected = {
        (n, torque, dt)
        for n in (1, 4)
        for torque in (1e-20, 2e-20, 2.5e-20, 3e-20, 3.5e-20)
        for dt in (1e-4, 1e-3)
    }
    observed = {
        (int(row["n_flagella"]), float(row["torque_Nm"]), float(row["dt_star"]))
        for row in rows
    }
    if observed != expected or len(rows) != len(expected):
        raise ValueError("Issue #244 torque-dt screen requires exactly 20 unique cells")
    return sorted(
        rows, key=lambda row: (row["n_flagella"], row["torque_Nm"], row["dt_star"])
    )


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _plot(rows: list[dict[str, Any]], *, n_flagella: int, output_path: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    torques = [1e-20, 2e-20, 2.5e-20, 3e-20, 3.5e-20]
    dts = [1e-4, 1e-3]
    by_cell = {(float(row["torque_Nm"]), float(row["dt_star"])): row for row in rows}
    panels = (("screen_status", "screen status"),) + METRICS
    figure, axes = plt.subplots(2, 4, figsize=(16, 8), constrained_layout=True)
    for axis, (metric, label) in zip(axes.flat, panels, strict=False):
        matrix = np.full((len(dts), len(torques)), np.nan)
        for row_index, dt in enumerate(dts):
            for column_index, torque in enumerate(torques):
                row = by_cell[(torque, dt)]
                value = row[metric]
                matrix[row_index, column_index] = (
                    {"fail": 0.0, "warning": 1.0, "pass": 2.0}[value]
                    if metric == "screen_status"
                    else float(value)
                )
        if metric == "screen_status":
            image = axis.imshow(
                matrix,
                aspect="auto",
                cmap=ListedColormap(["#c53d3d", "#e6b94f", "#278b6e"]),
                norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5], 3),
            )
            colorbar = figure.colorbar(image, ax=axis, shrink=0.8)
            colorbar.set_ticks([0, 1, 2], labels=["FAIL", "WARNING", "PASS"])
        else:
            image = axis.imshow(matrix, aspect="auto", cmap="viridis")
            figure.colorbar(image, ax=axis, shrink=0.8)
        axis.set_title(label, fontsize=10)
        axis.set_xticks(
            range(len(torques)), [f"{torque / 1e-20:g}" for torque in torques]
        )
        axis.set_yticks(range(len(dts)), [f"{dt:.0e}" for dt in dts])
        axis.set_xlabel("torque [1e-20 N m / flagellum]")
        axis.set_ylabel("dt_star")
        for (row_index, column_index), value in np.ndenumerate(matrix):
            text = (
                {0.0: "FAIL", 1.0: "WARN", 2.0: "PASS"}[value]
                if metric == "screen_status"
                else f"{value:.3g}"
            )
            axis.text(
                column_index, row_index, text, ha="center", va="center", fontsize=8
            )
    figure.suptitle(
        f"Issue #244: 2010 hex project 1tau torque–dt screen (n={n_flagella})"
    )
    figure.savefig(output_path, dpi=220)
    plt.close(figure)


def build_analysis(
    *, baseline_run_dir: Path, coarse_run_dir: Path, output_dir: Path
) -> dict[str, Path]:
    rows = feature_rows(
        baseline_run_dir=baseline_run_dir, coarse_run_dir=coarse_run_dir
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "issue244_torque_dt_summary.csv"
    _write_csv(csv_path, rows)
    outputs = {"summary_csv": csv_path}
    for n_flagella in (1, 4):
        path = output_dir / f"issue244_torque_dt_heatmap_nf{n_flagella:02d}.png"
        _plot(
            [row for row in rows if row["n_flagella"] == n_flagella],
            n_flagella=n_flagella,
            output_path=path,
        )
        outputs[f"heatmap_nf{n_flagella:02d}"] = path
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "kind": "issue244_torque_dt_analysis",
                "baseline_run_dir": str(baseline_run_dir),
                "coarse_run_dir": str(coarse_run_dir),
                "condition_count": len(rows),
                "outputs": {key: str(value) for key, value in outputs.items()},
                "warning_policy": "one-step initial hook-angle-only transient; strict simulator gates unchanged",
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    outputs["manifest"] = manifest_path
    return outputs


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-run-dir", type=Path, required=True)
    parser.add_argument("--coarse-run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    outputs = build_analysis(**vars(args))
    for path in outputs.values():
        print(path)


if __name__ == "__main__":
    main()
