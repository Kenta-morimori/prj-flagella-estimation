"""Assemble Issue #199 Task D 2015 tau-linked torque x dt diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from sim_swim.analysis.stage_a_2015_analysis import (
    _observed_metrics,
    load_threshold_contract,
)

DT_ORDER = (1.0e-3, 1.0e-4)
TORQUE_ORDER_NM = (1.0e-19, 2.5e-20, 1.0e-21)


def _float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"summary.csv is empty: {path}")
    return rows


def _truth(value: Any) -> bool:
    return str(value).strip().lower() == "true"


def _close(actual: Any, expected: float) -> bool:
    return math.isclose(_float(actual), expected, rel_tol=1.0e-12, abs_tol=1.0e-15)


def _grid_key(dt_star: Any, torque_Nm: Any) -> tuple[float, float]:
    """Canonical grid key, avoiding YAML/float representation noise."""
    return (round(_float(dt_star), 10), round(_float(torque_Nm), 30))


def _row_torque_Nm(row: dict[str, Any]) -> float:
    """Read new physical-torque summaries and pre-migration summaries alike."""
    return _float(row.get("motor_torque_Nm") or row.get("torque_Nm"))


def _load_campaign(
    run_root: Path, expected_dt: float
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    manifest = _read_json(run_root / "run_manifest.json")
    if manifest.get("kind") != "stage_a_2015" or manifest.get("stage") != "motor_on":
        raise ValueError(f"Not a 2015 Stage A motor-on root: {run_root}")
    for name, expected in (("dt_star", expected_dt), ("duration_tau", 1.0)):
        if not _close(manifest.get(name), expected):
            raise ValueError(
                f"{run_root}: {name}={manifest.get(name)!r}; expected {expected}"
            )
    if manifest.get("profiles") not in (None, ["project"]):
        raise ValueError(f"{run_root}: Task D requires only the project profile")
    if manifest.get("link_reference_torque") is not True:
        raise ValueError(f"{run_root}: link_reference_torque must be true")
    rows = _read_rows(run_root / "summary.csv")
    torques = {_float(row.get("motor_torque_Nm")) for row in rows}
    if torques != set(TORQUE_ORDER_NM):
        raise ValueError(
            f"{run_root}: expected torques {TORQUE_ORDER_NM}, got {sorted(torques)}"
        )
    return manifest, rows


def _condition_by_torque(manifest: dict[str, Any], torque_Nm: float) -> dict[str, Any]:
    matches = [
        item
        for item in manifest.get("conditions", [])
        if _close(item.get("motor_torque_Nm"), torque_Nm)
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one condition manifest at torque {torque_Nm}")
    condition = matches[0]
    overrides = condition.get("config_overrides", {})
    if overrides.get("time.scale_policy") != "reference_torque":
        raise ValueError(
            f"{condition.get('condition_id')}: missing reference_torque policy"
        )
    if not _close(
        overrides.get("motor.torque_Nm"), overrides.get("motor.reference_torque_Nm")
    ):
        raise ValueError(
            f"{condition.get('condition_id')}: motor/reference torque differ"
        )
    return condition


def _safety(row: dict[str, str], thresholds: dict[str, Any]) -> tuple[bool, list[str]]:
    failures: list[str] = []
    for field in ("status", "completion_pass", "finite_pass_all"):
        expected = "completed" if field == "status" else "true"
        if str(row.get(field, "")).strip().lower() != expected:
            failures.append(f"{field}={row.get(field, '')}")
    for metric, limit_raw in thresholds.items():
        value, limit = _float(row.get(metric)), _float(limit_raw)
        if not math.isfinite(value):
            failures.append(f"{metric}=nonfinite")
        elif value > limit:
            failures.append(f"{metric}={value:.4g}>{limit:.4g}")
    return not failures, failures


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _plot_heatmaps(rows: list[dict[str, Any]], path: Path) -> None:
    dt_order = tuple(sorted({_float(row["dt_star"]) for row in rows}, reverse=True))
    torque_order = tuple(sorted({_row_torque_Nm(row) for row in rows}, reverse=True))
    if len(dt_order) != 2 or len(torque_order) != 3:
        raise ValueError("Task D heatmap requires a 3-torque by 2-dt grid")
    lookup = {_grid_key(row["dt_star"], _row_torque_Nm(row)): row for row in rows}
    # Match the common Issue #199 heatmap layout: diagnostic QC first, then
    # shape/helix quantities, and a compact reading guide in the last panel.
    panels: list[tuple[str, str, float, bool]] = [
        ("safety_pass", "required safety QC (PASS=1)", 1.0, True),
        (
            "body_spring_max_stretch_ratio",
            "max body spring relative error",
            0.10,
            False,
        ),
        ("max_flag_bond_rel_err", "max flag bond relative error", 0.10, False),
        ("max_hook_len_rel_err", "max hook length relative error", 0.10, False),
        (
            "max_flag_helix_radius_abs_err_over_b",
            "max flag helix radius error / b",
            0.05,
            False,
        ),
        (
            "max_flag_helix_pitch_rel_err",
            "max flag helix pitch relative error",
            0.10,
            False,
        ),
        ("max_flag_torsion_err_deg", "max flag torsion error [deg]", 30.0, False),
        (
            "max_motor_torque_balance_residual_ratio",
            "max motor torque-balance residual ratio",
            0.10,
            False,
        ),
    ]
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    for ax, (field, title, cap, is_bool) in zip(axes.flat, panels):
        values = np.array(
            [
                [_float(lookup[_grid_key(dt, torque_Nm)].get(field)) for dt in dt_order]
                for torque_Nm in torque_order
            ]
        )
        if is_bool:
            values = values.astype(float)
            image = ax.imshow(values, vmin=0, vmax=1, cmap="RdYlGn")
            labels = [
                ["PASS" if value >= 0.5 else "FAIL" for value in line]
                for line in values
            ]
        else:
            image = ax.imshow(
                np.clip(values, 0.0, cap), vmin=0, vmax=cap, cmap="viridis"
            )
            labels = [
                [f"{value:.3g}" if math.isfinite(value) else "nan" for value in line]
                for line in values
            ]
        for i in range(3):
            for j in range(2):
                ax.text(j, i, labels[i][j], ha="center", va="center", fontsize=9)
        ax.set(
            title=title,
            xticks=range(2),
            xticklabels=[f"{dt:.0e}" for dt in dt_order],
            yticks=range(3),
            yticklabels=[f"{torque:.1e}" for torque in torque_order],
            xlabel="dt_star",
            ylabel="torque [N m / flagellum]",
        )
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    guide = axes.flat[8]
    guide.axis("off")
    guide.text(
        0.0,
        1.0,
        "How to read\n\n"
        "Rows: torque [N m / flagellum]\n"
        "Columns: dimensionless dt_star\n\n"
        "PASS is the conjunction of the locked Stage A safety checks. "
        "Other panels are full-run extrema; colors are clipped at the "
        "display cap while annotations retain the raw value.\n\n"
        "This is a τ-linked, 1τ diagnostic screen. It does not adopt "
        "dt_star=1e-3 or promote the 2015 project profile.",
        ha="left",
        va="top",
        fontsize=11,
        wrap=True,
    )
    fig.suptitle("Issue #199 Task D: 2015 project τ-linked torque × dt diagnostics")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def assemble(
    *,
    reference_run: Path,
    candidate_run: Path,
    threshold_contract: Path,
    output_dir: Path,
) -> Path:
    contract = load_threshold_contract(threshold_contract)
    thresholds = contract.get("thresholds", {})
    ref_manifest, ref_rows = _load_campaign(reference_run, 1.0e-4)
    candidate_manifest, candidate_rows = _load_campaign(candidate_run, 1.0e-3)
    output_dir.mkdir(parents=True, exist_ok=False)
    all_rows: list[dict[str, Any]] = []
    replay_conditions: list[dict[str, Any]] = []
    campaigns = (
        (1.0e-3, candidate_manifest, candidate_rows, "candidate"),
        (1.0e-4, ref_manifest, ref_rows, "reference"),
    )
    for grid_row, torque_Nm in enumerate(TORQUE_ORDER_NM):
        for grid_col, (dt, manifest, source_rows, role) in enumerate(campaigns):
            by_torque = {_float(row.get("motor_torque_Nm")): row for row in source_rows}
            source = dict(by_torque[torque_Nm])
            condition = _condition_by_torque(manifest, torque_Nm)
            # Stage A's summary has maxima useful for a compact display, but
            # the locked contract also includes body drifts.  Recompute those
            # from each condition's bounded diagnostics before screening.
            source.update(_observed_metrics(Path(str(condition["output_dir"]))))
            passed, failures = _safety(source, thresholds)
            replay_id = f"torque_{torque_Nm:.0e}_dt{dt:.0e}".replace("-", "m")
            source.update(
                {
                    "condition_id": replay_id,
                    "condition_label": f"torque={torque_Nm:.0e} N m / flagellum, dt*={dt:.0e}",
                    "grid_row_index": grid_row,
                    "grid_col_index": grid_col,
                    "dt_role": role,
                    "safety_pass": int(passed),
                    "safety_failures": "; ".join(failures),
                    "source_run": str(
                        reference_run if role == "reference" else candidate_run
                    ),
                }
            )
            all_rows.append(source)
            replay_conditions.append(
                {
                    **condition,
                    "condition_id": replay_id,
                    "condition_label": source["condition_label"],
                    "grid_row_index": grid_row,
                    "grid_col_index": grid_col,
                }
            )
    _write_csv(output_dir / "task_d_summary.csv", all_rows)
    _plot_heatmaps(all_rows, output_dir / "task_d_feature_heatmaps.png")
    replay_dir = output_dir / "replay_input"
    replay_dir.mkdir()
    _write_csv(replay_dir / "summary.csv", all_rows)
    (replay_dir / "run_manifest.json").write_text(
        json.dumps(
            {
                "kind": "issue199_task_d_2015_replay",
                "base_config": "conf/sim_swim_2015.yaml",
                "conditions": replay_conditions,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    decision = {
        "kind": "issue199_task_d_2015_safety_screen",
        "status": "diagnostic_only",
        "reference_run": str(reference_run),
        "candidate_run": str(candidate_run),
        "threshold_contract": str(threshold_contract),
        "conditions": len(all_rows),
        "safety_pass_count": sum(row["safety_pass"] for row in all_rows),
        "adoption": "No dt*=1e-3 adoption or 2015 supported-profile promotion is implied.",
    }
    (output_dir / "task_d_decision.json").write_text(
        json.dumps(decision, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return output_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-run", type=Path, required=True)
    parser.add_argument("--candidate-run", type=Path, required=True)
    parser.add_argument("--threshold-contract", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(
        assemble(
            reference_run=args.reference_run,
            candidate_run=args.candidate_run,
            threshold_contract=args.threshold_contract,
            output_dir=args.output_dir,
        )
    )
