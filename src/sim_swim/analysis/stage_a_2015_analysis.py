"""Issue #168 Stage A threshold proposal and result evaluation."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import yaml

THRESHOLD_FACTOR = 5.0
THRESHOLD_POLICY = {
    "body_spring_max_stretch_ratio": {"floor": 0.01, "cap": 0.08},
    "body_length_rel_drift_max": {"floor": 0.01, "cap": 0.10},
    "body_width_rel_drift_max": {"floor": 0.01, "cap": 0.10},
    "body_cross_section_area_rel_drift_max": {"floor": 0.02, "cap": 0.20},
    "max_flag_bond_rel_err": {"floor": 0.01, "cap": 0.08},
    "max_hook_len_rel_err": {"floor": 0.01, "cap": 0.08},
    "max_hook_angle_err_deg": {"floor": 5.0, "cap": 30.0},
    "max_flag_bend_err_deg": {"floor": 5.0, "cap": 30.0},
    "max_flag_torsion_err_deg": {"floor": 10.0, "cap": 60.0},
    "max_flag_helix_radius_abs_err_over_b": {"floor": 0.01, "cap": 0.035},
    "max_flag_helix_pitch_rel_err": {"floor": 0.02, "cap": 0.05},
}


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"CSV is empty: {path}")
    return rows


def _float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _max_field(rows: list[dict[str, str]], field: str) -> float:
    values = [_float(row.get(field)) for row in rows]
    finite = [value for value in values if math.isfinite(value)]
    return max(finite) if finite else float("nan")


def _relative_drift(rows: list[dict[str, str]], fields: tuple[str, ...]) -> float:
    drifts: list[float] = []
    for field in fields:
        values = [_float(row.get(field)) for row in rows]
        if not values or not all(math.isfinite(value) for value in values):
            return float("nan")
        initial = values[0]
        drifts.extend(
            abs(value - initial) / max(abs(initial), 1.0e-30) for value in values
        )
    return max(drifts) if drifts else float("nan")


def _observed_metrics(condition_dir: Path) -> dict[str, float]:
    step_rows = _read_csv(condition_dir / "step_summary.csv")
    body_rows = _read_csv(condition_dir / "body_constraint_diagnostics.csv")
    return {
        "body_spring_max_stretch_ratio": _max_field(
            body_rows, "body_spring_max_stretch_ratio"
        ),
        "body_length_rel_drift_max": _relative_drift(body_rows, ("body_length_um",)),
        "body_width_rel_drift_max": _relative_drift(
            body_rows,
            ("body_width_mean_um", "body_width_min_um", "body_width_max_um"),
        ),
        "body_cross_section_area_rel_drift_max": _relative_drift(
            body_rows,
            (
                "body_cross_section_area_min_um2",
                "body_cross_section_area_max_um2",
            ),
        ),
        "max_flag_bond_rel_err": _max_field(step_rows, "flag_bond_rel_err_max"),
        "max_hook_len_rel_err": _max_field(step_rows, "hook_len_rel_err_max"),
        "max_hook_angle_err_deg": _max_field(step_rows, "hook_angle_err_max_deg"),
        "max_flag_bend_err_deg": _max_field(step_rows, "flag_bend_err_max_deg"),
        "max_flag_torsion_err_deg": _max_field(step_rows, "flag_torsion_err_max_deg"),
        "max_flag_helix_radius_abs_err_over_b": _max_field(
            step_rows, "flag_helix_radius_abs_err_over_b_max"
        ),
        "max_flag_helix_pitch_rel_err": _max_field(
            step_rows, "flag_helix_pitch_rel_err_max"
        ),
        "max_net_force_residual_ratio": _max_field(
            step_rows, "net_force_residual_ratio"
        ),
        "max_net_torque_residual_ratio": _max_field(
            step_rows, "net_torque_residual_ratio"
        ),
        "max_motor_force_body_norm_N": _max_field(
            step_rows, "motor_net_force_body_norm_N"
        ),
        "max_motor_force_flag_norm_N": _max_field(
            step_rows, "motor_net_force_flag_norm_N"
        ),
    }


def _summary_by_profile(run_root: Path) -> dict[str, dict[str, str]]:
    rows = _read_csv(run_root / "summary.csv")
    by_profile = {str(row.get("profile", "")): row for row in rows}
    missing = sorted({"project", "paper"} - set(by_profile))
    if missing:
        raise ValueError("Missing profile rows: " + ", ".join(missing))
    return by_profile


def propose_thresholds(motor_off_run: Path) -> dict[str, Any]:
    manifest_path = motor_off_run / "run_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("smoke_steps") is not None:
            raise ValueError("Smoke run cannot be used to propose Stage A thresholds")
    summary = _summary_by_profile(motor_off_run)
    observed = {
        profile: _observed_metrics(motor_off_run / profile)
        for profile in ("project", "paper")
    }
    failures: list[str] = []
    for profile, row in summary.items():
        if str(row.get("status")) != "completed":
            failures.append(f"{profile}: status={row.get('status')}")
        if str(row.get("completion_pass", "")).lower() != "true":
            failures.append(f"{profile}: completion_pass is not true")
        if str(row.get("finite_pass_all", "")).lower() != "true":
            failures.append(f"{profile}: finite_pass_all is not true")
        for force_field in (
            "max_motor_force_body_norm_N",
            "max_motor_force_flag_norm_N",
        ):
            value = observed[profile][force_field]
            if not math.isfinite(value) or value > 1.0e-30:
                failures.append(f"{profile}: motor-off {force_field}={value:.6g}")

    thresholds: dict[str, float] = {}
    pilot_maxima: dict[str, float] = {}
    for metric, policy in THRESHOLD_POLICY.items():
        values = [observed[profile][metric] for profile in ("project", "paper")]
        if not all(math.isfinite(value) for value in values):
            failures.append(f"{metric}: pilot value is non-finite")
            continue
        pilot_max = max(values)
        pilot_maxima[metric] = pilot_max
        if pilot_max > policy["cap"]:
            failures.append(
                f"{metric}: pilot max {pilot_max:.6g} exceeds cap {policy['cap']:.6g}"
            )
        thresholds[metric] = min(
            policy["cap"], max(policy["floor"], THRESHOLD_FACTOR * pilot_max)
        )

    off_rotation_drift = max(
        _float(
            summary[profile].get("axis_center_body_relative_net_abs_revolutions_max")
        )
        for profile in ("project", "paper")
    )
    if not math.isfinite(off_rotation_drift):
        failures.append("motor-off body-relative rotation drift is non-finite")
    thresholds["min_axis_center_body_relative_revolutions"] = max(
        0.01,
        10.0 * off_rotation_drift if math.isfinite(off_rotation_drift) else 0.01,
    )
    thresholds["max_motor_force_balance_residual_ratio"] = 1.0e-8
    thresholds["max_motor_torque_balance_residual_ratio"] = 1.0e-8

    return {
        "kind": "stage_a_2015_threshold_contract",
        "issue": 168,
        "status": "proposed" if not failures else "pilot_failed",
        "source_motor_off_run": str(motor_off_run),
        "policy": {
            "factor": THRESHOLD_FACTOR,
            "metrics": THRESHOLD_POLICY,
        },
        "observed": observed,
        "pilot_maxima": pilot_maxima,
        "thresholds": thresholds,
        "failures": failures,
        "next_action": "lock_thresholds" if not failures else "diagnose_motor_off",
    }


def load_threshold_contract(path: Path) -> dict[str, Any]:
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise ValueError(f"Threshold contract must be a mapping: {path}")
    if raw.get("status") != "locked":
        raise ValueError(f"Threshold contract is not locked: {path}")
    return raw


def evaluate_motor_on(motor_on_run: Path, contract: dict[str, Any]) -> dict[str, Any]:
    summary = _summary_by_profile(motor_on_run)
    thresholds = dict(contract.get("thresholds", {}) or {})
    evaluations: dict[str, Any] = {}
    brace_profiles: list[str] = []
    for profile in ("project", "paper"):
        row = summary[profile]
        metrics = _observed_metrics(motor_on_run / profile)
        metrics.update(
            {
                "max_motor_force_balance_residual_ratio": _float(
                    row.get("max_motor_force_balance_residual_ratio")
                ),
                "max_motor_torque_balance_residual_ratio": _float(
                    row.get("max_motor_torque_balance_residual_ratio")
                ),
            }
        )
        reasons: list[str] = []
        if str(row.get("status")) != "completed":
            reasons.append(f"status={row.get('status')}")
        if str(row.get("completion_pass", "")).lower() != "true":
            reasons.append("completion_pass is not true")
        if str(row.get("finite_pass_all", "")).lower() != "true":
            reasons.append("finite_pass_all is not true")
        for metric, limit_raw in thresholds.items():
            if metric.startswith("min_") or metric not in metrics:
                continue
            value = metrics[metric]
            limit = float(limit_raw)
            if not math.isfinite(value) or value > limit:
                reasons.append(f"{metric}={value:.6g} > {limit:.6g}")
        rotation = _float(row.get("axis_center_body_relative_net_abs_revolutions_min"))
        rotation_limit = float(thresholds["min_axis_center_body_relative_revolutions"])
        if not math.isfinite(rotation) or rotation < rotation_limit:
            reasons.append(
                "axis_center_body_relative_net_abs_revolutions_min="
                f"{rotation:.6g} < {rotation_limit:.6g}"
            )
        body_failed = str(row.get("body_shape_pass", "")).lower() != "true"
        if body_failed:
            brace_profiles.append(profile)
        evaluations[profile] = {
            "automated_pass": not reasons,
            "reasons": reasons,
            "observed": metrics,
            "body_shape_pass": not body_failed,
            "visual_review_required": True,
        }
    return {
        "kind": "stage_a_2015_decision",
        "issue": 168,
        "source_motor_on_run": str(motor_on_run),
        "threshold_contract_source": contract.get("source_motor_off_run"),
        "profiles": evaluations,
        "brace_comparison_profiles": brace_profiles,
        "next_action": (
            "run_brace_comparison" if brace_profiles else "perform_visual_review"
        ),
    }


def write_markdown(path: Path, result: dict[str, Any]) -> None:
    lines = ["# Phase 2 Issue #168 Stage A analysis", ""]
    lines.append(f"- status: `{result.get('status', result.get('next_action'))}`")
    if "thresholds" in result:
        lines.extend(["", "## Proposed thresholds", ""])
        for metric, value in result["thresholds"].items():
            lines.append(f"- `{metric}`: `{float(value):.8g}`")
    if result.get("failures"):
        lines.extend(["", "## Failures", ""])
        lines.extend(f"- {reason}" for reason in result["failures"])
    if "profiles" in result:
        lines.extend(["", "## Motor-on evaluation", ""])
        for profile, evaluation in result["profiles"].items():
            lines.append(
                f"- `{profile}`: automated_pass=`{evaluation['automated_pass']}`, "
                "visual_review_required=`True`"
            )
            lines.extend(f"  - {reason}" for reason in evaluation["reasons"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_analysis(
    output_dir: Path,
    *,
    motor_off_run: Path,
    motor_on_run: Path | None = None,
    threshold_contract: Path | None = None,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    proposal = propose_thresholds(motor_off_run)
    proposal_json = output_dir / "threshold_proposal.json"
    proposal_md = output_dir / "threshold_proposal.md"
    proposal_json.write_text(
        json.dumps(proposal, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    write_markdown(proposal_md, proposal)
    outputs = {
        "threshold_proposal_json": proposal_json,
        "threshold_proposal_markdown": proposal_md,
    }
    if motor_on_run is not None:
        if threshold_contract is None:
            raise ValueError("threshold_contract is required for motor-on evaluation")
        locked = load_threshold_contract(threshold_contract)
        decision = evaluate_motor_on(motor_on_run, locked)
        decision_json = output_dir / "stage_a_decision.json"
        decision_md = output_dir / "stage_a_decision.md"
        decision_json.write_text(
            json.dumps(decision, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        write_markdown(decision_md, decision)
        outputs.update(
            {
                "stage_a_decision_json": decision_json,
                "stage_a_decision_markdown": decision_md,
            }
        )
    return outputs
