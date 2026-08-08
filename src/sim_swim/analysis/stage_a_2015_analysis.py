"""Issue #168 Stage A threshold proposal and result evaluation."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from sim_swim.analysis.flagella_count_behavior import load_state_archive

THRESHOLD_FACTOR = 5.0
DT_SENSITIVITY_LIMITS = {
    "body_relative_rotation_rel_difference_max": 0.10,
    "body_orientation_difference_deg_max": 5.0,
    "centered_bead_rms_difference_over_b_max": 0.10,
    "centered_bead_max_difference_over_b_max": 0.25,
}
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


def _load_run_manifest(run_root: Path) -> dict[str, Any]:
    path = run_root / "run_manifest.json"
    raw = json.loads(path.read_text(encoding="utf-8"))
    if raw.get("smoke_steps") is not None:
        raise ValueError(f"Smoke run cannot be used for dt comparison: {run_root}")
    return raw


def _validate_run_contract(
    run_root: Path,
    *,
    stage: str,
    dt_star: float,
    duration_tau: float,
) -> dict[str, Any]:
    manifest = _load_run_manifest(run_root)
    expected = {
        "stage": stage,
        "dt_star": dt_star,
        "duration_tau": duration_tau,
    }
    for field, expected_value in expected.items():
        actual = manifest.get(field)
        if isinstance(expected_value, float):
            valid = math.isclose(_float(actual), expected_value, rel_tol=1.0e-12)
        else:
            valid = actual == expected_value
        if not valid:
            raise ValueError(
                f"Unexpected {field} for {run_root}: {actual!r}; "
                f"expected {expected_value!r}"
            )
    return manifest


def _condition_manifest(manifest: dict[str, Any], profile: str) -> dict[str, Any]:
    for condition in manifest.get("conditions", []):
        if condition.get("condition_id") == profile:
            return condition
    raise ValueError(f"Missing condition manifest for profile={profile}")


def _quaternion_difference_deg(first: Any, second: Any) -> float:
    q1 = np.asarray(first, dtype=float)
    q2 = np.asarray(second, dtype=float)
    q1 /= max(float(np.linalg.norm(q1)), 1.0e-30)
    q2 /= max(float(np.linalg.norm(q2)), 1.0e-30)
    dot = float(np.clip(abs(np.dot(q1, q2)), 0.0, 1.0))
    return math.degrees(2.0 * math.acos(dot))


def _signed_body_relative_revolutions(
    path: Path, *, max_t_s: float | None = None
) -> dict[int, float]:
    rows = _read_csv(path)
    phases: dict[int, list[float]] = {}
    for row in rows:
        t_s = _float(row.get("t_s"))
        if max_t_s is not None and t_s > max_t_s + 1.0e-15:
            continue
        try:
            flag_id = int(row["flag_id"])
        except (KeyError, TypeError, ValueError):
            continue
        phase = _float(row.get("axis_center_body_relative_phase_deg"))
        if math.isfinite(phase):
            phases.setdefault(flag_id, []).append(phase)
    result: dict[int, float] = {}
    for flag_id, values in phases.items():
        if len(values) < 2:
            continue
        unwrapped = np.unwrap(np.radians(np.asarray(values, dtype=float)))
        result[flag_id] = float((unwrapped[-1] - unwrapped[0]) / (2.0 * math.pi))
    return result


def _paired_profile_comparison(
    *,
    profile: str,
    coarse_run: Path,
    fine_run: Path,
    coarse_manifest: dict[str, Any],
    fine_manifest: dict[str, Any],
) -> dict[str, Any]:
    coarse_states = load_state_archive(coarse_run / profile / "state_archive.npz")
    fine_states = load_state_archive(fine_run / profile / "state_archive.npz")
    if not coarse_states or not fine_states:
        raise ValueError(f"Missing states for profile={profile}")
    fine_final = fine_states[-1]
    coarse_final = min(coarse_states, key=lambda state: abs(state.t - fine_final.t))
    time_difference_s = abs(coarse_final.t - fine_final.t)
    coarse_interval_s = min(
        (
            abs(state.t - coarse_final.t)
            for state in coarse_states
            if state is not coarse_final
        ),
        default=0.0,
    )
    if time_difference_s > max(1.0e-15, 0.5 * coarse_interval_s):
        raise ValueError(
            f"No time-aligned coarse state for profile={profile}: "
            f"difference={time_difference_s:.6g} s"
        )

    coarse_condition = _condition_manifest(coarse_manifest, profile)
    fine_condition = _condition_manifest(fine_manifest, profile)
    coarse_scales = coarse_condition.get("comparison_scales", {})
    fine_scales = fine_condition.get("comparison_scales", {})
    if coarse_scales != fine_scales:
        raise ValueError(f"Comparison scales differ for profile={profile}")
    b_um = _float(coarse_scales.get("b_um"))
    body_beads = int(coarse_scales.get("body_beads", 0))
    if not math.isfinite(b_um) or b_um <= 0.0 or body_beads <= 0:
        raise ValueError(f"Invalid comparison scales for profile={profile}")

    coarse_positions = np.asarray(coarse_final.bead_positions_um, dtype=float)
    fine_positions = np.asarray(fine_final.bead_positions_um, dtype=float)
    if coarse_positions.shape != fine_positions.shape:
        raise ValueError(f"Bead shapes differ for profile={profile}")
    coarse_centered = coarse_positions - coarse_positions[:body_beads].mean(axis=0)
    fine_centered = fine_positions - fine_positions[:body_beads].mean(axis=0)
    bead_difference = np.linalg.norm(coarse_centered - fine_centered, axis=1) / b_um
    rms_difference = float(np.sqrt(np.mean(np.square(bead_difference))))
    max_difference = float(np.max(bead_difference))
    orientation_difference = _quaternion_difference_deg(
        coarse_final.quaternion, fine_final.quaternion
    )

    fine_axis_path = fine_run / profile / "flag_helix_axis_diagnostics.csv"
    coarse_axis_path = coarse_run / profile / "flag_helix_axis_diagnostics.csv"
    fine_rows = _read_csv(fine_axis_path)
    fine_last_t_s = max(_float(row.get("t_s")) for row in fine_rows)
    coarse_revolutions = _signed_body_relative_revolutions(
        coarse_axis_path, max_t_s=fine_last_t_s
    )
    fine_revolutions = _signed_body_relative_revolutions(fine_axis_path)
    if set(coarse_revolutions) != set(fine_revolutions) or not fine_revolutions:
        raise ValueError(f"Flag rotation series differ for profile={profile}")
    rotation_rows: list[dict[str, Any]] = []
    rotation_pass = True
    for flag_id in sorted(fine_revolutions):
        coarse_value = coarse_revolutions[flag_id]
        fine_value = fine_revolutions[flag_id]
        same_direction = (
            abs(coarse_value) <= 1.0e-12 and abs(fine_value) <= 1.0e-12
        ) or coarse_value * fine_value > 0.0
        relative_difference = abs(abs(coarse_value) - abs(fine_value)) / max(
            abs(coarse_value), abs(fine_value), 1.0e-12
        )
        flag_pass = (
            same_direction
            and relative_difference
            <= DT_SENSITIVITY_LIMITS["body_relative_rotation_rel_difference_max"]
        )
        rotation_pass = rotation_pass and flag_pass
        rotation_rows.append(
            {
                "flag_id": flag_id,
                "coarse_revolutions": coarse_value,
                "fine_revolutions": fine_value,
                "same_direction": same_direction,
                "relative_difference": relative_difference,
                "pass": flag_pass,
            }
        )

    reasons: list[str] = []
    if not rotation_pass:
        reasons.append("body-relative rotation differs by direction or more than 10%")
    if (
        orientation_difference
        > DT_SENSITIVITY_LIMITS["body_orientation_difference_deg_max"]
    ):
        reasons.append(f"body orientation difference={orientation_difference:.6g} deg")
    if (
        rms_difference
        > DT_SENSITIVITY_LIMITS["centered_bead_rms_difference_over_b_max"]
    ):
        reasons.append(f"centered bead RMS difference={rms_difference:.6g} b")
    if (
        max_difference
        > DT_SENSITIVITY_LIMITS["centered_bead_max_difference_over_b_max"]
    ):
        reasons.append(f"centered bead max difference={max_difference:.6g} b")
    return {
        "profile": profile,
        "comparison_t_s": fine_final.t,
        "time_alignment_difference_s": time_difference_s,
        "body_orientation_difference_deg": orientation_difference,
        "centered_bead_rms_difference_over_b": rms_difference,
        "centered_bead_max_difference_over_b": max_difference,
        "flag_rotation": rotation_rows,
        "paired_pass": not reasons,
        "reasons": reasons,
    }


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


def _evaluate_short_motor_on(
    motor_on_run: Path, contract: dict[str, Any]
) -> dict[str, Any]:
    summary = _summary_by_profile(motor_on_run)
    thresholds = dict(contract.get("thresholds", {}) or {})
    evaluations: dict[str, Any] = {}
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
        evaluations[profile] = {
            "automated_pass": not reasons,
            "reasons": reasons,
            "observed": metrics,
        }
    return evaluations


def evaluate_dt_sensitivity(
    *,
    baseline_motor_off_run: Path,
    coarse_motor_off_run: Path,
    coarse_motor_on_run: Path,
    fine_motor_on_short_run: Path,
    contract: dict[str, Any],
) -> dict[str, Any]:
    baseline_off_manifest = _validate_run_contract(
        baseline_motor_off_run, stage="motor_off", dt_star=1.0e-5, duration_tau=0.1
    )
    _validate_run_contract(
        coarse_motor_off_run, stage="motor_off", dt_star=1.0e-4, duration_tau=0.1
    )
    coarse_on_manifest = _validate_run_contract(
        coarse_motor_on_run, stage="motor_on", dt_star=1.0e-4, duration_tau=1.0
    )
    fine_on_manifest = _validate_run_contract(
        fine_motor_on_short_run,
        stage="motor_on",
        dt_star=1.0e-5,
        duration_tau=0.1,
    )
    del baseline_off_manifest  # Contract validation is the required baseline check.

    coarse_off_proposal = propose_thresholds(coarse_motor_off_run)
    coarse_on_evaluation = evaluate_motor_on(coarse_motor_on_run, contract)
    fine_short_evaluation = _evaluate_short_motor_on(fine_motor_on_short_run, contract)
    baseline_off_metrics = {
        profile: _observed_metrics(baseline_motor_off_run / profile)
        for profile in ("project", "paper")
    }
    coarse_off_metrics = {
        profile: _observed_metrics(coarse_motor_off_run / profile)
        for profile in ("project", "paper")
    }
    off_comparison: dict[str, Any] = {}
    for profile in ("project", "paper"):
        metrics: dict[str, Any] = {}
        for metric, fine_value in baseline_off_metrics[profile].items():
            coarse_value = coarse_off_metrics[profile][metric]
            metrics[metric] = {
                "dt1e5": fine_value,
                "dt1e4": coarse_value,
                "absolute_difference": abs(coarse_value - fine_value),
                "ratio_dt1e4_over_dt1e5": (
                    coarse_value / fine_value if abs(fine_value) > 1.0e-30 else None
                ),
            }
        off_comparison[profile] = metrics

    paired = {
        profile: _paired_profile_comparison(
            profile=profile,
            coarse_run=coarse_motor_on_run,
            fine_run=fine_motor_on_short_run,
            coarse_manifest=coarse_on_manifest,
            fine_manifest=fine_on_manifest,
        )
        for profile in ("project", "paper")
    }
    speedups: dict[str, Any] = {}
    baseline_summary = _summary_by_profile(baseline_motor_off_run)
    coarse_off_summary = _summary_by_profile(coarse_motor_off_run)
    coarse_on_summary = _summary_by_profile(coarse_motor_on_run)
    fine_on_summary = _summary_by_profile(fine_motor_on_short_run)
    for profile in ("project", "paper"):
        speedups[profile] = {
            "motor_off_steps_per_s_ratio": _float(
                coarse_off_summary[profile].get("steps_per_s")
            )
            / max(_float(baseline_summary[profile].get("steps_per_s")), 1.0e-30),
            "motor_on_steps_per_s_ratio": _float(
                coarse_on_summary[profile].get("steps_per_s")
            )
            / max(_float(fine_on_summary[profile].get("steps_per_s")), 1.0e-30),
        }

    reference_stable = (
        coarse_off_proposal["status"] == "proposed"
        and all(
            item["automated_pass"] for item in coarse_on_evaluation["profiles"].values()
        )
        and all(item["automated_pass"] for item in fine_short_evaluation.values())
    )
    adoption_candidate = reference_stable and all(
        item["paired_pass"] for item in paired.values()
    )
    return {
        "kind": "stage_a_2015_dt_sensitivity_decision",
        "issue": 168,
        "status": "adoption_candidate" if adoption_candidate else "reference_only",
        "limits": DT_SENSITIVITY_LIMITS,
        "runs": {
            "baseline_motor_off": str(baseline_motor_off_run),
            "coarse_motor_off": str(coarse_motor_off_run),
            "coarse_motor_on": str(coarse_motor_on_run),
            "fine_motor_on_short": str(fine_motor_on_short_run),
        },
        "reference_stable": reference_stable,
        "adoption_candidate": adoption_candidate,
        "coarse_motor_off": coarse_off_proposal,
        "coarse_motor_on": coarse_on_evaluation,
        "fine_motor_on_short": fine_short_evaluation,
        "motor_off_comparison": off_comparison,
        "paired_motor_on_comparison": paired,
        "performance_speedup": speedups,
        "visual_review_required": True,
        "next_action": (
            "perform_visual_review"
            if adoption_candidate
            else "retain_dt1e5_and_diagnose_dt_sensitivity"
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
    if result.get("kind") == "stage_a_2015_dt_sensitivity_decision":
        lines.extend(
            [
                "",
                "## dt sensitivity decision",
                "",
                f"- reference_stable: `{result['reference_stable']}`",
                f"- adoption_candidate: `{result['adoption_candidate']}`",
                f"- visual_review_required: `{result['visual_review_required']}`",
            ]
        )
        for profile, comparison in result["paired_motor_on_comparison"].items():
            lines.append(
                f"- `{profile}`: paired_pass=`{comparison['paired_pass']}`, "
                f"orientation=`{comparison['body_orientation_difference_deg']:.6g} deg`, "
                f"RMS=`{comparison['centered_bead_rms_difference_over_b']:.6g} b`, "
                f"max=`{comparison['centered_bead_max_difference_over_b']:.6g} b`"
            )
            lines.extend(f"  - {reason}" for reason in comparison["reasons"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_dt_comparison_csv(path: Path, result: dict[str, Any]) -> None:
    rows: list[dict[str, Any]] = []
    for profile, comparison in result["paired_motor_on_comparison"].items():
        speedup = result["performance_speedup"][profile]
        base = {
            "profile": profile,
            "paired_pass": comparison["paired_pass"],
            "comparison_t_s": comparison["comparison_t_s"],
            "body_orientation_difference_deg": comparison[
                "body_orientation_difference_deg"
            ],
            "centered_bead_rms_difference_over_b": comparison[
                "centered_bead_rms_difference_over_b"
            ],
            "centered_bead_max_difference_over_b": comparison[
                "centered_bead_max_difference_over_b"
            ],
            "motor_off_steps_per_s_ratio": speedup["motor_off_steps_per_s_ratio"],
            "motor_on_steps_per_s_ratio": speedup["motor_on_steps_per_s_ratio"],
        }
        for rotation in comparison["flag_rotation"]:
            row = dict(base)
            row.update(rotation)
            rows.append(row)
    if not rows:
        raise ValueError("dt comparison produced no rows")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_analysis(
    output_dir: Path,
    *,
    motor_off_run: Path,
    motor_on_run: Path | None = None,
    threshold_contract: Path | None = None,
    coarse_motor_off_run: Path | None = None,
    coarse_motor_on_run: Path | None = None,
    fine_motor_on_short_run: Path | None = None,
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
    dt_runs = (
        coarse_motor_off_run,
        coarse_motor_on_run,
        fine_motor_on_short_run,
    )
    if any(path is not None for path in dt_runs):
        if not all(path is not None for path in dt_runs):
            raise ValueError("All dt sensitivity run paths are required together")
        if threshold_contract is None:
            raise ValueError("threshold_contract is required for dt sensitivity")
        locked = load_threshold_contract(threshold_contract)
        dt_result = evaluate_dt_sensitivity(
            baseline_motor_off_run=motor_off_run,
            coarse_motor_off_run=coarse_motor_off_run,  # type: ignore[arg-type]
            coarse_motor_on_run=coarse_motor_on_run,  # type: ignore[arg-type]
            fine_motor_on_short_run=fine_motor_on_short_run,  # type: ignore[arg-type]
            contract=locked,
        )
        dt_json = output_dir / "dt_sensitivity_decision.json"
        dt_markdown = output_dir / "dt_sensitivity_decision.md"
        dt_csv = output_dir / "dt_sensitivity_comparison.csv"
        dt_json.write_text(
            json.dumps(dt_result, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )
        write_markdown(dt_markdown, dt_result)
        _write_dt_comparison_csv(dt_csv, dt_result)
        outputs.update(
            {
                "dt_sensitivity_decision_json": dt_json,
                "dt_sensitivity_decision_markdown": dt_markdown,
                "dt_sensitivity_comparison_csv": dt_csv,
            }
        )
    return outputs
