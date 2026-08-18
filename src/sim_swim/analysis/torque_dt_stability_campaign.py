"""Issue #61 runner and compact comparison for the 2010 torque-dt screen.

This is intentionally a separate opt-in campaign.  It must never reinterpret
the legacy 2010 ``tau_s=1`` baseline or the #183 tracking-reference outputs.
"""

from __future__ import annotations

import csv
from dataclasses import asdict
from datetime import datetime
import json
import math
from pathlib import Path
import time
from typing import Any
from zoneinfo import ZoneInfo

import numpy as np

from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.analysis.torque_dt_stability import build_plan
from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.core.run_context import init_run
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _now_jst() -> str:
    """Return a human-readable JST timestamp for the live campaign manifest."""

    return datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(timespec="seconds")


def _format_elapsed(seconds: float) -> str:
    """Format elapsed wall-clock time for terminal and run.log progress lines."""

    total_seconds = max(0, int(round(seconds)))
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}h{minutes:02d}m{seconds:02d}s"
    if minutes:
        return f"{minutes}m{seconds:02d}s"
    return f"{seconds}s"


def _write_run_manifest(
    root: Path,
    *,
    config_path: Path,
    raw: dict[str, Any],
    plan: dict[str, Any],
    records: list[dict[str, Any]],
    execution_status: str,
    started_at_jst: str,
    finished_at_jst: str | None = None,
    wall_seconds: float | None = None,
) -> None:
    """Persist campaign progress so it remains inspectable while the run is active."""

    execution: dict[str, Any] = {
        "status": execution_status,
        "started_at_jst": started_at_jst,
    }
    if finished_at_jst is not None:
        execution["finished_at_jst"] = finished_at_jst
    if wall_seconds is not None:
        execution["wall_seconds"] = wall_seconds
    _write_json(
        root / "run_manifest.json",
        {
            "kind": "phase2_2010_torque_linked_dt_stability_campaign",
            "contract_version": 1,
            "campaign_config": str(config_path),
            "base_config": raw["base_config"],
            "campaign_plan": str(root / "campaign_plan.json"),
            "execution": execution,
            "conditions": records,
        },
    )


def _assert_condition(cfg: SimulationConfig, condition: dict[str, Any]) -> None:
    profile = cfg.model_profile_manifest()
    if (
        profile.get("year") != 2010
        or profile.get("variant") != "project"
        or profile.get("resolution") != "legacy_project"
    ):
        raise ValueError("This campaign only supports the 2010 project profile")
    if cfg.time_scale_policy != "reference_torque":
        raise ValueError("This campaign requires time.scale_policy=reference_torque")
    if cfg.brownian.enabled:
        raise ValueError("This campaign requires brownian.enabled=false")
    if not cfg.motor.body_reaction_full_vector:
        raise ValueError("This campaign requires motor.body_reaction_full_vector=true")
    torque = float(condition["torque_Nm_per_flagellum"])
    if not (
        math.isclose(abs(cfg.motor_torque_Nm), torque, rel_tol=1e-12)
        and math.isclose(cfg.reference_torque_Nm, torque, rel_tol=1e-12)
        and math.isclose(cfg.torque_for_forces_Nm, torque, rel_tol=1e-12)
    ):
        raise ValueError(
            "motor/reference/force torques must all equal condition torque"
        )
    if not math.isclose(
        cfg.tau_s, cfg.viscosity_Pa_s * cfg.b_m**3 / torque, rel_tol=1e-12
    ):
        raise ValueError("tau_s does not follow eta*b^3/T")
    if not math.isclose(cfg.dt_s, cfg.dt_star * cfg.tau_s, rel_tol=1e-12):
        raise ValueError("dt_internal_s does not equal dt_star*tau_s")


def _validate_campaign_contract(raw: dict[str, Any]) -> None:
    contract = dict(raw.get("campaign_contract", {}) or {})
    name = contract.get("name")
    stage = contract.get("stage")
    if name not in {
        "2010_torque_linked_dt_stability",
        "2010_torque_linked_fixed_real_time_performance",
    }:
        raise ValueError("Missing supported 2010 torque-linked campaign contract")
    if stage == "fixed_real_time_performance":
        if name != "2010_torque_linked_fixed_real_time_performance":
            raise ValueError(
                "fixed_real_time_performance requires its dedicated contract"
            )
        if "duration_s" not in raw or "duration_tau" in raw:
            raise ValueError("fixed_real_time_performance requires duration_s only")
        if not math.isclose(float(raw["duration_s"]), 0.5, rel_tol=0.0):
            raise ValueError("fixed_real_time_performance requires duration_s=0.5")
        if {float(value) for value in raw.get("dt_stars", [])} != {1.0e-3}:
            raise ValueError("fixed_real_time_performance requires dt_star=1e-3")
        return
    if name != "2010_torque_linked_dt_stability" or stage != "initial_screen":
        raise ValueError("Only supported 2010 torque-linked stages are executable")
    active_dt = {float(value) for value in raw.get("dt_stars", [])}
    if active_dt != {1.0e-3, 1.0e-4}:
        raise ValueError("initial_screen must contain exactly dt_star=1e-3 and 1e-4")
    if {float(value) for value in raw.get("torques_Nm", [])} != {
        1.0e-21,
        2.5e-20,
        1.0e-19,
        1.2e-18,
    }:
        raise ValueError("initial_screen must contain its fixed four-torque grid")
    if not math.isclose(float(raw.get("duration_tau", 0.0)), 1.0, rel_tol=0.0):
        raise ValueError("initial_screen requires duration_tau=1")
    if int(raw.get("comparison_sample_count", 0)) != 201:
        raise ValueError("initial_screen requires comparison_sample_count=201")
    if float(contract.get("deferred_reference_dt_star", 0.0)) != 1.0e-5:
        raise ValueError(
            "initial_screen must defer dt_star=1e-5 as its formal reference"
        )


def _archive_arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as raw:
        return {
            key: np.asarray(raw[key])
            for key in (
                "t",
                "position_um",
                "quaternion",
                "velocity_um_s",
                "bead_positions_um",
            )
        }


def _validate_archive(
    path: Path, *, duration_s: float, count: int
) -> dict[str, np.ndarray]:
    values = _archive_arrays(path)
    expected = np.linspace(0.0, duration_s, count)
    if values["t"].shape != expected.shape or not np.allclose(
        values["t"], expected, rtol=0.0, atol=max(2e-14, duration_s * 2e-10)
    ):
        raise ValueError(
            "comparison archive must contain exactly the configured time grid"
        )
    return values


def _comparison_archive_states(states: list[Any], *, duration_s: float) -> list[Any]:
    """Exclude only the ceil-induced state after a requested final archive time."""

    tolerance = max(1e-15, duration_s * 1e-12)
    selected = [state for state in states if float(state.t) <= duration_s + tolerance]
    if not selected:
        raise ValueError("comparison archive has no states at or before duration")
    return selected


def _quat_angle_deg(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    dot = np.abs(np.sum(first * second, axis=1))
    return np.degrees(2.0 * np.arccos(np.clip(dot, -1.0, 1.0)))


def _body_frame_beads(
    beads: np.ndarray, positions: np.ndarray, quaternions: np.ndarray
) -> np.ndarray:
    """Return bead positions in the archived body frame (xyzw quaternion)."""
    x, y, z, w = (quaternions[:, index] for index in range(4))
    matrices = np.stack(
        (
            1 - 2 * (y * y + z * z),
            2 * (x * y - z * w),
            2 * (x * z + y * w),
            2 * (x * y + z * w),
            1 - 2 * (x * x + z * z),
            2 * (y * z - x * w),
            2 * (x * z - y * w),
            2 * (y * z + x * w),
            1 - 2 * (x * x + y * y),
        ),
        axis=1,
    ).reshape((-1, 3, 3))
    centered = beads - positions[:, None, :]
    return np.einsum("tji,tbj->tbi", matrices, centered)


def _rotation_degrees(
    beads: np.ndarray,
    positions: np.ndarray,
    quaternions: np.ndarray,
    cfg: SimulationConfig,
) -> list[float]:
    """Estimate each flagellum's body-frame helix phase displacement.

    The definition is deliberately comparison-only: its fixed initial
    body-frame root-to-tip axis removes body translation/attitude before the
    transverse phase is unwrapped.
    """
    rig = Simulator(cfg).model
    beads = _body_frame_beads(beads, positions, quaternions)
    result: list[float] = []
    for indices in rig.flagella_indices:
        root, probe, tip = (
            int(indices[0]),
            int(indices[min(3, len(indices) - 1)]),
            int(indices[-1]),
        )
        initial = beads[0]
        axis = initial[tip] - initial[root]
        axis /= max(float(np.linalg.norm(axis)), 1e-30)
        vector0 = initial[probe] - initial[root]
        vector0 -= axis * float(np.dot(vector0, axis))
        vector0 /= max(float(np.linalg.norm(vector0)), 1e-30)
        basis2 = np.cross(axis, vector0)
        phases: list[float] = []
        for frame in beads:
            vector = frame[probe] - frame[root]
            vector -= axis * float(np.dot(vector, axis))
            phases.append(
                math.atan2(
                    float(np.dot(vector, basis2)), float(np.dot(vector, vector0))
                )
            )
        result.append(
            float(
                np.degrees(
                    np.unwrap(np.asarray(phases))[-1] - np.unwrap(np.asarray(phases))[0]
                )
            )
        )
    return result


def _condition_record(
    root: Path, condition: dict[str, Any], cfg: SimulationConfig, qc: dict[str, Any]
) -> dict[str, Any]:
    output = root / condition["condition_id"]
    record: dict[str, Any] = {
        "condition_id": condition["condition_id"],
        "torque_Nm_per_flagellum": condition["torque_Nm_per_flagellum"],
        "dt_star": condition["dt_star"],
        "comparison_role": condition["comparison_role"],
        "output_dir": str(output),
        "config_overrides": condition["config_overrides"],
        "time": cfg.time_manifest(),
        "comparison_archive_interval_s": condition["comparison_archive_interval_s"],
        "comparison_sample_count": condition["comparison_sample_count"],
        "qc": qc,
    }
    summary_path = output / "run_summary.json"
    performance_path = output / "performance.json"
    if summary_path.is_file():
        record["run_summary_json"] = str(summary_path)
    if performance_path.is_file():
        record["performance_json"] = str(performance_path)
    return record


def _run_condition(
    root: Path,
    base: dict[str, Any],
    condition: dict[str, Any],
    qc: dict[str, Any],
) -> dict[str, Any]:
    cfg = SimulationConfig.from_dict(base).with_overrides(condition["config_overrides"])
    _assert_condition(cfg, condition)
    directory = root / condition["condition_id"]
    directory.mkdir(parents=True, exist_ok=False)
    _write_json(
        directory / "effective_config.json",
        {
            "config_overrides": condition["config_overrides"],
            "effective_config": asdict(cfg),
            "time": cfg.time_manifest(),
        },
    )
    try:
        simulator = Simulator(cfg)
        states = simulator.run(
            cfg.time.duration_s,
            step_summary_dir=directory,
            stop_on_shape_fail=False,
            record_body_diagnostics=True,
        )
        save_state_archive(
            directory / "state_archive.npz",
            _comparison_archive_states(states, duration_s=cfg.time.duration_s),
        )
        _validate_archive(
            directory / "state_archive.npz",
            duration_s=cfg.time.duration_s,
            count=int(condition["comparison_sample_count"]),
        )
        return _condition_record(root, condition, cfg, qc)
    except Exception as exc:
        _write_json(
            directory / "failure.json",
            {"error_type": type(exc).__name__, "error_message": str(exc)},
        )
        record = _condition_record(root, condition, cfg, qc)
        record["error_type"] = type(exc).__name__
        record["error_message"] = str(exc)
        return record


def _max_metric(summary: dict[str, Any], key: str) -> float | None:
    item = dict(summary.get("all_step_metrics", {}).get(key, {}) or {})
    value = item.get("max")
    return float(value) if value is not None and math.isfinite(float(value)) else None


def _safety(record: dict[str, Any], qc: dict[str, Any]) -> dict[str, Any]:
    if record.get("error_type"):
        return {"status": "fail", "reason": "simulation exception"}
    output = Path(record["output_dir"])
    summary = _read_json(output / "run_summary.json")
    gates = dict(summary.get("gates", {}) or {})
    complete = dict(summary.get("execution", {}) or {}).get("status") == "completed"
    required = ("finite", "shape_nonbody", "shape_body")
    gate_ok = all(
        not bool(dict(gates.get(name, {}) or {}).get("any_fail", True))
        for name in required
    )
    force = _max_metric(summary, "motor_force_balance_residual_ratio")
    torque = _max_metric(summary, "motor_torque_balance_residual_ratio")
    action_ok = (
        force is not None
        and torque is not None
        and force <= float(qc["max_motor_force_balance_residual_ratio"])
        and torque <= float(qc["max_motor_torque_balance_residual_ratio"])
    )
    return {
        "status": "pass" if complete and gate_ok and action_ok else "fail",
        "completed": complete,
        "gates_pass": gate_ok,
        "action_reaction_pass": action_ok,
        "max_motor_force_balance_residual_ratio": force,
        "max_motor_torque_balance_residual_ratio": torque,
    }


def _performance_rows(
    records: list[dict[str, Any]], safety_rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Combine simulator loop performance with mandatory safety results."""

    safety_by_id = {str(row["condition_id"]): row for row in safety_rows}
    rows: list[dict[str, Any]] = []
    for record in records:
        performance_path = Path(str(record["output_dir"])) / "performance.json"
        performance = _read_json(performance_path) if performance_path.is_file() else {}
        time_info = dict(record.get("time", {}) or {})
        rows.append(
            {
                "condition_id": record["condition_id"],
                "torque_Nm_per_flagellum": record["torque_Nm_per_flagellum"],
                "dt_star": record["dt_star"],
                "duration_s": time_info.get("duration_s"),
                "duration_tau": time_info.get("duration_tau"),
                "tau_s": time_info.get("tau_s"),
                "dt_internal_s": time_info.get("dt_internal_s"),
                "total_steps": performance.get("total_steps"),
                "completed_steps": performance.get("completed_steps"),
                "simulation_wall_time_s": performance.get("wall_time_s"),
                "steps_per_s": performance.get("steps_per_s"),
                "safety_status": safety_by_id[str(record["condition_id"])]["status"],
            }
        )
    return rows


def _completed_result_record(record: dict[str, Any]) -> bool:
    """Return whether a record has the artifacts needed for offline analysis."""

    output = Path(str(record["output_dir"]))
    return (
        dict(record.get("execution", {}) or {}).get("status") == "completed"
        and (output / "run_summary.json").is_file()
        and (output / "state_archive.npz").is_file()
        and _performance_artifact_error(output) is None
    )


def _performance_artifact_error(output: Path) -> str | None:
    """Return why a performance artifact cannot support a completed result."""

    path = output / "performance.json"
    if not path.is_file():
        return "missing performance.json"
    try:
        performance = _read_json(path)
    except (json.JSONDecodeError, OSError):
        return "invalid performance.json"
    required = ("wall_time_s", "steps_per_s", "total_steps", "completed_steps")
    for key in required:
        value = performance.get(key)
        try:
            if not math.isfinite(float(value)):
                return f"invalid performance.json ({key})"
        except (TypeError, ValueError):
            return f"invalid performance.json ({key})"
    return None


def _incomplete_reason(record: dict[str, Any]) -> str:
    execution = dict(record.get("execution", {}) or {})
    status = str(execution.get("status") or "unknown")
    output = Path(str(record["output_dir"]))
    missing = [
        name
        for name in ("run_summary.json", "state_archive.npz")
        if not (output / name).is_file()
    ]
    if missing:
        return f"execution={status}; missing {', '.join(missing)}"
    performance_error = _performance_artifact_error(output)
    if performance_error is not None:
        return f"execution={status}; {performance_error}"
    return f"execution={status}"


def _replay_summary_row(
    record: dict[str, Any], safety: dict[str, Any]
) -> dict[str, Any]:
    """Make the generic replay's small CSV row from the compact QC summary."""

    summary = _read_json(Path(str(record["output_dir"])) / "run_summary.json")
    gates = dict(summary.get("gates", {}) or {})
    first_fail: tuple[float, str] | None = None
    for gate in gates.values():
        item = dict(gate or {})
        observed = item.get("first_observed_fail_t_s")
        if observed is None:
            continue
        candidate = (float(observed), str(item.get("first_failure_category") or ""))
        if first_fail is None or candidate[0] < first_fail[0]:
            first_fail = candidate
    time_info = dict(record.get("time", {}) or {})
    return {
        "condition_id": record["condition_id"],
        "condition_label": (
            f"T={float(record['torque_Nm_per_flagellum']):.2e} N m\n"
            f"dt*={float(record['dt_star']):.0e}"
        ),
        "torque_Nm_per_flagellum": record["torque_Nm_per_flagellum"],
        "dt_star": record["dt_star"],
        "duration_s": time_info.get("duration_s"),
        "duration_tau": time_info.get("duration_tau"),
        "status": safety["status"],
        "final_shape_pass_nonbody": safety["gates_pass"],
        # phase2_replay historically calls this nonbody, but a body-gate failure
        # must be visible in a qualitative replay as well.
        "first_fail_t_s": first_fail[0] if first_fail is not None else "",
        "first_fail_category_nonbody": first_fail[1] if first_fail else "",
    }


def _observables(record: dict[str, Any], cfg: SimulationConfig) -> dict[str, Any]:
    archive = _validate_archive(
        Path(record["output_dir"]) / "state_archive.npz",
        duration_s=cfg.time.duration_s,
        count=int(record["comparison_sample_count"]),
    )
    displacement = archive["position_um"] - archive["position_um"][0]
    speed = np.linalg.norm(archive["velocity_um_s"], axis=1)
    return {
        "archive": archive,
        "rotation_deg_per_flagellum": _rotation_degrees(
            archive["bead_positions_um"],
            archive["position_um"],
            archive["quaternion"],
            cfg,
        ),
        "com_displacement_over_b": float(
            np.linalg.norm(displacement[-1]) / cfg.scale.b_um
        ),
        "mean_speed_um_s": float(np.mean(speed)),
        "max_speed_um_s": float(np.max(speed)),
        "body_orientation_change_deg": float(
            _quat_angle_deg(archive["quaternion"][:1], archive["quaternion"][-1:])[0]
        ),
    }


def _effective_config(record: dict[str, Any], base: dict[str, Any]) -> SimulationConfig:
    """Load the run-time config snapshot, with legacy-run fallback only."""

    snapshot_path = Path(str(record["output_dir"])) / "effective_config.json"
    if snapshot_path.is_file():
        snapshot = _read_json(snapshot_path).get("effective_config")
        if isinstance(snapshot, dict):
            return SimulationConfig.from_dict(snapshot)
    return SimulationConfig.from_dict(base).with_overrides(record["config_overrides"])


def _compare_same_torque(
    reference: dict[str, Any],
    candidate: dict[str, Any],
    ref_cfg: SimulationConfig,
    cand_cfg: SimulationConfig,
    qc: dict[str, Any],
) -> dict[str, Any]:
    if (
        _safety(reference, qc)["status"] != "pass"
        or _safety(candidate, qc)["status"] != "pass"
    ):
        return {
            "status": "not_comparable",
            "reason": "reference or candidate safety failure",
        }
    ref, cand = _observables(reference, ref_cfg), _observables(candidate, cand_cfg)
    ref_archive, cand_archive = ref["archive"], cand["archive"]
    angle = _quat_angle_deg(ref_archive["quaternion"], cand_archive["quaternion"])
    ref_beads = (
        ref_archive["bead_positions_um"] - ref_archive["position_um"][:, None, :]
    )
    cand_beads = (
        cand_archive["bead_positions_um"] - cand_archive["position_um"][:, None, :]
    )
    delta = np.linalg.norm(ref_beads - cand_beads, axis=2) / ref_cfg.scale.b_um
    rotation_errors = [
        abs(a - b) / max(abs(a), 1e-12)
        for a, b in zip(
            ref["rotation_deg_per_flagellum"], cand["rotation_deg_per_flagellum"]
        )
    ]
    direction_match = all(
        (a == 0.0 and b == 0.0) or a * b > 0.0
        for a, b in zip(
            ref["rotation_deg_per_flagellum"], cand["rotation_deg_per_flagellum"]
        )
    )
    values = {
        "max_body_orientation_difference_deg": float(np.max(angle)),
        "bead_rms_over_b": float(np.sqrt(np.mean(delta**2))),
        "max_bead_displacement_over_b": float(np.max(delta)),
        "max_rotation_relative_error": float(max(rotation_errors, default=0.0)),
        "rotation_direction_match": direction_match,
        "reference_rotation_deg_per_flagellum": ref["rotation_deg_per_flagellum"],
        "candidate_rotation_deg_per_flagellum": cand["rotation_deg_per_flagellum"],
        "reference_com_displacement_over_b": ref["com_displacement_over_b"],
        "candidate_com_displacement_over_b": cand["com_displacement_over_b"],
        "reference_mean_speed_um_s": ref["mean_speed_um_s"],
        "candidate_mean_speed_um_s": cand["mean_speed_um_s"],
        "swimming_stability_policy": qc["swimming_stability_policy"],
    }
    passed = (
        direction_match
        and values["max_rotation_relative_error"]
        <= float(qc["max_rotation_relative_error"])
        and values["max_body_orientation_difference_deg"]
        <= float(qc["max_body_orientation_difference_deg"])
        and values["bead_rms_over_b"] <= float(qc["max_bead_rms_over_b"])
        and values["max_bead_displacement_over_b"]
        <= float(qc["max_bead_displacement_over_b"])
    )
    return {"status": "pass" if passed else "fail", **values}


def summarize_campaign(
    root: Path, *, config_path: Path, allow_incomplete: bool = False
) -> dict[str, Any]:
    """Rebuild campaign summaries without simulation.

    ``allow_incomplete`` is only available for the fixed-real-time performance
    stage.  It is an explicit qualitative-review mode: missing conditions are
    recorded as exclusions and are never silently omitted from a conclusion.
    """

    raw = load_yaml(config_path)
    _validate_campaign_contract(raw)
    plan = build_plan(config_path)
    qc = dict(raw.get("qc", {}) or {})
    base = load_yaml(Path(str(raw["base_config"])))
    manifest = _read_json(root / "run_manifest.json")
    records = list(manifest["conditions"])
    stage = dict(raw.get("campaign_contract", {}) or {}).get("stage")
    if stage == "fixed_real_time_performance":
        completed_records = [row for row in records if _completed_result_record(row)]
        incomplete_records = [row for row in records if row not in completed_records]
    else:
        # Preserve the #61 initial-screen contract and its validation behavior.
        completed_records = records
        incomplete_records: list[dict[str, Any]] = []
    if incomplete_records and not allow_incomplete:
        details = "; ".join(
            f"{row['condition_id']} ({_incomplete_reason(row)})"
            for row in incomplete_records
        )
        raise ValueError(
            "Campaign is incomplete; pass allow_incomplete=True only for "
            f"fixed-real-time qualitative review. Incomplete: {details}"
        )
    if not completed_records:
        raise ValueError("No completed condition has replayable artifacts")
    configs: dict[str, SimulationConfig] = {}
    safety_rows: list[dict[str, Any]] = []
    replay_rows: list[dict[str, Any]] = []
    for row in completed_records:
        cfg = _effective_config(row, base)
        configs[row["condition_id"]] = cfg
        safety = _safety(row, qc)
        row["safety"] = safety
        replay_rows.append(_replay_summary_row(row, safety))
        safety_rows.append(
            {
                "condition_id": row["condition_id"],
                "torque_Nm_per_flagellum": row["torque_Nm_per_flagellum"],
                "dt_star": row["dt_star"],
                **safety,
            }
        )
    if stage == "fixed_real_time_performance":
        performance_rows = _performance_rows(completed_records, safety_rows)
        exclusions = [
            {
                "condition_id": row["condition_id"],
                "reason": _incomplete_reason(row),
                "included_in_qualitative_replay": False,
            }
            for row in incomplete_records
        ]
        _write_csv(root / "summary.csv", replay_rows)
        _write_csv(root / "performance_summary.csv", performance_rows)
        payload = {
            "kind": "phase2_2010_torque_linked_fixed_real_time_performance_qc",
            "contract_version": 1,
            "source_config": str(config_path),
            "campaign_plan_kind": plan["kind"],
            "qc_thresholds": qc,
            "conditions": safety_rows,
            "performance": performance_rows,
            "exclusions": exclusions,
            "interpretation": {
                "scope": "fixed_real_time_performance_only",
                "duration_s": float(raw["duration_s"]),
                "dt_star": 1.0e-3,
                "completion_status": (
                    "partial_qualitative_review_only" if exclusions else "complete"
                ),
                "prohibitions": [
                    "Do not interpret this as a same-normalized-time similarity test.",
                    "Do not use this campaign to adopt dt_star or torque.",
                    "Do not treat excluded conditions as pass, fail, or replay evidence.",
                ],
            },
        }
        _write_json(root / "qc_summary.json", payload)
        manifest["conditions"] = records
        manifest["outputs"] = {
            "summary_csv": str(root / "summary.csv"),
            "performance_summary_csv": str(root / "performance_summary.csv"),
            "qc_summary_json": str(root / "qc_summary.json"),
        }
        _write_json(root / "run_manifest.json", manifest)
        return payload
    comparisons: list[dict[str, Any]] = []
    torques = sorted(
        {float(row["torque_Nm_per_flagellum"]) for row in completed_records}
    )
    refs = {
        float(row["torque_Nm_per_flagellum"]): row
        for row in completed_records
        if row["comparison_role"] == "formal_reference"
    }
    screen_comparators = {
        float(row["torque_Nm_per_flagellum"]): row
        for row in completed_records
        if row["comparison_role"] == "screen_comparator"
    }
    for torque in torques:
        reference = refs.get(torque)
        if reference is None:
            reference = screen_comparators.get(torque)
            comparison_stage = "provisional_screen"
        else:
            comparison_stage = "formal_reference"
        if reference is None:
            continue
        for candidate in completed_records:
            if (
                float(candidate["torque_Nm_per_flagellum"]) == torque
                and candidate["condition_id"] != reference["condition_id"]
            ):
                comparison = {
                    "torque_Nm_per_flagellum": torque,
                    "reference_condition_id": reference["condition_id"],
                    "candidate_condition_id": candidate["condition_id"],
                    "candidate_dt_star": candidate["dt_star"],
                    "comparison_stage": comparison_stage,
                    **_compare_same_torque(
                        reference,
                        candidate,
                        configs[reference["condition_id"]],
                        configs[candidate["condition_id"]],
                        qc,
                    ),
                }
                if comparison_stage == "provisional_screen" and comparison[
                    "status"
                ] in {
                    "pass",
                    "fail",
                }:
                    comparison["status"] = "diagnostic_only"
                    comparison["formal_adoption_status"] = "awaiting_1e-5_reference"
                comparisons.append(comparison)
    # Similarity is deliberately descriptive: compare reference trajectories after COM removal.
    similarity: list[dict[str, Any]] = []
    reference_rows = [
        (refs.get(value) or screen_comparators[value])
        for value in torques
        if value in refs or value in screen_comparators
    ]
    if reference_rows:
        anchor = reference_rows[0]
        for row in reference_rows[1:]:
            if (
                _safety(anchor, qc)["status"] != "pass"
                or _safety(row, qc)["status"] != "pass"
            ):
                similarity.append(
                    {
                        "reference_condition_id": anchor["condition_id"],
                        "candidate_condition_id": row["condition_id"],
                        "status": "not_comparable",
                    }
                )
                continue
            a, b = (
                _observables(anchor, configs[anchor["condition_id"]]),
                _observables(row, configs[row["condition_id"]]),
            )
            beads_a = (
                a["archive"]["bead_positions_um"]
                - a["archive"]["position_um"][:, None, :]
            )
            beads_b = (
                b["archive"]["bead_positions_um"]
                - b["archive"]["position_um"][:, None, :]
            )
            similarity.append(
                {
                    "reference_condition_id": anchor["condition_id"],
                    "candidate_condition_id": row["condition_id"],
                    "status": "recorded",
                    "reference_torque_Nm_per_flagellum": anchor[
                        "torque_Nm_per_flagellum"
                    ],
                    "candidate_torque_Nm_per_flagellum": row["torque_Nm_per_flagellum"],
                    "bead_rms_over_b": float(
                        np.sqrt(np.mean((beads_a - beads_b) ** 2))
                        / configs[anchor["condition_id"]].scale.b_um
                    ),
                    "max_body_orientation_difference_deg": float(
                        np.max(
                            _quat_angle_deg(
                                a["archive"]["quaternion"], b["archive"]["quaternion"]
                            )
                        )
                    ),
                    "reference_rotation_deg_per_flagellum": a[
                        "rotation_deg_per_flagellum"
                    ],
                    "candidate_rotation_deg_per_flagellum": b[
                        "rotation_deg_per_flagellum"
                    ],
                }
            )
    _write_csv(root / "summary.csv", safety_rows)
    _write_csv(root / "dt_comparison.csv", comparisons)
    _write_csv(root / "torque_similarity.csv", similarity)
    payload = {
        "kind": "phase2_2010_torque_linked_dt_stability_qc",
        "contract_version": 1,
        "source_config": str(config_path),
        "campaign_plan_kind": plan["kind"],
        "qc_thresholds": qc,
        "conditions": safety_rows,
        "dt_comparisons": comparisons,
        "torque_similarity": similarity,
        "interpretation": {
            "dt_candidate_default": 1.0e-4,
            "dt_reference": 1.0e-5,
            "formal_dt_adoption_status": (
                "awaiting_1e-5_reference"
                if not refs
                else "reference_comparison_available"
            ),
            "dt_boundary_screen": 1.0e-3,
            "swimming_stability": "recorded; review required; not a standalone 1-tau adoption gate",
            "torque_similarity": "diagnostic only; not a dt adoption gate",
        },
    }
    _write_json(root / "qc_summary.json", payload)
    manifest["conditions"] = records
    manifest["outputs"] = {
        "summary_csv": str(root / "summary.csv"),
        "dt_comparison_csv": str(root / "dt_comparison.csv"),
        "torque_similarity_csv": str(root / "torque_similarity.csv"),
        "qc_summary_json": str(root / "qc_summary.json"),
    }
    _write_json(root / "run_manifest.json", manifest)
    return payload


def render_fixed_real_time_qualitative_replay(
    root: Path,
    *,
    config_path: Path,
    allow_incomplete: bool,
    fps_out_3d: float = 20.0,
    target_frame_count: int = 101,
) -> Path:
    """Render completed fixed-real-time conditions without restarting simulation.

    This deliberately delegates drawing to the common Phase 2 replay renderer.
    The summary generated immediately before rendering contains only completed
    artifacts, while the root QC and replay manifests retain every excluded
    condition and its reason.
    """

    if not allow_incomplete:
        raise ValueError(
            "Qualitative replay requires explicit allow_incomplete=True; "
            "it must not silently omit planned conditions"
        )
    raw = load_yaml(config_path)
    _validate_campaign_contract(raw)
    stage = dict(raw.get("campaign_contract", {}) or {}).get("stage")
    if stage != "fixed_real_time_performance":
        raise ValueError(
            "Qualitative partial replay is only for fixed_real_time_performance"
        )
    payload = summarize_campaign(
        root, config_path=config_path, allow_incomplete=allow_incomplete
    )
    output_dir = root / "analysis" / "fixed_real_time_qualitative_replay"
    exclusions = list(payload.get("exclusions", []) or [])
    figure_note = (
        "Partial qualitative review: "
        + "; ".join(f"{row['condition_id']} not run (excluded)" for row in exclusions)
        if exclusions
        else ""
    )
    from sim_swim.analysis.phase2_replay import main as render_phase2_replay

    render_phase2_replay(
        [
            "--input-dir",
            str(root),
            "--output-dir",
            str(output_dir),
            "--mode",
            "render-only",
            "--fps-out-3d",
            str(fps_out_3d),
            "--target-frame-count",
            str(target_frame_count),
            "--overwrite",
        ]
        + (["--figure-note", figure_note] if figure_note else [])
    )
    replay_manifest_path = output_dir / "manifest.json"
    replay_manifest = _read_json(replay_manifest_path)
    replay_manifest["qualitative_review"] = {
        "scope": "same-real-time completed conditions only",
        "source_duration_s": float(raw["duration_s"]),
        "excluded_conditions": exclusions,
        "not_a_same_normalized_time_comparison": True,
    }
    _write_json(replay_manifest_path, replay_manifest)
    root_manifest = _read_json(root / "run_manifest.json")
    root_manifest.setdefault("outputs", {})["qualitative_replay_manifest_json"] = str(
        replay_manifest_path
    )
    _write_json(root / "run_manifest.json", root_manifest)
    return output_dir


def run_torque_linked_campaign(
    config_path: Path, *, dry_run: bool = False, output_root: Path | None = None
) -> Path:
    raw = load_yaml(config_path)
    _validate_campaign_contract(raw)
    plan = build_plan(config_path)
    if dry_run:
        for condition in plan["conditions"]:
            print(condition["condition_id"])
        return Path()
    output = dict(raw.get("output", {}) or {})
    ctx = init_run(
        base_dir=output_root or Path(str(output.get("base_dir"))),
        timestamp_subdir=bool(output.get("timestamp_subdir", True)),
        run_name=(str(output["run_name"]) if output.get("run_name") else None),
        input_info={"campaign_config": str(config_path), "campaign_plan": plan},
        source_config_path=raw["base_config"],
        model_profile=SimulationConfig.from_dict(
            load_yaml(Path(str(raw["base_config"])))
        ).model_profile_manifest(),
    )
    root = ctx.out.root
    _write_json(root / "campaign_plan.json", plan)
    qc = dict(raw["qc"])
    base = load_yaml(Path(str(raw["base_config"])))
    started_at_jst = _now_jst()
    campaign_started = time.perf_counter()
    records: list[dict[str, Any]] = [
        {
            "condition_id": condition["condition_id"],
            "torque_Nm_per_flagellum": condition["torque_Nm_per_flagellum"],
            "dt_star": condition["dt_star"],
            "comparison_role": condition["comparison_role"],
            "output_dir": str(root / condition["condition_id"]),
            "execution": {"status": "planned"},
        }
        for condition in plan["conditions"]
    ]
    _write_run_manifest(
        root,
        config_path=config_path,
        raw=raw,
        plan=plan,
        records=records,
        execution_status="running",
        started_at_jst=started_at_jst,
    )
    ctx.logger.info("Torque-dt campaign start: conditions=%d", len(records))
    for index, condition in enumerate(plan["conditions"]):
        condition_started_at_jst = _now_jst()
        condition_started = time.perf_counter()
        records[index]["execution"] = {
            "status": "running",
            "started_at_jst": condition_started_at_jst,
        }
        _write_run_manifest(
            root,
            config_path=config_path,
            raw=raw,
            plan=plan,
            records=records,
            execution_status="running",
            started_at_jst=started_at_jst,
        )
        ctx.logger.info(
            "Condition start [%d/%d]: %s",
            index + 1,
            len(records),
            condition["condition_id"],
        )
        record = _run_condition(root, base, condition, qc)
        condition_wall_seconds = time.perf_counter() - condition_started
        condition_finished_at_jst = _now_jst()
        condition_status = "failed" if record.get("error_type") else "completed"
        record["execution"] = {
            "status": condition_status,
            "started_at_jst": condition_started_at_jst,
            "finished_at_jst": condition_finished_at_jst,
            "wall_seconds": condition_wall_seconds,
        }
        records[index] = record
        _write_run_manifest(
            root,
            config_path=config_path,
            raw=raw,
            plan=plan,
            records=records,
            execution_status="running",
            started_at_jst=started_at_jst,
        )
        ctx.logger.info(
            "Condition %s [%d/%d]: %s (elapsed=%s, completed=%d/%d)",
            condition_status,
            index + 1,
            len(records),
            condition["condition_id"],
            _format_elapsed(condition_wall_seconds),
            sum(
                1
                for item in records
                if item.get("execution", {}).get("status") in {"completed", "failed"}
            ),
            len(records),
        )
    campaign_wall_seconds = time.perf_counter() - campaign_started
    campaign_finished_at_jst = _now_jst()
    _write_run_manifest(
        root,
        config_path=config_path,
        raw=raw,
        plan=plan,
        records=records,
        execution_status="completed",
        started_at_jst=started_at_jst,
        finished_at_jst=campaign_finished_at_jst,
        wall_seconds=campaign_wall_seconds,
    )
    failed_count = sum(1 for record in records if record.get("error_type"))
    ctx.logger.info(
        "Torque-dt campaign completed: completed=%d, failed=%d, elapsed=%s",
        len(records) - failed_count,
        failed_count,
        _format_elapsed(campaign_wall_seconds),
    )
    summarize_campaign(root, config_path=config_path)
    root_manifest_path = root / "manifest.json"
    root_manifest = _read_json(root_manifest_path)
    stage = dict(raw.get("campaign_contract", {}) or {}).get("stage")
    campaign_outputs = {
        "campaign_plan_json": str(root / "campaign_plan.json"),
        "run_manifest_json": str(root / "run_manifest.json"),
        "summary_csv": str(root / "summary.csv"),
        "qc_summary_json": str(root / "qc_summary.json"),
    }
    if stage == "fixed_real_time_performance":
        campaign_outputs["performance_summary_csv"] = str(
            root / "performance_summary.csv"
        )
    else:
        campaign_outputs.update(
            {
                "dt_comparison_csv": str(root / "dt_comparison.csv"),
                "torque_similarity_csv": str(root / "torque_similarity.csv"),
            }
        )
    root_manifest.setdefault("input", {})["campaign"] = {
        "kind": "phase2_2010_torque_linked_dt_stability_campaign",
        "config": str(config_path),
        "plan": str(root / "campaign_plan.json"),
        "condition_count": len(records),
    }
    root_manifest.setdefault("outputs", {}).update(campaign_outputs)
    _write_json(root_manifest_path, root_manifest)
    return root
