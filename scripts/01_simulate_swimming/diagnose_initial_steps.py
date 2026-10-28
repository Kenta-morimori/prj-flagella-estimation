"""Record a bounded, reproducible first-step force and shape diagnosis."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np

from sim_swim.analysis.bundle_rotation import motor_root_lateral_arms_um
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.model.builder import ModelBuilder
from sim_swim.model.initial_bundle import initial_bundle_force_metrics
from sim_swim.sim.debug_summary import StepSummaryRecorder
from sim_swim.sim.params import SimulationConfig


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _minimum_torsion_plane_sine(model) -> float:
    quads = model.torsion_quads[model.torsion_flag_ids >= 0]
    points = model.positions_m[quads]
    b0 = points[:, 0] - points[:, 1]
    b1 = points[:, 2] - points[:, 1]
    b2 = points[:, 3] - points[:, 2]
    direction = b1 / np.linalg.norm(b1, axis=1)[:, None]
    first = b0 - np.sum(b0 * direction, axis=1)[:, None] * direction
    last = b2 - np.sum(b2 * direction, axis=1)[:, None] * direction
    return float(
        min(
            np.min(np.linalg.norm(first, axis=1) / np.linalg.norm(b0, axis=1)),
            np.min(np.linalg.norm(last, axis=1) / np.linalg.norm(b2, axis=1)),
        )
    )


def diagnose(
    campaign_path: Path,
    condition_id: str,
    *,
    steps: int,
    dt_star: float | None,
    motor_off: bool,
    torsion_off: bool,
    recorder_dir: Path,
) -> dict:
    campaign = normalize_campaign_config(load_yaml(campaign_path))
    condition = next(
        (
            row
            for row in build_campaign_conditions(campaign)
            if row["condition_id"] == condition_id
        ),
        None,
    )
    if condition is None:
        raise ValueError(f"Unknown campaign condition: {condition_id}")
    base_path = Path(campaign["base_config"])
    cfg = SimulationConfig.from_dict(load_yaml(base_path)).with_overrides(
        condition["config_overrides"]
    )
    if dt_star is not None:
        cfg = cfg.with_overrides({"time": {"integration": {"dt_star": dt_star}}})
    if torsion_off:
        cfg = cfg.with_overrides({"potentials": {"torsion": {"kt_over_T": 0.0}}})
    if motor_off:
        target_dt_s = cfg.dt_s
        cfg = cfg.with_overrides(
            {
                "motor": {
                    "enabled": False,
                    "allow_reference_torque_mismatch": True,
                    "torque_for_forces_override_Nm": cfg.torque_for_forces_Nm,
                }
            }
        )
        cfg = cfg.with_overrides(
            {"time": {"integration": {"dt_star": target_dt_s / cfg.tau_s}}}
        )
    model = ModelBuilder(cfg).build()
    engine = DynamicsEngine(model, cfg)
    recorder = StepSummaryRecorder(model, cfg, recorder_dir, emit_csv=False)
    initial = initial_bundle_force_metrics(cfg, model)
    rows = []
    for step in range(steps):
        diag = engine.step(cfg.dt_star)
        recorder.record(step=step, t_star=engine.t_star, diag=diag)
        summary = recorder.last_row
        rows.append(
            {
                "step": step + 1,
                "t_s": float(engine.t_star * cfg.tau_s),
                "shape_pass_nonbody": bool(summary["shape_pass_nonbody"]),
                "flag_bond_rel_err_max": float(summary["flag_bond_rel_err_max"]),
                "flag_bend_err_max_deg": float(summary["flag_bend_err_max_deg"]),
                "flag_torsion_err_max_deg": float(summary["flag_torsion_err_max_deg"]),
                "minimum_torsion_plane_sine": _minimum_torsion_plane_sine(model),
                "max_displacement_um": float(
                    np.max(
                        np.linalg.norm(
                            diag.positions_after_m - diag.positions_before_m, axis=1
                        )
                    )
                    * 1e6
                ),
                "motor_root_lateral_arms_um": motor_root_lateral_arms_um(
                    diag.positions_before_m, model.flagella_indices
                ).tolist(),
                "force_max_N": {
                    name: float(np.max(np.linalg.norm(getattr(diag, name), axis=1)))
                    for name in (
                        "spring_forces",
                        "bend_forces",
                        "torsion_forces",
                        "hook_forces",
                        "repulsion_forces",
                        "motor_forces",
                        "total_forces",
                    )
                },
            }
        )
    return {
        "condition_id": condition_id,
        "campaign": str(campaign_path),
        "campaign_sha256": _sha(campaign_path),
        "base_config_sha256": _sha(base_path),
        "state_sha256": cfg.flagella.initial_bundle_state_sha256,
        "dt_star": cfg.dt_star,
        "dt_s": cfg.dt_s,
        "motor_off": motor_off,
        "torsion_off_diagnostic_only": torsion_off,
        "torque_for_forces_Nm": cfg.torque_for_forces_Nm,
        "initial": initial,
        "steps": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--steps", type=int, default=2)
    parser.add_argument("--dt-star", type=float)
    parser.add_argument("--motor-off", action="store_true")
    parser.add_argument("--torsion-off", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.steps <= 20:
        parser.error("--steps must be in [1, 20]")
    output = diagnose(
        args.campaign,
        args.condition,
        steps=args.steps,
        dt_star=args.dt_star,
        motor_off=args.motor_off,
        torsion_off=args.torsion_off,
        recorder_dir=args.output_dir,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    data_path = args.output_dir / "initial_step_diagnostic.json"
    data_path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    (args.output_dir / "run.log").write_text(
        f"{args.condition}: {args.steps} internal steps; "
        f"motor_off={args.motor_off}; torsion_off={args.torsion_off}\n",
        encoding="utf-8",
    )
    (args.output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "kind": "initial_step_diagnostic",
                "created_at_jst": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
                "diagnostic_sha256": _sha(data_path),
                "status": "completed",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
