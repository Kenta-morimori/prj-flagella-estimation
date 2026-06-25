#!/usr/bin/env python3
"""Evaluate Phase 2.6 root-torque segment-couple conditions."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
import json
import logging
from pathlib import Path
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.sim.core import SimulationState, Simulator
from sim_swim.sim.helix_retention_gate import (
    summarize_single_flagellum_helix_retention,
)
from sim_swim.sim.params import SimulationConfig

LOCAL_SCALE_TARGETS = (
    "local_spring_scale",
    "local_bend_scale",
    "local_torsion_scale",
    "local_hook_scale",
)
DEFAULT_TORQUES = (0.5e-20, 1.0e-20, 2.0e-20, 4.0e-20, 6.0e-20, 8.0e-20, 1.0e-19)
DEFAULT_SCALE_VALUES = (1.0, 1.1, 1.2, 1.5, 2.0)

SUMMARY_FIELDS = (
    "condition_id",
    "stage",
    "target",
    "value",
    "torque_Nm",
    "force_distribution",
    "n_flagella",
    "stub_mode",
    "duration_s",
    "dt_star",
    "dt_internal_s",
    "local_hook_scale",
    "local_spring_scale",
    "local_bend_scale",
    "local_torsion_scale",
    "output_dir",
    "shape_pass",
    "helix_retention_pass",
    "first_fail_category",
    "first_fail_reason",
    "step_count",
    "net_abs_flag_helix_spin_revolutions",
    "signed_flag_helix_spin_rate_hz",
    "flag_helix_spin_direction_consistency",
    "helix_to_root_net_rotation_ratio",
    "median_abs_flag_helix_spin_rate_hz",
    "min_flag_helix_spin_fit_r2",
    "max_hook_len_rel_err",
    "max_local_attach_first_rel_err",
    "max_flag_bond_rel_err",
    "max_flag_bend_err_deg",
    "max_flag_torsion_err_deg",
    "final_shape_pass_nonbody",
    "final_first_fail_category_nonbody",
    "final_flag_helix_spin_rate_hz",
    "final_local_twist_root_orientation_deg",
    "final_local_twist_tip_orientation_deg",
    "final_local_twist_tip_activity_ratio",
    "body_displacement_um",
    "body_path_length_um",
    "body_mean_speed_um_s",
    "body_path_mean_speed_um_s",
    "body_axis_angle_change_deg",
)


@dataclass(frozen=True)
class ConditionSpec:
    stage: str
    target: str
    value: float
    torque_Nm: float
    local_scales: dict[str, float]


def _parse_float_list(text: str) -> list[float]:
    values: list[float] = []
    for raw_item in text.split(","):
        item = raw_item.strip()
        if item:
            values.append(float(item))
    if not values:
        raise argparse.ArgumentTypeError("at least one numeric value is required")
    return values


def _parse_positive_float(text: str) -> float:
    value = float(text)
    if value <= 0.0:
        raise argparse.ArgumentTypeError("value must be positive")
    return value


def _default_output_dir() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / ("phase2_6_torque_model_evaluation")
    )


def _load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _read_step_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _axis_from_quaternion(q: tuple[float, float, float, float]) -> np.ndarray:
    x, y, z, w = q
    return np.array(
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y + z * w),
            2.0 * (x * z - y * w),
        ],
        dtype=float,
    )


def _summarize_body_motion(states: list[SimulationState]) -> dict[str, float]:
    if len(states) < 2:
        return {
            "body_displacement_um": float("nan"),
            "body_path_length_um": float("nan"),
            "body_mean_speed_um_s": float("nan"),
            "body_path_mean_speed_um_s": float("nan"),
            "body_axis_angle_change_deg": float("nan"),
        }

    positions = np.asarray([state.position_um for state in states], dtype=float)
    displacement = float(np.linalg.norm(positions[-1] - positions[0]))
    path_length = float(np.sum(np.linalg.norm(np.diff(positions, axis=0), axis=1)))
    duration_s = max(float(states[-1].t - states[0].t), 1e-30)

    axis0 = _axis_from_quaternion(states[0].quaternion)
    axis1 = _axis_from_quaternion(states[-1].quaternion)
    axis0 = axis0 / max(float(np.linalg.norm(axis0)), 1e-18)
    axis1 = axis1 / max(float(np.linalg.norm(axis1)), 1e-18)
    axis_angle = float(
        np.rad2deg(np.arccos(np.clip(float(np.dot(axis0, axis1)), -1.0, 1.0)))
    )

    return {
        "body_displacement_um": displacement,
        "body_path_length_um": path_length,
        "body_mean_speed_um_s": displacement / duration_s,
        "body_path_mean_speed_um_s": path_length / duration_s,
        "body_axis_angle_change_deg": axis_angle,
    }


def _base_overrides(
    *,
    duration_s: float,
    dt_star: float,
    torque_Nm: float,
    local_scales: dict[str, float],
) -> dict[str, Any]:
    motor = {
        "torque_Nm": float(torque_Nm),
        "force_distribution": "root_torque_segment_couples",
        "local_hook_scale": float(local_scales["local_hook_scale"]),
        "local_spring_scale": float(local_scales["local_spring_scale"]),
        "local_bend_scale": float(local_scales["local_bend_scale"]),
        "local_torsion_scale": float(local_scales["local_torsion_scale"]),
        "enable_switching": False,
    }
    return {
        "flagella": {"n_flagella": 1, "stub_mode": "full_flagella"},
        "motor": motor,
        "time": {"duration_s": float(duration_s), "dt_star": float(dt_star)},
        "brownian": {"enabled": False},
        "render": {
            "render_flagella": False,
            "save_frames_3d": False,
            "save_frames_2d": False,
        },
    }


def _build_config(
    base_cfg: SimulationConfig,
    *,
    duration_s: float,
    dt_star: float,
    condition: ConditionSpec,
) -> SimulationConfig:
    return base_cfg.with_overrides(
        _base_overrides(
            duration_s=duration_s,
            dt_star=dt_star,
            torque_Nm=condition.torque_Nm,
            local_scales=condition.local_scales,
        )
    )


def _run_condition(
    *,
    base_cfg: SimulationConfig,
    condition: ConditionSpec,
    condition_id: int,
    run_dir: Path,
    duration_s: float,
    dt_star: float,
    min_net_spin_revolutions: float,
    logger: logging.Logger,
) -> dict[str, Any]:
    cfg = _build_config(
        base_cfg,
        duration_s=duration_s,
        dt_star=dt_star,
        condition=condition,
    )
    logger.info(
        "condition %03d: stage=%s target=%s value=%g torque=%.3e",
        condition_id,
        condition.stage,
        condition.target,
        condition.value,
        condition.torque_Nm,
    )
    for key, value in cfg.motor_local_scale_deviations().items():
        logger.warning(
            (
                "Paper model extension: motor.%s=%.6g differs from the "
                "paper-aligned default 1.0. Treat this condition as a local "
                "stiffness scaling diagnostic."
            ),
            key,
            value,
        )
    run_dir.mkdir(parents=True, exist_ok=True)
    sim = Simulator(cfg)
    states = sim.run(cfg.time.duration_s, step_summary_dir=run_dir)
    rows = _read_step_rows(run_dir / "step_summary.csv")
    last = rows[-1] if rows else {}
    helix = summarize_single_flagellum_helix_retention(
        rows,
        min_net_abs_spin_revolutions=min_net_spin_revolutions,
    )
    motion = _summarize_body_motion(states)
    shape_pass = bool(helix["helix_retention_pass"])

    return {
        "condition_id": condition_id,
        "stage": condition.stage,
        "target": condition.target,
        "value": condition.value,
        "torque_Nm": condition.torque_Nm,
        "force_distribution": cfg.motor.force_distribution,
        "n_flagella": cfg.flagella.n_flagella,
        "stub_mode": cfg.flagella.stub_mode,
        "duration_s": cfg.time.duration_s,
        "dt_star": cfg.dt_star,
        "dt_internal_s": cfg.dt_s,
        "local_hook_scale": cfg.motor.local_hook_scale,
        "local_spring_scale": cfg.motor.local_spring_scale,
        "local_bend_scale": cfg.motor.local_bend_scale,
        "local_torsion_scale": cfg.motor.local_torsion_scale,
        "output_dir": str(run_dir),
        "shape_pass": shape_pass,
        "helix_retention_pass": bool(helix["helix_retention_pass"]),
        "first_fail_category": helix["first_fail_category"],
        "first_fail_reason": helix["first_fail_reason"],
        "step_count": helix["step_count"],
        "net_abs_flag_helix_spin_revolutions": helix[
            "net_abs_flag_helix_spin_revolutions"
        ],
        "signed_flag_helix_spin_rate_hz": helix["signed_flag_helix_spin_rate_hz"],
        "flag_helix_spin_direction_consistency": helix[
            "flag_helix_spin_direction_consistency"
        ],
        "helix_to_root_net_rotation_ratio": helix["helix_to_root_net_rotation_ratio"],
        "median_abs_flag_helix_spin_rate_hz": helix[
            "median_abs_flag_helix_spin_rate_hz"
        ],
        "min_flag_helix_spin_fit_r2": helix["min_flag_helix_spin_fit_r2"],
        "max_hook_len_rel_err": helix["max_hook_len_rel_err"],
        "max_local_attach_first_rel_err": helix["max_local_attach_first_rel_err"],
        "max_flag_bond_rel_err": helix["max_flag_bond_rel_err"],
        "max_flag_bend_err_deg": helix["max_flag_bend_err_deg"],
        "max_flag_torsion_err_deg": helix["max_flag_torsion_err_deg"],
        "final_shape_pass_nonbody": last.get("shape_pass_nonbody", ""),
        "final_first_fail_category_nonbody": last.get(
            "first_fail_category_nonbody", ""
        ),
        "final_flag_helix_spin_rate_hz": last.get("flag_helix_spin_rate_hz", ""),
        "final_local_twist_root_orientation_deg": last.get(
            "local_twist_root_orientation_deg", ""
        ),
        "final_local_twist_tip_orientation_deg": last.get(
            "local_twist_tip_orientation_deg", ""
        ),
        "final_local_twist_tip_activity_ratio": last.get(
            "local_twist_tip_activity_ratio", ""
        ),
        **motion,
    }


def _baseline_conditions(torques: list[float]) -> list[ConditionSpec]:
    scales = {target: 1.0 for target in LOCAL_SCALE_TARGETS}
    return [
        ConditionSpec(
            stage="baseline_all_local_scales_1",
            target="all_local_scales",
            value=1.0,
            torque_Nm=float(torque),
            local_scales=dict(scales),
        )
        for torque in torques
    ]


def _scale_conditions(
    *,
    torques: list[float],
    values: list[float],
) -> list[ConditionSpec]:
    conditions: list[ConditionSpec] = []
    for torque in torques:
        for target in LOCAL_SCALE_TARGETS:
            for value in values:
                scales = {scale_target: 1.0 for scale_target in LOCAL_SCALE_TARGETS}
                scales[target] = float(value)
                conditions.append(
                    ConditionSpec(
                        stage="one_factor_scale",
                        target=target,
                        value=float(value),
                        torque_Nm=float(torque),
                        local_scales=scales,
                    )
                )
    return conditions


def _select_auto_scale_torques(
    baseline_rows: list[dict[str, Any]],
    *,
    torques: list[float],
) -> list[float]:
    shape_failing = [
        float(row["torque_Nm"])
        for row in baseline_rows
        if str(row["helix_retention_pass"]) != "True"
        and str(row.get("first_fail_category", "")) != "motor_no_rotation"
    ]
    selected = {2.0e-20, max(torques)}
    if shape_failing:
        selected.add(min(shape_failing))
    return sorted(torque for torque in selected if torque in set(torques))


def _write_manifest(
    path: Path,
    *,
    args: argparse.Namespace,
    summary_path: Path,
    condition_count: int,
    status: str,
) -> None:
    git_commit = "unknown"
    try:
        git_commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        pass

    manifest = {
        "task": "P2-6-009 / issue #54 torque transmission model evaluation",
        "config": str(args.config),
        "torques_Nm": args.torques,
        "scale_values": args.scale_values,
        "scale_torques": args.scale_torques,
        "duration_s": args.duration,
        "dt_star": args.dt_star,
        "min_net_spin_revolutions": args.min_net_spin_revolutions,
        "condition_count": condition_count,
        "summary_csv": str(summary_path),
        "git_commit": git_commit,
        "status": status,
    }
    path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )


def _setup_logger(output_dir: Path) -> logging.Logger:
    logger = logging.getLogger("phase2_6_torque_model_evaluation")
    for handler in list(logger.handlers):
        handler.close()
    logger.handlers.clear()
    logger.propagate = False
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    file_handler = logging.FileHandler(output_dir / "run.log", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run Phase 2.6 torque transmission model evaluation."
    )
    parser.add_argument("--config", type=Path, default=Path("conf/sim_swim.yaml"))
    parser.add_argument(
        "--torques",
        type=_parse_float_list,
        default=list(DEFAULT_TORQUES),
        help="Comma-separated torque_Nm values. Default includes 1.0e-19.",
    )
    parser.add_argument(
        "--scale-values",
        type=_parse_float_list,
        default=list(DEFAULT_SCALE_VALUES),
        help="Comma-separated local scale values for one-factor sweep.",
    )
    parser.add_argument(
        "--scale-torques",
        default="auto",
        help=(
            "'auto', 'all', 'none', or comma-separated torque_Nm values for "
            "one-factor scale sweeps."
        ),
    )
    parser.add_argument("--duration", type=_parse_positive_float, default=0.5)
    parser.add_argument("--dt-star", type=_parse_positive_float, default=1.0e-4)
    parser.add_argument(
        "--min-net-spin-revolutions",
        type=_parse_positive_float,
        default=1.0,
        help="Minimum net helix spin revolutions for PASS.",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--max-conditions",
        type=int,
        default=None,
        help="Optional limit for smoke tests and quick diagnostics.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    output_dir = args.output_dir or _default_output_dir()
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = _setup_logger(output_dir)

    torques = [float(torque) for torque in args.torques]
    if any(torque <= 0.0 for torque in torques):
        raise SystemExit("all torques must be positive")
    if max(torques) > 1.0e-19:
        raise SystemExit("P2-6-009 torque upper bound is 1.0e-19 N m")

    base_cfg = SimulationConfig.from_dict(_load_config(args.config))
    summary_path = output_dir / "phase2_6_torque_model_evaluation_summary.csv"
    _write_manifest(
        output_dir / "manifest.json",
        args=args,
        summary_path=summary_path,
        condition_count=0,
        status="started",
    )

    rows: list[dict[str, Any]] = []
    condition_id = 0
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()

        baseline_conditions = _baseline_conditions(torques)
        if args.max_conditions is not None:
            baseline_conditions = baseline_conditions[: max(0, args.max_conditions)]

        for condition in baseline_conditions:
            condition_id += 1
            row = _run_condition(
                base_cfg=base_cfg,
                condition=condition,
                condition_id=condition_id,
                run_dir=output_dir / f"condition_{condition_id:03d}",
                duration_s=float(args.duration),
                dt_star=float(args.dt_star),
                min_net_spin_revolutions=float(args.min_net_spin_revolutions),
                logger=logger,
            )
            rows.append(row)
            writer.writerow(row)
            handle.flush()

        if args.scale_torques == "none":
            scale_torques: list[float] = []
        elif args.scale_torques == "all":
            scale_torques = torques
        elif args.scale_torques == "auto":
            scale_torques = _select_auto_scale_torques(rows, torques=torques)
        else:
            scale_torques = _parse_float_list(str(args.scale_torques))

        scale_conditions = _scale_conditions(
            torques=scale_torques,
            values=[float(value) for value in args.scale_values],
        )
        if args.max_conditions is not None:
            remaining = max(0, int(args.max_conditions) - len(rows))
            scale_conditions = scale_conditions[:remaining]

        for condition in scale_conditions:
            condition_id += 1
            row = _run_condition(
                base_cfg=base_cfg,
                condition=condition,
                condition_id=condition_id,
                run_dir=output_dir / f"condition_{condition_id:03d}",
                duration_s=float(args.duration),
                dt_star=float(args.dt_star),
                min_net_spin_revolutions=float(args.min_net_spin_revolutions),
                logger=logger,
            )
            rows.append(row)
            writer.writerow(row)
            handle.flush()

    _write_manifest(
        output_dir / "manifest.json",
        args=args,
        summary_path=summary_path,
        condition_count=len(rows),
        status="completed",
    )
    logger.info("summary saved to %s", summary_path)


if __name__ == "__main__":
    main()
