#!/usr/bin/env python3
"""Run Issue #82 hook-overstretch mitigation diagnostics."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
import math
from pathlib import Path
import sys
from typing import Any
from zoneinfo import ZoneInfo

import yaml

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.sim.core import Simulator
from sim_swim.sim.debug_summary import PROXIMAL_FLAG_BOND_REL_ERR_FIELDS
from sim_swim.sim.helix_retention_gate import summarize_single_flagellum_helix_retention
from sim_swim.sim.params import SimulationConfig

FIRST_FAIL_PROXIMAL_FLAG_BOND_REL_ERR_FIELDS = tuple(
    f"first_fail_{field}" for field in PROXIMAL_FLAG_BOND_REL_ERR_FIELDS
)
MAX_PROXIMAL_FLAG_BOND_REL_ERR_FIELDS = tuple(
    f"max_{field}" for field in PROXIMAL_FLAG_BOND_REL_ERR_FIELDS
)

SUMMARY_FIELDS = (
    "condition_id",
    "mode",
    "description",
    "output_dir",
    "duration_s",
    "dt_star",
    "torque_Nm",
    "n_flagella",
    "local_attach_first_spring_scale",
    "local_attach_first_body_axis_angle_scale",
    "local_first_second_spring_scale",
    "local_attach_frame_position_scale",
    "local_attach_frame_tangent_scale",
    "final_t_s",
    "final_shape_pass_nonbody",
    "final_first_fail_category_nonbody",
    "hook_len_rel_err_max",
    "local_attach_first_rel_err",
    "local_first_second_rel_err",
    "local_attach_first_vs_body_axis_angle_deg",
    "local_attach_first_vs_body_axis_err_deg",
    "local_attach_frame_position_rel_err",
    "local_attach_frame_position_angle_err_deg",
    "local_attach_frame_tangent_angle_err_deg",
    "hook_len_rel_err_max_flag_id",
    "hook_len_rel_err_max_attach_body_bead_index",
    "hook_len_rel_err_max_flag_first_bead_index",
    "hook_len_rel_err_max_len_over_b",
    "hook_len_rel_err_per_flag",
    "local_attach_first_rel_err_per_flag",
    "local_first_second_rel_err_per_flag",
    "local_attach_first_vs_body_axis_err_deg_per_flag",
    "local_attach_frame_position_rel_err_per_flag",
    "local_attach_frame_position_angle_err_deg_per_flag",
    "local_attach_frame_tangent_angle_err_deg_per_flag",
    "first_fail_t_s",
    "first_fail_category_nonbody",
    "first_fail_hook_len_rel_err_max",
    "first_fail_hook_len_rel_err_max_flag_id",
    "first_fail_local_attach_first_rel_err",
    "first_fail_local_first_second_rel_err",
    "first_fail_local_attach_frame_position_rel_err",
    "first_fail_local_attach_frame_position_angle_err_deg",
    "first_fail_local_attach_frame_tangent_angle_err_deg",
    "max_hook_len_rel_err_t_s",
    "max_hook_len_rel_err",
    "max_hook_len_rel_err_flag_id",
    "max_hook_len_rel_err_attach_body_bead_index",
    "max_hook_len_rel_err_flag_first_bead_index",
    "max_hook_len_rel_err_len_over_b",
    "max_hook_local_attach_frame_position_rel_err",
    "max_hook_local_attach_frame_position_angle_err_deg",
    "max_hook_local_attach_frame_tangent_angle_err_deg",
    "flag_bond_rel_err_max",
    "flag_bond_rel_err_max_flag_id",
    "flag_bond_rel_err_max_bead_i",
    "flag_bond_rel_err_max_bead_j",
    "flag_bond_rel_err_max_local_bead_i",
    "flag_bond_rel_err_max_local_bead_j",
    "flag_bond_rel_err_max_len_over_b",
    "flag_bond_rel_err_per_flag",
    *PROXIMAL_FLAG_BOND_REL_ERR_FIELDS,
    "first_fail_flag_bond_rel_err_max",
    "first_fail_flag_bond_rel_err_max_flag_id",
    "first_fail_flag_bond_rel_err_max_bead_i",
    "first_fail_flag_bond_rel_err_max_bead_j",
    "first_fail_flag_bond_rel_err_max_local_bead_i",
    "first_fail_flag_bond_rel_err_max_local_bead_j",
    "first_fail_flag_bond_rel_err_max_len_over_b",
    *FIRST_FAIL_PROXIMAL_FLAG_BOND_REL_ERR_FIELDS,
    "max_flag_bond_rel_err_t_s",
    "max_flag_bond_rel_err",
    "max_flag_bond_rel_err_flag_id",
    "max_flag_bond_rel_err_bead_i",
    "max_flag_bond_rel_err_bead_j",
    "max_flag_bond_rel_err_local_bead_i",
    "max_flag_bond_rel_err_local_bead_j",
    "max_flag_bond_rel_err_len_over_b",
    *MAX_PROXIMAL_FLAG_BOND_REL_ERR_FIELDS,
    "flag_bend_err_max_deg",
    "flag_torsion_err_max_deg",
    "net_abs_flag_helix_spin_revolutions",
    "flag_helix_spin_direction_consistency",
)


@dataclass(frozen=True)
class Condition:
    condition_id: str
    mode: str
    description: str
    scales: dict[str, float]


CONDITIONS = (
    Condition("baseline", "preset", "no extra local mitigation", {}),
    Condition(
        "attach_first_spring_2",
        "preset",
        "strengthen body attach to first bead distance",
        {"local_attach_first_spring_scale": 2.0},
    ),
    Condition(
        "body_axis_angle_2",
        "preset",
        "keep attach-first perpendicular to body long axis",
        {"local_attach_first_body_axis_angle_scale": 2.0},
    ),
    Condition(
        "first_second_spring_2",
        "preset",
        "strengthen first to second bead distance",
        {"local_first_second_spring_scale": 2.0},
    ),
    Condition(
        "attach_angle_2",
        "preset",
        "combine attach-first distance and body-axis 90 degree mitigation",
        {
            "local_attach_first_spring_scale": 2.0,
            "local_attach_first_body_axis_angle_scale": 2.0,
        },
    ),
    Condition(
        "attach_angle_first_second_2",
        "preset",
        "combine all Issue #82 local mitigations",
        {
            "local_attach_first_spring_scale": 2.0,
            "local_attach_first_body_axis_angle_scale": 2.0,
            "local_first_second_spring_scale": 2.0,
        },
    ),
    Condition(
        "attach_angle_first_second_4",
        "preset",
        "strong combined local mitigation",
        {
            "local_attach_first_spring_scale": 4.0,
            "local_attach_first_body_axis_angle_scale": 4.0,
            "local_first_second_spring_scale": 4.0,
        },
    ),
)


def _default_output_dir() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / "phase2_82_hook_overstretch"
    )


def _load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _read_step_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _parse_bool(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _parse_float(value: str | float | None) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")


def _parse_int(value: str | int | None) -> int | None:
    try:
        return int(float(value))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _first_fail_row(rows: list[dict[str, str]]) -> dict[str, str] | None:
    for row in rows[1:]:
        if not _parse_bool(row.get("shape_pass_nonbody")):
            return row
    return None


def _max_hook_row(rows: list[dict[str, str]]) -> dict[str, str] | None:
    finite_rows = [
        row
        for row in rows
        if math.isfinite(_parse_float(row.get("hook_len_rel_err_max")))
    ]
    if not finite_rows:
        return None
    return max(
        finite_rows, key=lambda row: _parse_float(row.get("hook_len_rel_err_max"))
    )


def _max_flag_bond_row(rows: list[dict[str, str]]) -> dict[str, str] | None:
    finite_rows = [
        row
        for row in rows
        if math.isfinite(_parse_float(row.get("flag_bond_rel_err_max")))
    ]
    if not finite_rows:
        return None
    return max(
        finite_rows, key=lambda row: _parse_float(row.get("flag_bond_rel_err_max"))
    )


def _flag_bond_local_bead_index(
    cfg: SimulationConfig,
    row: dict[str, str] | None,
    bead_field: str,
    flag_field: str = "flag_bond_rel_err_max_flag_id",
) -> str:
    if row is None:
        return ""
    flag_id = _parse_int(row.get(flag_field))
    bead_index = _parse_int(row.get(bead_field))
    if flag_id is None or bead_index is None:
        return ""

    n_body = cfg.compute_body_n_layers() * int(cfg.body.prism.n_prism)
    n_flag = int(cfg.flagella.n_beads_per_flagellum)
    local_index = bead_index - n_body - flag_id * n_flag
    if local_index < 0 or local_index >= n_flag:
        return ""
    return str(local_index)


def _copy_prefixed_fields(
    target: dict[str, str | float],
    source: dict[str, str] | None,
    fields: tuple[str, ...],
    prefix: str,
) -> None:
    for field in fields:
        target[f"{prefix}{field}"] = "" if source is None else source.get(field, "")


def _parse_float_values(text: str) -> list[float]:
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise argparse.ArgumentTypeError("at least one float value is required")
    return values


def _format_scale(value: float) -> str:
    return f"{value:g}".replace(".", "p").replace("-", "m")


def build_conditions(args: argparse.Namespace) -> tuple[Condition, ...]:
    if args.mode == "preset":
        return CONDITIONS

    conditions: list[Condition] = []
    if args.mode == "body-first-grid":
        first_second_scale = float(args.fixed_first_second_spring_scale)
        for attach_scale in args.attach_first_spring_scales:
            for angle_scale in args.body_axis_angle_scales:
                condition_id = (
                    f"af{_format_scale(attach_scale)}"
                    f"_axis{_format_scale(angle_scale)}"
                    f"_fs{_format_scale(first_second_scale)}"
                )
                conditions.append(
                    Condition(
                        condition_id=condition_id,
                        mode=args.mode,
                        description=(
                            "body-first distance and body-axis 90 degree grid"
                        ),
                        scales={
                            "local_attach_first_spring_scale": float(attach_scale),
                            "local_attach_first_body_axis_angle_scale": float(
                                angle_scale
                            ),
                            "local_first_second_spring_scale": first_second_scale,
                        },
                    )
                )
        return tuple(conditions)

    if args.mode == "first-second-grid":
        attach_scale = float(args.fixed_attach_first_spring_scale or 2.0)
        angle_scale = float(args.fixed_body_axis_angle_scale or 2.0)
        frame_position_scale = getattr(args, "fixed_attach_frame_position_scale", None)
        frame_tangent_scale = getattr(args, "fixed_attach_frame_tangent_scale", None)
        for first_second_scale in args.first_second_spring_scales:
            condition_id = (
                f"af{_format_scale(attach_scale)}"
                f"_axis{_format_scale(angle_scale)}"
                f"_fs{_format_scale(first_second_scale)}"
            )
            scales = {
                "local_attach_first_spring_scale": attach_scale,
                "local_attach_first_body_axis_angle_scale": angle_scale,
                "local_first_second_spring_scale": float(first_second_scale),
            }
            if frame_position_scale is not None:
                position_scale = float(frame_position_scale)
                scales["local_attach_frame_position_scale"] = position_scale
                condition_id += f"_fp{_format_scale(position_scale)}"
            if frame_tangent_scale is not None:
                tangent_scale = float(frame_tangent_scale)
                scales["local_attach_frame_tangent_scale"] = tangent_scale
                condition_id += f"_ft{_format_scale(tangent_scale)}"
            conditions.append(
                Condition(
                    condition_id=condition_id,
                    mode=args.mode,
                    description="first-second distance grid after body-first fix",
                    scales=scales,
                )
            )
        return tuple(conditions)

    if args.mode == "attach-frame-grid":
        attach_scale = float(args.fixed_attach_first_spring_scale or 1.0)
        angle_scale = float(args.fixed_body_axis_angle_scale or 1.0)
        first_second_scale = float(args.fixed_first_second_spring_scale)
        for position_scale in args.attach_frame_position_scales:
            for tangent_scale in args.attach_frame_tangent_scales:
                condition_id = (
                    f"af{_format_scale(attach_scale)}"
                    f"_axis{_format_scale(angle_scale)}"
                    f"_fs{_format_scale(first_second_scale)}"
                    f"_fp{_format_scale(position_scale)}"
                    f"_ft{_format_scale(tangent_scale)}"
                )
                conditions.append(
                    Condition(
                        condition_id=condition_id,
                        mode=args.mode,
                        description="attach-frame position and tangent grid",
                        scales={
                            "local_attach_first_spring_scale": attach_scale,
                            "local_attach_first_body_axis_angle_scale": angle_scale,
                            "local_first_second_spring_scale": first_second_scale,
                            "local_attach_frame_position_scale": float(position_scale),
                            "local_attach_frame_tangent_scale": float(tangent_scale),
                        },
                    )
                )
        return tuple(conditions)

    raise ValueError(f"Unsupported mode: {args.mode}")


def _overrides_for_condition(
    args: argparse.Namespace, condition: Condition
) -> dict[str, Any]:
    motor_overrides = {
        "torque_Nm": args.torque_nm,
        "enable_switching": False,
        "force_distribution": "root_torque_segment_couples",
    }
    motor_overrides.update(condition.scales)
    return {
        "flagella": {
            "n_flagella": args.n_flagella,
            "placement_mode": "seeded_surface",
            "initial_phase_mode": "seeded",
            "initial_helix_axis_from_rear_deg": 0.0,
        },
        "motor": motor_overrides,
        "seed": {"attach_seed": args.attach_seed, "phase_seed": args.phase_seed},
        "time": {"duration_s": args.duration_s, "dt_star": args.dt_star},
        "output_sampling": {"out_all_steps_3d": False},
        "render": {"save_frames_3d": False, "save_frames_2d": False},
    }


def _summary_row(
    cfg: SimulationConfig,
    condition: Condition,
    output_dir: Path,
    last: dict[str, str],
    rows: list[dict[str, str]],
    helix_summary: dict[str, bool | int | float | str],
) -> dict[str, str | float]:
    first_fail = _first_fail_row(rows)
    max_hook = _max_hook_row(rows)
    max_flag_bond = _max_flag_bond_row(rows)
    row: dict[str, str | float] = {
        "condition_id": condition.condition_id,
        "mode": condition.mode,
        "description": condition.description,
        "output_dir": str(output_dir),
        "duration_s": cfg.time.duration_s,
        "dt_star": cfg.dt_star,
        "torque_Nm": cfg.motor_torque_Nm,
        "n_flagella": cfg.flagella.n_flagella,
        "local_attach_first_spring_scale": cfg.motor.local_attach_first_spring_scale,
        "local_attach_first_body_axis_angle_scale": (
            cfg.motor.local_attach_first_body_axis_angle_scale
        ),
        "local_first_second_spring_scale": cfg.motor.local_first_second_spring_scale,
        "local_attach_frame_position_scale": (
            cfg.motor.local_attach_frame_position_scale
        ),
        "local_attach_frame_tangent_scale": (
            cfg.motor.local_attach_frame_tangent_scale
        ),
        "first_fail_t_s": "" if first_fail is None else first_fail.get("t_s", ""),
        "first_fail_category_nonbody": (
            "none"
            if first_fail is None
            else first_fail.get("first_fail_category_nonbody", "")
        ),
        "first_fail_hook_len_rel_err_max": (
            "" if first_fail is None else first_fail.get("hook_len_rel_err_max", "")
        ),
        "first_fail_hook_len_rel_err_max_flag_id": (
            ""
            if first_fail is None
            else first_fail.get("hook_len_rel_err_max_flag_id", "")
        ),
        "first_fail_local_attach_first_rel_err": (
            ""
            if first_fail is None
            else first_fail.get("local_attach_first_rel_err", "")
        ),
        "first_fail_local_first_second_rel_err": (
            ""
            if first_fail is None
            else first_fail.get("local_first_second_rel_err", "")
        ),
        "first_fail_local_attach_frame_position_rel_err": (
            ""
            if first_fail is None
            else first_fail.get("local_attach_frame_position_rel_err", "")
        ),
        "first_fail_local_attach_frame_position_angle_err_deg": (
            ""
            if first_fail is None
            else first_fail.get("local_attach_frame_position_angle_err_deg", "")
        ),
        "first_fail_local_attach_frame_tangent_angle_err_deg": (
            ""
            if first_fail is None
            else first_fail.get("local_attach_frame_tangent_angle_err_deg", "")
        ),
        "first_fail_flag_bond_rel_err_max": (
            "" if first_fail is None else first_fail.get("flag_bond_rel_err_max", "")
        ),
        "first_fail_flag_bond_rel_err_max_flag_id": (
            ""
            if first_fail is None
            else first_fail.get("flag_bond_rel_err_max_flag_id", "")
        ),
        "first_fail_flag_bond_rel_err_max_bead_i": (
            ""
            if first_fail is None
            else first_fail.get("flag_bond_rel_err_max_bead_i", "")
        ),
        "first_fail_flag_bond_rel_err_max_bead_j": (
            ""
            if first_fail is None
            else first_fail.get("flag_bond_rel_err_max_bead_j", "")
        ),
        "first_fail_flag_bond_rel_err_max_local_bead_i": (
            _flag_bond_local_bead_index(cfg, first_fail, "flag_bond_rel_err_max_bead_i")
        ),
        "first_fail_flag_bond_rel_err_max_local_bead_j": (
            _flag_bond_local_bead_index(cfg, first_fail, "flag_bond_rel_err_max_bead_j")
        ),
        "first_fail_flag_bond_rel_err_max_len_over_b": (
            ""
            if first_fail is None
            else first_fail.get("flag_bond_rel_err_max_len_over_b", "")
        ),
        "max_hook_len_rel_err_t_s": "" if max_hook is None else max_hook.get("t_s", ""),
        "max_hook_len_rel_err": (
            "" if max_hook is None else max_hook.get("hook_len_rel_err_max", "")
        ),
        "max_hook_len_rel_err_flag_id": (
            "" if max_hook is None else max_hook.get("hook_len_rel_err_max_flag_id", "")
        ),
        "max_hook_len_rel_err_attach_body_bead_index": (
            ""
            if max_hook is None
            else max_hook.get("hook_len_rel_err_max_attach_body_bead_index", "")
        ),
        "max_hook_len_rel_err_flag_first_bead_index": (
            ""
            if max_hook is None
            else max_hook.get("hook_len_rel_err_max_flag_first_bead_index", "")
        ),
        "max_hook_len_rel_err_len_over_b": (
            ""
            if max_hook is None
            else max_hook.get("hook_len_rel_err_max_len_over_b", "")
        ),
        "max_hook_local_attach_frame_position_rel_err": (
            ""
            if max_hook is None
            else max_hook.get("local_attach_frame_position_rel_err", "")
        ),
        "max_hook_local_attach_frame_position_angle_err_deg": (
            ""
            if max_hook is None
            else max_hook.get("local_attach_frame_position_angle_err_deg", "")
        ),
        "max_hook_local_attach_frame_tangent_angle_err_deg": (
            ""
            if max_hook is None
            else max_hook.get("local_attach_frame_tangent_angle_err_deg", "")
        ),
        "max_flag_bond_rel_err_t_s": (
            "" if max_flag_bond is None else max_flag_bond.get("t_s", "")
        ),
        "max_flag_bond_rel_err": (
            ""
            if max_flag_bond is None
            else max_flag_bond.get("flag_bond_rel_err_max", "")
        ),
        "max_flag_bond_rel_err_flag_id": (
            ""
            if max_flag_bond is None
            else max_flag_bond.get("flag_bond_rel_err_max_flag_id", "")
        ),
        "max_flag_bond_rel_err_bead_i": (
            ""
            if max_flag_bond is None
            else max_flag_bond.get("flag_bond_rel_err_max_bead_i", "")
        ),
        "max_flag_bond_rel_err_bead_j": (
            ""
            if max_flag_bond is None
            else max_flag_bond.get("flag_bond_rel_err_max_bead_j", "")
        ),
        "max_flag_bond_rel_err_local_bead_i": (
            _flag_bond_local_bead_index(
                cfg, max_flag_bond, "flag_bond_rel_err_max_bead_i"
            )
        ),
        "max_flag_bond_rel_err_local_bead_j": (
            _flag_bond_local_bead_index(
                cfg, max_flag_bond, "flag_bond_rel_err_max_bead_j"
            )
        ),
        "max_flag_bond_rel_err_len_over_b": (
            ""
            if max_flag_bond is None
            else max_flag_bond.get("flag_bond_rel_err_max_len_over_b", "")
        ),
    }
    _copy_prefixed_fields(
        row,
        first_fail,
        PROXIMAL_FLAG_BOND_REL_ERR_FIELDS,
        prefix="first_fail_",
    )
    _copy_prefixed_fields(
        row,
        max_flag_bond,
        PROXIMAL_FLAG_BOND_REL_ERR_FIELDS,
        prefix="max_",
    )
    for field in SUMMARY_FIELDS:
        if field in row:
            continue
        source_key = field.removeprefix("final_")
        row[field] = helix_summary.get(field, last.get(source_key, last.get(field, "")))
    row["flag_bond_rel_err_max_local_bead_i"] = _flag_bond_local_bead_index(
        cfg, last, "flag_bond_rel_err_max_bead_i"
    )
    row["flag_bond_rel_err_max_local_bead_j"] = _flag_bond_local_bead_index(
        cfg, last, "flag_bond_rel_err_max_bead_j"
    )
    return row


def _run_condition(
    base_cfg: dict[str, Any],
    args: argparse.Namespace,
    condition: Condition,
) -> dict[str, str | float]:
    cfg = SimulationConfig.from_dict(base_cfg).with_overrides(
        _overrides_for_condition(args, condition)
    )
    condition_dir = args.output_dir / condition.condition_id
    condition_dir.mkdir(parents=True, exist_ok=args.overwrite)
    simulator = Simulator(cfg)
    simulator.run(
        cfg.time.duration_s,
        step_summary_dir=condition_dir,
        stop_on_shape_fail=False,
        progress_interval=args.progress_interval,
    )
    step_summary_path = condition_dir / "step_summary.csv"
    rows = _read_step_rows(step_summary_path)
    helix_summary = summarize_single_flagellum_helix_retention(
        rows,
        min_steps=min(50, max(len(rows) - 1, 1)),
        min_net_abs_spin_revolutions=0.0,
    )
    if not rows:
        raise RuntimeError(f"No rows found in {step_summary_path}")
    last = rows[-1]
    return _summary_row(cfg, condition, condition_dir, last, rows, helix_summary)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=(
            "preset",
            "body-first-grid",
            "first-second-grid",
            "attach-frame-grid",
        ),
        default="preset",
        help="Condition generation mode.",
    )
    parser.add_argument("--config", type=Path, default=Path("conf/sim_swim.yaml"))
    parser.add_argument("--output-dir", type=Path, default=_default_output_dir())
    parser.add_argument("--duration-s", type=float, default=0.5)
    parser.add_argument("--dt-star", type=float, default=1.0e-4)
    parser.add_argument("--torque-nm", type=float, default=2.5e-20)
    parser.add_argument("--n-flagella", type=int, default=3)
    parser.add_argument("--attach-seed", type=int, default=0)
    parser.add_argument("--phase-seed", type=int, default=0)
    parser.add_argument(
        "--attach-first-spring-scales",
        type=_parse_float_values,
        default=_parse_float_values("1,1.25,1.5,2,3"),
        help="Comma-separated scale values for --mode body-first-grid.",
    )
    parser.add_argument(
        "--body-axis-angle-scales",
        type=_parse_float_values,
        default=_parse_float_values("1,1.25,1.5,2,3"),
        help="Comma-separated scale values for --mode body-first-grid.",
    )
    parser.add_argument(
        "--first-second-spring-scales",
        type=_parse_float_values,
        default=_parse_float_values("1,1.25,1.5,2,3"),
        help="Comma-separated scale values for --mode first-second-grid.",
    )
    parser.add_argument(
        "--attach-frame-position-scales",
        type=_parse_float_values,
        default=_parse_float_values("1,1.25,1.5,2,3"),
        help="Comma-separated position scale values for --mode attach-frame-grid.",
    )
    parser.add_argument(
        "--attach-frame-tangent-scales",
        type=_parse_float_values,
        default=_parse_float_values("1,1.25,1.5,2"),
        help="Comma-separated tangent scale values for --mode attach-frame-grid.",
    )
    parser.add_argument("--fixed-attach-first-spring-scale", type=float, default=None)
    parser.add_argument("--fixed-body-axis-angle-scale", type=float, default=None)
    parser.add_argument("--fixed-first-second-spring-scale", type=float, default=1.0)
    parser.add_argument("--fixed-attach-frame-position-scale", type=float, default=None)
    parser.add_argument("--fixed-attach-frame-tangent-scale", type=float, default=None)
    parser.add_argument("--sample-limit", type=int, default=None)
    parser.add_argument("--progress-interval", type=int, default=1000)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    base_cfg = _load_yaml(args.config)
    all_conditions = build_conditions(args)
    conditions = (
        all_conditions[: args.sample_limit] if args.sample_limit else all_conditions
    )

    if args.dry_run:
        for condition in conditions:
            print(condition.condition_id)
        return

    args.output_dir.mkdir(parents=True, exist_ok=args.overwrite)
    rows = [_run_condition(base_cfg, args, condition) for condition in conditions]
    summary_path = args.output_dir / "phase2_82_hook_scale_sweep_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(summary_path)


if __name__ == "__main__":
    main()
