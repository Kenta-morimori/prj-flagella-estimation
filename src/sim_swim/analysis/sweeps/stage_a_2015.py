"""Run the Issue #168 Stage A validation campaign."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
import math
from pathlib import Path
import time
import traceback
from typing import Any
from zoneinfo import ZoneInfo

from sim_swim.analysis.flagella_count_behavior import (
    save_state_archive,
    write_trajectory_csv,
)
from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.analysis.sweeps.shape_stability_grid import (
    Condition,
    _axis_center_phase_summary,
    _read_step_rows,
    _summary_row,
)
from sim_swim.core.run_context import init_run
from sim_swim.sim.body_shape_gate import summarize_body_shape_diagnostics_csv
from sim_swim.sim.core import Simulator
from sim_swim.sim.helix_retention_gate import summarize_single_flagellum_helix_retention
from sim_swim.sim.params import SimulationConfig

STAGE_CONTRACT = {
    "motor_off": {"duration_tau": 0.1, "motor_enabled": False},
    "motor_on": {"duration_tau": 1.0, "motor_enabled": True},
}
PROFILE_DEFAULTS = {
    "project": Path("conf/sim_swim_2015.yaml"),
    "paper": Path("conf/sim_swim_2015_paper.yaml"),
}


def _parse_csv_names(raw: str) -> list[str]:
    values = [item.strip() for item in raw.split(",") if item.strip()]
    unknown = sorted(set(values) - set(PROFILE_DEFAULTS))
    if not values or unknown:
        raise argparse.ArgumentTypeError(
            "profiles must contain project and/or paper; unknown=" + ",".join(unknown)
        )
    return values


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=sorted(STAGE_CONTRACT), required=True)
    parser.add_argument(
        "--project-config", type=Path, default=PROFILE_DEFAULTS["project"]
    )
    parser.add_argument("--paper-config", type=Path, default=PROFILE_DEFAULTS["paper"])
    parser.add_argument(
        "--profiles", type=_parse_csv_names, default=["project", "paper"]
    )
    parser.add_argument("--output-base-dir", type=Path, default=Path("outputs"))
    parser.add_argument("--state-sample-count", type=int, default=201)
    parser.add_argument("--progress-interval", type=int, default=1000)
    parser.add_argument("--dt-star", type=float, default=1.0e-5)
    parser.add_argument("--duration-tau", type=float, default=None)
    parser.add_argument("--comparison-role", default="canonical_stage_a")
    parser.add_argument("--diagonal-braces-enabled", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--smoke-steps", type=int, default=None)
    return parser.parse_args(argv)


def _max_numeric(rows: list[dict[str, str]], field: str) -> float:
    values: list[float] = []
    for row in rows:
        try:
            value = float(row.get(field, "nan"))
        except (TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
    return max(values) if values else float("nan")


def _successful_row(
    cfg: SimulationConfig,
    *,
    profile_name: str,
    condition_dir: Path,
    rows: list[dict[str, str]],
    expected_steps: int,
    elapsed_s: float,
    sampled_state_count: int,
    sample_interval_steps: int,
) -> dict[str, Any]:
    helix_summary = summarize_single_flagellum_helix_retention(
        rows,
        min_steps=min(50, max(len(rows) - 1, 1)),
        min_median_abs_spin_rate_hz=0.0,
        min_net_abs_spin_revolutions=0.0,
        min_direction_consistency=0.0,
        min_helix_spin_fit_r2=0.0,
    )
    axis_summary = _axis_center_phase_summary(
        _read_step_rows(condition_dir / "flag_helix_axis_diagnostics.csv")
    )
    base = _summary_row(
        cfg,
        Condition(
            condition_id=profile_name,
            mode="2015-stage-a",
            description=f"2015 {profile_name} profile",
            scales={},
        ),
        condition_dir,
        rows[-1],
        rows,
        helix_summary,
        axis_summary,
    )
    base.update(
        summarize_body_shape_diagnostics_csv(
            condition_dir / "body_constraint_diagnostics.csv"
        )
    )
    completed_steps = len(rows)
    base.update(
        {
            "status": "completed"
            if completed_steps == expected_steps
            else "incomplete",
            "profile": profile_name,
            "expected_steps": expected_steps,
            "completed_steps": completed_steps,
            "completion_pass": completed_steps == expected_steps,
            "finite_pass_all": all(
                str(row.get("finite_pass", "")).lower() == "true" for row in rows
            ),
            "wall_time_s": elapsed_s,
            "steps_per_s": completed_steps / max(elapsed_s, 1.0e-30),
            "sampled_state_count": sampled_state_count,
            "state_sample_interval_steps": sample_interval_steps,
            "dt_star": cfg.dt_star,
            "duration_tau": cfg.duration_star,
            "max_hook_angle_err_deg": _max_numeric(rows, "hook_angle_err_max_deg"),
            "max_flag_helix_radius_abs_err_over_b": _max_numeric(
                rows, "flag_helix_radius_abs_err_over_b_max"
            ),
            "max_flag_helix_pitch_rel_err": _max_numeric(
                rows, "flag_helix_pitch_rel_err_max"
            ),
            "max_net_force_residual_ratio": _max_numeric(
                rows, "net_force_residual_ratio"
            ),
            "max_net_torque_residual_ratio": _max_numeric(
                rows, "net_torque_residual_ratio"
            ),
            "max_motor_force_balance_residual_ratio": _max_numeric(
                rows, "motor_force_balance_residual_ratio"
            ),
            "max_motor_torque_balance_residual_ratio": _max_numeric(
                rows, "motor_torque_balance_residual_ratio"
            ),
            "error_type": "",
            "error_message": "",
        }
    )
    return base


def _write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_stage_a(argv: list[str] | None = None) -> Path:
    args = _parse_args(argv)
    contract = STAGE_CONTRACT[args.stage]
    if not math.isfinite(args.dt_star) or args.dt_star <= 0.0:
        raise ValueError("dt_star must be finite and positive")
    if args.duration_tau is not None and (
        not math.isfinite(args.duration_tau) or args.duration_tau <= 0.0
    ):
        raise ValueError("duration_tau must be finite and positive")
    if args.smoke_steps is not None and args.smoke_steps <= 0:
        raise ValueError("smoke_steps must be positive")
    duration_tau = (
        float(args.smoke_steps) * args.dt_star
        if args.smoke_steps is not None
        else float(
            contract["duration_tau"] if args.duration_tau is None else args.duration_tau
        )
    )
    if args.state_sample_count < 2:
        raise ValueError("state_sample_count must be >= 2")

    config_paths = {
        "project": args.project_config,
        "paper": args.paper_config,
    }
    if args.dry_run:
        for profile_name in args.profiles:
            print(
                f"{profile_name}\tstage={args.stage}\t"
                f"duration_tau={duration_tau}\t"
                f"dt_star={args.dt_star}\t"
                f"motor_enabled={contract['motor_enabled']}\t"
                f"comparison_role={args.comparison_role}\t"
                f"config={config_paths[profile_name]}"
            )
        return Path()
    first_cfg = SimulationConfig.from_dict(load_yaml(config_paths[args.profiles[0]]))
    ctx = init_run(
        base_dir=args.output_base_dir,
        timestamp_subdir=True,
        overwrite=args.overwrite,
        input_info={
            "kind": "stage_a_2015",
            "stage": args.stage,
            "profiles": args.profiles,
            "diagonal_braces_enabled": args.diagonal_braces_enabled,
            "smoke_steps": args.smoke_steps,
            "dt_star": args.dt_star,
            "duration_tau": duration_tau,
            "comparison_role": args.comparison_role,
        },
        source_config_path=config_paths[args.profiles[0]],
        model_profile=first_cfg.model_profile_manifest(),
    )
    logger = ctx.logger
    logger.info("Issue #168 Stage A start: stage=%s", args.stage)

    summary_rows: list[dict[str, Any]] = []
    condition_records: list[dict[str, Any]] = []
    summary_path = ctx.out.root / "summary.csv"
    for profile_name in args.profiles:
        config_path = config_paths[profile_name]
        condition_dir = ctx.out.root / profile_name
        condition_dir.mkdir(parents=True, exist_ok=False)
        raw = load_yaml(config_path)
        cfg = SimulationConfig.from_dict(raw).with_overrides(
            {
                "time": {
                    "duration": {
                        "value": duration_tau,
                        "unit": "tau",
                    },
                    "integration": {"dt_star": args.dt_star},
                },
                "motor": {
                    "enabled": contract["motor_enabled"],
                    "enable_switching": False,
                },
                "body": {
                    "prism": {"diagonal_braces_enabled": args.diagonal_braces_enabled}
                },
            }
        )
        expected_steps = cfg.total_steps
        sample_interval = max(
            1,
            int(math.ceil(expected_steps / max(args.state_sample_count - 1, 1))),
        )
        started = time.perf_counter()
        row: dict[str, Any]
        implementation_manifest: dict[str, Any] | None = None
        try:
            simulator = Simulator(cfg)
            states = simulator.run(
                cfg.time.duration_s,
                logger=logger,
                progress_interval=args.progress_interval,
                step_summary_dir=condition_dir,
                stop_on_shape_fail=False,
                flush_interval_steps=100,
                record_body_diagnostics=True,
                record_body_local_diagnostics=False,
                state_sample_interval_steps=sample_interval,
            )
            elapsed_s = time.perf_counter() - started
            save_state_archive(condition_dir / "state_archive.npz", states)
            write_trajectory_csv(condition_dir / "trajectory.csv", states)
            step_rows = _read_step_rows(condition_dir / "step_summary.csv")
            if not step_rows:
                raise RuntimeError("step_summary.csv has no rows")
            row = _successful_row(
                cfg,
                profile_name=profile_name,
                condition_dir=condition_dir,
                rows=step_rows,
                expected_steps=expected_steps,
                elapsed_s=elapsed_s,
                sampled_state_count=len(states),
                sample_interval_steps=sample_interval,
            )
            implementation_manifest = simulator.implementation_manifest()
        except Exception as exc:  # noqa: BLE001 - campaign must preserve other conditions
            elapsed_s = time.perf_counter() - started
            failure = {
                "status": "exception",
                "profile": profile_name,
                "error_type": type(exc).__name__,
                "error_message": str(exc),
                "traceback": traceback.format_exc(),
            }
            (condition_dir / "failure.json").write_text(
                json.dumps(failure, ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )
            partial_rows = _read_step_rows(condition_dir / "step_summary.csv")
            row = {
                "condition_id": profile_name,
                "status": "exception",
                "profile": profile_name,
                "output_dir": str(condition_dir),
                "expected_steps": expected_steps,
                "completed_steps": len(partial_rows),
                "completion_pass": False,
                "finite_pass_all": False,
                "wall_time_s": elapsed_s,
                "steps_per_s": len(partial_rows) / max(elapsed_s, 1.0e-30),
                "sampled_state_count": 0,
                "state_sample_interval_steps": sample_interval,
                "dt_star": cfg.dt_star,
                "duration_tau": cfg.duration_star,
                "error_type": type(exc).__name__,
                "error_message": str(exc),
            }
            logger.exception("Stage A condition failed: profile=%s", profile_name)

        summary_rows.append(row)
        _write_summary(summary_path, summary_rows)
        condition_record = {
            "condition_id": profile_name,
            "source_config_path": str(config_path),
            "config_overrides": {
                "time.duration": {"value": duration_tau, "unit": "tau"},
                "time.integration.dt_star": args.dt_star,
                "motor.enabled": contract["motor_enabled"],
                "motor.enable_switching": False,
                "body.prism.diagonal_braces_enabled": (args.diagonal_braces_enabled),
            },
            "output_dir": str(condition_dir),
            "status": row["status"],
            "time": cfg.time_manifest(),
            "performance": {
                "expected_steps": expected_steps,
                "completed_steps": row["completed_steps"],
                "wall_time_s": row["wall_time_s"],
                "steps_per_s": row["steps_per_s"],
            },
            "comparison_scales": {
                "b_um": cfg.scale.b_um,
                "body_beads": cfg.model_profile.body_beads,
            },
        }
        if implementation_manifest is not None:
            condition_record.update(implementation_manifest)
        condition_records.append(condition_record)

    performance = {
        "kind": "stage_a_2015_performance",
        "stage": args.stage,
        "dt_star": args.dt_star,
        "duration_tau": duration_tau,
        "comparison_role": args.comparison_role,
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "short_profile_reference": {
            "condition": "2015 paper motor-on, 1 step, cProfile",
            "date": "2026-08-03",
            "seconds_per_step_with_profiler": 0.500,
            "primary_bottleneck": "compute_segment_repulsion_forces",
            "cumulative_time_fraction": {
                "segment_repulsion": 0.706,
                "torsion": 0.144,
                "rpy_mobility": 0.110,
            },
            "note": "Profiler overhead is included; long-run wall time is recorded per condition.",
        },
        "conditions": [
            record["performance"] | {"profile": record["condition_id"]}
            for record in condition_records
        ],
    }
    performance_path = ctx.out.root / "performance.json"
    performance_path.write_text(
        json.dumps(performance, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    run_manifest = {
        "kind": "stage_a_2015",
        "issue": 168,
        "stage": args.stage,
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "duration_tau": duration_tau,
        "dt_star": args.dt_star,
        "comparison_role": args.comparison_role,
        "motor_enabled": contract["motor_enabled"],
        "diagonal_braces_enabled": args.diagonal_braces_enabled,
        "state_sample_count_target": args.state_sample_count,
        "smoke_steps": args.smoke_steps,
        "base_config": str(config_paths[args.profiles[0]]),
        "output_root": str(ctx.out.root),
        "summary_csv": str(summary_path),
        "performance_json": str(performance_path),
        "git": {
            "commit": ctx.git.commit,
            "commit_short": ctx.git.commit_short,
            "branch": ctx.git.branch,
            "is_clean": ctx.git.is_clean,
        },
        "condition_order": list(args.profiles),
        "conditions": condition_records,
    }
    run_manifest_path = ctx.out.root / "run_manifest.json"
    run_manifest_path.write_text(
        json.dumps(run_manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    manifest_path = ctx.out.root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["kind"] = "stage_a_2015"
    manifest["stage"] = args.stage
    manifest["outputs"].update(
        {
            "summary_csv": str(summary_path),
            "performance_json": str(performance_path),
            "run_manifest_json": str(run_manifest_path),
        }
    )
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    logger.info("Wrote Stage A summary: %s", summary_path)
    print(summary_path)
    return summary_path


def main(argv: list[str] | None = None) -> None:
    run_stage_a(argv)


if __name__ == "__main__":
    main()
