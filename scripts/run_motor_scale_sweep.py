#!/usr/bin/env python3
"""Sweep motor-local scale factors and summarize observables.

This helper is intended for observation-first experiments where the user
weakens one motor-local scale at a time and checks which shape or rotation
signal degrades first.

Example:
  uv run python -m scripts.run_motor_scale_sweep \
    --target local_hook_scale --values 8,4,2,1,0.5,0.25,0.0
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _parse_values(text: str) -> list[float]:
    values: list[float] = []
    for raw_item in text.split(","):
        item = raw_item.strip()
        if not item:
            continue
        values.append(float(item))
    if not values:
        raise ValueError("At least one value is required.")
    return values


def _base_cfg() -> dict[str, Any]:
    return {
        "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
        "body": {
            "prism": {
                "n_prism": 3,
                "dz_over_b": 0.5,
                "radius_over_b": 0.5,
                "axis": "x",
            },
            "length_total_um": 2.0,
        },
        "flagella": {
            "n_flagella": 1,
            "placement_mode": "uniform",
            "init_mode": "paper_table1",
            "stub_mode": "minimal_basal_stub",
            "n_beads_per_flagellum": 11,
            "discretization": {"ds_over_b": 0.58},
            "bond_L_over_b": 0.58,
            "length_over_b": 5.8,
            "helix_init": {"radius_over_b": 0.25, "pitch_over_b": 2.5},
        },
        "fluid": {"viscosity_Pa_s": 1.0e-3},
        "motor": {
            "torque_Nm": 4.0e-21,
            "reverse_n_flagella": 1,
            "enable_switching": False,
            "torque_ramp_enabled": False,
            "torque_ramp_duration_s": 0.0,
            "torque_for_forces_override_Nm": 0.0,
            "local_hook_scale": 8.0,
            "local_spring_scale": 5.0,
            "local_bend_scale": 4.0,
            "local_torsion_scale": 4.0,
        },
        "potentials": {
            "spring": {"H_over_T_over_b": 20.0, "s": 0.2},
            "bend": {
                "kb_over_T": 20.0,
                "theta0_deg": {
                    "normal": 142.0,
                    "semicoiled": 90.0,
                    "curly1": 105.0,
                },
            },
            "torsion": {
                "kt_over_T": 10.0,
                "fd_eps_over_b": 1.0e-3,
                "phi0_deg": {
                    "normal": -60.0,
                    "semicoiled": 65.0,
                    "curly1": 120.0,
                },
            },
            "spring_spring_repulsion": {
                "A_ss_over_T": 1.0,
                "a_ss_over_b": 0.2,
                "cutoff_over_b": 0.2,
            },
        },
        "hook": {"enabled": True, "threshold_deg": 90.0, "kb_over_T": 20.0},
        "run_tumble": {
            "run_tau": 20.0,
            "tumble_tau": 8.0,
            "semicoiled_tau": 4.0,
            "curly1_tau": 4.0,
        },
        "time": {"duration_s": 0.01, "dt_s": 1.0e-3},
        "output_sampling": {"out_all_steps_3d": True, "fps_out_2d": 25.0},
        "brownian": {
            "enabled": False,
            "temperature_K": 298.0,
            "method": "cholesky",
            "jitter": 1.0e-20,
        },
        "render": {
            "image_size_px": 128,
            "pixel_size_um": 0.203,
            "flagella_linewidth_px": 2.0,
            "render_flagella": False,
            "save_frames_3d": False,
            "follow_camera_3d": True,
            "view_range_um": 8.0,
            "timestamp_3d": False,
            "label_flagella": True,
            "follow_camera_2d": False,
            "save_frames_2d": False,
        },
        "seed": {"global_seed": 0},
        "output": {"base_dir": "outputs"},
    }


def _collect_last_step_metrics(csv_path: Path) -> dict[str, str]:
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"No step rows found in {csv_path}")
    return rows[-1]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a sweep over motor-local scale factors."
    )
    parser.add_argument(
        "--target",
        choices=(
            "local_hook_scale",
            "local_spring_scale",
            "local_bend_scale",
            "local_torsion_scale",
        ),
        default="local_hook_scale",
        help="Motor-local scale to sweep.",
    )
    parser.add_argument(
        "--values",
        type=_parse_values,
        default=_parse_values("8,4,2,1,0.5,0.25,0.0"),
        help="Comma-separated sweep values.",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=0.01,
        help="Simulation duration in seconds.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="outputs/motor_scale_sweep",
        help="Directory for sweep results.",
    )
    args = parser.parse_args()

    sweep_dir = Path(args.output_dir)
    sweep_dir.mkdir(parents=True, exist_ok=True)

    summary_path = sweep_dir / f"{args.target}_sweep_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "target",
                "value",
                "output_dir",
                "pos_all_finite",
                "any_nan",
                "any_inf",
                "flag_root_azimuth_deg",
                "flag_phase_deg",
                "flag_phase_rate_hz",
                "flag_body_phase_diff_deg",
                "local_attach_first_rel_err",
                "local_first_second_rel_err",
                "hook_angle_err_max_deg",
                "hook_len_rel_err_max",
                "flag_bond_rel_err_max",
                "flag_bend_err_max_deg",
                "flag_torsion_err_max_deg",
                "motor_split_residual_norm",
                "motor_degenerate_axis_count",
                "motor_split_rank_deficient_count",
                "motor_bond_length_clipped_count",
            ],
        )
        writer.writeheader()

        for value in args.values:
            cfg_dict = _base_cfg()
            cfg_dict["motor"][args.target] = float(value)
            cfg_dict["time"]["duration_s"] = float(args.duration)
            cfg = SimulationConfig.from_dict(cfg_dict)

            run_dir = sweep_dir / f"{args.target}_{value:g}"
            run_dir.mkdir(parents=True, exist_ok=True)
            sim = Simulator(cfg)
            sim.run(cfg.time.duration_s, step_summary_dir=run_dir)
            last = _collect_last_step_metrics(run_dir / "step_summary.csv")

            writer.writerow(
                {
                    "target": args.target,
                    "value": float(value),
                    "output_dir": str(run_dir),
                    "pos_all_finite": last.get("pos_all_finite", ""),
                    "any_nan": last.get("any_nan", ""),
                    "any_inf": last.get("any_inf", ""),
                    "flag_root_azimuth_deg": last.get("flag_root_azimuth_deg", ""),
                    "flag_phase_deg": last.get("flag_phase_deg", ""),
                    "flag_phase_rate_hz": last.get("flag_phase_rate_hz", ""),
                    "flag_body_phase_diff_deg": last.get(
                        "flag_body_phase_diff_deg", ""
                    ),
                    "local_attach_first_rel_err": last.get(
                        "local_attach_first_rel_err", ""
                    ),
                    "local_first_second_rel_err": last.get(
                        "local_first_second_rel_err", ""
                    ),
                    "hook_angle_err_max_deg": last.get("hook_angle_err_max_deg", ""),
                    "hook_len_rel_err_max": last.get("hook_len_rel_err_max", ""),
                    "flag_bond_rel_err_max": last.get("flag_bond_rel_err_max", ""),
                    "flag_bend_err_max_deg": last.get("flag_bend_err_max_deg", ""),
                    "flag_torsion_err_max_deg": last.get(
                        "flag_torsion_err_max_deg", ""
                    ),
                    "motor_split_residual_norm": last.get(
                        "motor_split_residual_norm", ""
                    ),
                    "motor_degenerate_axis_count": last.get(
                        "motor_degenerate_axis_count", ""
                    ),
                    "motor_split_rank_deficient_count": last.get(
                        "motor_split_rank_deficient_count", ""
                    ),
                    "motor_bond_length_clipped_count": last.get(
                        "motor_bond_length_clipped_count", ""
                    ),
                }
            )

    print(f"Sweep summary saved to {summary_path}")


if __name__ == "__main__":
    main()
