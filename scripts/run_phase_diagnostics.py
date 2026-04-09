#!/usr/bin/env python3
"""Phase別の短時間 diagnostics レポート実行。

使用例:
  uv run python scripts/run_phase_diagnostics.py phase0a
  uv run python scripts/run_phase_diagnostics.py phase0b
  uv run python scripts/run_phase_diagnostics.py phase1
  uv run python scripts/run_phase_diagnostics.py phase2
  uv run python scripts/run_phase_diagnostics.py phase3
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _make_phase_cfg(phase: str, duration_s: float = 0.01) -> SimulationConfig:
    """各フェーズの canonical config を生成する。"""
    base_cfg = {
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
            "n_flagella": 3,
            "placement_mode": "uniform",
            "init_mode": "paper_table1",
            "stub_mode": "full_flagella",
            "n_beads_per_flagellum": 11,
            "discretization": {"ds_over_b": 0.58},
            "bond_L_over_b": 0.58,
            "length_over_b": 2.32,
            "helix_init": {"radius_over_b": 0.25, "pitch_over_b": 2.5},
        },
        "fluid": {"viscosity_Pa_s": 1.0e-3},
        "motor": {"torque_Nm": -1.0, "reverse_n_flagella": 1},
        "potentials": {
            "spring": {
                "H_over_T_over_b": 20.0,
                "s": 0.2,
            },  # Trial: Phase B stabilization (H:10->20, s:0.1->0.2)
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
        "time": {"duration_s": duration_s, "dt_s": 1.0e-3},
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
            "render_flagella": True,
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

    if phase == "phase0a":
        # body-only static
        base_cfg["flagella"]["n_flagella"] = 0
        base_cfg["motor"]["torque_Nm"] = 0.0
    elif phase == "phase0b":
        # body + hook + minimal_basal_stub static
        base_cfg["flagella"]["n_flagella"] = 1
        base_cfg["flagella"]["stub_mode"] = "minimal_basal_stub"
        base_cfg["motor"]["torque_Nm"] = 0.0
    elif phase == "phase1":
        # body-only + surrogate torque
        base_cfg["flagella"]["n_flagella"] = 0
        base_cfg["motor"]["torque_Nm"] = 0.0
        base_cfg["body_equiv_load"] = {
            "enabled": True,
            "mode": "pure_couple",
            "target_torque_Nm": 1.0e-20,
            "target_force_N": 0.0,
            "attach_region_id": 0,
        }
    elif phase == "phase2":
        # minimal_basal_stub + actual motor
        base_cfg["flagella"]["n_flagella"] = 1
        base_cfg["flagella"]["stub_mode"] = "minimal_basal_stub"
        base_cfg["motor"]["torque_Nm"] = 4.0e-21
    elif phase == "phase3":
        # full_flagella + actual motor
        base_cfg["flagella"]["n_flagella"] = 1
        base_cfg["flagella"]["stub_mode"] = "full_flagella"
        base_cfg["motor"]["torque_Nm"] = 4.0e-21
    else:
        raise ValueError(f"Unknown phase: {phase}")

    return SimulationConfig.from_dict(base_cfg)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run phase-specific diagnostics report."
    )
    parser.add_argument(
        "phase",
        choices=["phase0a", "phase0b", "phase1", "phase2", "phase3"],
        help="Phase to run",
    )
    parser.add_argument(
        "--duration",
        type=float,
        default=0.01,
        help="Simulation duration in seconds (default: 0.01)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="outputs/phase_diagnostics",
        help="Output directory",
    )
    args = parser.parse_args()

    print(f"Running {args.phase} diagnostics (duration={args.duration:.3f} s)...")

    cfg = _make_phase_cfg(args.phase, duration_s=args.duration)
    sim = Simulator(cfg)

    output_dir = Path(args.output_dir) / args.phase
    output_dir.mkdir(parents=True, exist_ok=True)

    print("  Simulation starting...")
    sim.run(cfg.time.duration_s, step_summary_dir=output_dir)

    # Load and summarize
    csv_path = output_dir / "step_summary.csv"
    if not csv_path.exists():
        print("ERROR: step_summary.csv not found")
        return

    with csv_path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))

    if not rows:
        print("ERROR: No data in CSV")
        return

    print(f"\n{args.phase} Report:")
    print(f"  Duration: {args.duration:.3f} s")
    print(f"  Steps: {len(rows)}")
    print("  Config settings:")
    print(f"    n_flagella: {cfg.flagella.n_flagella}")
    if cfg.flagella.n_flagella > 0:
        print(f"    stub_mode: {cfg.flagella.stub_mode}")
    print(f"    motor.torque_Nm: {cfg.motor.torque_Nm}")
    if cfg.body_equiv_load.enabled:
        print(f"    body_equiv_load.mode: {cfg.body_equiv_load.mode}")

    first = rows[0]
    last = rows[-1]

    print("\n  Status (first -> last):")
    print(
        f"    pos_all_finite: {first.get('pos_all_finite', 'N/A')} -> "
        f"{last.get('pos_all_finite', 'N/A')}"
    )
    print(f"    any_nan: {first.get('any_nan', 'N/A')} -> {last.get('any_nan', 'N/A')}")
    print(f"    any_inf: {first.get('any_inf', 'N/A')} -> {last.get('any_inf', 'N/A')}")

    print("\n  Key diagnostics (last step):")

    if cfg.flagella.n_flagella > 0:
        try:
            attach_err = float(last.get("local_attach_first_rel_err", 0))
            first_second_err = float(last.get("local_first_second_rel_err", 0))
            print(f"    local_attach_first_rel_err: {attach_err:.3%}")
            print(f"    local_first_second_rel_err: {first_second_err:.3%}")
        except (ValueError, TypeError):
            pass

        if cfg.motor.torque_Nm != 0.0:
            try:
                bond_err = float(last.get("flag_bond_rel_err_max", 0))
                print(f"    flag_bond_rel_err_max: {bond_err:.3%}")
            except (ValueError, TypeError):
                pass

    try:
        f_total = float(last.get("F_total_mean_all", 0))
        print(f"    F_total_mean_all: {f_total:.3e} N")
    except (ValueError, TypeError):
        pass

    print(f"\n  Output: {output_dir}")
    print(f"  CSV: {csv_path}")


if __name__ == "__main__":
    main()
