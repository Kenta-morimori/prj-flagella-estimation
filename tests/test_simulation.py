from __future__ import annotations

import csv
from pathlib import Path
import numpy as np

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _make_cfg() -> SimulationConfig:
    return SimulationConfig.from_dict(
        {
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
                "discretization": {"ds_over_b": 0.58},
                "bond_L_over_b": 0.58,
                "length_over_b": 2.32,
                "helix_init": {"radius_over_b": 0.2, "pitch_over_b": 1.0},
            },
            "fluid": {"viscosity_Pa_s": 1.0e-3},
            "motor": {"torque_Nm": -1.0, "reverse_n_flagella": 1},
            "potentials": {
                "spring": {"H_over_T_over_b": 10.0, "s": 0.1},
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
                    "phi0_deg": {"normal": -60.0, "semicoiled": 65.0, "curly1": 120.0},
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
            "time": {"duration_s": 5.0e-5, "dt_s": 1.0e-3},
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
    )


def test_short_run_no_nan_inf() -> None:
    cfg = _make_cfg()
    sim = Simulator(cfg)
    states = sim.run(cfg.time.duration_s)

    assert len(states) >= 2

    arr_pos = np.array([s.position_um for s in states], dtype=float)
    arr_q = np.array([s.quaternion for s in states], dtype=float)
    arr_v = np.array([s.velocity_um_s for s in states], dtype=float)
    arr_w = np.array([s.omega_rad_s for s in states], dtype=float)

    assert np.isfinite(arr_pos).all()
    assert np.isfinite(arr_q).all()
    assert np.isfinite(arr_v).all()
    assert np.isfinite(arr_w).all()


def test_run_writes_step_summary_csv(tmp_path: Path) -> None:
    cfg = _make_cfg()
    sim = Simulator(cfg)
    sim.run(cfg.time.duration_s, sim_debug_dir=tmp_path / "sim_debug")

    step_csv = tmp_path / "sim_debug" / "step_summary.csv"
    step_full_csv = tmp_path / "sim_debug" / "step_summary_full.csv"
    assert step_csv.is_file()
    assert step_full_csv.is_file()

    with step_csv.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows
    first = rows[0]
    assert "step" in first
    assert "F_total_mean_all" in first
    assert first["brownian_enabled"] in {"False", "false", "0"}

    with step_full_csv.open("r", encoding="utf-8", newline="") as f:
        rows_full = list(csv.DictReader(f))
    assert rows_full
    first_full = rows_full[0]
    assert "F_motor_mean_flag" in first_full
    assert "F_repulsion_mean_body" in first_full
