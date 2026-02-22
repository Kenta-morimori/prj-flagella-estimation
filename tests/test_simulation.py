from __future__ import annotations

import numpy as np

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _make_cfg() -> SimulationConfig:
    return SimulationConfig.from_dict(
        {
            "discretization": {"ds_um": 0.8},
            "body": {"length_total_um": 3.0, "diameter_um": 0.8, "bond_L_um": None},
            "flagella": {
                "n_flagella": 3,
                "length_um": 5.6,
                "bond_L_um": None,
                "pitch_um": 2.3,
                "radius_um": 0.2,
                "filament_diameter_um": 0.02,
                "placement_mode": "uniform",
            },
            "scale": {"bead_radius_a_um": 0.15},
            "fluid": {"viscosity_Pa_s": 1.0e-3},
            "motor": {"torque_Nm": 4.0e-18, "reverse_n_flagella": 1},
            "potentials": {
                "spring": {"H": 1.0e-4, "s_um": 0.15},
                "bend": {
                    "kb": 2.0e-19,
                    "theta0_deg": {"normal": 25.0, "semicoiled": 55.0, "curly1": 75.0},
                },
                "torsion": {
                    "kt": 2.0e-19,
                    "phi0_deg": {"normal": 15.0, "semicoiled": 95.0, "curly1": 145.0},
                },
                "spring_spring_repulsion": {
                    "A_ss": 2.0e-19,
                    "a_ss_um": 0.1,
                    "cutoff_um": 0.6,
                },
            },
            "hook": {"enabled": True, "kb": 8.0e-20, "threshold_deg": 90.0},
            "run_tumble": {
                "run_tau": 0.2,
                "tumble_tau": 0.08,
                "semicoiled_tau": 0.03,
                "curly1_tau": 0.03,
            },
            "brownian": {
                "enabled": False,
                "temperature_K": 298.0,
                "method": "cholesky",
                "jitter": 1.0e-20,
            },
            "time": {"dt": 0.002, "fps_out": 20.0, "duration_s": 0.04},
            "render": {
                "image_size_px": 128,
                "pixel_size_um": 0.203,
                "render_flagella": False,
                "flagella_linewidth_px": 2.0,
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
