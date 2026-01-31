from __future__ import annotations

import numpy as np

from flagella_sim.sim.core import Simulator
from flagella_sim.sim.params import SimulationConfig


def _make_cfg(**overrides):
    base = {
        "body": {"length_total_um": 3.0, "diameter_um": 0.8},
        "flagella": {
            "n_flagella": 5,
            "length_um": 12.0,
            "pitch_um": 2.3,
            "radius_um": 0.2,
            "filament_diameter_um": 0.02,
            "placement_mode": "uniform",
            "motor_freq_hz": 100.0,
        },
        "environment": {
            "viscosity_mpas": 1.0,
            "temperature_k": 298,
            "include_brownian": True,
            "include_gravity": False,
        },
        "time": {
            "fps_out": 50,
            "dt_sim": 0.0005,
            "duration_s": 0.05,
        },
        "render": {
            "image_size_px": 256,
            "pixel_size_um": 0.203,
            "render_flagella": False,
            "flagella_linewidth_px": 3.0,
        },
        "seed": {"global_seed": 0},
    }
    for k, v in overrides.items():
        base[k] = {**base.get(k, {}), **v} if isinstance(v, dict) else v
    return SimulationConfig.from_dict(base)


def test_smoke_runs_and_outputs_positions():
    cfg = _make_cfg()
    sim = Simulator(cfg)
    states = sim.run(0.05)
    assert len(states) >= 2
    arr = np.array([s.position_um for s in states], dtype=float)
    assert not np.isnan(arr).any()


def test_reproducibility_same_seed_matches():
    cfg = _make_cfg(seed={"global_seed": 123})
    sim1 = Simulator(cfg)
    sim2 = Simulator(cfg)
    s1 = sim1.run(0.05)
    s2 = sim2.run(0.05)
    a1 = np.array([s.position_um for s in s1])
    a2 = np.array([s.position_um for s in s2])
    assert np.allclose(a1, a2)


def test_motor_freq_increases_speed():
    cfg_low = _make_cfg(flagella={"motor_freq_hz": 50.0})
    cfg_high = _make_cfg(flagella={"motor_freq_hz": 200.0})
    sim_low = Simulator(cfg_low)
    sim_high = Simulator(cfg_high)
    s_low = sim_low.run(0.05)
    s_high = sim_high.run(0.05)
    disp_low = np.linalg.norm(np.array(s_low[-1].position_um))
    disp_high = np.linalg.norm(np.array(s_high[-1].position_um))
    assert disp_high > disp_low * 1.2


def test_brownian_msd_nonzero_when_no_flagella():
    cfg = _make_cfg(flagella={"n_flagella": 0}, seed={"global_seed": 7})
    sim = Simulator(cfg)
    states = sim.run(0.05)
    positions = np.array([s.position_um for s in states])
    diffs = positions - positions[0]
    msd = np.mean(np.sum(diffs**2, axis=1))
    assert msd > 1e-4
