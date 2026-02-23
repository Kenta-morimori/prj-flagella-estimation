from __future__ import annotations

import pytest

from sim_swim.sim.params import SimulationConfig


def _base_cfg() -> dict:
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
        "fluid": {"viscosity_Pa_s": 1.0e-3},
        "motor": {"torque_Nm": 4.0e-18, "reverse_n_flagella": 1},
    }


def test_dt_over_tau_input_is_rejected() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_over_tau": 0.001}

    with pytest.raises(ValueError, match="time.dt_over_tau"):
        SimulationConfig.from_dict(cfg)


def test_validate_time_scaling_rejects_non_paper_dt_star() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-7}
    sim_cfg = SimulationConfig.from_dict(cfg)

    with pytest.raises(ValueError, match="time.dt_s"):
        sim_cfg.validate_time_scaling()


def test_body_n_layers_is_derived_from_length_and_spacing() -> None:
    cfg = _base_cfg()
    cfg["body"]["length_total_um"] = 2.5
    cfg["time"] = {"duration_s": 0.1, "dt_s": 2.5e-7}
    sim_cfg = SimulationConfig.from_dict(cfg)
    assert sim_cfg.compute_body_n_layers() == 6


def test_body_n_layers_requires_integer_multiple() -> None:
    cfg = _base_cfg()
    cfg["body"]["length_total_um"] = 2.3
    cfg["time"] = {"duration_s": 0.1, "dt_s": 2.5e-7}
    sim_cfg = SimulationConfig.from_dict(cfg)

    with pytest.raises(ValueError, match="整数倍"):
        sim_cfg.compute_body_n_layers()
