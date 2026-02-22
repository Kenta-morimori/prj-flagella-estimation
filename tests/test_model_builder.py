from __future__ import annotations

from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig


def test_ds_generates_expected_bead_counts() -> None:
    cfg = SimulationConfig.from_dict(
        {
            "discretization": {"ds_um": 0.5},
            "body": {"length_total_um": 3.0, "diameter_um": 0.8},
            "flagella": {"n_flagella": 2, "length_um": 2.2},
            "time": {"dt": 0.002, "fps_out": 10, "duration_s": 0.01},
            "brownian": {"enabled": False},
        }
    )

    model = ModelBuilder(cfg).build()

    n_body_expected = int(3.0 // 0.5) + 1
    n_flag_expected = int(2.2 // 0.5) + 1

    assert model.body_indices.shape[0] == n_body_expected
    assert all(idx.shape[0] == n_flag_expected for idx in model.flagella_indices)
