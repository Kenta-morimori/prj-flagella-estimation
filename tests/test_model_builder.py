from __future__ import annotations

import numpy as np

from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig


def test_body_prism_bead_count_and_layout() -> None:
    cfg = SimulationConfig.from_dict(
        {
            "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
            "body": {
                "prism": {
                    "n_prism": 4,
                    "dz_over_b": 0.5,
                    "radius_over_b": 0.4,
                    "axis": "x",
                },
                "length_total_um": 2.5,
            },
            "flagella": {
                "n_flagella": 2,
                "placement_mode": "uniform",
                "discretization": {"ds_over_b": 0.5},
                "bond_L_over_b": 0.5,
                "length_over_b": 2.0,
                "helix_init": {"radius_over_b": 0.2, "pitch_over_b": 1.0},
            },
            "time": {"duration_s": 0.02, "dt_s": 1.0e-3},
            "brownian": {"enabled": False},
        }
    )

    model = ModelBuilder(cfg).build()

    assert model.body_indices.shape[0] == 4 * 6

    first_layer = model.body_layer_indices[0]
    pts = model.positions_m[first_layer] / 1e-6
    radii = np.linalg.norm(pts[:, 1:3], axis=1)
    assert np.allclose(radii, 0.4, atol=1e-6)

    x_coords = [
        np.mean(model.positions_m[layer][:, 0] / 1e-6)
        for layer in model.body_layer_indices
    ]
    diffs = np.diff(x_coords)
    assert np.allclose(diffs, 0.5, atol=1e-6)


def test_body_flag_spring_connection_exists() -> None:
    cfg = SimulationConfig.from_dict(
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
            "time": {"duration_s": 0.02, "dt_s": 1.0e-3},
            "brownian": {"enabled": False},
        }
    )
    model = ModelBuilder(cfg).build()

    body_flag_pairs = 0
    for i, j in model.spring_pairs:
        if bool(model.bead_is_body[int(i)]) ^ bool(model.bead_is_body[int(j)]):
            body_flag_pairs += 1
    assert body_flag_pairs == cfg.flagella.n_flagella
