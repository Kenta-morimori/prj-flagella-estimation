from __future__ import annotations

import numpy as np
import pytest

from sim_swim.model.builder import ModelBuilder
from sim_swim.model.types import SimModel
from sim_swim.sim.params import SimulationConfig


def _make_cfg(
    *,
    n_prism: int = 3,
    n_flagella: int = 3,
    dz_over_b: float = 0.5,
    radius_over_b: float = 0.5,
    body_length_um: float = 2.0,
    ds_over_b: float = 0.58,
) -> SimulationConfig:
    return SimulationConfig.from_dict(
        {
            "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
            "body": {
                "prism": {
                    "n_prism": n_prism,
                    "dz_over_b": dz_over_b,
                    "radius_over_b": radius_over_b,
                    "axis": "x",
                },
                "length_total_um": body_length_um,
            },
            "flagella": {
                "n_flagella": n_flagella,
                "placement_mode": "uniform",
                "discretization": {"ds_over_b": ds_over_b},
                "bond_L_over_b": ds_over_b,
                "length_over_b": 2.32,
                "helix_init": {"radius_over_b": 0.2, "pitch_over_b": 1.0},
            },
            "time": {"duration_s": 0.02, "dt_s": 1.0e-3},
            "brownian": {"enabled": False},
        }
    )


def _body_flag_pairs(model: SimModel) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for i, j in model.spring_pairs:
        ii = int(i)
        jj = int(j)
        i_is_body = bool(model.bead_is_body[ii])
        j_is_body = bool(model.bead_is_body[jj])
        if i_is_body ^ j_is_body:
            body_idx = ii if i_is_body else jj
            flag_idx = jj if i_is_body else ii
            pairs.append((body_idx, flag_idx))
    return pairs


def test_body_prism_bead_count_and_layout() -> None:
    cfg = _make_cfg(
        n_prism=3,
        n_flagella=2,
        dz_over_b=0.5,
        radius_over_b=0.4,
        body_length_um=2.5,
        ds_over_b=0.5,
    )

    model = ModelBuilder(cfg).build()

    assert model.body_indices.shape[0] == 3 * 6

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


def test_mvp_requires_n_prism_eq_3() -> None:
    cfg = _make_cfg(n_prism=4)
    with pytest.raises(ValueError, match="MVP: body.prism.n_prism must be 3"):
        ModelBuilder(cfg).build()


@pytest.mark.parametrize("n_flagella", [-1, 4])
def test_mvp_requires_n_flagella_in_range(n_flagella: int) -> None:
    cfg = _make_cfg(n_flagella=n_flagella)
    with pytest.raises(
        ValueError, match=r"MVP: flagella.n_flagella must be in \[0,3\]"
    ):
        ModelBuilder(cfg).build()


def test_body_flag_spring_connection_exists() -> None:
    cfg = _make_cfg(n_flagella=3)
    model = ModelBuilder(cfg).build()

    assert len(_body_flag_pairs(model)) == cfg.flagella.n_flagella


def test_flagella_attach_to_center_layer() -> None:
    cfg = _make_cfg(n_flagella=3)
    model = ModelBuilder(cfg).build()

    center_layer = set(
        model.body_layer_indices[len(model.body_layer_indices) // 2].tolist()
    )
    last_layer = set(model.body_layer_indices[-1].tolist())

    attached_body_indices: set[int] = set()
    for body_idx, _ in _body_flag_pairs(model):
        attached_body_indices.add(body_idx)

    assert attached_body_indices
    assert attached_body_indices.issubset(center_layer)
    assert attached_body_indices.isdisjoint(last_layer)


def test_body_flag_hook_length_is_initialized_to_0p25b() -> None:
    cfg = _make_cfg(n_flagella=3)
    model = ModelBuilder(cfg).build()

    body_flag_pairs = _body_flag_pairs(model)
    dists_m = np.asarray(
        [
            np.linalg.norm(model.positions_m[flag_idx] - model.positions_m[body_idx])
            for body_idx, flag_idx in body_flag_pairs
        ],
        dtype=float,
    )
    assert dists_m.shape[0] == cfg.flagella.n_flagella
    assert np.allclose(dists_m, 0.25 * cfg.b_m, atol=1e-12)


def test_spring_rest_lengths_match_initial_distances() -> None:
    cfg = _make_cfg(n_flagella=3)
    model = ModelBuilder(cfg).build()

    spring_pairs = model.spring_pairs
    actual = np.linalg.norm(
        model.positions_m[spring_pairs[:, 1]] - model.positions_m[spring_pairs[:, 0]],
        axis=1,
    )
    assert np.allclose(model.spring_rest_lengths_m, actual, atol=1e-15)
