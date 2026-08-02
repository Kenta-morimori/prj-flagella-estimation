from __future__ import annotations

import math

import numpy as np
import pytest

from sim_swim.dynamics.forces import (
    compute_bending_forces,
    compute_hook_forces,
    compute_segment_repulsion_forces,
    compute_spring_forces,
    compute_torsion_forces,
)
from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.model.builder import ModelBuilder, _segment_pairs_without_neighbors
from sim_swim.sim.params import SimulationConfig

pytestmark = pytest.mark.light


def _spring_force(
    distance: float,
    *,
    rest_length: float = 0.58,
    h_const: float = 2.0,
    s_limit: float = 0.1,
    formulation: str,
) -> np.ndarray:
    return compute_spring_forces(
        positions_m=np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]]),
        spring_pairs=np.array([[0, 1]], dtype=int),
        spring_rest_lengths_m=np.array([rest_length]),
        h_const=h_const,
        s=s_limit,
        b_m=1.0,
        formulation=formulation,
    )


@pytest.mark.parametrize("formulation", ["legacy", "fene_fraenkel"])
def test_spring_force_is_zero_at_rest_length(formulation: str) -> None:
    forces = _spring_force(0.58, formulation=formulation)

    assert np.allclose(forces, 0.0, atol=1.0e-14)


@pytest.mark.parametrize(
    ("distance", "expected_force_on_left"),
    [
        (0.55, -2.0 * 0.03 / (1.0 - (0.03 / 0.1) ** 2) ** 2),
        (0.61, 2.0 * 0.03 / (1.0 - (0.03 / 0.1) ** 2) ** 2),
    ],
)
def test_legacy_spring_matches_existing_analytic_force(
    distance: float, expected_force_on_left: float
) -> None:
    forces = _spring_force(distance, formulation="legacy")

    assert forces[0, 0] == pytest.approx(expected_force_on_left)
    assert forces[1, 0] == pytest.approx(-expected_force_on_left)


@pytest.mark.parametrize("relative_extension", [-0.03, 0.03])
def test_fene_fraenkel_matches_relative_extension_analytic_force(
    relative_extension: float,
) -> None:
    rest_length = 0.58
    h_const = 2.0
    s_limit = 0.1
    distance = rest_length * (1.0 + relative_extension)
    expected = (
        h_const * relative_extension / (1.0 - (relative_extension / s_limit) ** 2)
    )

    forces = _spring_force(
        distance,
        rest_length=rest_length,
        h_const=h_const,
        s_limit=s_limit,
        formulation="fene_fraenkel",
    )

    assert forces[0, 0] == pytest.approx(expected)
    assert forces[1, 0] == pytest.approx(-expected)


def test_fene_fraenkel_limit_scales_with_each_bond_rest_length() -> None:
    relative_extension = 0.05
    short = _spring_force(
        0.4 * (1.0 + relative_extension),
        rest_length=0.4,
        formulation="fene_fraenkel",
    )
    long = _spring_force(
        0.8 * (1.0 + relative_extension),
        rest_length=0.8,
        formulation="fene_fraenkel",
    )

    assert short[0, 0] == pytest.approx(long[0, 0])


def _angle_triplet(theta: float) -> np.ndarray:
    return np.array(
        [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [math.cos(theta), math.sin(theta), 0.0]]
    )


@pytest.mark.parametrize("offset", [-0.15, 0.15])
def test_bending_force_restores_both_sides_of_equilibrium(offset: float) -> None:
    theta0 = 1.4
    displaced = _angle_triplet(theta0 + offset)
    equilibrium = _angle_triplet(theta0)

    forces = compute_bending_forces(
        positions_m=displaced,
        triplets=np.array([[0, 1, 2]], dtype=int),
        theta0_rad=np.array([theta0]),
        kb=3.0,
    )

    assert float(np.sum(forces * (equilibrium - displaced))) > 0.0
    assert np.allclose(forces.sum(axis=0), 0.0, atol=1.0e-12)


def _dihedral_quad(phi: float) -> np.ndarray:
    return np.array(
        [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, math.cos(phi), math.sin(phi)],
        ]
    )


@pytest.mark.parametrize("offset", [-0.15, 0.15])
def test_torsion_force_restores_both_sides_of_equilibrium(offset: float) -> None:
    phi0 = 0.7
    displaced = _dihedral_quad(phi0 + offset)
    equilibrium = _dihedral_quad(phi0)

    forces = compute_torsion_forces(
        positions_m=displaced,
        quads=np.array([[0, 1, 2, 3]], dtype=int),
        phi0_rad=np.array([phi0]),
        kt=2.0,
        fd_eps_m=1.0e-5,
    )

    assert float(np.sum(forces * (equilibrium - displaced))) > 0.0
    assert np.allclose(forces.sum(axis=0), 0.0, atol=1.0e-8)


def test_torsion_uses_shortest_error_across_pi_wrap() -> None:
    phi0 = -math.pi + 0.02
    displaced = _dihedral_quad(math.pi - 0.02)
    equilibrium = _dihedral_quad(phi0)

    forces = compute_torsion_forces(
        positions_m=displaced,
        quads=np.array([[0, 1, 2, 3]], dtype=int),
        phi0_rad=np.array([phi0]),
        kt=2.0,
        fd_eps_m=1.0e-5,
    )

    assert float(np.sum(forces * (equilibrium - displaced))) > 0.0


def test_hook_is_active_at_or_below_threshold_only() -> None:
    triplets = np.array([[0, 1, 2]], dtype=int)
    below = compute_hook_forces(_angle_triplet(math.radians(80.0)), triplets, 2.0, 90.0)
    boundary = compute_hook_forces(
        _angle_triplet(math.radians(90.0)), triplets, 2.0, 90.0
    )
    above = compute_hook_forces(
        _angle_triplet(math.radians(100.0)), triplets, 2.0, 90.0
    )

    assert np.linalg.norm(below) > 0.0
    assert np.allclose(boundary, 0.0, atol=1.0e-12)
    assert np.allclose(above, 0.0, atol=1.0e-12)


def test_segment_repulsion_balances_force_and_distributes_at_midpoints() -> None:
    positions = np.array(
        [
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.1, -1.0],
            [0.0, 0.1, 1.0],
        ]
    )
    forces = compute_segment_repulsion_forces(
        positions_m=positions,
        spring_pairs=np.array([[0, 1], [2, 3]], dtype=int),
        segment_pair_indices=np.array([[0, 1]], dtype=int),
        a_ss=2.0,
        cutoff=0.2,
        a_length=0.5,
    )
    expected_magnitude = 2.0 / 0.5 * math.exp(-0.1 / 0.5)

    assert np.allclose(forces.sum(axis=0), 0.0, atol=1.0e-12)
    assert forces[0, 1] == pytest.approx(-0.5 * expected_magnitude)
    assert forces[1] == pytest.approx(forces[0])
    assert forces[2, 1] == pytest.approx(0.5 * expected_magnitude)
    assert forces[3] == pytest.approx(forces[2])


def test_segment_repulsion_is_zero_at_cutoff() -> None:
    positions = np.array(
        [[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.2, -1.0], [0.0, 0.2, 1.0]]
    )
    forces = compute_segment_repulsion_forces(
        positions_m=positions,
        spring_pairs=np.array([[0, 1], [2, 3]], dtype=int),
        segment_pair_indices=np.array([[0, 1]], dtype=int),
        a_ss=2.0,
        cutoff=0.2,
        a_length=0.5,
    )

    assert np.allclose(forces, 0.0)


def test_repulsion_candidates_exclude_segments_with_shared_endpoint() -> None:
    spring_pairs = np.array([[0, 1], [1, 2], [3, 4]], dtype=int)

    candidates = _segment_pairs_without_neighbors(spring_pairs)

    assert candidates.tolist() == [[0, 2], [1, 2]]


def test_repulsion_excludes_body_body_and_keeps_other_segment_pairs() -> None:
    cfg = SimulationConfig.from_dict(_minimal_config_with_spring(None))
    engine = DynamicsEngine(ModelBuilder(cfg).build(), cfg)
    candidates = engine.segment_pair_indices_for_repulsion
    first_is_body = engine.body_spring_mask[candidates[:, 0]]
    second_is_body = engine.body_spring_mask[candidates[:, 1]]

    assert not np.any(first_is_body & second_is_body)
    assert np.any(first_is_body ^ second_is_body)
    assert np.any(~first_is_body & ~second_is_body)


def _minimal_config_with_spring(spring: dict[str, object] | None) -> dict[str, object]:
    config: dict[str, object] = {
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
        "motor": {"torque_Nm": 0.0, "reverse_n_flagella": 0},
        "time": {"duration_s": 0.1, "dt_s": 1.0e-3},
    }
    if spring is not None:
        config["potentials"] = {"spring": spring}
    return config


def test_missing_spring_formulation_falls_back_to_legacy() -> None:
    cfg = SimulationConfig.from_dict(
        _minimal_config_with_spring({"H_over_T_over_b": 10.0, "s": 0.1})
    )

    assert cfg.potentials.spring.formulation == "legacy"


def test_spring_formulation_accepts_fene_fraenkel() -> None:
    cfg = SimulationConfig.from_dict(
        _minimal_config_with_spring(
            {"formulation": "fene_fraenkel", "H_over_T_over_b": 10.0, "s": 0.1}
        )
    )

    assert cfg.potentials.spring.formulation == "fene_fraenkel"


@pytest.mark.parametrize(
    "spring",
    [
        {"formulation": "unknown", "H_over_T_over_b": 10.0, "s": 0.1},
        {"formulation": "legacy", "H_over_T_over_b": 0.0, "s": 0.1},
        {"formulation": "legacy", "H_over_T_over_b": 10.0, "s": 0.0},
        {"formulation": "fene_fraenkel", "H_over_T_over_b": 10.0, "s": 1.0},
    ],
)
def test_invalid_spring_parameters_are_rejected(spring: dict[str, object]) -> None:
    with pytest.raises(ValueError):
        SimulationConfig.from_dict(_minimal_config_with_spring(spring))
