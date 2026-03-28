from __future__ import annotations

import numpy as np

from sim_swim.dynamics.forces import compute_bead_steric_exclusion_forces


def test_wca_steric_force_nonzero_when_within_cutoff() -> None:
    sigma = 1.0e-7
    cutoff = (2.0 ** (1.0 / 6.0)) * sigma
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.9 * sigma, 0.0, 0.0],
        ],
        dtype=float,
    )
    pairs = np.array([[0, 1]], dtype=int)

    forces = compute_bead_steric_exclusion_forces(
        positions_m=positions,
        bead_pair_indices=pairs,
        epsilon=1.0e-18,
        sigma=sigma,
        cutoff=cutoff,
    )

    assert float(np.linalg.norm(forces[0])) > 0.0
    assert float(np.linalg.norm(forces[1])) > 0.0


def test_wca_steric_force_zero_when_outside_cutoff() -> None:
    sigma = 1.0e-7
    cutoff = (2.0 ** (1.0 / 6.0)) * sigma
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.5 * cutoff, 0.0, 0.0],
        ],
        dtype=float,
    )
    pairs = np.array([[0, 1]], dtype=int)

    forces = compute_bead_steric_exclusion_forces(
        positions_m=positions,
        bead_pair_indices=pairs,
        epsilon=1.0e-18,
        sigma=sigma,
        cutoff=cutoff,
    )

    assert np.allclose(forces, np.zeros_like(forces))


def test_wca_steric_force_is_action_reaction_symmetric() -> None:
    sigma = 1.0e-7
    cutoff = (2.0 ** (1.0 / 6.0)) * sigma
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.95 * sigma, 0.0, 0.0],
        ],
        dtype=float,
    )
    pairs = np.array([[0, 1]], dtype=int)

    forces = compute_bead_steric_exclusion_forces(
        positions_m=positions,
        bead_pair_indices=pairs,
        epsilon=2.0e-18,
        sigma=sigma,
        cutoff=cutoff,
    )

    assert np.allclose(forces[0] + forces[1], np.zeros(3), atol=1e-18)
