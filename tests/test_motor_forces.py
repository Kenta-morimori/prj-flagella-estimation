from __future__ import annotations

import numpy as np

from sim_swim.dynamics.forces import compute_motor_forces


def test_motor_force_couple_matches_target_torque() -> None:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],  # b0
            [1.0, 0.0, 0.0],  # f1
            [1.0, 1.0, 0.0],  # f2
        ],
        dtype=float,
    )
    motor_triplets = np.array([[0, 1, 2]], dtype=int)
    torque_per_flag = np.array([2.0], dtype=float)

    forces, diag = compute_motor_forces(
        positions_m=positions,
        motor_triplets=motor_triplets,
        torque_per_flag=torque_per_flag,
    )

    assert np.allclose(forces.sum(axis=0), np.zeros(3), atol=1e-12)

    total_torque = np.sum(np.cross(positions, forces), axis=0)
    expected = np.array([0.0, 2.0, 0.0], dtype=float)
    assert np.allclose(total_torque, expected, atol=1e-8)
    assert diag.degenerate_axis_count == 0
    assert diag.bond_length_clipped_count == 0


def test_motor_force_skips_degenerate_axis() -> None:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],  # r_b == 0
        ],
        dtype=float,
    )
    motor_triplets = np.array([[0, 1, 2]], dtype=int)
    torque_per_flag = np.array([1.0], dtype=float)

    forces, diag = compute_motor_forces(
        positions_m=positions,
        motor_triplets=motor_triplets,
        torque_per_flag=torque_per_flag,
    )

    assert np.allclose(forces, np.zeros_like(forces))
    assert diag.degenerate_axis_count == 1
