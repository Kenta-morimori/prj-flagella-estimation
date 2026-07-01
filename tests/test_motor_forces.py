from __future__ import annotations

import math

import numpy as np
import pytest

from sim_swim.dynamics.forces import (
    compute_attach_first_body_axis_angle_forces,
    compute_attach_frame_target_forces,
    compute_motor_forces,
    compute_root_torque_axis_projection_forces,
    compute_root_torque_segment_couples_forces,
)


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


def test_attach_frame_target_forces_zero_at_targets() -> None:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    hook_triplets = np.array([[0, 1, 2]], dtype=int)
    forces = compute_attach_frame_target_forces(
        positions_m=positions,
        hook_triplets=hook_triplets,
        attach_first_target_vectors_m=np.array([[1.0, 0.0, 0.0]], dtype=float),
        first_second_target_vectors_m=np.array([[0.0, 1.0, 0.0]], dtype=float),
        attach_first_rest_lengths_m=np.array([1.0], dtype=float),
        first_second_rest_lengths_m=np.array([1.0], dtype=float),
        k_position=2.0,
        k_tangent=3.0,
    )

    assert np.allclose(forces, np.zeros_like(forces))


def test_attach_frame_target_forces_are_pair_balanced() -> None:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],
            [1.2, 1.3, 0.0],
        ],
        dtype=float,
    )
    hook_triplets = np.array([[0, 1, 2]], dtype=int)
    forces = compute_attach_frame_target_forces(
        positions_m=positions,
        hook_triplets=hook_triplets,
        attach_first_target_vectors_m=np.array([[1.0, 0.0, 0.0]], dtype=float),
        first_second_target_vectors_m=np.array([[0.0, 1.0, 0.0]], dtype=float),
        attach_first_rest_lengths_m=np.array([1.0], dtype=float),
        first_second_rest_lengths_m=np.array([1.0], dtype=float),
        k_position=2.0,
        k_tangent=3.0,
    )

    assert np.allclose(forces.sum(axis=0), np.zeros(3), atol=1e-12)
    assert forces[1, 0] < 0.0
    assert forces[2, 1] < 0.0


def test_motor_force_split_limits_short_basal_link_force() -> None:
    # Short basal link (attach-first) and longer first-second link.
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.25, 0.0, 0.0],
            [0.25, 0.58, 0.0],
        ],
        dtype=float,
    )
    motor_triplets = np.array([[0, 1, 2]], dtype=int)
    torque_per_flag = np.array([5.0e-20], dtype=float)

    forces, _ = compute_motor_forces(
        positions_m=positions,
        motor_triplets=motor_triplets,
        torque_per_flag=torque_per_flag,
    )

    total_torque = np.sum(np.cross(positions, forces), axis=0)
    expected = np.array([0.0, 5.0e-20, 0.0], dtype=float)
    assert np.allclose(total_torque, expected, atol=1e-28)

    # Phase E objective: avoid unrealistically large basal-link force from short ra.
    attach_force_norm = float(np.linalg.norm(forces[0]))
    expected_scale = float(torque_per_flag[0] / 0.25)
    # Keep the attach force at the same order as T/ra (not >10x explosion).
    assert attach_force_norm <= 10.0 * expected_scale


def test_root_torque_axis_projection_applies_balanced_torque() -> None:
    body = np.array(
        [
            [0.0, -1.0, -1.0],
            [0.0, 1.0, -1.0],
            [0.0, 1.0, 1.0],
            [0.0, -1.0, 1.0],
        ],
        dtype=float,
    )
    theta = np.linspace(0.0, 2.0 * math.pi, 8, endpoint=False)
    flag = np.column_stack(
        [
            np.linspace(0.0, 3.0, theta.size),
            0.25 * np.cos(theta),
            0.25 * np.sin(theta),
        ]
    )
    positions = np.vstack([body, flag])
    body_indices = np.arange(body.shape[0], dtype=int)
    flag_indices = [np.arange(body.shape[0], positions.shape[0], dtype=int)]
    torque = 2.0e-20

    forces, diag = compute_root_torque_axis_projection_forces(
        positions_m=positions,
        flagella_indices=flag_indices,
        body_indices=body_indices,
        torque_per_flag=np.array([torque], dtype=float),
    )

    assert diag.degenerate_axis_count == 0
    assert np.allclose(forces.sum(axis=0), np.zeros(3), atol=1e-30)

    flag_forces = forces[flag_indices[0]]
    body_forces = forces[body_indices]
    flag_torque = np.sum(np.cross(flag, flag_forces), axis=0)
    body_torque = np.sum(np.cross(body, body_forces), axis=0)
    centered = flag - np.mean(flag, axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    axis = vh[0]
    if float(np.dot(axis, flag[-1] - flag[0])) < 0.0:
        axis = -axis
    assert np.isclose(abs(float(flag_torque @ axis)), torque, rtol=1.0e-8)
    assert np.isclose(abs(float(body_torque @ axis)), torque, rtol=1.0e-8)
    assert np.isclose(float((flag_torque + body_torque) @ axis), 0.0, atol=1e-28)


def test_root_torque_segment_couples_applies_local_balanced_torque() -> None:
    body = np.array(
        [
            [0.0, -1.0, -1.0],
            [0.0, 1.0, -1.0],
            [0.0, 1.0, 1.0],
            [0.0, -1.0, 1.0],
        ],
        dtype=float,
    )
    theta = np.linspace(0.0, 2.0 * math.pi, 8, endpoint=False)
    flag = np.column_stack(
        [
            np.linspace(0.0, 3.0, theta.size),
            0.25 * np.cos(theta),
            0.25 * np.sin(theta),
        ]
    )
    positions = np.vstack([body, flag])
    body_indices = np.arange(body.shape[0], dtype=int)
    flag_indices = [np.arange(body.shape[0], positions.shape[0], dtype=int)]
    torque = 2.0e-20
    segment_weights = [np.linspace(1.0, 0.2, flag.shape[0] - 1)]

    forces, diag = compute_root_torque_segment_couples_forces(
        positions_m=positions,
        flagella_indices=flag_indices,
        body_indices=body_indices,
        torque_per_flag=np.array([torque], dtype=float),
        segment_weights=segment_weights,
    )

    assert diag.degenerate_axis_count == 0
    assert np.allclose(forces.sum(axis=0), np.zeros(3), atol=1e-30)

    flag_forces = forces[flag_indices[0]]
    body_forces = forces[body_indices]
    origin = flag[0]
    centered = flag - np.mean(flag, axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    axis = vh[0]
    if float(np.dot(axis, flag[-1] - flag[0])) < 0.0:
        axis = -axis
    flag_torque = np.sum(np.cross(flag - origin, flag_forces), axis=0)
    body_torque = np.sum(np.cross(body - origin, body_forces), axis=0)

    assert np.isclose(abs(float(flag_torque @ axis)), torque, rtol=1.0e-8)
    assert np.isclose(abs(float(body_torque @ axis)), torque, rtol=1.0e-8)
    assert np.isclose(float((flag_torque + body_torque) @ axis), 0.0, atol=1e-28)


def test_root_torque_segment_couples_reacts_to_applied_flag_torque() -> None:
    body = np.array(
        [
            [0.0, -1.0, -1.0],
            [0.0, 1.0, -1.0],
            [0.0, 1.0, 1.0],
            [0.0, -1.0, 1.0],
        ],
        dtype=float,
    )
    flag = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.25, 0.0],
            [3.0, 0.25, 0.0],
            [4.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    positions = np.vstack([body, flag])
    body_indices = np.arange(body.shape[0], dtype=int)
    flag_indices = [np.arange(body.shape[0], positions.shape[0], dtype=int)]
    torque = 2.0e-20
    segment_weights = [np.ones(flag.shape[0] - 1, dtype=float)]

    forces, diag = compute_root_torque_segment_couples_forces(
        positions_m=positions,
        flagella_indices=flag_indices,
        body_indices=body_indices,
        torque_per_flag=np.array([torque], dtype=float),
        segment_weights=segment_weights,
    )

    assert diag.degenerate_axis_count > 0
    assert np.allclose(forces.sum(axis=0), np.zeros(3), atol=1e-30)

    flag_forces = forces[flag_indices[0]]
    body_forces = forces[body_indices]
    origin = flag[0]
    centered = flag - np.mean(flag, axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    axis = vh[0]
    if float(np.dot(axis, flag[-1] - flag[0])) < 0.0:
        axis = -axis
    flag_torque = float(np.sum(np.cross(flag - origin, flag_forces), axis=0) @ axis)
    body_torque = float(np.sum(np.cross(body - origin, body_forces), axis=0) @ axis)

    assert abs(flag_torque) < torque
    assert np.isclose(flag_torque + body_torque, 0.0, atol=1e-28)


def test_root_torque_axis_projection_skips_nonfinite_flag_positions() -> None:
    body = np.array(
        [
            [0.0, -1.0, -1.0],
            [0.0, 1.0, -1.0],
            [0.0, 1.0, 1.0],
            [0.0, -1.0, 1.0],
        ],
        dtype=float,
    )
    flag = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [float("nan"), 0.25, 0.0],
            [3.0, 0.25, 0.0],
            [4.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    positions = np.vstack([body, flag])
    body_indices = np.arange(body.shape[0], dtype=int)
    flag_indices = [np.arange(body.shape[0], positions.shape[0], dtype=int)]

    forces, diag = compute_root_torque_axis_projection_forces(
        positions_m=positions,
        flagella_indices=flag_indices,
        body_indices=body_indices,
        torque_per_flag=np.array([2.0e-20], dtype=float),
    )

    assert diag.degenerate_axis_count == 1
    assert np.allclose(forces, np.zeros_like(positions), equal_nan=False)


def test_root_torque_segment_couples_skips_when_axis_svd_fails(monkeypatch) -> None:
    body = np.array(
        [
            [0.0, -1.0, -1.0],
            [0.0, 1.0, -1.0],
            [0.0, 1.0, 1.0],
            [0.0, -1.0, 1.0],
        ],
        dtype=float,
    )
    flag = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.25, 0.0],
            [3.0, 0.25, 0.0],
            [4.0, 0.0, 0.0],
        ],
        dtype=float,
    )
    positions = np.vstack([body, flag])
    body_indices = np.arange(body.shape[0], dtype=int)
    flag_indices = [np.arange(body.shape[0], positions.shape[0], dtype=int)]

    def raise_linalg_error(*_args, **_kwargs):
        raise np.linalg.LinAlgError("SVD did not converge")

    monkeypatch.setattr(np.linalg, "svd", raise_linalg_error)

    forces, diag = compute_root_torque_segment_couples_forces(
        positions_m=positions,
        flagella_indices=flag_indices,
        body_indices=body_indices,
        torque_per_flag=np.array([2.0e-20], dtype=float),
        segment_weights=[np.ones(flag.shape[0] - 1, dtype=float)],
    )

    assert diag.degenerate_axis_count == 1
    assert np.allclose(forces, np.zeros_like(positions))


def test_attach_first_body_axis_angle_force_penalizes_body_axis_component() -> None:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],  # body attach
            [0.1, 0.25, 0.0],  # first bead: has body-axis component
            [0.1, 0.25, 0.58],  # second bead: not used by this force
        ],
        dtype=float,
    )
    forces = compute_attach_first_body_axis_angle_forces(
        positions_m=positions,
        hook_triplets=np.array([[0, 1, 2]], dtype=int),
        attach_first_rest_lengths_m=np.array([0.25], dtype=float),
        body_axis_unit=np.array([1.0, 0.0, 0.0], dtype=float),
        k_angle=2.0,
    )

    assert np.allclose(forces.sum(axis=0), np.zeros(3), atol=1e-12)
    assert forces[1, 0] < 0.0
    assert forces[0, 0] > 0.0
    assert forces[1, 1] == pytest.approx(0.0)
    assert forces[1, 2] == pytest.approx(0.0)


def test_attach_first_body_axis_angle_force_is_zero_when_perpendicular() -> None:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.25, 0.0],
            [0.0, 0.25, 0.58],
        ],
        dtype=float,
    )
    forces = compute_attach_first_body_axis_angle_forces(
        positions_m=positions,
        hook_triplets=np.array([[0, 1, 2]], dtype=int),
        attach_first_rest_lengths_m=np.array([0.25], dtype=float),
        body_axis_unit=np.array([1.0, 0.0, 0.0], dtype=float),
        k_angle=2.0,
    )

    assert np.allclose(forces, np.zeros_like(forces), atol=1e-12)
