"""Closed-segment distances must preserve geometry at simulation length scales."""

from __future__ import annotations

import math

import numpy as np
import pytest

from sim_swim.dynamics.forces import (
    _closest_points_on_segments,
    compute_segment_repulsion_forces,
)

pytestmark = pytest.mark.light


@pytest.mark.parametrize("scale", [1e-9, 1e-6, 1.0, 1e6])
@pytest.mark.parametrize(
    ("points", "expected_a", "expected_b"),
    [
        # Interior/interior closest points, including the observed SI regression.
        ([[-1, 0, 0], [1, 0, 0], [0, -1, 0.1], [0, 1, 0.1]], [0, 0, 0], [0, 0, 0.1]),
        # Infinite-line solution lies outside the second segment; both parameters
        # cannot simply be clipped independently.
        ([[0, 0, 0], [1, 0, 0], [1, 1, 0], [2, 2, 0]], [1, 0, 0], [1, 1, 0]),
        # Parallel, disjoint segments.
        ([[0, 0, 0], [1, 0, 0], [2, 0.1, 0], [3, 0.1, 0]], [1, 0, 0], [2, 0.1, 0]),
        # Point/segment and point/point degeneracy.
        ([[0, 0, 0], [0, 0, 0], [-1, 0.1, 0], [1, 0.1, 0]], [0, 0, 0], [0, 0.1, 0]),
        ([[0, 0, 0], [0, 0, 0], [1, 2, 0], [1, 2, 0]], [0, 0, 0], [1, 2, 0]),
    ],
)
def test_closed_segment_closest_points_across_length_scales(
    scale: float,
    points: list[list[float]],
    expected_a: list[float],
    expected_b: list[float],
) -> None:
    p = np.array(points, dtype=float) * scale
    a, b, s, t, distance = _closest_points_on_segments(*p)
    np.testing.assert_allclose(a / scale, expected_a, atol=1e-12)
    np.testing.assert_allclose(b / scale, expected_b, atol=1e-12)
    assert 0 <= s <= 1 and 0 <= t <= 1
    assert distance / scale == pytest.approx(
        np.linalg.norm(np.array(expected_a) - expected_b), abs=1e-12
    )


@pytest.mark.parametrize("scale", [1e-9, 1e-6, 1.0])
def test_nearly_parallel_interior_closest_points(scale: float) -> None:
    points = np.array([[-1, 0, 0], [1, 0, 0], [-1, -1e-8, 0.1], [1, 1e-8, 0.1]]) * scale
    a, b, s, t, distance = _closest_points_on_segments(*points)
    np.testing.assert_allclose(a / scale, [0, 0, 0], atol=1e-10)
    np.testing.assert_allclose(b / scale, [0, 0, 0.1], atol=1e-10)
    assert s == pytest.approx(0.5)
    assert t == pytest.approx(0.5)
    assert distance / scale == pytest.approx(0.1)


def test_segment_distance_is_covariant_under_rotation_translation_and_order() -> None:
    p = np.array([[-1, 0, 0], [1, 0, 0], [0, -1, 0.1], [0, 1, 0.1]]) * 1e-6
    rotation, _ = np.linalg.qr(np.random.default_rng(263).normal(size=(3, 3)))
    shift = np.array([2, -3, 5]) * 1e-6
    moved = p @ rotation.T + shift
    a, b, _, _, distance = _closest_points_on_segments(*moved)
    np.testing.assert_allclose(a, shift, atol=1e-18)
    np.testing.assert_allclose(
        b, np.array([0, 0, 0.1e-6]) @ rotation.T + shift, atol=1e-18
    )
    for order in [[1, 0, 2, 3], [0, 1, 3, 2], [2, 3, 0, 1]]:
        aa, bb, _, _, dd = _closest_points_on_segments(*moved[order])
        assert dd == pytest.approx(distance, abs=1e-18)
        expected = (b, a) if order[0] == 2 else (a, b)
        np.testing.assert_allclose(aa, expected[0], atol=1e-18)
        np.testing.assert_allclose(bb, expected[1], atol=1e-18)


@pytest.mark.parametrize("offset", [0.0, 0.1])
def test_parallel_overlapping_segments_have_correct_distance(offset: float) -> None:
    p = np.array([[0, 0, 0], [2, 0, 0], [1, offset, 0], [3, offset, 0]]) * 1e-6
    a, b, _, _, distance = _closest_points_on_segments(*p)
    assert distance == pytest.approx(offset * 1e-6, abs=1e-18)
    assert 1e-6 <= a[0] <= 2e-6
    assert 1e-6 <= b[0] <= 2e-6


def test_repulsion_detects_si_interior_contact_and_balances_force_and_torque() -> None:
    positions = np.array([[-1, 0, 0], [1, 0, 0], [0, -1, 0.1], [0, 1, 0.1]]) * 1e-6
    amplitude, length = 2.5e-20, 0.2e-6
    force = compute_segment_repulsion_forces(
        positions_m=positions,
        spring_pairs=np.array([[0, 1], [2, 3]]),
        segment_pair_indices=np.array([[0, 1]]),
        a_ss=amplitude,
        cutoff=0.2e-6,
        a_length=length,
    )
    expected = amplitude / length * math.exp(-0.1e-6 / length)
    np.testing.assert_allclose(force[:2, 2], -expected / 2, rtol=1e-12, atol=1e-28)
    np.testing.assert_allclose(force[2:, 2], expected / 2, rtol=1e-12, atol=1e-28)
    np.testing.assert_allclose(force[:, :2], 0, atol=1e-28)
    np.testing.assert_allclose(force.sum(0), 0, atol=1e-28)
    np.testing.assert_allclose(np.cross(positions, force).sum(0), 0, atol=1e-33)
