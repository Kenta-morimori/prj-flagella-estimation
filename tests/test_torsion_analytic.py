"""Physical checks for the opt-in gradient; historical FD remains the default."""

import numpy as np
import pytest

from sim_swim.dynamics.forces import (
    _torsion_energy_quad,
    compute_torsion_forces,
)
from sim_swim.sim.params import SimulationConfig


def geometry():
    points = (
        np.array(
            [
                [0.0, 0.2, 0.1],
                [0.5, 0.0, 0.0],
                [1.0, 0.1, 0.0],
                [1.5, 0.3, 0.4],
                [2.0, -0.1, 0.6],
            ]
        )
        * 1e-6
    )
    return points, np.array([[0, 1, 2, 3], [1, 2, 3, 4]]), np.array([-0.7, 0.2])


def test_torsion_gradient_conservation_covariance():
    points, quads, rest = geometry()
    kt = 2.5e-19
    force = compute_torsion_forces(points, quads, rest, kt, 1e-9, "analytic")
    numerical = np.zeros_like(points)

    def energy(p):
        return sum(_torsion_energy_quad(p[q], phi, kt) for q, phi in zip(quads, rest))

    eps = 1e-12
    for i in range(len(points)):
        for axis in range(3):
            high, low = points.copy(), points.copy()
            high[i, axis] += eps
            low[i, axis] -= eps
            numerical[i, axis] = -(energy(high) - energy(low)) / (2 * eps)
    np.testing.assert_allclose(
        force, numerical, rtol=2e-7, atol=np.linalg.norm(force) * 1e-9
    )
    assert np.linalg.norm(force.sum(0)) < np.linalg.norm(force) * 1e-14
    assert (
        np.linalg.norm(np.cross(points - points.mean(0), force).sum(0))
        < np.linalg.norm(force) * 1e-6 * 1e-13
    )
    rotation, _ = np.linalg.qr(np.random.default_rng(12).normal(size=(3, 3)))
    if np.linalg.det(rotation) < 0:
        rotation[:, 0] *= -1
    moved = points @ rotation.T + np.array([3.0, -2.0, 4.0]) * 1e-6
    rotated = compute_torsion_forces(moved, quads, rest, kt, 1e-9, "analytic")
    np.testing.assert_allclose(
        rotated, force @ rotation.T, rtol=1e-12, atol=np.linalg.norm(force) * 1e-13
    )
    assert energy(moved) == pytest.approx(energy(points), rel=1e-12)


def test_fd_converges_to_analytic_and_legacy_default_is_exact():
    p, q, r = geometry()
    analytic = compute_torsion_forces(p, q, r, 2.5e-19, 1e-9, "analytic")
    errors = []
    for eps in (1e-8, 1e-9, 1e-10):
        old = compute_torsion_forces(p, q, r, 2.5e-19, eps)
        np.testing.assert_array_equal(
            old, compute_torsion_forces(p, q, r, 2.5e-19, eps, "finite_difference")
        )
        errors.append(np.linalg.norm(old - analytic))
    assert errors[2] < errors[1] < errors[0]
    assert errors[2] < np.linalg.norm(analytic) * 1e-5


def test_analytic_degenerate_and_cusp_fail_explicitly():
    p = (
        np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
        * 1e-6
    )
    q = np.array([[0, 1, 2, 3]])
    with pytest.raises(RuntimeError, match="cusp"):
        compute_torsion_forces(p, q, np.array([np.pi]), 1.0, 1e-9, "analytic")
    p[0] = p[1]
    with pytest.raises(RuntimeError, match="degenerate"):
        compute_torsion_forces(p, q, np.array([0.0]), 1.0, 1e-9, "analytic")


def test_torsion_method_config_validation():
    assert (
        SimulationConfig.from_dict(
            {"time": {"dt_s": 0.001}}
        ).potentials.torsion.force_method
        == "finite_difference"
    )
    with pytest.raises(ValueError, match="force_method"):
        SimulationConfig.from_dict(
            {"potentials": {"torsion": {"force_method": "fallback"}}}
        )
