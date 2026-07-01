from __future__ import annotations

import numpy as np

from sim_swim.sim.helix_axis import (
    HelixAxisEstimate,
    angle_deg_between,
    estimate_flag_helix_axis,
    helix_axis_centered_metrics,
    helix_axis_alignment_metrics,
)


def test_helix_axis_uses_beads_after_hook_bead() -> None:
    positions = np.array(
        [
            [0.0, 10.0, 0.0],
            [0.0, 0.2, 0.0],
            [1.0, -0.2, 0.1],
            [2.0, 0.2, -0.1],
            [3.0, -0.2, 0.0],
        ],
        dtype=float,
    )
    flag_indices = np.array([0, 1, 2, 3, 4], dtype=int)

    estimate = estimate_flag_helix_axis(positions, flag_indices, flag_id=0)

    assert not estimate.degenerate
    assert angle_deg_between(estimate.axis, np.array([1.0, 0.0, 0.0])) < 10.0
    assert abs(float(estimate.origin[1])) < 0.1


def test_helix_axis_reports_degenerate_when_only_one_helix_bead_exists() -> None:
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=float)

    estimate = estimate_flag_helix_axis(positions, np.array([0, 1]), flag_id=0)

    assert estimate.degenerate
    assert np.isnan(estimate.fit_r2)
    assert np.isnan(estimate.axis).all()


def test_helix_axis_centered_metrics_detects_clean_helix() -> None:
    theta = np.linspace(0.0, 6.0 * np.pi, 36, endpoint=False)
    helix = np.column_stack(
        [
            np.linspace(0.0, 20.0, theta.size),
            0.25 * np.cos(theta),
            0.25 * np.sin(theta),
        ]
    )
    positions = np.vstack([[0.0, 0.25, 0.0], helix])
    flag_indices = np.arange(positions.shape[0], dtype=int)
    estimate = HelixAxisEstimate(
        flag_id=0,
        origin=np.zeros(3, dtype=float),
        axis=np.array([1.0, 0.0, 0.0], dtype=float),
        line_start=np.array([0.0, 0.0, 0.0], dtype=float),
        line_end=np.array([20.0, 0.0, 0.0], dtype=float),
        fit_r2=1.0,
        degenerate=False,
    )

    metrics = helix_axis_centered_metrics(
        positions,
        flag_indices,
        estimate,
        reference_direction=np.array([0.0, 1.0, 0.0]),
    )

    assert not metrics.degenerate
    assert metrics.fit_r2 > 0.95
    assert metrics.radius_cv < 0.1
    assert metrics.radius_mean_m > 0.0


def test_helix_axis_reports_degenerate_for_nonfinite_positions() -> None:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [float("nan"), 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ],
        dtype=float,
    )

    estimate = estimate_flag_helix_axis(positions, np.array([0, 1, 2, 3]), flag_id=2)

    assert estimate.flag_id == 2
    assert estimate.degenerate
    assert np.isnan(estimate.fit_r2)
    assert np.isnan(estimate.axis).all()


def test_helix_axis_centered_metrics_reports_degenerate_for_nonfinite() -> None:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.1, 0.0],
            [2.0, float("nan"), 0.0],
            [3.0, 0.0, 0.1],
        ],
        dtype=float,
    )
    flag_indices = np.array([0, 1, 2, 3], dtype=int)
    estimate = estimate_flag_helix_axis(positions, flag_indices, flag_id=0)

    metrics = helix_axis_centered_metrics(
        positions,
        flag_indices,
        estimate,
        reference_direction=np.array([0.0, 1.0, 0.0]),
    )

    assert metrics.degenerate
    assert np.isnan(metrics.fit_r2)
    assert np.isnan(metrics.radius_cv)


def test_helix_axis_reports_degenerate_when_svd_fails(monkeypatch) -> None:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ],
        dtype=float,
    )

    def raise_linalg_error(*_args, **_kwargs):
        raise np.linalg.LinAlgError("SVD did not converge")

    monkeypatch.setattr(np.linalg, "svd", raise_linalg_error)

    estimate = estimate_flag_helix_axis(positions, np.array([0, 1, 2, 3]), flag_id=3)

    assert estimate.flag_id == 3
    assert estimate.degenerate
    assert np.isnan(estimate.fit_r2)
    assert np.isnan(estimate.axis).all()


def test_helix_axis_alignment_metrics_detects_aligned_axes() -> None:
    metrics = helix_axis_alignment_metrics(
        [
            np.array([1.0, 0.0, 0.0]),
            np.array([0.99, 0.05, 0.0]),
            np.array([0.98, -0.05, 0.0]),
        ]
    )

    assert metrics.pair_angle_deg_max < 10.0
    assert metrics.mean_deviation_deg_max < 5.0
    assert metrics.alignment_order > 0.99


def test_helix_axis_alignment_metrics_detects_one_outlier() -> None:
    metrics = helix_axis_alignment_metrics(
        [
            np.array([1.0, 0.0, 0.0]),
            np.array([0.99, 0.05, 0.0]),
            np.array([0.0, 1.0, 0.0]),
        ]
    )

    assert metrics.pair_angle_deg_max > 80.0
    assert metrics.mean_deviation_deg_max > 40.0
    assert metrics.alignment_order < 0.9
