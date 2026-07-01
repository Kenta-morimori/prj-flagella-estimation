from __future__ import annotations

import numpy as np

from sim_swim.sim.helix_axis import (
    angle_deg_between,
    estimate_flag_helix_axis,
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
