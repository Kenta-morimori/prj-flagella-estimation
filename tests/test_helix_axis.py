from __future__ import annotations

import numpy as np

from sim_swim.sim.helix_axis import (
    angle_deg_between,
    estimate_flag_helix_axis,
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
