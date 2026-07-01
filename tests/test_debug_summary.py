from __future__ import annotations

import numpy as np

from sim_swim.sim.debug_summary import (
    _estimate_helix_spin_offset_deg,
    _flag_helix_bundle_radius_stats_um,
)


def test_helix_spin_offset_returns_nan_for_nonfinite_points() -> None:
    phase_deg, radius_um = _estimate_helix_spin_offset_deg(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.1, 0.0],
                [2.0, 0.0, 0.1],
                [3.0, -0.1, 0.0],
                [float("nan"), 0.0, -0.1],
            ],
            dtype=float,
        )
    )

    assert np.isnan(phase_deg)
    assert np.isnan(radius_um)


def test_bundle_radius_returns_nan_when_svd_fails(monkeypatch) -> None:
    positions_m = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.1, 0.0],
            [3.0, 0.0, 0.1],
        ],
        dtype=float,
    )
    flagella_indices = [np.array([0, 1, 2, 3], dtype=int)]

    def raise_linalg_error(*_args, **_kwargs):
        raise np.linalg.LinAlgError("SVD did not converge")

    monkeypatch.setattr(np.linalg, "svd", raise_linalg_error)

    mean_um, max_um = _flag_helix_bundle_radius_stats_um(
        positions_m,
        flagella_indices,
    )

    assert np.isnan(mean_um)
    assert np.isnan(max_um)
