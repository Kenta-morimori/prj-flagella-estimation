from __future__ import annotations

import numpy as np

from sim_swim.dynamics.hydro_rpy import compute_rpy_mobility


def test_rpy_is_symmetric_and_almost_psd() -> None:
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.2e-6, 0.0, 0.0],
            [0.0, 1.1e-6, 0.0],
            [0.0, 0.0, 1.3e-6],
        ],
        dtype=float,
    )
    h = compute_rpy_mobility(
        positions_m=positions,
        bead_radius_m=0.15e-6,
        viscosity_Pa_s=1.0e-3,
    )

    assert np.allclose(h, h.T, atol=1e-12)

    evals = np.linalg.eigvalsh(0.5 * (h + h.T))
    assert evals.min() > -1e-12
