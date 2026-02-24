from __future__ import annotations

import numpy as np

from sim_swim.dynamics.brownian import K_B, sample_brownian_displacement
from sim_swim.dynamics.hydro_rpy import compute_rpy_mobility


def test_brownian_covariance_matches_target_scale() -> None:
    rng = np.random.default_rng(123)
    positions = np.array([[0.0, 0.0, 0.0], [1.0e-6, 0.0, 0.0]], dtype=float)
    h = compute_rpy_mobility(
        positions_m=positions, bead_radius_m=0.12e-6, viscosity_Pa_s=1e-3
    )

    dt = 2.0e-3
    temp = 300.0
    target = 2.0 * K_B * temp * h * dt

    n_sample = 3000
    samples = np.array(
        [
            sample_brownian_displacement(
                mobility=h,
                dt=dt,
                temperature_K=temp,
                rng=rng,
                method="eigh",
                jitter=1e-20,
            )
            for _ in range(n_sample)
        ],
        dtype=float,
    )

    cov = np.cov(samples, rowvar=False, bias=True)
    diag_target = np.diag(target)
    diag_cov = np.diag(cov)

    rel = np.abs(diag_cov - diag_target) / np.maximum(np.abs(diag_target), 1e-30)
    assert np.percentile(rel, 75) < 0.35
