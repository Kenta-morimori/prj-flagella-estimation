"""自由空間RPY mobility。"""

from __future__ import annotations

import math

import numpy as np


def compute_rpy_pair_mobility(
    displacement_m: np.ndarray,
    bead_radius_m: float,
    viscosity_Pa_s: float,
) -> np.ndarray:
    """Return the RPY mobility block from a source bead to an observation point.

    ``displacement_m`` is the observation position minus the source position.
    The finite-size branch also makes this well defined at a bead centre.
    """

    eye3 = np.eye(3, dtype=float)
    a = max(float(bead_radius_m), 1e-12)
    eta = max(float(viscosity_Pa_s), 1e-12)
    self_block = (1.0 / (6.0 * math.pi * eta * a)) * eye3
    displacement = np.asarray(displacement_m, dtype=float)
    dist = float(np.linalg.norm(displacement))
    if dist < 1e-15:
        return self_block
    rhat = displacement / dist
    rr = np.outer(rhat, rhat)
    if dist >= 2.0 * a:
        c1 = 1.0 / (8.0 * math.pi * eta * dist)
        return c1 * (
            (1.0 + (2.0 * a * a) / (3.0 * dist * dist)) * eye3
            + (1.0 - (2.0 * a * a) / (dist * dist)) * rr
        )
    c2 = 1.0 / (6.0 * math.pi * eta * a)
    return c2 * (
        (1.0 - (9.0 * dist) / (32.0 * a)) * eye3 + (3.0 * dist / (32.0 * a)) * rr
    )


def compute_rpy_mobility(
    positions_m: np.ndarray,
    bead_radius_m: float,
    viscosity_Pa_s: float,
) -> np.ndarray:
    """RPY mobility tensorを構築する。

    Args:
        positions_m: ビーズ座標 [m]。shape=(N,3)
        bead_radius_m: ビーズ半径 [m]
        viscosity_Pa_s: 粘性係数 [Pa*s]

    Returns:
        mobility行列 H。shape=(3N,3N)
    """

    n = int(positions_m.shape[0])
    h = np.zeros((3 * n, 3 * n), dtype=float)
    eye3 = np.eye(3, dtype=float)

    a = max(bead_radius_m, 1e-12)
    eta = max(viscosity_Pa_s, 1e-12)
    self_block = (1.0 / (6.0 * math.pi * eta * a)) * eye3

    for i in range(n):
        h[3 * i : 3 * i + 3, 3 * i : 3 * i + 3] = self_block

    for i in range(n):
        ri = positions_m[i]
        for j in range(i + 1, n):
            rj = positions_m[j]
            rij = ri - rj
            block = compute_rpy_pair_mobility(rij, a, eta)

            i_slice = slice(3 * i, 3 * i + 3)
            j_slice = slice(3 * j, 3 * j + 3)
            h[i_slice, j_slice] = block
            h[j_slice, i_slice] = block.T

    return 0.5 * (h + h.T)
