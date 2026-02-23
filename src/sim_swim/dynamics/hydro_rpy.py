"""自由空間RPY mobility。"""

from __future__ import annotations

import math

import numpy as np


def compute_rpy_mobility(
    positions_m: np.ndarray,
    bead_radius_m: float,
    viscosity_Pa_s: float,
) -> np.ndarray:
    """RPY mobility tensorを構築する。

    Args:
        positions_m: ビーズ座標 [m]。shape=(N,3)
        bead_radius_m: ビーズ半径 a [m]
        viscosity_Pa_s: 粘性係数 [Pa*s]

    Returns:
        mobility行列 H。shape=(3N,3N), 単位は [m/(N*s)]。
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
            dist = float(np.linalg.norm(rij))
            if dist < 1e-15:
                block = self_block.copy()
            else:
                rhat = rij / dist
                rr = np.outer(rhat, rhat)
                if dist >= 2.0 * a:
                    c1 = 1.0 / (8.0 * math.pi * eta * dist)
                    block = c1 * (
                        (1.0 + (2.0 * a * a) / (3.0 * dist * dist)) * eye3
                        + (1.0 - (2.0 * a * a) / (dist * dist)) * rr
                    )
                else:
                    c2 = 1.0 / (6.0 * math.pi * eta * a)
                    block = c2 * (
                        (1.0 - (9.0 * dist) / (32.0 * a)) * eye3
                        + (3.0 * dist / (32.0 * a)) * rr
                    )

            i_slice = slice(3 * i, 3 * i + 3)
            j_slice = slice(3 * j, 3 * j + 3)
            h[i_slice, j_slice] = block
            h[j_slice, i_slice] = block.T

    return 0.5 * (h + h.T)
