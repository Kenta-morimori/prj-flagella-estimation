"""Brownian変位の生成。"""

from __future__ import annotations

import numpy as np

K_B = 1.380649e-23


def sample_brownian_displacement(
    mobility: np.ndarray,
    dt: float,
    temperature_K: float,
    rng: np.random.Generator,
    method: str = "cholesky",
    jitter: float = 1.0e-20,
) -> np.ndarray:
    """共分散 `2 kB T H dt` を満たすBrownian変位を生成する。

    Args:
        mobility: RPY mobility行列 H。
        dt: 時間刻み [s]。
        temperature_K: 温度 [K]。
        rng: 乱数生成器。
        method: `cholesky` または `eigh`。
        jitter: 数値安定化用の対角jitter。

    Returns:
        変位ベクトル xi（shape=(3N,)）。
    """

    sym_h = 0.5 * (mobility + mobility.T)
    cov = 2.0 * K_B * max(temperature_K, 0.0) * sym_h * max(dt, 0.0)
    n = cov.shape[0]
    z = rng.normal(0.0, 1.0, n)

    if method == "cholesky":
        eye = np.eye(n, dtype=float)
        last_exc: Exception | None = None
        for scale in (1.0, 10.0, 100.0, 1000.0):
            try:
                lmat = np.linalg.cholesky(cov + eye * max(jitter, 0.0) * scale)
                return lmat @ z
            except np.linalg.LinAlgError as exc:
                last_exc = exc
        raise RuntimeError("Brownian cholesky 分解に失敗しました。") from last_exc

    if method == "eigh":
        evals, evecs = np.linalg.eigh(cov)
        evals = np.clip(evals, 0.0, None)
        sqrt_cov = (evecs * np.sqrt(evals)) @ evecs.T
        return sqrt_cov @ z

    raise ValueError(f"Unknown brownian.method: {method}")
