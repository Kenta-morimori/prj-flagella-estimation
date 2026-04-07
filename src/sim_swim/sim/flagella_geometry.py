"""レンダリング用トポロジ定義。"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FlagellaRig:
    """描画に必要なビーズ接続情報。"""

    body_layer_indices: list[np.ndarray]
    body_ring_edges: np.ndarray
    body_vertical_edges: np.ndarray
    body_spring_edges: np.ndarray
    flagella_indices: list[np.ndarray]
    hook_triplets: np.ndarray
