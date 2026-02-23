"""bead-spring モデルのデータ構造。"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np


class PolymorphState(IntEnum):
    """べん毛多形状態。"""

    NORMAL = 0
    SEMICOILED = 1
    CURLY1 = 2


@dataclass
class SimModel:
    """シミュレーションに必要な幾何・トポロジを保持する。"""

    positions_m: np.ndarray

    body_indices: np.ndarray
    body_layer_indices: list[np.ndarray]
    body_ring_edges: np.ndarray
    body_vertical_edges: np.ndarray

    flagella_indices: list[np.ndarray]

    spring_pairs: np.ndarray
    spring_rest_lengths_m: np.ndarray
    bending_triplets: np.ndarray
    bending_flag_ids: np.ndarray
    torsion_quads: np.ndarray
    torsion_flag_ids: np.ndarray

    hook_triplets: np.ndarray
    motor_triplets: np.ndarray
    segment_pair_indices: np.ndarray

    bead_radius_m: float
    b_m: float

    reverse_flagella: np.ndarray
    flag_states: np.ndarray
    torque_signs: np.ndarray
