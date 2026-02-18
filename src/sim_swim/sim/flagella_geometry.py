"""べん毛のらせん中心線を生成するユーティリティ。"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np

from dataclasses import dataclass

from sim_swim.sim.params import FlagellumParams


@dataclass(frozen=True)
class FlagellaRig:
    base_offsets_body: np.ndarray  # (N,3)
    helix_local: np.ndarray  # (M,3)


def generate_helix_points(
    params: FlagellumParams, step_um: float | None = None
) -> np.ndarray:
    """らせん中心線の3D点列を生成する。

    Args:
        params: べん毛パラメータ。
        step_um: サンプリング間隔（Noneなら params.helix_step_um）。

    Returns:
        (N, 3) のnumpy配列。z方向に長さ L_f、x-y 平面に半径 R のらせん。
    """
    s = step_um or params.helix_step_um
    n_steps = max(2, int(math.ceil(params.length_um / max(s, 1e-6))))
    z = np.linspace(0.0, params.length_um, n_steps)
    theta = 2.0 * math.pi * z / params.pitch_um
    x = params.radius_um * np.cos(theta)
    y = params.radius_um * np.sin(theta)
    return np.stack([x, y, z], axis=1)


def rotate_and_translate(
    points: np.ndarray, rot_mat: np.ndarray, origin: np.ndarray
) -> np.ndarray:
    """点列に回転・並進を適用する。"""
    return (points @ rot_mat.T) + origin


def sample_base_directions(n: int, rng: np.random.Generator) -> Iterable[np.ndarray]:
    """菌体中心から等方的に基部方向を与える。"""
    for _ in range(n):
        v = rng.normal(0.0, 1.0, 3)
        v /= max(np.linalg.norm(v), 1e-9)
        yield v
