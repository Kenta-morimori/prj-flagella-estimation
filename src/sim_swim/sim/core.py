"""bead-spring + RPY のシミュレーションラッパ。"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import List, Tuple

import numpy as np

from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.flagella_geometry import FlagellaRig
from sim_swim.sim.params import SimulationConfig

M_TO_UM = 1e6


def _quat_normalize(q: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(q)
    if norm == 0:
        return np.array([0.0, 0.0, 0.0, 1.0], dtype=float)
    return q / norm


def _quat_multiply(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    """Hamilton積でクォータニオンを合成する。"""
    x1, y1, z1, w1 = q1
    x2, y2, z2, w2 = q2
    return np.array(
        [
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        ],
        dtype=float,
    )


def _rotate_vec(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """クォータニオンでベクトルを回転する。"""
    v_q = np.array([v[0], v[1], v[2], 0.0], dtype=float)
    q_conj = np.array([-q[0], -q[1], -q[2], q[3]], dtype=float)
    return _quat_multiply(_quat_multiply(q, v_q), q_conj)[:3]


def _quat_to_rotmat(q: np.ndarray) -> np.ndarray:
    """クォータニオンから回転行列を計算する。"""
    x, y, z, w = q
    return np.array(
        [
            [1 - 2 * (y**2 + z**2), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x**2 + z**2), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x**2 + y**2)],
        ],
        dtype=float,
    )


def _quat_from_two_vectors(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a_n = a / max(np.linalg.norm(a), 1e-12)
    b_n = b / max(np.linalg.norm(b), 1e-12)
    dot = float(np.clip(np.dot(a_n, b_n), -1.0, 1.0))

    if dot < -0.999999:
        axis = np.array([0.0, 1.0, 0.0], dtype=float)
        if abs(a_n[1]) > 0.9:
            axis = np.array([0.0, 0.0, 1.0], dtype=float)
        axis = axis - np.dot(axis, a_n) * a_n
        axis = axis / max(np.linalg.norm(axis), 1e-12)
        return np.array([axis[0], axis[1], axis[2], 0.0], dtype=float)

    cross = np.cross(a_n, b_n)
    q = np.array([cross[0], cross[1], cross[2], 1.0 + dot], dtype=float)
    return _quat_normalize(q)


def _omega_from_quats(q_prev: np.ndarray, q_curr: np.ndarray, dt: float) -> np.ndarray:
    if dt <= 0:
        return np.zeros(3, dtype=float)
    q_prev_conj = np.array([-q_prev[0], -q_prev[1], -q_prev[2], q_prev[3]], dtype=float)
    dq = _quat_multiply(q_curr, q_prev_conj)
    dq = _quat_normalize(dq)

    if dq[3] < 0.0:
        dq = -dq

    v = dq[:3]
    w = float(np.clip(dq[3], -1.0, 1.0))
    nv = np.linalg.norm(v)
    if nv < 1e-12:
        return np.zeros(3, dtype=float)

    axis = v / nv
    angle = 2.0 * math.atan2(nv, w)
    return axis * (angle / dt)


@dataclass
class SimulationState:
    """1時刻分の観測互換データ。"""

    t: float
    position_um: Tuple[float, float, float]
    quaternion: Tuple[float, float, float, float]
    velocity_um_s: Tuple[float, float, float]
    omega_rad_s: Tuple[float, float, float]


class Simulator:
    """時間ループと出力変換を担う薄いラッパ。"""

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.model = ModelBuilder(config).build()
        self.engine = DynamicsEngine(self.model, config)
        self.rig = FlagellaRig(
            base_offsets_body=self.model.base_offsets_body_um,
            helix_local=self.model.helix_local_um,
        )

    def _observe(self, t: float, prev: SimulationState | None) -> SimulationState:
        body_pts = self.model.positions_m[self.model.body_indices]
        center_m = body_pts.mean(axis=0)
        axis = body_pts[-1] - body_pts[0]
        if np.linalg.norm(axis) < 1e-15:
            axis = np.array([1.0, 0.0, 0.0], dtype=float)
        q = _quat_from_two_vectors(np.array([1.0, 0.0, 0.0], dtype=float), axis)

        if prev is None:
            vel_um_s = np.zeros(3, dtype=float)
            omega = np.zeros(3, dtype=float)
        else:
            dt = max(t - prev.t, 1e-12)
            prev_center_m = np.array(prev.position_um, dtype=float) / M_TO_UM
            vel_um_s = ((center_m - prev_center_m) / dt) * M_TO_UM
            prev_q = np.array(prev.quaternion, dtype=float)
            omega = _omega_from_quats(prev_q, q, dt)

        return SimulationState(
            t=round(t, 9),
            position_um=tuple((center_m * M_TO_UM).tolist()),
            quaternion=tuple(q.tolist()),
            velocity_um_s=tuple(vel_um_s.tolist()),
            omega_rad_s=tuple(omega.tolist()),
        )

    def run(self, duration_s: float) -> List[SimulationState]:
        """与えた時間だけシミュレーションして状態列を返す。"""

        dt = max(self.config.time.dt, 1e-9)
        dt_out = max(self.config.time.dt_out, dt)
        total_steps = max(1, int(math.ceil(duration_s / dt)))

        states: List[SimulationState] = []
        next_sample = 0.0

        for step in range(total_steps + 1):
            t_now = self.engine.t
            if t_now + 1e-12 >= next_sample or step == total_steps:
                prev = states[-1] if states else None
                obs = self._observe(t_now, prev)
                states.append(obs)
                next_sample += dt_out

            if step < total_steps:
                self.engine.step(dt)

        return states
