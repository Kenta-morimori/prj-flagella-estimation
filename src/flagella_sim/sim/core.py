"""剛体菌体＋剛体べん毛による簡易遊泳シミュレーション（MVP版）。"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import List, Tuple

import numpy as np

from flagella_sim.sim.params import SimulationConfig


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


def _quat_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """軸と角度からクォータニオンを生成する（右手系）。"""
    axis_norm = np.linalg.norm(axis)
    if axis_norm == 0.0 or angle == 0.0:
        return np.array([0.0, 0.0, 0.0, 1.0], dtype=float)
    a = axis / axis_norm
    s = math.sin(angle / 2.0)
    c = math.cos(angle / 2.0)
    return np.array([a[0] * s, a[1] * s, a[2] * s, c], dtype=float)


def _rotate_vec(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """クォータニオンでベクトル v を回転させる。"""
    v_q = np.array([v[0], v[1], v[2], 0.0], dtype=float)
    q_conj = np.array([-q[0], -q[1], -q[2], q[3]], dtype=float)
    return _quat_multiply(_quat_multiply(q, v_q), q_conj)[:3]


def _brownian_sigma_um2_per_s(mu_mpas: float, temp_k: float) -> float:
    """簡易な拡散係数モデル（近似）。"""
    # 目安：水中で 0.3 um^2/s 程度を基準に、粘度に反比例させる。
    _ = temp_k  # 現段階では温度は係数へ反映しない簡略モデル
    return 0.3 / max(mu_mpas, 1e-6)


@dataclass
class SimulationState:
    """1時刻分の状態を保持する簡易データ構造。"""

    t: float
    position_um: Tuple[float, float, float]
    quaternion: Tuple[float, float, float, float]
    velocity_um_s: Tuple[float, float, float]
    omega_rad_s: Tuple[float, float, float]


class Simulator:
    """べん毛駆動による簡易遊泳シミュレータ。"""

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.rng = np.random.default_rng(config.seed.global_seed)

    def run(self, duration_s: float) -> List[SimulationState]:
        """与えられた時間だけ遊泳させ、状態リストを返す。"""

        cfg = self.config
        dt_sim = float(cfg.time.dt_sim)
        t = 0.0

        # 初期姿勢を一様乱数で設定
        q = _quat_normalize(
            np.array(
                self.rng.normal(0.0, 1.0, 4),
                dtype=float,
            )
        )
        pos = np.zeros(3, dtype=float)

        swim_coeff = 0.0025  # [um/s] per (flagella * Hz) の経験的スケール
        spin_coeff = 0.05  # [rad/s] per Hz
        base_diff = _brownian_sigma_um2_per_s(
            cfg.env.viscosity_mpas, cfg.env.temperature_k
        )
        # 推進がある場合はブラウンを弱め、モータ依存性が見えやすいようにする
        trans_diff = base_diff * (1.0 if cfg.flagella.n_flagella == 0 else 0.05)
        rot_diff = 0.05 * (1.0 if cfg.flagella.n_flagella == 0 else 0.2)  # [rad^2/s]

        states: List[SimulationState] = []
        total_steps = max(1, int(math.ceil(duration_s / dt_sim)))

        for _ in range(total_steps + 1):
            # 現在の軸（菌体 x 軸）方向
            body_axis = _rotate_vec(q, np.array([1.0, 0.0, 0.0]))

            # 推進速度
            v_prop_mag = (
                swim_coeff * cfg.flagella.n_flagella * cfg.flagella.motor_freq_hz
            )
            v_prop = v_prop_mag * body_axis  # um/s

            # 角速度（counter-rotation）
            omega = -spin_coeff * cfg.flagella.motor_freq_hz * body_axis  # rad/s

            # ブラウン運動
            disp_brown = np.zeros(3, dtype=float)
            omega_brown = np.zeros(3, dtype=float)
            if cfg.env.include_brownian or cfg.flagella.n_flagella == 0:
                sigma = math.sqrt(max(0.0, 2.0 * trans_diff * dt_sim))
                disp_brown = self.rng.normal(0.0, sigma, 3)
                rot_sigma = math.sqrt(max(0.0, 2.0 * rot_diff * dt_sim))
                omega_brown = self.rng.normal(0.0, rot_sigma, 3)

            # 並進更新
            pos = pos + v_prop * dt_sim + disp_brown

            # 回転更新
            total_omega = omega + omega_brown
            omega_norm = np.linalg.norm(total_omega)
            dq = _quat_from_axis_angle(
                total_omega
                if omega_norm == 0
                else total_omega / max(omega_norm, 1e-12),
                omega_norm * dt_sim,
            )
            q = _quat_normalize(_quat_multiply(dq, q))

            # 全ステップ記録
            states.append(
                SimulationState(
                    t=round(t, 9),
                    position_um=tuple(pos.tolist()),
                    quaternion=tuple(q.tolist()),
                    velocity_um_s=tuple(v_prop.tolist()),
                    omega_rad_s=tuple(omega.tolist()),
                )
            )

            t += dt_sim

        return states
