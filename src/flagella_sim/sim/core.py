"""シミュレーションコアの骨組み（計算ロジックは後続実装）。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

from flagella_sim.sim.params import SimulationConfig


@dataclass
class SimulationState:
    """1時刻分の状態を保持する簡易データ構造。"""

    t: float
    position_um: Tuple[float, float, float]
    quaternion: Tuple[float, float, float, float]
    velocity_um_s: Tuple[float, float, float]
    omega_rad_s: Tuple[float, float, float]


class Simulator:
    """べん毛駆動による遊泳シミュレータ（暫定スタブ）。"""

    def __init__(self, config: SimulationConfig):
        self.config = config

    def run(self, duration_s: float) -> List[SimulationState]:
        """与えられた時間だけ遊泳させ、状態のリストを返す（スタブ）。

        現段階では骨組みのみで、返り値は原点に静止した単一状態となる。
        """

        state = SimulationState(
            t=0.0,
            position_um=(0.0, 0.0, 0.0),
            quaternion=(0.0, 0.0, 0.0, 1.0),
            velocity_um_s=(0.0, 0.0, 0.0),
            omega_rad_s=(0.0, 0.0, 0.0),
        )
        _ = duration_s  # 将来の実装で使用
        return [state]
