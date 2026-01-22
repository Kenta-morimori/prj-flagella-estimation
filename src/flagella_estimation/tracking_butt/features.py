from __future__ import annotations

import logging
import math
from typing import Callable, Dict, List

from flagella_estimation.tracking_butt.types import TrackUpdate


def _speed(update: TrackUpdate) -> float:
    """速度の大きさを計算する。

    Args:
        update: トラック更新情報。

    Returns:
        float: 速度ノルム。
    """
    return math.hypot(update.vx, update.vy)


def _vel_axis_dot(update: TrackUpdate) -> float:
    """速度ベクトルと主軸方向の内積を計算する。

    Args:
        update: トラック更新情報。

    Returns:
        float: 内積。theta が無い場合は NaN。
    """
    if update.detection.theta is None:
        return math.nan
    return update.vx * math.cos(update.detection.theta) + update.vy * math.sin(
        update.detection.theta
    )


def _heading_change(update: TrackUpdate) -> float:
    """前フレームと比べた速度方向の角度差を計算する。

    Args:
        update: トラック更新情報。

    Returns:
        float: 角度差（ラジアン）。計算できない場合は NaN。
    """
    if update.prev_velocity is None:
        return math.nan
    prev_vx, prev_vy = update.prev_velocity
    if math.hypot(prev_vx, prev_vy) < 1e-6 or math.hypot(update.vx, update.vy) < 1e-6:
        return math.nan
    prev_angle = math.atan2(prev_vy, prev_vx)
    curr_angle = math.atan2(update.vy, update.vx)
    diff = curr_angle - prev_angle
    if diff > math.pi:
        diff -= 2 * math.pi
    elif diff < -math.pi:
        diff += 2 * math.pi
    return diff


SUPPORTED_FEATURES: Dict[str, Callable[[TrackUpdate], float]] = {
    "speed": _speed,
    "vel_axis_dot": _vel_axis_dot,
    "heading_change": _heading_change,
}


class FeatureComputer:
    def __init__(self, requested: List[str], logger: logging.Logger) -> None:
        """特徴量計算器を初期化する。

        Args:
            requested: 計算したい特徴量名のリスト。
            logger: ロガー。
        """
        self.requested = requested
        self.logger = logger
        self._warned: set[str] = set()

    def compute(self, update: TrackUpdate) -> Dict[str, float]:
        """要求された特徴量を計算する。未対応は警告してNaNを埋める。

        Args:
            update: トラック更新情報。

        Returns:
            dict: 特徴量名をキー、値を値とする辞書。
        """
        values: Dict[str, float] = {}
        for name in self.requested:
            func = SUPPORTED_FEATURES.get(name)
            if func is None:
                if name not in self._warned:
                    self.logger.warning(
                        "Feature '%s' not implemented; filling with NaN.", name
                    )
                    self._warned.add(name)
                values[name] = math.nan
            else:
                values[name] = func(update)
        return values
