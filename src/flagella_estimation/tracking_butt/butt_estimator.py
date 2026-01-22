from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field
from typing import Deque, Dict

from flagella_estimation.tracking_butt.types import ButtEstimate, TrackUpdate


@dataclass
class ButtState:
    point: tuple[float, float] | None = None
    sign_history: Deque[int] = field(default_factory=deque)
    last_sign: int = 1


def _unit(vec: tuple[float, float]) -> tuple[float, float]:
    """Return unit vector; fallback to upward direction when norm is tiny."""
    norm = math.hypot(vec[0], vec[1])
    if norm < 1e-8:
        return (0.0, -1.0)
    return (vec[0] / norm, vec[1] / norm)


class ButtEstimator:
    def __init__(self, smooth_window: int, freeze_speed_thresh: float) -> None:
        """Initialize butt estimator with smoothing window and freeze threshold."""
        self.smooth_window = max(1, smooth_window)
        self.freeze_speed_thresh = freeze_speed_thresh
        self.states: Dict[int, ButtState] = {}

    def estimate(self, update: TrackUpdate) -> ButtEstimate:
        """Estimate butt point and flagella direction for a track update."""
        state = self.states.setdefault(
            update.track_id,
            ButtState(sign_history=deque(maxlen=self.smooth_window)),
        )
        det = update.detection
        cx, cy = det.cx, det.cy

        if (
            det.theta is None
            or det.major is None
            or det.minor is None
            or det.major <= 0
            or det.minor <= 0
            or not math.isfinite(det.major)
            or not math.isfinite(det.minor)
        ):
            point = state.point or (cx, cy)
            flagella_dir = _unit((point[0] - cx, point[1] - cy))
            state.point = point
            return ButtEstimate(
                point=point, flagella_dir=flagella_dir, conf=0.0, frozen=True
            )

        u = (math.cos(det.theta), math.sin(det.theta))
        half_major = det.major / 2.0
        p1 = (cx + half_major * u[0], cy + half_major * u[1])
        p2 = (cx - half_major * u[0], cy - half_major * u[1])

        speed = math.hypot(update.vx, update.vy)
        if speed < self.freeze_speed_thresh and state.point is not None:
            point = state.point
            frozen = True
            chosen_sign = state.last_sign
        else:
            dot1 = update.vx * (p1[0] - cx) + update.vy * (p1[1] - cy)
            dot2 = update.vx * (p2[0] - cx) + update.vy * (p2[1] - cy)
            sign = 1 if dot1 <= dot2 else -1
            state.sign_history.append(sign)
            chosen_sign = 1 if sum(state.sign_history) >= 0 else -1
            point = p1 if chosen_sign >= 0 else p2
            frozen = False

        state.last_sign = chosen_sign
        state.point = point

        conf = abs(update.vx * u[0] + update.vy * u[1]) / (speed + 1e-8)
        flagella_dir = _unit((point[0] - cx, point[1] - cy))

        return ButtEstimate(
            point=point,
            flagella_dir=flagella_dir,
            conf=float(conf),
            frozen=frozen,
        )
