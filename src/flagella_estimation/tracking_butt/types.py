from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray


@dataclass
class Detection:
    frame_idx: int
    contour: NDArray[np.int32]
    area: float
    bbox: tuple[int, int, int, int]
    cx: float
    cy: float
    theta: float | None
    major: float | None
    minor: float | None
    angle_deg: float | None
    is_valid: bool


@dataclass
class TrackUpdate:
    frame_idx: int
    track_id: int
    detection: Detection
    vx: float
    vy: float
    prev_velocity: Optional[Tuple[float, float]]


@dataclass
class ButtEstimate:
    point: tuple[float, float]
    flagella_dir: tuple[float, float]
    conf: float
    frozen: bool
