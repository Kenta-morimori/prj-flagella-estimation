from __future__ import annotations

import numpy as np
import pytest

from flagella_estimation.tracking_butt.butt_estimator import ButtEstimator
from flagella_estimation.tracking_butt.types import Detection, TrackUpdate


def _make_detection() -> Detection:
    contour = np.array(
        [[[0, 0]], [[1, 0]], [[1, 1]], [[0, 1]], [[0, 0]]], dtype=np.int32
    )
    return Detection(
        frame_idx=0,
        contour=contour,
        area=4.0,
        bbox=(0, 0, 2, 2),
        cx=0.0,
        cy=0.0,
        theta=0.0,
        major=10.0,
        minor=5.0,
        angle_deg=0.0,
        is_valid=True,
    )


def test_butt_estimator_chooses_opposite_to_velocity_and_freezes_when_slow() -> None:
    det = _make_detection()
    estimator = ButtEstimator(smooth_window=3, freeze_speed_thresh=0.5)

    update = TrackUpdate(
        frame_idx=0,
        track_id=0,
        detection=det,
        vx=1.0,
        vy=0.0,
        prev_velocity=None,
    )
    butt1 = estimator.estimate(update)
    assert butt1.point[0] < det.cx  # butt should be opposite to +x motion
    assert butt1.frozen is False

    frozen_update = TrackUpdate(
        frame_idx=1,
        track_id=0,
        detection=det,
        vx=0.0,
        vy=0.0,
        prev_velocity=(1.0, 0.0),
    )
    butt2 = estimator.estimate(frozen_update)
    assert butt2.frozen is True
    assert butt2.point == pytest.approx(butt1.point)
