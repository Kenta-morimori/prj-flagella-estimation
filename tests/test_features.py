from __future__ import annotations

import logging
import math

import numpy as np
import pytest

from flagella_estimation.tracking_butt.features import FeatureComputer
from flagella_estimation.tracking_butt.types import Detection, TrackUpdate


def _make_detection(theta: float = 0.0) -> Detection:
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
        theta=theta,
        major=10.0,
        minor=5.0,
        angle_deg=0.0,
        is_valid=True,
    )


def test_feature_computation_handles_supported_and_heading_change() -> None:
    logger = logging.getLogger("feature-test")
    computer = FeatureComputer(
        ["speed", "vel_axis_dot", "heading_change"], logger=logger
    )

    det = _make_detection(theta=0.0)
    update = TrackUpdate(
        frame_idx=0,
        track_id=0,
        detection=det,
        vx=3.0,
        vy=4.0,
        prev_velocity=None,
    )
    features = computer.compute(update)
    assert features["speed"] == pytest.approx(5.0)
    assert features["vel_axis_dot"] == pytest.approx(3.0)
    assert math.isnan(features["heading_change"])

    update2 = TrackUpdate(
        frame_idx=1,
        track_id=0,
        detection=det,
        vx=0.0,
        vy=1.0,
        prev_velocity=(1.0, 0.0),
    )
    features2 = computer.compute(update2)
    assert features2["heading_change"] == pytest.approx(math.pi / 2)
