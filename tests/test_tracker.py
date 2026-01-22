from __future__ import annotations

import numpy as np

from flagella_estimation.tracking_butt.tracker import Tracker
from flagella_estimation.tracking_butt.types import Detection


def _make_detection(frame_idx: int, cx: float, cy: float) -> Detection:
    contour = np.array(
        [[[0, 0]], [[1, 0]], [[1, 1]], [[0, 1]], [[0, 0]]], dtype=np.int32
    )
    return Detection(
        frame_idx=frame_idx,
        contour=contour,
        area=4.0,
        bbox=(0, 0, 2, 2),
        cx=cx,
        cy=cy,
        theta=0.0,
        major=10.0,
        minor=5.0,
        angle_deg=0.0,
        is_valid=True,
    )


def test_tracker_links_nearest_and_counts_missing_frames() -> None:
    tracker = Tracker(max_link_distance=15.0, max_inactive=5)

    det0 = _make_detection(0, 0.0, 0.0)
    updates0 = tracker.step(0, [det0])
    assert updates0[0].track_id == 0
    assert updates0[0].vx == 0

    det1 = _make_detection(1, 3.0, 4.0)
    updates1 = tracker.step(1, [det1])
    assert updates1[0].track_id == 0
    assert updates1[0].vx == 3
    assert updates1[0].vy == 4
    assert tracker.tracks[0].missing_frames == 0

    det2 = _make_detection(3, 3.0, 4.0)
    updates2 = tracker.step(3, [det2])
    assert updates2[0].track_id == 0
    assert tracker.tracks[0].missing_frames == 1  # frame 2 missing
