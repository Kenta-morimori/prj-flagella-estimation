from __future__ import annotations

import logging

import cv2
import numpy as np

from flagella_estimation.tracking_butt.config import (
    DetectionConfig,
    FilterConfig,
    PreprocessConfig,
    ThresholdConfig,
)
from flagella_estimation.tracking_butt.detector import detect_frame


def _make_cfg(
    min_area: float = 10.0,
    max_area_frac: float | None = 0.2,
    reject_border: bool = True,
    invert: bool = False,
) -> DetectionConfig:
    return DetectionConfig(
        preprocess=PreprocessConfig(method="none", kernel_size=5),
        threshold=ThresholdConfig(method="otsu", invert=invert, block_size=35),
        filter=FilterConfig(
            min_area_px=min_area,
            max_area_px=None,
            max_area_frac=max_area_frac,
            reject_border_touch=reject_border,
        ),
    )


def test_detect_frame_filters_large_and_border_touch() -> None:
    logger = logging.getLogger("detector-test")

    cfg = _make_cfg(max_area_frac=0.2)
    frame = np.zeros((100, 100), dtype=np.uint8)
    cv2.rectangle(frame, (30, 30), (70, 70), 255, -1)
    detections = detect_frame(frame, 0, cfg, logger=logger)
    assert len(detections) == 1

    cfg_large = _make_cfg(max_area_frac=0.05)
    full = np.full((100, 100), 255, dtype=np.uint8)
    large_dets = detect_frame(full, 1, cfg_large, logger=logger)
    assert len(large_dets) == 0

    cfg_border = _make_cfg(max_area_frac=0.2, reject_border=True)
    border_frame = np.zeros((100, 100), dtype=np.uint8)
    cv2.rectangle(border_frame, (0, 0), (90, 10), 255, -1)  # touches top/left
    border_dets = detect_frame(border_frame, 2, cfg_border, logger=logger)
    assert len(border_dets) == 0
