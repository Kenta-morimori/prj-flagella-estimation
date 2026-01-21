from __future__ import annotations

from typing import Dict, List

import cv2
import numpy as np

from flagella_estimation.tracking_butt.config import DetectionConfig
from flagella_estimation.tracking_butt.types import Detection


def _fit_ellipse(
    contour: np.ndarray,
) -> tuple[float, float, float, float, float] | None:
    """Fit an ellipse to a contour, returning center, axes, angle or None."""
    if len(contour) < 5:
        return None
    ellipse = cv2.fitEllipse(contour)
    (cx, cy), (major, minor), angle_deg = ellipse
    # Ensure major >= minor for consistency
    if major < minor:
        major, minor = minor, major
        angle_deg = angle_deg + 90.0
    return float(cx), float(cy), float(major), float(minor), float(angle_deg)


def _apply_preprocess(gray: np.ndarray, cfg: DetectionConfig) -> np.ndarray:
    """Apply background correction according to config."""
    if cfg.preprocess.method == "tophat":
        k = max(3, cfg.preprocess.kernel_size | 1)
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
        return cv2.morphologyEx(gray, cv2.MORPH_TOPHAT, kernel)
    if cfg.preprocess.method == "bg_subtract":
        k = max(3, cfg.preprocess.kernel_size | 1)
        blur = cv2.GaussianBlur(gray, (k, k), 0)
        return cv2.subtract(gray, blur)
    return gray


def _threshold_image(processed: np.ndarray, cfg: DetectionConfig, invert: bool) -> np.ndarray:
    """Threshold preprocessed image using Otsu or adaptive; allows inversion."""
    if cfg.threshold.method == "adaptive":
        block_size = max(3, cfg.threshold.block_size | 1)
        return cv2.adaptiveThreshold(
            processed,
            255,
            cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
            cv2.THRESH_BINARY_INV if invert else cv2.THRESH_BINARY,
            block_size,
            2,
        )

    flag = cv2.THRESH_BINARY_INV if invert else cv2.THRESH_BINARY
    _, thr = cv2.threshold(processed, 0, 255, flag + cv2.THRESH_OTSU)
    return thr


def _touches_border(bbox: tuple[int, int, int, int], width: int, height: int) -> bool:
    """Return True if bbox touches two or more image borders."""
    x, y, w, h = bbox
    margin = 1
    touch_left = x <= margin
    touch_right = x + w >= width - margin
    touch_top = y <= margin
    touch_bottom = y + h >= height - margin
    touches = sum([touch_left, touch_right, touch_top, touch_bottom])
    return touches >= 2


def _max_area_limit(cfg: DetectionConfig, image_area: float) -> float | None:
    """Compute maximum allowed area from px and frac limits."""
    candidates = []
    if cfg.filter.max_area_px:
        candidates.append(cfg.filter.max_area_px)
    if cfg.filter.max_area_frac and cfg.filter.max_area_frac > 0:
        candidates.append(cfg.filter.max_area_frac * image_area)
    if not candidates:
        return None
    return min(candidates)


def _filter_and_build(
    contours: list[np.ndarray],
    frame_idx: int,
    cfg: DetectionConfig,
    logger,
    width: int,
    height: int,
) -> tuple[List[Detection], Dict[str, int]]:
    """Filter contours by size/border and build Detection objects."""
    stats: Dict[str, int] = {"total": len(contours), "kept": 0}
    max_area = _max_area_limit(cfg, image_area=float(width * height))
    detections: List[Detection] = []

    for cnt in contours:
        area = float(cv2.contourArea(cnt))
        x, y, w, h = cv2.boundingRect(cnt)
        reason = None
        if area <= 0 or area < cfg.filter.min_area_px:
            reason = "too_small"
        elif max_area is not None and area > max_area:
            reason = "too_large"
        elif cfg.filter.reject_border_touch and _touches_border((x, y, w, h), width, height):
            reason = "border_touch"

        if reason:
            stats[reason] = stats.get(reason, 0) + 1
            logger.info(
                "frame %d contour skipped (%s): area=%.1f bbox=(%d,%d,%d,%d)",
                frame_idx,
                reason,
                area,
                x,
                y,
                w,
                h,
            )
            continue

        ellipse = _fit_ellipse(cnt)
        if ellipse is not None:
            cx, cy, major, minor, angle_deg = ellipse
            if not (
                np.isfinite(cx)
                and np.isfinite(cy)
                and np.isfinite(major)
                and np.isfinite(minor)
                and major > 0
                and minor > 0
            ):
                ellipse = None

        if ellipse is None:
            cx = x + w / 2.0
            cy = y + h / 2.0
            major = float(max(w, h))
            minor = float(min(w, h))
            theta = None
            angle_deg = None
            is_valid = False
        else:
            theta = np.deg2rad(angle_deg)
            is_valid = True

        detections.append(
            Detection(
                frame_idx=frame_idx,
                contour=cnt,
                area=area,
                bbox=(int(x), int(y), int(w), int(h)),
                cx=float(cx),
                cy=float(cy),
                theta=float(theta) if theta is not None else None,
                major=float(major) if major is not None else None,
                minor=float(minor) if minor is not None else None,
                angle_deg=float(angle_deg) if angle_deg is not None else None,
                is_valid=is_valid,
            )
        )
    stats["kept"] = len(detections)
    return detections, stats


def _detect_once(
    gray: np.ndarray,
    frame_idx: int,
    cfg: DetectionConfig,
    invert: bool,
    logger,
) -> tuple[List[Detection], Dict[str, int]]:
    """Run preprocessing, thresholding, contour extraction once."""
    processed = _apply_preprocess(gray, cfg)
    thr = _threshold_image(processed, cfg, invert=invert)
    contours, _ = cv2.findContours(thr, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    detections, stats = _filter_and_build(
        contours=contours,
        frame_idx=frame_idx,
        cfg=cfg,
        logger=logger,
        width=gray.shape[1],
        height=gray.shape[0],
    )
    return detections, stats


def detect_frame(
    frame: np.ndarray,
    frame_idx: int,
    cfg: DetectionConfig,
    logger,
) -> List[Detection]:
    """Detect multiple cells in a frame with optional invert fallback."""
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) if frame.ndim == 3 else frame
    invert_flag = bool(cfg.threshold.invert)

    detections, stats = _detect_once(gray, frame_idx, cfg, invert=invert_flag, logger=logger)
    fallback_used = False
    if not detections:
        fallback_used = True
        detections, stats = _detect_once(
            gray, frame_idx, cfg, invert=not invert_flag, logger=logger
        )
        logger.info("frame %d invert retry executed (kept=%d)", frame_idx, len(detections))

    logger.info(
        "frame %d detections kept=%d/%d%s",
        frame_idx,
        len(detections),
        stats.get("total", 0),
        " (fallback invert)" if fallback_used else "",
    )
    return detections
