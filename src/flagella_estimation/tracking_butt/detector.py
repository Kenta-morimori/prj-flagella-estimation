from __future__ import annotations

from typing import List

import cv2
import numpy as np

from flagella_estimation.tracking_butt.types import Detection


def _fit_ellipse(
    contour: np.ndarray,
) -> tuple[float, float, float, float, float] | None:
    if len(contour) < 5:
        return None
    ellipse = cv2.fitEllipse(contour)
    (cx, cy), (major, minor), angle_deg = ellipse
    # Ensure major >= minor for consistency
    if major < minor:
        major, minor = minor, major
        angle_deg = angle_deg + 90.0
    return float(cx), float(cy), float(major), float(minor), float(angle_deg)


def detect_frame(frame: np.ndarray, frame_idx: int, min_area: float = 20.0) -> List[Detection]:
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) if frame.ndim == 3 else frame
    blur = cv2.GaussianBlur(gray, (5, 5), 0)

    _, thr = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    contours, _ = cv2.findContours(thr, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        _, thr_inv = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
        contours, _ = cv2.findContours(thr_inv, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    detections: List[Detection] = []
    for cnt in contours:
        area = float(cv2.contourArea(cnt))
        if area < min_area:
            continue
        x, y, w, h = cv2.boundingRect(cnt)
        ellipse = _fit_ellipse(cnt)

        if ellipse is not None:
            cx, cy, major, minor, angle_deg = ellipse
            theta = np.deg2rad(angle_deg)
            is_valid = True
        else:
            cx = x + w / 2.0
            cy = y + h / 2.0
            major = float(max(w, h))
            minor = float(min(w, h))
            theta = None
            angle_deg = None
            is_valid = False

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
    return detections
