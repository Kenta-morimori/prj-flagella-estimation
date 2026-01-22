from __future__ import annotations

from typing import Dict, List

import cv2
import numpy as np

from flagella_estimation.tracking_butt.config import DetectionConfig
from flagella_estimation.tracking_butt.types import Detection


def _fit_ellipse(
    contour: np.ndarray,
) -> tuple[float, float, float, float, float] | None:
    """輪郭に楕円をフィットし、中心・軸長・角度を返す。

    Args:
        contour: OpenCV輪郭配列。

    Returns:
        tuple | None: (cx, cy, major, minor, angle_deg)。フィットできない場合はNone。
    """
    if len(contour) < 5:
        return None
    ellipse = cv2.fitEllipse(contour)
    (cx, cy), (major, minor), angle_deg = ellipse
    # Ensure major >= minor for consistency
    if major < minor:
        major, minor = minor, major
        angle_deg = (angle_deg + 90.0) % 180.0
    return float(cx), float(cy), float(major), float(minor), float(angle_deg)


def _apply_preprocess(gray: np.ndarray, cfg: DetectionConfig) -> np.ndarray:
    """前処理を設定に従って適用する。

    Args:
        gray: 入力グレースケール画像。
        cfg: 検知設定。

    Returns:
        np.ndarray: 前処理後の画像。
    """
    if cfg.preprocess.method == "tophat":
        k = max(3, cfg.preprocess.kernel_size | 1)
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
        return cv2.morphologyEx(gray, cv2.MORPH_TOPHAT, kernel)
    if cfg.preprocess.method == "bg_subtract":
        k = max(3, cfg.preprocess.kernel_size | 1)
        blur = cv2.GaussianBlur(gray, (k, k), 0)
        return cv2.subtract(gray, blur)
    return gray


def _threshold_image(
    processed: np.ndarray, cfg: DetectionConfig, invert: bool
) -> np.ndarray:
    """前処理済み画像をOtsuまたは適応的閾値で二値化する。

    Args:
        processed: 前処理済み画像。
        cfg: 検知設定。
        invert: 反転するかどうか。

    Returns:
        np.ndarray: 二値化画像。
    """
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
    """バウンディングボックスが2辺以上の枠に接しているか判定する。

    Args:
        bbox: (x, y, w, h) 形式のバウンディングボックス。
        width: 画像幅。
        height: 画像高さ。

    Returns:
        bool: 2辺以上に接していればTrue。
    """
    x, y, w, h = bbox
    margin = 1
    touch_left = x <= margin
    touch_right = x + w >= width - margin
    touch_top = y <= margin
    touch_bottom = y + h >= height - margin
    touches = sum([touch_left, touch_right, touch_top, touch_bottom])
    return touches >= 2


def _max_area_limit(cfg: DetectionConfig, image_area: float) -> float | None:
    """面積上限をpx指定と面積比指定から計算する。

    Args:
        cfg: 検知設定。
        image_area: 画像全体の面積。

    Returns:
        float | None: 上限面積。指定が無ければNone。
    """
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
    expected_minor_px: float,
    width: int,
    height: int,
) -> tuple[List[Detection], Dict[str, int]]:
    """輪郭をフィルタし Detection オブジェクトへ変換する。

    Args:
        contours: 抽出した輪郭リスト。
        frame_idx: フレーム番号。
        cfg: 検知設定。
        expected_minor_px: 期待する短径[px]。
        width: 画像幅。
        height: 画像高さ。

    Returns:
        tuple: (detections, stats)。stats はフィルタ理由などの統計。
    """
    stats: Dict[str, int] = {"total": len(contours), "kept": 0}
    max_area = _max_area_limit(cfg, image_area=float(width * height))
    detections: List[Detection] = []

    minor_limit = (
        expected_minor_px * cfg.filter.max_minor_factor
        if expected_minor_px > 0 and cfg.filter.max_minor_factor
        else None
    )

    for cnt in contours:
        area = float(cv2.contourArea(cnt))
        x, y, w, h = cv2.boundingRect(cnt)
        reason = None
        if area <= 0 or area < cfg.filter.min_area_px:
            reason = "too_small"
        elif max_area is not None and area > max_area:
            reason = "too_large"
        elif cfg.filter.reject_border_touch and _touches_border(
            (x, y, w, h), width, height
        ):
            reason = "border_touch"

        if reason:
            stats[reason] = stats.get(reason, 0) + 1
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
            if minor_limit is not None and minor > minor_limit:
                stats["too_thick"] = stats.get("too_thick", 0) + 1
                continue
        else:
            theta = np.deg2rad(angle_deg)
            is_valid = True
            if minor_limit is not None and minor > minor_limit:
                stats["too_thick"] = stats.get("too_thick", 0) + 1
                continue

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
    stats["mean_area"] = (
        float(np.mean([d.area for d in detections])) if detections else 0.0
    )
    return detections, stats


def _detect_once(
    gray: np.ndarray,
    frame_idx: int,
    cfg: DetectionConfig,
    invert: bool,
    expected_minor_px: float,
) -> tuple[List[Detection], Dict[str, int]]:
    """1回分の前処理・二値化・輪郭抽出を行う。

    Args:
        gray: 入力グレースケール画像。
        frame_idx: フレーム番号。
        cfg: 検知設定。
        invert: 反転するかどうか。
        expected_minor_px: 期待する短径[px]。

    Returns:
        tuple: (detections, stats) 抽出結果と統計。
    """
    processed = _apply_preprocess(gray, cfg)
    thr = _threshold_image(processed, cfg, invert=invert)
    # 小さな切れ目を繋ぐため軽いクロージング＋膨張を適用
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3))
    thr = cv2.morphologyEx(thr, cv2.MORPH_CLOSE, kernel)
    thr = cv2.dilate(thr, kernel, iterations=1)
    contours, _ = cv2.findContours(thr, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    detections, stats = _filter_and_build(
        contours=contours,
        frame_idx=frame_idx,
        cfg=cfg,
        expected_minor_px=expected_minor_px,
        width=gray.shape[1],
        height=gray.shape[0],
    )
    return detections, stats


def detect_frame(
    frame: np.ndarray,
    frame_idx: int,
    cfg: DetectionConfig,
    expected_minor_px: float,
    logger,
) -> List[Detection]:
    """設定された反転有無で1フレームの菌体検知を行う。

    Args:
        frame: 入力フレーム（BGR またはグレースケール）。
        frame_idx: フレーム番号。
        cfg: 検知設定。
        expected_minor_px: 期待する短径[px]。
        logger: ロガー。

    Returns:
        list[Detection]: 検知結果リスト。
    """
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY) if frame.ndim == 3 else frame
    invert_flag = bool(cfg.threshold.invert)

    detections, stats = _detect_once(
        gray, frame_idx, cfg, invert=invert_flag, expected_minor_px=expected_minor_px
    )

    if frame_idx % 10 == 0:
        logger.info(
            "Processed %d frames. kept=%d/%d",
            frame_idx + 1,
            len(detections),
            stats.get("total", 0),
        )
    return detections


def choose_invert_flag(
    gray: np.ndarray, cfg: DetectionConfig, expected_minor_px: float
) -> bool:
    """初期フレームから極性（反転有無）を1回だけ選択する。

    極性反転しない前提のため、ここで決めた極性を全フレームに固定する。

    Args:
        gray: 初期フレームのグレースケール画像。
        cfg: 検知設定。
        expected_minor_px: 期待する短径[px]。

    Returns:
        bool: 反転する場合はTrue。
    """
    options = {bool(cfg.threshold.invert), not bool(cfg.threshold.invert)}
    scored: list[tuple[float, float, bool]] = []
    for inv in options:
        dets, stats = _detect_once(
            gray, frame_idx=0, cfg=cfg, invert=inv, expected_minor_px=expected_minor_px
        )
        mean_area = stats.get("mean_area", float("inf"))
        score = (mean_area if dets else float("inf"), -len(dets))
        scored.append((score[0], score[1], inv))
    scored.sort(key=lambda x: (x[0], x[1]))
    return scored[0][2]
