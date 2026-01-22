from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
from typing import Any, Dict, List

import cv2
import numpy as np
import pandas as pd

from flagella_estimation.core.run_context import init_run
from flagella_estimation.tracking_butt.butt_estimator import ButtEstimator
from flagella_estimation.tracking_butt.config import (
    DetectionConfig,
    apply_overrides,
    load_config,
    with_save_contour,
)

from flagella_estimation.tracking_butt.detector import choose_invert_flag, detect_frame
from flagella_estimation.tracking_butt.features import FeatureComputer
from flagella_estimation.tracking_butt.overlay import OverlayRenderer
from flagella_estimation.tracking_butt.tracker import Tracker
from flagella_estimation.tracking_butt.types import ButtEstimate, TrackUpdate


def _write_json(path: Path, data: Dict[str, Any]) -> None:
    """JSONをUTF-8・インデント付きで保存する。

    Args:
        path: 保存先パス。
        data: JSONシリアライズ対象。

    Returns:
        None
    """
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def _init_video_writer(
    path: Path, width: int, height: int, fps: float
) -> cv2.VideoWriter:
    """mp4出力用のVideoWriterを初期化する。

    Args:
        path: 出力ファイルパス。
        width: フレーム幅[px]。
        height: フレーム高さ[px]。
        fps: フレームレート。

    Returns:
        cv2.VideoWriter: mp4を書き出すライター。
    """
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    return cv2.VideoWriter(str(path), fourcc, fps, (width, height))


def _butt_record(butt: ButtEstimate) -> Dict[str, float | bool]:
    """ButtEstimate を JSON に保存しやすい辞書へ変換する。

    Args:
        butt: お尻推定結果。

    Returns:
        dict: butt.json 用の辞書。
    """
    return {
        "x": float(butt.point[0]),
        "y": float(butt.point[1]),
        "conf": float(butt.conf),
        "frozen": bool(butt.frozen),
        "flagella_dir_x": float(butt.flagella_dir[0]),
        "flagella_dir_y": float(butt.flagella_dir[1]),
    }


def _save_contour(
    contour_dir: Path, frame_idx: int, track_id: int, contour: np.ndarray
) -> None:
    """輪郭を指定ディレクトリに npy 形式で保存する。

    Args:
        contour_dir: 保存先ディレクトリ。
        frame_idx: フレーム番号。
        track_id: トラックID。
        contour: OpenCVの輪郭配列。

    Returns:
        None
    """
    contour_dir.mkdir(parents=True, exist_ok=True)
    fname = contour_dir / f"frame_{frame_idx:06d}_track_{track_id:04d}.npy"
    np.save(fname, contour[:, 0, :])


def _prepare_track_row(
    update: TrackUpdate, features: Dict[str, float]
) -> Dict[str, Any]:
    """track.csv の1行分の辞書を生成する。

    Args:
        update: トラッキング更新情報。
        features: 追加特徴量の辞書。

    Returns:
        dict: track.csv 用の1行データ。
    """
    row = {
        "frame": update.frame_idx,
        "track_id": update.track_id,
        "cx": update.detection.cx,
        "cy": update.detection.cy,
        "theta": update.detection.theta,
        "major": update.detection.major,
        "minor": update.detection.minor,
        "vx": update.vx,
        "vy": update.vy,
        "is_valid": update.detection.is_valid,
    }
    row.update(features)
    return row


def _process_frames(
    cap: cv2.VideoCapture,
    detection_cfg: DetectionConfig,
    expected_minor_px: float,
    tracker: Tracker,
    butt_estimator: ButtEstimator,
    feature_comp: FeatureComputer,
    overlay: OverlayRenderer,
    writer: cv2.VideoWriter,
    contour_dir: Path | None,
    logger,
) -> tuple[List[Dict[str, Any]], Dict[int, Dict[int, Dict[str, float | bool]]]]:
    """全フレームを処理して検知・トラッキング・お尻推定・描画を行う。

    Args:
        cap: 入力動画キャプチャ。
        detection_cfg: 検知設定。
        expected_minor_px: 期待する短径[px]。
        tracker: トラッカー。
        butt_estimator: お尻推定器。
        feature_comp: 特徴量計算器。
        overlay: オーバーレイ描画器。
        writer: 出力動画ライター。
        contour_dir: 輪郭保存ディレクトリ（無効ならNone）。
        logger: ロガー。

    Returns:
        tuple: (track_rows, butt_store) track_rowsはtrack.csv用の辞書リスト、
            butt_storeは butt.json 用のネスト辞書。
    """
    track_rows: List[Dict[str, Any]] = []
    butt_store: Dict[int, Dict[int, Dict[str, float | bool]]] = {}

    frame_idx = 0
    while True:
        ret, frame = cap.read()
        if not ret:
            break

        detections = detect_frame(
            frame,
            frame_idx,
            detection_cfg,
            expected_minor_px=expected_minor_px,
            logger=logger,
        )
        updates = tracker.step(frame_idx, detections)

        for upd in updates:
            butt = butt_estimator.estimate(upd)
            features = feature_comp.compute(upd)
            track_rows.append(_prepare_track_row(upd, features))

            butt_store.setdefault(upd.track_id, {})[upd.frame_idx] = _butt_record(butt)

            if contour_dir is not None:
                _save_contour(
                    contour_dir, upd.frame_idx, upd.track_id, upd.detection.contour
                )

            overlay.draw(frame, upd.detection, upd.track_id, butt, frame_idx=frame_idx)

        overlay.draw_scale_bar(frame)
        writer.write(frame)
        frame_idx += 1

        if frame_idx % 10 == 0:
            logger.info("Processed %d frames", frame_idx)

    logger.info("Total processed frames: %d", frame_idx)
    return track_rows, butt_store


def run_tracking_butt(
    config_path: Path,
    save_contour_flag: bool = False,
    overrides: Dict[str, Any] | None = None,
) -> None:
    """トラッキング＋お尻推定のエントリーポイント。

    Args:
        config_path: 設定ファイルパス。
        save_contour_flag: 輪郭保存を上書きするかどうか。
        overrides: key=value 形式の設定上書き辞書。

    Returns:
        None
    """
    cfg = load_config(config_path)
    if overrides:
        cfg = apply_overrides(cfg, overrides)
    if save_contour_flag:
        cfg = with_save_contour(cfg, True)

    ctx = init_run(
        base_dir=cfg.output.base_dir,
        input_info={
            "config": str(config_path),
            "video_path": str(cfg.data.video_path),
            "save_contour": cfg.tracking_butt.save.contour,
        },
    )
    logger = ctx.logger
    logger.info("Loaded config: %s", cfg)

    cap = cv2.VideoCapture(str(cfg.data.video_path))
    if not cap.isOpened():
        raise RuntimeError(f"Failed to open video: {cfg.data.video_path}")

    fps = cfg.data.fps if cfg.data.fps > 0 else cap.get(cv2.CAP_PROP_FPS)
    if not fps or fps <= 1e-6:
        fps = 30.0

    ret, first_frame = cap.read()
    if not ret:
        raise RuntimeError("Video contains no frames.")
    height, width = first_frame.shape[:2]

    expected_minor_px = (
        cfg.data.bac_short_axis_length_um * cfg.data.px_per_um
        if cfg.data.bac_short_axis_length_um > 0 and cfg.data.px_per_um > 0
        else 0.0
    )

    detection_cfg = cfg.tracking_butt.detection
    chosen_invert = choose_invert_flag(
        first_frame
        if first_frame.ndim == 2
        else cv2.cvtColor(first_frame, cv2.COLOR_BGR2GRAY),
        detection_cfg,
        expected_minor_px=expected_minor_px,
    )
    detection_cfg = replace(
        detection_cfg, threshold=replace(detection_cfg.threshold, invert=chosen_invert)
    )
    logger.info("Invert selected globally: %s", chosen_invert)

    writer = _init_video_writer(
        ctx.out.tracking_dir / "overlay.mp4", width, height, fps
    )
    cap.set(cv2.CAP_PROP_POS_FRAMES, 0)

    tracker = Tracker(max_link_distance=cfg.tracking_butt.tracking.max_link_distance)
    butt_estimator = ButtEstimator(
        smooth_window=cfg.tracking_butt.butt_estimation.smooth_window,
        freeze_speed_thresh=cfg.tracking_butt.butt_estimation.freeze_speed_thresh,
    )
    feature_comp = FeatureComputer(
        cfg.tracking_butt.butt_estimation.features, logger=logger
    )
    overlay = OverlayRenderer(
        scale_bar_px=expected_minor_px,
        scale_bar_um=cfg.data.bac_short_axis_length_um,
    )
    initial_overlay_path = ctx.out.tracking_dir / "initial_detection.png"
    initial_frame = first_frame.copy()
    init_detections = detect_frame(
        first_frame,
        0,
        detection_cfg,
        expected_minor_px=expected_minor_px,
        logger=logger,
    )
    overlay.draw_scale_bar(initial_frame)
    for idx, det in enumerate(init_detections):
        center_pt = (int(round(det.cx)), int(round(det.cy)))
        if det.is_valid and det.angle_deg is not None and det.major and det.minor:
            ellipse = (
                (float(det.cx), float(det.cy)),
                (float(det.major), float(det.minor)),
                float(det.angle_deg),
            )
            cv2.ellipse(initial_frame, ellipse, (0, 200, 0), 2)
        else:
            x, y, w, h = det.bbox
            cv2.rectangle(initial_frame, (x, y), (x + w, y + h), (0, 150, 0), 1)
        cv2.circle(initial_frame, center_pt, 3, (0, 255, 255), -1)
        cv2.putText(
            initial_frame,
            f"det{idx}",
            (center_pt[0] + 5, center_pt[1] - 5),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.5,
            (255, 255, 255),
            1,
            cv2.LINE_AA,
        )
    cv2.imwrite(str(initial_overlay_path), initial_frame)

    contour_dir = (
        ctx.out.tracking_dir / "contours" if cfg.tracking_butt.save.contour else None
    )

    try:
        track_rows, butt_store = _process_frames(
            cap=cap,
            detection_cfg=detection_cfg,
            expected_minor_px=expected_minor_px,
            tracker=tracker,
            butt_estimator=butt_estimator,
            feature_comp=feature_comp,
            overlay=overlay,
            writer=writer,
            contour_dir=contour_dir,
            logger=logger,
        )
    finally:
        cap.release()
        writer.release()

    base_columns = [
        "frame",
        "track_id",
        "cx",
        "cy",
        "theta",
        "major",
        "minor",
        "vx",
        "vy",
        "is_valid",
    ]
    columns = base_columns + cfg.tracking_butt.butt_estimation.features
    track_df = pd.DataFrame(track_rows, columns=columns)
    track_path = ctx.out.tracking_dir / "track.csv"
    track_df.to_csv(track_path, index=False)

    butt_json = {
        str(k): {str(f): v for f, v in frames.items()}
        for k, frames in butt_store.items()
    }
    _write_json(ctx.out.tracking_dir / "butt.json", butt_json)

    qc_summary = {str(k): v for k, v in tracker.qc_summary().items()}
    _write_json(ctx.out.tracking_dir / "qc.json", qc_summary)

    logger.info("Saved track.csv to %s", track_path)
    logger.info("Saved butt.json and qc.json to %s", ctx.out.tracking_dir)
