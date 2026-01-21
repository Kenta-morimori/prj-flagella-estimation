from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List

import cv2
import numpy as np
import pandas as pd

from flagella_estimation.core.run_context import init_run
from flagella_estimation.tracking_butt.butt_estimator import ButtEstimator
from flagella_estimation.tracking_butt.config import load_config, with_save_contour
from flagella_estimation.tracking_butt.detector import detect_frame
from flagella_estimation.tracking_butt.features import FeatureComputer
from flagella_estimation.tracking_butt.overlay import OverlayRenderer
from flagella_estimation.tracking_butt.tracker import Tracker
from flagella_estimation.tracking_butt.types import ButtEstimate, TrackUpdate


def _write_json(path: Path, data: Dict[str, Any]) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def _init_video_writer(path: Path, width: int, height: int, fps: float) -> cv2.VideoWriter:
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    return cv2.VideoWriter(str(path), fourcc, fps, (width, height))


def _butt_record(butt: ButtEstimate) -> Dict[str, float | bool]:
    return {
        "x": float(butt.point[0]),
        "y": float(butt.point[1]),
        "conf": float(butt.conf),
        "frozen": bool(butt.frozen),
        "flagella_dir_x": float(butt.flagella_dir[0]),
        "flagella_dir_y": float(butt.flagella_dir[1]),
    }


def _save_contour(contour_dir: Path, frame_idx: int, track_id: int, contour: np.ndarray) -> None:
    contour_dir.mkdir(parents=True, exist_ok=True)
    fname = contour_dir / f"frame_{frame_idx:06d}_track_{track_id:04d}.npy"
    np.save(fname, contour[:, 0, :])


def _prepare_track_row(update: TrackUpdate, features: Dict[str, float]) -> Dict[str, Any]:
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
    tracker: Tracker,
    butt_estimator: ButtEstimator,
    feature_comp: FeatureComputer,
    overlay: OverlayRenderer,
    writer: cv2.VideoWriter,
    contour_dir: Path | None,
    logger,
) -> tuple[List[Dict[str, Any]], Dict[int, Dict[int, Dict[str, float | bool]]]]:
    track_rows: List[Dict[str, Any]] = []
    butt_store: Dict[int, Dict[int, Dict[str, float | bool]]] = {}

    frame_idx = 0
    while True:
        ret, frame = cap.read()
        if not ret:
            break

        detections = detect_frame(frame, frame_idx)
        updates = tracker.step(frame_idx, detections)

        for upd in updates:
            butt = butt_estimator.estimate(upd)
            features = feature_comp.compute(upd)
            track_rows.append(_prepare_track_row(upd, features))

            butt_store.setdefault(upd.track_id, {})[upd.frame_idx] = _butt_record(butt)

            if contour_dir is not None:
                _save_contour(contour_dir, upd.frame_idx, upd.track_id, upd.detection.contour)

            overlay.draw(frame, upd.detection, upd.track_id, butt)

        writer.write(frame)
        frame_idx += 1

        if frame_idx % 50 == 0:
            logger.info("Processed %d frames", frame_idx)

    logger.info("Total processed frames: %d", frame_idx)
    return track_rows, butt_store


def run_tracking_butt(config_path: Path, save_contour_flag: bool = False) -> None:
    cfg = load_config(config_path)
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

    writer = _init_video_writer(ctx.out.tracking_dir / "overlay.mp4", width, height, fps)
    cap.set(cv2.CAP_PROP_POS_FRAMES, 0)

    tracker = Tracker(max_link_distance=cfg.tracking_butt.tracking.max_link_distance)
    butt_estimator = ButtEstimator(
        smooth_window=cfg.tracking_butt.butt_estimation.smooth_window,
        freeze_speed_thresh=cfg.tracking_butt.butt_estimation.freeze_speed_thresh,
    )
    feature_comp = FeatureComputer(
        cfg.tracking_butt.butt_estimation.features, logger=logger
    )
    overlay = OverlayRenderer()

    contour_dir = (
        ctx.out.tracking_dir / "contours" if cfg.tracking_butt.save.contour else None
    )

    try:
        track_rows, butt_store = _process_frames(
            cap=cap,
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

    butt_json = {str(k): {str(f): v for f, v in frames.items()} for k, frames in butt_store.items()}
    _write_json(ctx.out.tracking_dir / "butt.json", butt_json)

    qc_summary = {str(k): v for k, v in tracker.qc_summary().items()}
    _write_json(ctx.out.tracking_dir / "qc.json", qc_summary)

    logger.info("Saved track.csv to %s", track_path)
    logger.info("Saved butt.json and qc.json to %s", ctx.out.tracking_dir)
