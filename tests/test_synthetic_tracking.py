from __future__ import annotations

import itertools
from dataclasses import replace
from pathlib import Path
from typing import List, Tuple

import cv2
import numpy as np

from flagella_estimation.tracking_butt.config import (
    DetectionConfig,
    FilterConfig,
    PreprocessConfig,
    ThresholdConfig,
)
from flagella_estimation.tracking_butt.detector import choose_invert_flag, detect_frame
from flagella_estimation.tracking_butt.tracker import Tracker


class _NullLogger:
    """テスト用のダミーロガー。info/warningを捨てる。"""

    def info(self, *_args, **_kwargs) -> None:
        return

    def warning(self, *_args, **_kwargs) -> None:
        return


def _generate_synthetic_frames(
    num_frames: int,
    image_size: int,
    major_px: float,
    minor_px: float,
    speeds: List[Tuple[float, float]],
    bounds: List[Tuple[Tuple[float, float], Tuple[float, float]]],
) -> tuple[List[np.ndarray], np.ndarray]:
    """重ならない楕円を描いた合成フレームを生成する。

    Args:
        num_frames: 総フレーム数。
        image_size: 画像の一辺[px]。
        major_px: 楕円の長軸長[px]。
        minor_px: 楕円の短軸長[px]。
        speeds: 各オブジェクトの速度ベクトル初期値。
        bounds: 各オブジェクトの移動範囲[(x_low, x_high), (y_low, y_high)]。

    Returns:
        tuple: (frames, gt_positions)。framesはBGR画像リスト、gt_positionsは
        (frame, obj, 2) の中心座標配列。
    """
    num_objects = len(speeds)
    positions = np.array(
        [
            [bounds[i][0][0] + 20.0, bounds[i][1][0] + 20.0]
            for i in range(num_objects)
        ],
        dtype=float,
    )
    velocities = np.array(speeds, dtype=float)

    frames: list[np.ndarray] = []
    gt_positions = np.zeros((num_frames, num_objects, 2), dtype=float)

    for f in range(num_frames):
        frame = np.full((image_size, image_size, 3), 180, dtype=np.uint8)
        for i in range(num_objects):
            center = (int(round(positions[i, 0])), int(round(positions[i, 1])))
            axes = (
                max(1, int(np.ceil(major_px / 2))),
                max(1, int(np.ceil(minor_px / 2))),
            )
            cv2.ellipse(
                frame,
                center,
                axes,
                0.0,
                0.0,
                360.0,
                (60, 60, 60),
                -1,
            )
            gt_positions[f, i] = positions[i]

        frames.append(frame)

        for i in range(num_objects):
            positions[i] += velocities[i]
            # x方向の壁反射
            if positions[i, 0] < bounds[i][0][0]:
                positions[i, 0] = bounds[i][0][0] + (
                    bounds[i][0][0] - positions[i, 0]
                )
                velocities[i, 0] *= -1
            elif positions[i, 0] > bounds[i][0][1]:
                positions[i, 0] = bounds[i][0][1] - (
                    positions[i, 0] - bounds[i][0][1]
                )
                velocities[i, 0] *= -1

            # y方向の壁反射
            if positions[i, 1] < bounds[i][1][0]:
                positions[i, 1] = bounds[i][1][0] + (
                    bounds[i][1][0] - positions[i, 1]
                )
                velocities[i, 1] *= -1
            elif positions[i, 1] > bounds[i][1][1]:
                positions[i, 1] = bounds[i][1][1] - (
                    positions[i, 1] - bounds[i][1][1]
                )
                velocities[i, 1] *= -1

    return frames, gt_positions


def test_synthetic_tracking_stable_ids(tmp_path: Path) -> None:
    """合成映像で重心トラッキングが安定することを確認する。"""
    num_frames = 200
    image_size = 512
    um_per_px = 0.2
    bac_short_axis_length_um = 1.0
    expected_minor_px = bac_short_axis_length_um / um_per_px
    major_px = expected_minor_px * 3.0

    # オブジェクトごとに重ならないY範囲を設定
    bounds = [
        ((60.0, image_size - 60.0), (80.0, 220.0)),
        ((60.0, image_size - 60.0), (300.0, 460.0)),
    ]
    speeds = [(6.0, 7.0), (-7.0, -5.0)]

    frames, gt_positions = _generate_synthetic_frames(
        num_frames=num_frames,
        image_size=image_size,
        major_px=major_px,
        minor_px=expected_minor_px,
        speeds=speeds,
        bounds=bounds,
    )

    # AVI を生成して VideoCapture 経由でも読めることを確認
    video_path = tmp_path / "synthetic.avi"
    writer = cv2.VideoWriter(
        str(video_path),
        cv2.VideoWriter_fourcc(*"MJPG"),
        30.0,
        (image_size, image_size),
    )
    for frame in frames:
        writer.write(frame)
    writer.release()

    detection_cfg = DetectionConfig(
        preprocess=PreprocessConfig(method="none", kernel_size=31),
        threshold=ThresholdConfig(method="otsu", invert=False, block_size=35),
        filter=FilterConfig(
            min_area_px=20.0,
            max_area_px=None,
            max_area_frac=0.05,
            reject_border_touch=True,
        ),
    )

    cap = cv2.VideoCapture(str(video_path))
    assert cap.isOpened()
    ret, first_frame = cap.read()
    assert ret
    first_gray = cv2.cvtColor(first_frame, cv2.COLOR_BGR2GRAY)
    invert_flag = choose_invert_flag(
        first_gray, detection_cfg, expected_minor_px=expected_minor_px
    )
    detection_cfg = replace(
        detection_cfg,
        threshold=replace(detection_cfg.threshold, invert=invert_flag),
    )
    cap.set(cv2.CAP_PROP_POS_FRAMES, 0)

    logger = _NullLogger()
    tracker = Tracker(max_link_distance=25.0)
    track_positions: dict[int, list[np.ndarray]] = {}
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
        assert len(detections) == gt_positions.shape[1]
        updates = tracker.step(frame_idx, detections)
        for upd in updates:
            track_positions.setdefault(upd.track_id, []).append(
                np.array([upd.detection.cx, upd.detection.cy], dtype=float)
            )
        frame_idx += 1

    cap.release()

    assert frame_idx == num_frames
    assert len(track_positions) == gt_positions.shape[1]
    for traj in track_positions.values():
        assert len(traj) == num_frames

    # IDの整合性を平均距離で確認（2オブジェクトなので全順列を探索）
    track_ids = sorted(track_positions.keys())
    traj_arr = [np.stack(track_positions[tid]) for tid in track_ids]
    errors = np.zeros((len(traj_arr), gt_positions.shape[1]))
    for i, traj in enumerate(traj_arr):
        for j in range(gt_positions.shape[1]):
            errors[i, j] = float(
                np.mean(np.linalg.norm(traj - gt_positions[:, j, :], axis=1))
            )

    best_total = float("inf")
    for perm in itertools.permutations(range(gt_positions.shape[1])):
        total = sum(errors[i, perm[i]] for i in range(len(traj_arr)))
        best_total = min(best_total, total)

    assert best_total < 1.0
