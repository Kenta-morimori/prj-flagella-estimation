from __future__ import annotations

import itertools
import csv
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from typing import List, Tuple
from zoneinfo import ZoneInfo

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
    minors_px: List[float],
    major_multipliers: List[float],
    speeds: List[Tuple[float, float]],
    bounds: List[Tuple[Tuple[float, float], Tuple[float, float]]],
) -> tuple[List[np.ndarray], np.ndarray]:
    """重ならない楕円を描いた合成フレームを生成する。

    Args:
        num_frames: 総フレーム数。
        image_size: 画像の一辺[px]。
        minors_px: 各オブジェクトの短軸長[px]。
        major_multipliers: 短軸に対する長軸倍率。
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
            minor_px = minors_px[i]
            major_px = minor_px * major_multipliers[i]
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
    """合成映像で重心トラッキングが安定することを確認し、成果物も保存する。"""
    _ = tmp_path  # fixture 使用を明示（出力は outputs/ に保存する）
    num_frames = 200
    image_size = 512
    um_per_px = 0.2
    bac_short_axis_length_um = 1.0
    expected_minor_px = bac_short_axis_length_um / um_per_px
    # 長短軸比を 2〜6倍でばらつかせ、短軸も±20%程度ばらつかせる
    num_objects = 10
    rng = np.random.default_rng(0)
    minors = expected_minor_px * rng.uniform(0.8, 1.2, size=num_objects)
    major_mults = rng.uniform(2.0, 6.0, size=num_objects)
    try:
        now = datetime.now(ZoneInfo("Asia/Tokyo"))
    except Exception:
        now = datetime.now()
    date_dir = now.strftime("%Y-%m-%d")
    time_dir = now.strftime("%H%M%S") + "_test"
    out_dir = Path("outputs") / date_dir / time_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # オブジェクトごとに重ならないY範囲を設定
    # 縦方向に分けた帯域内でジグザグ移動（上下反射＋横方向も反射）
    band_height = (image_size - 120.0) / num_objects
    bounds: list[Tuple[Tuple[float, float], Tuple[float, float]]] = []
    speeds: list[Tuple[float, float]] = []
    for i in range(num_objects):
        y_low = 60.0 + band_height * i
        y_high = y_low + band_height - 10.0
        bounds.append(((60.0, image_size - 60.0), (y_low, y_high)))
        # 速度は 5〜10px/frame 相当で符号も混ぜる
        vx = rng.uniform(5.0, 10.0) * (-1 if i % 2 else 1)
        vy = rng.uniform(5.0, 10.0) * (1 if i % 3 else -1)
        speeds.append((vx, vy))

    frames, gt_positions = _generate_synthetic_frames(
        num_frames=num_frames,
        image_size=image_size,
        minors_px=minors.tolist(),
        major_multipliers=major_mults.tolist(),
        speeds=speeds,
        bounds=bounds,
    )

    # AVI を生成して VideoCapture 経由でも読めることを確認
    video_path = out_dir / "synthetic_input.avi"
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
            max_minor_factor=3.0,
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

    # IDの整合性を平均距離で確認（複数オブジェクトなので全順列を探索）
    track_ids = sorted(track_positions.keys())
    traj_arr = [np.stack(track_positions[tid]) for tid in track_ids]
    errors = np.zeros((len(traj_arr), gt_positions.shape[1]))
    for i, traj in enumerate(traj_arr):
        for j in range(gt_positions.shape[1]):
            errors[i, j] = float(
                np.mean(np.linalg.norm(traj - gt_positions[:, j, :], axis=1))
            )

    best_total = float("inf")
    best_perm: tuple[int, ...] | None = None
    for perm in itertools.permutations(range(gt_positions.shape[1])):
        total = sum(errors[i, perm[i]] for i in range(len(traj_arr)))
        if total < best_total:
            best_total = total
            best_perm = perm

    # 合成ノイズや量子化のぶれを考慮し、平均誤差は <5px を許容
    assert best_total < 5.0
    assert best_perm is not None

    track_to_gt = {track_ids[i]: best_perm[i] for i in range(len(track_ids))}

    # 重心CSVを保存
    csv_path = out_dir / "synthetic_centroids.csv"
    with csv_path.open("w", newline="") as f:
        writer_csv = csv.writer(f)
        writer_csv.writerow(
            ["frame", "track_id", "gt_id", "pred_cx", "pred_cy", "gt_cx", "gt_cy", "error"]
        )
        for tid in track_ids:
            gt_idx = track_to_gt[tid]
            pred_traj = np.stack(track_positions[tid])
            gt_traj = gt_positions[:, gt_idx, :]
            for frame_idx, (pred_pt, gt_pt) in enumerate(zip(pred_traj, gt_traj)):
                err = float(np.linalg.norm(pred_pt - gt_pt))
                writer_csv.writerow(
                    [
                        frame_idx,
                        tid,
                        gt_idx,
                        float(pred_pt[0]),
                        float(pred_pt[1]),
                        float(gt_pt[0]),
                        float(gt_pt[1]),
                        err,
                    ]
                )

    # 予測とGT重心の重ね合わせ映像を保存
    overlay_path = out_dir / "synthetic_overlay.avi"
    overlay_writer = cv2.VideoWriter(
        str(overlay_path),
        cv2.VideoWriter_fourcc(*"MJPG"),
        30.0,
        (image_size, image_size),
    )
    for idx in range(num_frames):
        frame = frames[idx].copy()
        # GTを青
        for gt_idx in range(gt_positions.shape[1]):
            gt_pt = gt_positions[idx, gt_idx]
            cv2.circle(
                frame,
                (int(round(gt_pt[0])), int(round(gt_pt[1]))),
                4,
                (255, 0, 0),
                -1,
            )
            cv2.putText(
                frame,
                f"G{gt_idx}",
                (int(round(gt_pt[0])) + 4, int(round(gt_pt[1])) - 4),
                cv2.FONT_HERSHEY_SIMPLEX,
                0.45,
                (255, 0, 0),
                1,
                cv2.LINE_AA,
            )

        # 予測を緑
        for tid in track_ids:
            pred_pt = track_positions[tid][idx]
            gt_idx = track_to_gt[tid]
            cv2.circle(
                frame,
                (int(round(pred_pt[0])), int(round(pred_pt[1]))),
                4,
                (0, 255, 0),
                -1,
            )
            cv2.putText(
                frame,
                f"T{tid}->G{gt_idx}",
                (int(round(pred_pt[0])) + 4, int(round(pred_pt[1])) + 12),
                cv2.FONT_HERSHEY_SIMPLEX,
                0.45,
                (0, 200, 0),
                1,
                cv2.LINE_AA,
            )

        overlay_writer.write(frame)

    overlay_writer.release()
