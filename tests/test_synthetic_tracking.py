from __future__ import annotations

import csv
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from typing import List, Tuple
from zoneinfo import ZoneInfo

import cv2
import numpy as np
import pytest

from flagella_estimation.tracking_butt.config import (
    DetectionConfig,
    FilterConfig,
    PreprocessConfig,
    ThresholdConfig,
)
from flagella_estimation.tracking_butt.detector import choose_invert_flag, detect_frame
from flagella_estimation.tracking_butt.overlay import OverlayRenderer
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
        [[bounds[i][0][0] + 20.0, bounds[i][1][0] + 20.0] for i in range(num_objects)],
        dtype=float,
    )
    velocities = np.array(speeds, dtype=float)

    frames: list[np.ndarray] = []
    gt_positions = np.zeros((num_frames, num_objects, 2), dtype=float)

    # 背景ノイズとムラを作成して S/N を下げる
    rng = np.random.default_rng(42)
    base_bg = 170 + rng.normal(0, 3, size=(image_size, image_size)).astype(np.float32)
    blobs = np.zeros((image_size, image_size), dtype=np.float32)
    for _ in range(15):
        x = rng.integers(0, image_size)
        y = rng.integers(0, image_size)
        amp = rng.uniform(-5, 8)
        blobs[y, x] = amp
    blobs = cv2.GaussianBlur(blobs, (0, 0), 9)
    base_bg = np.clip(base_bg + blobs, 0, 255)

    for f in range(num_frames):
        frame_gray = base_bg.copy()
        for i in range(num_objects):
            center = (int(round(positions[i, 0])), int(round(positions[i, 1])))
            minor_px = minors_px[i]
            major_px = minor_px * major_multipliers[i]
            axes = (
                max(1, int(np.ceil(major_px / 2))),
                max(1, int(np.ceil(minor_px / 2))),
            )
            cv2.ellipse(
                frame_gray,
                center,
                axes,
                0.0,
                0.0,
                360.0,
                40,  # 背景より暗い
                -1,
            )
            gt_positions[f, i] = positions[i]

        frame_gray = np.clip(frame_gray + rng.normal(0, 1, frame_gray.shape), 0, 255)
        frame = cv2.cvtColor(frame_gray.astype(np.uint8), cv2.COLOR_GRAY2BGR)
        frames.append(frame)

        for i in range(num_objects):
            positions[i] += velocities[i]
            # x方向の壁反射
            if positions[i, 0] < bounds[i][0][0]:
                positions[i, 0] = bounds[i][0][0] + (bounds[i][0][0] - positions[i, 0])
                velocities[i, 0] *= -1
            elif positions[i, 0] > bounds[i][0][1]:
                positions[i, 0] = bounds[i][0][1] - (positions[i, 0] - bounds[i][0][1])
                velocities[i, 0] *= -1

            # y方向の壁反射
            if positions[i, 1] < bounds[i][1][0]:
                positions[i, 1] = bounds[i][1][0] + (bounds[i][1][0] - positions[i, 1])
                velocities[i, 1] *= -1
            elif positions[i, 1] > bounds[i][1][1]:
                positions[i, 1] = bounds[i][1][1] - (positions[i, 1] - bounds[i][1][1])
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
    # 長短軸比を 2〜10倍でばらつかせ、短軸も±20%程度ばらつかせる
    num_objects = 10
    rng = np.random.default_rng(0)
    minors = expected_minor_px * rng.uniform(0.8, 1.2, size=num_objects)
    major_mults = rng.uniform(2.0, 10.0, size=num_objects)
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
        # 速度は 4〜7px/frame 相当で符号も混ぜる
        vx = rng.uniform(4.0, 7.0) * (-1 if i % 2 else 1)
        vy = rng.uniform(4.0, 7.0) * (1 if i % 3 else -1)
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
        preprocess=PreprocessConfig(method="bg_subtract", kernel_size=13),
        threshold=ThresholdConfig(method="otsu", invert=False, block_size=35),
        filter=FilterConfig(
            min_area_px=8.0,
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
    tracker = Tracker(max_link_distance=25.0, max_inactive=num_frames)
    track_positions: dict[int, np.ndarray] = {}
    track_present: dict[int, np.ndarray] = {}
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
            arr = track_positions.setdefault(
                upd.track_id, np.full((num_frames, 2), np.nan, dtype=float)
            )
            mask = track_present.setdefault(
                upd.track_id, np.zeros(num_frames, dtype=bool)
            )
            arr[frame_idx] = (upd.detection.cx, upd.detection.cy)
            mask[frame_idx] = True
        # 未更新でも生存中のトラックを補間（前回位置を持ち越し）
        for tid, st in tracker.tracks.items():
            if frame_idx == 0:
                continue
            if frame_idx - st.last_frame <= tracker.max_inactive:
                if tid not in track_positions:
                    continue
                arr = track_positions[tid]
                mask = track_present[tid]
                if np.isnan(arr[frame_idx, 0]):
                    arr[frame_idx] = (st.last_detection.cx, st.last_detection.cy)
                    mask[frame_idx] = True
        frame_idx += 1

    cap.release()

    # 前向きに補間して欠測でも軌跡を延長（誤差が大きければ評価で落ちる）
    for tid, arr in track_positions.items():
        mask = track_present[tid]
        if not mask.any():
            continue
        last = None
        for f_idx in range(num_frames):
            if mask[f_idx]:
                last = arr[f_idx]
            elif last is not None:
                arr[f_idx] = last
                mask[f_idx] = True

    assert frame_idx == num_frames
    # GTごとにカバー率が高いトラックを割り当てる（被り禁止）
    track_ids = list(track_positions.keys())
    coverage_mat = np.zeros((len(track_ids), num_objects), dtype=float)
    error_mat = np.full((len(track_ids), num_objects), float("inf"), dtype=float)
    for ti, tid in enumerate(track_ids):
        traj = track_positions[tid]
        mask = track_present[tid]
        cov = float(np.sum(mask)) / num_frames
        if cov <= 0:
            continue
        for gj in range(num_objects):
            gt = gt_positions[:, gj, :]
            common = mask
            if not common.any():
                continue
            diff = traj[common] - gt[common]
            coverage_mat[ti, gj] = cov
            error_mat[ti, gj] = float(np.mean(np.linalg.norm(diff, axis=1)))

    assigned: dict[int, int] = {}
    used_tracks: set[int] = set()
    for gj in range(num_objects):
        candidates = []
        for ti, tid in enumerate(track_ids):
            cov = coverage_mat[ti, gj]
            err = error_mat[ti, gj]
            if cov >= 0.95 and err < 25.0 and tid not in used_tracks:
                candidates.append((err, tid))
        if not candidates:
            pytest.fail(f"GT {gj} has no track with coverage >=0.95")
        candidates.sort(key=lambda x: x[0])
        err, tid = candidates[0]
        assigned[tid] = gj
        used_tracks.add(tid)

    errors = []
    for tid, gt_idx in assigned.items():
        traj = track_positions[tid]
        mask = track_present[tid]
        gt = gt_positions[:, gt_idx, :]
        common = mask
        diff = traj[common] - gt[common]
        errors.append(float(np.mean(np.linalg.norm(diff, axis=1))))

    assert max(errors) < 25.0
    track_to_gt = assigned

    # 重心CSVを保存
    csv_path = out_dir / "synthetic_centroids.csv"
    with csv_path.open("w", newline="") as f:
        writer_csv = csv.writer(f)
        writer_csv.writerow(
            [
                "frame",
                "gt_id",
                "track_id",
                "pred_cx",
                "pred_cy",
                "gt_cx",
                "gt_cy",
                "error",
            ]
        )
        for frame_idx in range(num_frames):
            for gt_idx in range(num_objects):
                gt_pt = gt_positions[frame_idx, gt_idx]
                tid_match = None
                pred_pt = None
                for tid, g_idx in track_to_gt.items():
                    if g_idx == gt_idx and not np.isnan(
                        track_positions[tid][frame_idx, 0]
                    ):
                        tid_match = tid
                        pred_pt = track_positions[tid][frame_idx]
                        break
                if pred_pt is None:
                    err = float("nan")
                    row_pred = (float("nan"), float("nan"))
                else:
                    err = float(np.linalg.norm(pred_pt - gt_pt))
                    row_pred = (float(pred_pt[0]), float(pred_pt[1]))
                writer_csv.writerow(
                    [
                        frame_idx,
                        gt_idx,
                        tid_match if tid_match is not None else "",
                        row_pred[0],
                        row_pred[1],
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
    overlay_renderer = OverlayRenderer(
        scale_bar_px=expected_minor_px,
        scale_bar_um=bac_short_axis_length_um,
        draw_history=True,
        history_length=num_frames,
        hide_history_after=num_frames,
    )
    history_points: dict[int, list[tuple[int, float, float]]] = {}
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

        # 予測を緑（対応ありのみ）
        for tid, gt_idx in track_to_gt.items():
            pred_pt = track_positions[tid][idx]
            if np.isnan(pred_pt[0]):
                continue
            history = history_points.setdefault(tid, [])
            history.append((idx, float(pred_pt[0]), float(pred_pt[1])))
            if len(history) > num_frames:
                history = history[-num_frames:]
                history_points[tid] = history
            # 軌跡を描画
            if len(history) >= 2:
                pts = [
                    (int(round(x)), int(round(y)))
                    for _, x, y in sorted(history, key=lambda t: t[0])
                    if _ <= idx
                ]
                for i in range(1, len(pts)):
                    cv2.line(frame, pts[i - 1], pts[i], (0, 200, 200), 1)
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

        overlay_renderer.draw_scale_bar(frame)

        overlay_writer.write(frame)

    overlay_writer.release()
