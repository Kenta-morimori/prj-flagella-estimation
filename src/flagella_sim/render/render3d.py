"""3D軌跡の簡易可視化を行う。"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

import cv2
import numpy as np

from flagella_sim.sim.core import SimulationState


def _normalize_coords(points: np.ndarray, img_size: int) -> np.ndarray:
    """XY平面を画像座標へスケーリングして中央に配置する。"""
    if points.size == 0:
        return points
    # 正規化のためのスケールを求める（最大半径が画像サイズの 0.4 以内）
    max_span = max(np.ptp(points[:, 0]), np.ptp(points[:, 1]), 1e-6)
    scale = (img_size * 0.4) / max_span
    pts = points[:, :2] * scale
    pts[:, 0] += img_size / 2
    pts[:, 1] += img_size / 2
    return pts


def save_swim_movie(states: Iterable[SimulationState], out_dir: Path) -> None:
    """3D軌跡の可視化動画/最終フレームを保存する。"""

    states_list: List[SimulationState] = list(states)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not states_list:
        (out_dir / "swim3d_final.png").write_text("no states", encoding="utf-8")
        return

    coords = np.array([st.position_um for st in states_list], dtype=float)
    xy = _normalize_coords(coords.copy(), img_size=512)

    frames: List[np.ndarray] = []
    for i in range(len(states_list)):
        img = np.zeros((512, 512, 3), dtype=np.uint8)
        pts = xy[: i + 1]
        for j in range(1, pts.shape[0]):
            p0 = tuple(np.round(pts[j - 1]).astype(int))
            p1 = tuple(np.round(pts[j]).astype(int))
            cv2.line(img, p0, p1, (0, 255, 180), 2, cv2.LINE_AA)
        if pts.shape[0] > 0:
            p_last = tuple(np.round(pts[-1]).astype(int))
            cv2.circle(img, p_last, 4, (255, 200, 50), -1, cv2.LINE_AA)
        frames.append(img)

    final_path = out_dir / "swim3d_final.png"
    cv2.imwrite(str(final_path), frames[-1])

    video_path = out_dir / "swim3d.mp4"
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    writer = cv2.VideoWriter(str(video_path), fourcc, 15.0, (512, 512))
    for frame in frames:
        writer.write(frame)
    writer.release()
