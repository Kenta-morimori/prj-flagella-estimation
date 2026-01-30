"""3D軌跡の簡易可視化を行う。"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

import cv2
import numpy as np

from flagella_sim.sim.core import SimulationState, _quat_to_rotmat
from flagella_sim.sim.flagella_geometry import FlagellaRig
from flagella_sim.sim.params import SimulationConfig


def _normalize_coords(points: np.ndarray, img_size: int) -> np.ndarray:
    """XY平面を画像座標へスケーリングして中央に配置する。"""
    if points.size == 0:
        return points
    max_span = max(np.ptp(points[:, 0]), np.ptp(points[:, 1]), 1e-6)
    scale = (img_size * 0.4) / max_span
    pts = points[:, :2] * scale
    pts[:, 0] += img_size / 2
    pts[:, 1] += img_size / 2
    return pts


def _flagella_colors(n: int) -> list[tuple[int, int, int]]:
    colors = []
    for i in range(max(1, n)):
        hsv = np.uint8([[[int((i * 40) % 180), 200, 230]]])
        bgr = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)[0, 0]
        colors.append(tuple(int(c) for c in bgr))
    return colors


def save_swim_movie(
    states: Iterable[SimulationState],
    cfg: SimulationConfig,
    rig: FlagellaRig,
    out_dir: Path,
) -> None:
    """3D軌跡の可視化動画/最終フレームを保存する（白背景＋色分けべん毛）。"""

    states_list: List[SimulationState] = list(states)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not states_list:
        (out_dir / "swim3d_final.png").write_text("no states", encoding="utf-8")
        return

    img_size = 720
    coords = np.array([st.position_um for st in states_list], dtype=float)
    xy = _normalize_coords(coords.copy(), img_size=img_size)
    colors = _flagella_colors(cfg.flagella.n_flagella)

    frames: List[np.ndarray] = []
    for i, st in enumerate(states_list):
        img = np.full((img_size, img_size, 3), 255, dtype=np.uint8)
        # 軌跡
        pts = xy[: i + 1]
        for j in range(1, pts.shape[0]):
            p0 = tuple(np.round(pts[j - 1]).astype(int))
            p1 = tuple(np.round(pts[j]).astype(int))
            cv2.line(img, p0, p1, (100, 180, 180), 2, cv2.LINE_AA)

        if pts.shape[0] > 0:
            p_last = tuple(np.round(pts[-1]).astype(int))
            cv2.circle(img, p_last, 6, (0, 0, 0), -1, cv2.LINE_AA)

        # べん毛ヘリックスをXYに投影
        rot = _quat_to_rotmat(np.array(st.quaternion, dtype=float))
        for idx_f, base_off in enumerate(rig.base_offsets_body):
            color = colors[idx_f % len(colors)]
            base_world = rot @ base_off + np.array(st.position_um)
            helix_world = rig.helix_local @ rot.T + base_world
            pts = helix_world[:, :2]
            pts = _normalize_coords(pts, img_size=img_size)
            pts_int = np.round(pts).astype(int)
            cv2.polylines(img, [pts_int], False, color, 2, cv2.LINE_AA)
            cv2.putText(
                img,
                f"F{idx_f}",
                tuple(pts_int[-1]),
                cv2.FONT_HERSHEY_SIMPLEX,
                0.4,
                color,
                1,
                cv2.LINE_AA,
            )

        frames.append(img)

    final_path = out_dir / "swim3d_final.png"
    cv2.imwrite(str(final_path), frames[-1])

    video_path = out_dir / "swim3d.mp4"
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    writer = cv2.VideoWriter(str(video_path), fourcc, 30.0, (img_size, img_size))
    for frame in frames:
        writer.write(frame)
    writer.release()
