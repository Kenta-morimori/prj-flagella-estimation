"""3D軌跡の2D orthographic 投影を描画する。"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, List

import cv2
import numpy as np

from flagella_sim.sim.core import SimulationState, _rotate_vec, _quat_to_rotmat
from flagella_sim.sim.flagella_geometry import FlagellaRig
from flagella_sim.sim.params import SimulationConfig


def heading_from_quat(q: tuple[float, float, float, float]) -> float:
    """クォータニオンからXY平面上のヘディング角度[rad]を得る。"""
    v = _rotate_vec(np.array(q, dtype=float), np.array([1.0, 0.0, 0.0]))
    return math.atan2(v[1], v[0])


def _legend(img: np.ndarray, colors: list[tuple[int, int, int]]) -> None:
    x0, y0 = 10, 10
    for i, c in enumerate(colors):
        cv2.rectangle(img, (x0, y0 + i * 18), (x0 + 14, y0 + 14 + i * 18), c, -1)
        cv2.putText(
            img,
            f"F{i}",
            (x0 + 20, y0 + 12 + i * 18),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.45,
            (0, 0, 0),
            1,
            cv2.LINE_AA,
        )


def project_states(
    states: Iterable[SimulationState],
    cfg: SimulationConfig,
    rig: FlagellaRig,
    out_dir: Path,
) -> None:
    """軌跡を2Dへ投影しPNG連番とmp4を出力する。"""

    states_list: List[SimulationState] = list(states)
    if not states_list:
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    frames_dir = out_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    img_size = cfg.render.image_size_px
    px_per_um = 1.0 / cfg.render.pixel_size_um
    body_major_px = int(cfg.body.length_total_um * px_per_um)
    body_minor_px = int(cfg.body.diameter_um * px_per_um)
    thickness = max(1, int(round(body_minor_px / 6)))

    colors = [
        tuple(
            int(c)
            for c in cv2.cvtColor(
                np.uint8([[[int((i * 40) % 180), 200, 230]]]), cv2.COLOR_HSV2BGR
            )[0, 0]
        )
        for i in range(max(1, cfg.flagella.n_flagella))
    ]
    frames: List[np.ndarray] = []

    for idx, st in enumerate(states_list):
        img = np.full((img_size, img_size, 3), 255, dtype=np.uint8)  # 白背景

        cx = int(img_size // 2 + st.position_um[0] * px_per_um)
        cy = int(img_size // 2 + st.position_um[1] * px_per_um)
        heading = heading_from_quat(st.quaternion)
        rot = _quat_to_rotmat(np.array(st.quaternion, dtype=float))

        # べん毛ヘリックス
        if cfg.render.render_flagella and rig.base_offsets_body.size > 0:
            for idx_f, base_off in enumerate(rig.base_offsets_body):
                color = colors[idx_f % len(colors)]
                base_world = rot @ base_off + np.array(st.position_um)
                helix_world = rig.helix_local @ rot.T + base_world
                pts = helix_world[:, :2] * px_per_um + np.array(
                    [img_size / 2, img_size / 2]
                )
                pts_int = np.round(pts).astype(int)
                cv2.polylines(img, [pts_int], False, color, thickness, cv2.LINE_AA)
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

        axes = (max(1, body_major_px // 2), max(1, body_minor_px // 2))
        cv2.ellipse(
            img,
            (cx, cy),
            axes,
            math.degrees(heading),
            0,
            360,
            (80, 80, 80),
            thickness,
            cv2.LINE_AA,
        )

        cv2.circle(img, (cx, cy), max(1, thickness), (0, 0, 0), -1, cv2.LINE_AA)
        if cfg.render.render_flagella:
            _legend(img, colors)

        frame_path = frames_dir / f"frame_{idx:06d}.png"
        cv2.imwrite(str(frame_path), img)
        frames.append(img)

    # mp4書き出し
    video_path = out_dir / "projection.mp4"
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    writer = cv2.VideoWriter(
        str(video_path), fourcc, cfg.time.fps_out, (img_size, img_size)
    )
    for frame in frames:
        writer.write(frame)
    writer.release()
