"""2D固定カメラ投影を描画する。"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import cv2
import numpy as np

from sim_swim.sim.core import SimulationState
from sim_swim.sim.flagella_geometry import FlagellaRig
from sim_swim.sim.params import SimulationConfig


def _flagella_colors(n: int) -> list[tuple[int, int, int]]:
    if n <= 0:
        return []
    colors = []
    for i in range(n):
        hsv = np.uint8([[[int((i * 40) % 180), 200, 230]]])
        bgr = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)[0, 0]
        colors.append(tuple(int(c) for c in bgr))
    return colors


def _select_2d_frames(
    states: list[SimulationState], fps: float
) -> list[SimulationState]:
    if not states:
        return []
    interval = 1.0 / max(fps, 1e-9)
    selected: list[SimulationState] = []
    next_t = states[0].t
    for st in states:
        if st.t + 1e-12 >= next_t:
            selected.append(st)
            next_t += interval
    if selected and selected[-1] is not states[-1]:
        selected.append(states[-1])
    return selected


def project_states(
    states: Iterable[SimulationState],
    cfg: SimulationConfig,
    rig: FlagellaRig,
    out_dir: Path,
) -> None:
    """固定カメラ2D投影を25Hzで出力する。"""

    states_list = list(states)
    if not states_list:
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    frames_dir = out_dir / "frames"
    if cfg.render.save_frames_2d:
        frames_dir.mkdir(parents=True, exist_ok=True)

    fps_2d = cfg.output_sampling.fps_out_2d
    sampled = _select_2d_frames(states_list, fps_2d)

    img_size = cfg.render.image_size_px
    px_per_um = 1.0 / max(cfg.render.pixel_size_um, 1e-9)
    line_w = max(1, int(round(cfg.render.flagella_linewidth_px)))

    body_color = (90, 90, 90)
    colors = _flagella_colors(len(rig.flagella_indices))

    video_path = out_dir / "projection.mp4"
    writer = cv2.VideoWriter(
        str(video_path),
        cv2.VideoWriter_fourcc(*"mp4v"),
        fps_2d,
        (img_size, img_size),
    )

    for idx, st in enumerate(sampled):
        img = np.full((img_size, img_size, 3), 255, dtype=np.uint8)
        beads = st.bead_positions_um

        if cfg.render.follow_camera_2d:
            cam_center = np.array(st.position_um[:2], dtype=float)
        else:
            cam_center = np.zeros(2, dtype=float)

        def to_px(p: np.ndarray) -> tuple[int, int]:
            x = int(round((p[0] - cam_center[0]) * px_per_um + img_size / 2.0))
            y = int(round((p[1] - cam_center[1]) * px_per_um + img_size / 2.0))
            return x, y

        for i, j in rig.body_ring_edges:
            p = to_px(beads[int(i)])
            q = to_px(beads[int(j)])
            cv2.line(img, p, q, body_color, line_w, cv2.LINE_AA)

        for i, j in rig.body_vertical_edges:
            p = to_px(beads[int(i)])
            q = to_px(beads[int(j)])
            cv2.line(img, p, q, body_color, line_w, cv2.LINE_AA)

        if cfg.render.render_flagella_2d:
            for f_id, idxs in enumerate(rig.flagella_indices):
                color = colors[f_id % len(colors)] if colors else (30, 120, 220)
                pts = [to_px(beads[int(i)]) for i in idxs]
                for p, q in zip(pts[:-1], pts[1:]):
                    cv2.line(img, p, q, color, line_w, cv2.LINE_AA)
                if cfg.render.label_flagella and pts:
                    cv2.putText(
                        img,
                        f"F{f_id}",
                        pts[-1],
                        cv2.FONT_HERSHEY_SIMPLEX,
                        0.4,
                        color,
                        1,
                        cv2.LINE_AA,
                    )

        if cfg.render.save_frames_2d:
            cv2.imwrite(str(frames_dir / f"frame_{idx:06d}.png"), img)
        writer.write(img)

    writer.release()
