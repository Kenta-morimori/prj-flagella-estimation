"""3D軌跡の2D orthographic 投影を描画する。"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable, List

import cv2
import numpy as np

from flagella_sim.sim.core import SimulationState, _rotate_vec
from flagella_sim.sim.params import SimulationConfig


def heading_from_quat(q: tuple[float, float, float, float]) -> float:
    """クォータニオンからXY平面上のヘディング角度[rad]を得る。"""
    v = _rotate_vec(np.array(q, dtype=float), np.array([1.0, 0.0, 0.0]))
    return math.atan2(v[1], v[0])


def _draw_flagella(
    img: np.ndarray,
    center_px: tuple[int, int],
    heading: float,
    cfg: SimulationConfig,
    rng: np.random.Generator,
) -> None:
    """デバッグ用に太線のべん毛を描画する（白背景向け配色）。"""
    n = max(0, cfg.flagella.n_flagella)
    if n == 0:
        return
    length_px = min(
        int(cfg.flagella.length_um / cfg.render.pixel_size_um),
        img.shape[0],
    )
    base_angles = np.linspace(0, 2 * math.pi, n, endpoint=False)
    base_angles += rng.normal(0.0, 0.1, size=n)
    base_angles += heading
    thickness = int(max(1, round(cfg.render.flagella_linewidth_px)))
    cx, cy = center_px
    for i, ang in enumerate(base_angles):
        dx = int(math.cos(ang) * length_px * 0.6)
        dy = int(math.sin(ang) * length_px * 0.6)
        hsv = np.uint8([[[int((i * 40) % 180), 200, 230]]])
        color = tuple(int(c) for c in cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)[0, 0])
        cv2.line(img, (cx, cy), (cx + dx, cy + dy), color, thickness, cv2.LINE_AA)
        cv2.putText(
            img,
            f"F{i}",
            (cx + dx, cy + dy),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.4,
            color,
            1,
            cv2.LINE_AA,
        )


def project_states(
    states: Iterable[SimulationState], cfg: SimulationConfig, out_dir: Path
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

    rng = np.random.default_rng(cfg.seed.global_seed)
    frames: List[np.ndarray] = []

    for idx, st in enumerate(states_list):
        img = np.full((img_size, img_size, 3), 255, dtype=np.uint8)  # 白背景

        cx = int(img_size // 2 + st.position_um[0] * px_per_um)
        cy = int(img_size // 2 + st.position_um[1] * px_per_um)
        heading = heading_from_quat(st.quaternion)

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

        if cfg.render.render_flagella:
            _draw_flagella(img, (cx, cy), heading, cfg, rng)

        cv2.circle(img, (cx, cy), max(1, thickness), (0, 0, 0), -1, cv2.LINE_AA)

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
