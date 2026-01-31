"""3D軌跡をMatplotlibの3Dグリッドで可視化する。"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List

import cv2
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from flagella_sim.sim.core import SimulationState, _quat_to_rotmat
from flagella_sim.sim.flagella_geometry import FlagellaRig
from flagella_sim.sim.params import SimulationConfig


def _flagella_colors(n: int) -> list[tuple[int, int, int]]:
    if n <= 0:
        return []
    colors = []
    for i in range(n):
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
    """Matplotlibで3Dグリッド付きの可視化を生成する。"""

    states_list: List[SimulationState] = list(states)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not states_list:
        (out_dir / "swim3d_final.png").write_text("no states", encoding="utf-8")
        return

    coords = np.array([st.position_um for st in states_list], dtype=float)
    colors = _flagella_colors(cfg.flagella.n_flagella)

    # 初期位置を原点に寄せ、FOVは image_size_px * pixel_size_um の立方体
    render_offset = -coords[0]
    fov_um = cfg.render.image_size_px * cfg.render.pixel_size_um
    half = fov_um / 2.0

    frames: List[np.ndarray] = []
    for i, st in enumerate(states_list):
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(111, projection="3d")
        ax.set_facecolor("white")

        # グリッドと立方体範囲
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)
        ax.set_zlim(-half, half)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("x [µm]")
        ax.set_ylabel("y [µm]")
        ax.set_zlabel("z [µm]")
        ax.grid(True, which="both")

        # 軌跡
        path = coords[: i + 1] + render_offset
        ax.plot(path[:, 0], path[:, 1], path[:, 2], color="#64b6b6", linewidth=2)
        ax.scatter(
            path[-1, 0], path[-1, 1], path[-1, 2], color="k", s=20, depthshade=False
        )

        # 菌体（中心軸のみ簡易描画）
        rot = _quat_to_rotmat(np.array(st.quaternion, dtype=float))
        half_len = cfg.body.length_total_um / 2.0
        body_axis = np.array([half_len, 0.0, 0.0])
        p1 = (rot @ body_axis) + np.array(st.position_um) + render_offset
        p2 = (rot @ (-body_axis)) + np.array(st.position_um) + render_offset
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            [p1[2], p2[2]],
            color=(0.4, 0.4, 0.4),
            linewidth=4,
        )
        ax.scatter(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            [p1[2], p2[2]],
            color=(0.2, 0.2, 0.2),
            s=40,
            depthshade=False,
        )

        # べん毛
        for idx_f, base_off in enumerate(rig.base_offsets_body):
            color = np.array(colors[idx_f % len(colors)]) / 255.0
            base_world = rot @ base_off + np.array(st.position_um) + render_offset
            helix_world = rig.helix_local @ rot.T + base_world
            ax.plot(
                helix_world[:, 0],
                helix_world[:, 1],
                helix_world[:, 2],
                color=color,
                linewidth=2,
            )
            ax.text(
                helix_world[-1, 0],
                helix_world[-1, 1],
                helix_world[-1, 2],
                f"F{idx_f}",
                color=color,
                fontsize=8,
            )

        # 凡例
        handles = [
            plt.Line2D([0], [0], color=np.array(c) / 255.0, lw=2, label=f"F{i}")
            for i, c in enumerate(colors)
        ]
        if handles:
            ax.legend(handles=handles, loc="upper right", fontsize=8)

        fig.tight_layout()
        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())
        img_bgr = cv2.cvtColor(buf, cv2.COLOR_RGBA2BGR)
        frames.append(img_bgr)
        plt.close(fig)

    final_path = out_dir / "swim3d_final.png"
    cv2.imwrite(str(final_path), frames[-1])

    video_path = out_dir / "swim3d.mp4"
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    writer = cv2.VideoWriter(
        str(video_path),
        fourcc,
        cfg.time.fps_out,
        (frames[0].shape[1], frames[0].shape[0]),
    )
    for frame in frames:
        writer.write(frame)
    writer.release()
