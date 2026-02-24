"""3Dの実ビーズ連結表示を生成する。"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import cv2
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sim_swim.sim.core import SimulationState
from sim_swim.sim.flagella_geometry import FlagellaRig
from sim_swim.sim.params import SimulationConfig


def _flagella_colors(n: int) -> list[tuple[float, float, float]]:
    if n <= 0:
        return []
    colors: list[tuple[float, float, float]] = []
    for i in range(n):
        hsv = np.uint8([[[int((i * 40) % 180), 200, 230]]])
        bgr = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)[0, 0]
        colors.append(
            (float(bgr[2]) / 255.0, float(bgr[1]) / 255.0, float(bgr[0]) / 255.0)
        )
    return colors


def _select_frames(
    states: list[SimulationState], out_all_steps_3d: bool, fps_hint: float
) -> list[SimulationState]:
    if not states:
        return []
    if out_all_steps_3d:
        return states

    interval = 1.0 / max(fps_hint, 1e-9)
    selected: list[SimulationState] = []
    next_t = states[0].t
    for st in states:
        if st.t + 1e-12 >= next_t:
            selected.append(st)
            next_t += interval
    if selected and selected[-1] is not states[-1]:
        selected.append(states[-1])
    return selected


def save_swim_movie(
    states: Iterable[SimulationState],
    cfg: SimulationConfig,
    rig: FlagellaRig,
    out_dir: Path,
) -> None:
    """3Dの連結ビーズ可視化をPNG連番と動画で保存する。"""

    states_list = list(states)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not states_list:
        (out_dir / "swim3d_final.png").write_text("no states", encoding="utf-8")
        return

    render_states = _select_frames(
        states_list,
        out_all_steps_3d=cfg.output_sampling.out_all_steps_3d,
        fps_hint=cfg.output_sampling.fps_out_2d,
    )

    frames_dir = out_dir / "frames_3d"
    if cfg.render.save_frames_3d:
        frames_dir.mkdir(parents=True, exist_ok=True)

    colors = _flagella_colors(len(rig.flagella_indices))
    view_range = max(cfg.render.view_range_um, 1e-6)

    movie_path = out_dir / "swim3d.mp4"
    writer: cv2.VideoWriter | None = None
    last_frame: np.ndarray | None = None

    fps_3d = min(60.0, max(1.0, 1.0 / max(cfg.dt_s, 1e-9)))

    for idx, st in enumerate(render_states):
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(111, projection="3d")
        ax.set_facecolor("white")

        beads = st.bead_positions_um
        center = np.array(
            st.position_um if cfg.render.follow_camera_3d else (0.0, 0.0, 0.0),
            dtype=float,
        )

        ax.set_xlim(center[0] - view_range, center[0] + view_range)
        ax.set_ylim(center[1] - view_range, center[1] + view_range)
        ax.set_zlim(center[2] - view_range, center[2] + view_range)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlabel("x [um]")
        ax.set_ylabel("y [um]")
        ax.set_zlabel("z [um]")
        ax.grid(True)

        for i, j in rig.body_ring_edges:
            p = beads[int(i)]
            q = beads[int(j)]
            ax.plot(
                [p[0], q[0]],
                [p[1], q[1]],
                [p[2], q[2]],
                color=(0.35, 0.35, 0.35),
                linewidth=1.8,
            )

        for i, j in rig.body_vertical_edges:
            p = beads[int(i)]
            q = beads[int(j)]
            ax.plot(
                [p[0], q[0]],
                [p[1], q[1]],
                [p[2], q[2]],
                color=(0.35, 0.35, 0.35),
                linewidth=1.8,
            )

        body_pts = beads[np.concatenate(rig.body_layer_indices)]
        ax.scatter(
            body_pts[:, 0],
            body_pts[:, 1],
            body_pts[:, 2],
            color="k",
            s=8,
            depthshade=False,
        )

        handles = []
        if cfg.render.render_flagella:
            for f_id, idxs in enumerate(rig.flagella_indices):
                color = colors[f_id % len(colors)] if colors else (0.1, 0.4, 0.7)
                pts = beads[idxs]
                (line,) = ax.plot(
                    pts[:, 0], pts[:, 1], pts[:, 2], color=color, linewidth=2.0
                )
                handles.append((line, f"F{f_id}"))
                ax.scatter(
                    pts[:, 0],
                    pts[:, 1],
                    pts[:, 2],
                    color=[color],
                    s=6,
                    depthshade=False,
                )
                if cfg.render.label_flagella:
                    end = pts[-1]
                    ax.text(end[0], end[1], end[2], f"F{f_id}", color=color, fontsize=8)

        if handles:
            ax.legend(
                [h[0] for h in handles],
                [h[1] for h in handles],
                loc="upper right",
                fontsize=8,
            )

        if cfg.render.timestamp_3d:
            label = cfg.render.timestamp_fmt.format(t=st.t)
            ax.text2D(0.02, 0.96, label, transform=ax.transAxes)

        fig.tight_layout()
        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())
        frame = cv2.cvtColor(buf, cv2.COLOR_RGBA2BGR)
        plt.close(fig)

        if cfg.render.save_frames_3d:
            cv2.imwrite(str(frames_dir / f"frame_{idx:06d}.png"), frame)

        if writer is None:
            writer = cv2.VideoWriter(
                str(movie_path),
                cv2.VideoWriter_fourcc(*"mp4v"),
                fps_3d,
                (frame.shape[1], frame.shape[0]),
            )
        writer.write(frame)
        last_frame = frame

    if writer is not None:
        writer.release()

    if last_frame is not None:
        cv2.imwrite(str(out_dir / "swim3d_final.png"), last_frame)
