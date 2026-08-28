"""3Dの実ビーズ連結表示を生成する。"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import cv2
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg

from sim_swim.model.types import PolymorphState
from sim_swim.sim.core import SimulationState
from sim_swim.sim.flagella_geometry import FlagellaRig
from sim_swim.sim.helix_axis import estimate_flag_helix_axis
from sim_swim.sim.params import SimulationConfig
from sim_swim.render.video_writer import (
    VideoRenderResult,
    VideoWriterSelection,
    open_mp4_writer,
)


def _flagella_colors(n: int) -> list[tuple[float, float, float]]:
    if n <= 0:
        return []
    colors: list[tuple[float, float, float]] = []
    for i in range(n):
        hsv = np.array([[[int((i * 40) % 180), 200, 230]]], dtype=np.uint8)
        bgr = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)[0, 0]
        colors.append(
            (
                float(bgr[2]) / 255.0,
                float(bgr[1]) / 255.0,
                float(bgr[0]) / 255.0,
            )
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


def _run_tumble_label(st: SimulationState, cfg: SimulationConfig) -> str:
    if not cfg.motor.enable_switching:
        return "RUN"

    reverse_ids = np.asarray(st.reverse_flagella, dtype=int)
    if reverse_ids.size == 0:
        return "RUN"

    states = np.asarray(st.flag_states, dtype=int)
    # If reverse_flagella contains out-of-bounds indices, treat it as a safe
    # RUN.
    # Avoid falling back to a time-based heuristic that may mislabel frames.
    max_id = int(np.max(reverse_ids))
    if max_id >= states.size:
        return "RUN"

    rev_states = states[reverse_ids]
    is_tumble = np.any(
        (rev_states == int(PolymorphState.SEMICOILED))
        | (rev_states == int(PolymorphState.CURLY1))
    )
    return "TUMBLE" if bool(is_tumble) else "RUN"


def format_simulation_time_label(t_s: float, cfg: SimulationConfig) -> str:
    """Return the required time label for every 3D render path."""

    tau_s = max(float(cfg.tau_s), 1.0e-30)
    dt_internal_s = max(float(cfg.dt_star) * tau_s, 1.0e-30)
    return (
        f"t = {float(t_s) / tau_s:.3f} τ "
        f"({float(t_s):.6f} s, {int(round(float(t_s) / dt_internal_s)):,} steps)"
    )


def _frame_status_lines(
    st: SimulationState,
    cfg: SimulationConfig,
    *,
    extra_lines: Iterable[str] = (),
) -> list[str]:
    lines = [_run_tumble_label(st, cfg)]
    lines.append(format_simulation_time_label(st.t, cfg))
    lines.append(f"motor_torque_Nm = {cfg.motor_torque_Nm:.3e}")
    lines.append(f"follow_camera_3d = {cfg.render.follow_camera_3d}")
    lines.extend(str(line) for line in extra_lines if str(line).strip())
    return lines


def _plot_segments_3d(
    ax: plt.Axes,
    beads: np.ndarray,
    edges: np.ndarray,
    *,
    color: tuple[float, float, float] | str,
    linewidth: float,
) -> None:
    bead_count = int(beads.shape[0])
    for i, j in edges:
        if int(i) >= bead_count or int(j) >= bead_count:
            continue
        p = beads[int(i)]
        q = beads[int(j)]
        ax.plot(
            [p[0], q[0]],
            [p[1], q[1]],
            [p[2], q[2]],
            color=color,
            linewidth=linewidth,
        )


def _hook_edges(hook_triplets: np.ndarray) -> np.ndarray:
    if hook_triplets.size == 0:
        return np.empty((0, 2), dtype=int)

    triplets = np.asarray(hook_triplets, dtype=int)
    first = triplets[:, [0, 1]]
    second = triplets[:, [1, 2]]
    return np.vstack((first, second))


def _resolve_view_range_um(cfg: SimulationConfig, rig: FlagellaRig) -> float:
    if cfg.render.view_range_um != 5.0:
        return max(cfg.render.view_range_um, 1e-6)

    default_view_range = 3.0 if len(rig.flagella_indices) == 0 else 5.0
    return max(default_view_range, 1e-6)


def _plot_flagella_helix_axis_3d(
    ax: plt.Axes,
    beads_um: np.ndarray,
    flag_indices: np.ndarray,
    flag_id: int,
    color: tuple[float, float, float],
) -> None:
    axis = estimate_flag_helix_axis(beads_um, flag_indices, flag_id)
    if axis.degenerate or not np.isfinite(axis.line_start).all():
        return
    ax.plot(
        [axis.line_start[0], axis.line_end[0]],
        [axis.line_start[1], axis.line_end[1]],
        [axis.line_start[2], axis.line_end[2]],
        color=color,
        linewidth=1.1,
        linestyle="--",
        alpha=0.85,
    )


def plot_swim_frame_3d(
    ax: plt.Axes,
    st: SimulationState,
    cfg: SimulationConfig,
    rig: FlagellaRig,
    *,
    title: str | None = None,
    extra_status_lines: Iterable[str] = (),
    hide_ticks: bool = False,
    show_status: bool = True,
    show_legend: bool = True,
    body_linewidth: float = 1.6,
    flagella_linewidth: float = 2.0,
    hook_linewidth: float = 2.6,
    body_marker_size: float = 8.0,
    flagella_marker_size: float = 6.0,
    flow_vectors: tuple[np.ndarray, np.ndarray] | None = None,
) -> None:
    """Plot one 3D swim frame using the same visual contract as swim movies."""

    ax.set_facecolor("white")
    beads = st.bead_positions_um
    center = np.array(
        st.position_um if cfg.render.follow_camera_3d else (0.0, 0.0, 0.0),
        dtype=float,
    )
    view_range = _resolve_view_range_um(cfg, rig)
    ax.set_xlim(center[0] - view_range, center[0] + view_range)
    ax.set_ylim(center[1] - view_range, center[1] + view_range)
    ax.set_zlim(center[2] - view_range, center[2] + view_range)
    ax.set_box_aspect((1, 1, 1))
    if hide_ticks:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
    else:
        ax.set_xlabel("x [um]")
        ax.set_ylabel("y [um]")
        ax.set_zlabel("z [um]")
    ax.grid(True)
    if title:
        ax.set_title(title, fontsize=9, pad=6)

    if flow_vectors is not None:
        flow_points_um, flow_velocity_um_s = flow_vectors
        points = np.asarray(flow_points_um, dtype=float)
        vectors = np.asarray(flow_velocity_um_s, dtype=float)
        if points.shape != vectors.shape or points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("flow_vectors must be a pair of (M, 3) arrays")
        ax.quiver(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            vectors[:, 0],
            vectors[:, 1],
            vectors[:, 2],
            length=0.35,
            normalize=True,
            color="tab:blue",
            linewidth=0.7,
            alpha=0.7,
        )

    _plot_segments_3d(
        ax,
        beads,
        rig.body_spring_edges,
        color=(0.35, 0.35, 0.35),
        linewidth=body_linewidth,
    )

    body_indices = np.concatenate(rig.body_layer_indices)
    body_indices = body_indices[body_indices < beads.shape[0]]
    if body_indices.size == 0:
        body_indices = np.arange(beads.shape[0], dtype=int)
    body_pts = beads[body_indices]
    if body_pts.size > 0:
        ax.scatter(
            body_pts[:, 0],
            body_pts[:, 1],
            body_pts[:, 2],
            color="k",
            s=body_marker_size,
            depthshade=False,
        )

    handles = []
    colors = _flagella_colors(len(rig.flagella_indices))
    if cfg.render.render_flagella:
        for f_id, idxs in enumerate(rig.flagella_indices):
            idxs = idxs[idxs < beads.shape[0]]
            if idxs.size == 0:
                continue
            color = colors[f_id % len(colors)] if colors else (0.1, 0.4, 0.7)
            pts = beads[idxs]
            (line,) = ax.plot(
                pts[:, 0],
                pts[:, 1],
                pts[:, 2],
                color=color,
                linewidth=flagella_linewidth,
            )
            handles.append((line, f"F{f_id}"))
            ax.scatter(
                pts[:, 0],
                pts[:, 1],
                pts[:, 2],
                color=[color],
                s=flagella_marker_size,
                depthshade=False,
            )
            if cfg.render.label_flagella:
                end = pts[-1]
                ax.text(
                    end[0],
                    end[1],
                    end[2],
                    f"F{f_id}",
                    color=color,
                    fontsize=8,
                )
            if cfg.render.show_flagella_helix_axis_3d:
                _plot_flagella_helix_axis_3d(ax, beads, idxs, f_id, color)

    if rig.hook_triplets.size > 0:
        hook_edges = _hook_edges(rig.hook_triplets)
        _plot_segments_3d(
            ax,
            beads,
            hook_edges,
            color=(1.0, 0.85, 0.05),
            linewidth=hook_linewidth,
        )

    if handles and show_legend:
        ax.legend(
            [h[0] for h in handles],
            [h[1] for h in handles],
            loc="upper right",
            fontsize=8,
        )

    if show_status:
        ax.text2D(
            0.08,
            0.96,
            "\n".join(_frame_status_lines(st, cfg, extra_lines=extra_status_lines)),
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
        )


def render_swim_frame_3d(
    st: SimulationState,
    cfg: SimulationConfig,
    rig: FlagellaRig,
    *,
    width_px: int = 500,
    height_px: int = 500,
    title: str | None = None,
    extra_status_lines: Iterable[str] = (),
    hide_ticks: bool = False,
) -> np.ndarray:
    """Render one 3D swim frame to a BGR image."""

    fig = plt.figure(figsize=(width_px / 100.0, height_px / 100.0), dpi=100)
    ax = fig.add_subplot(111, projection="3d")
    plot_swim_frame_3d(
        ax,
        st,
        cfg,
        rig,
        title=title,
        extra_status_lines=extra_status_lines,
        hide_ticks=hide_ticks,
    )
    fig.tight_layout()
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())
    frame = cv2.cvtColor(buf, cv2.COLOR_RGBA2BGR)
    plt.close(fig)
    if frame.shape[:2] != (height_px, width_px):
        frame = cv2.resize(frame, (width_px, height_px), interpolation=cv2.INTER_AREA)
    return frame


def save_swim_movie(
    states: Iterable[SimulationState],
    cfg: SimulationConfig,
    rig: FlagellaRig,
    out_dir: Path,
) -> VideoRenderResult | None:
    """3Dの連結ビーズ可視化をPNG連番と動画で保存する。"""

    states_list = list(states)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not states_list:
        (out_dir / "swim3d_final.png").write_text("no states", encoding="utf-8")
        return None

    render_states = _select_frames(
        states_list,
        out_all_steps_3d=cfg.output_sampling.out_all_steps_3d,
        fps_hint=cfg.output_sampling.fps_out_3d,
    )

    frames_dir = out_dir / "frames_3d"
    if cfg.render.save_frames_3d:
        frames_dir.mkdir(parents=True, exist_ok=True)

    movie_path = out_dir / "swim3d.mp4"
    writer: cv2.VideoWriter | None = None
    writer_selection: VideoWriterSelection | None = None
    last_frame: np.ndarray | None = None
    frame_count = 0

    if cfg.output_sampling.out_all_steps_3d:
        fps_3d = min(60.0, max(1.0, 1.0 / max(cfg.output_dt_s, 1e-9)))
    else:
        fps_3d = max(1.0, float(cfg.output_sampling.fps_out_3d))

    for idx, st in enumerate(render_states):
        frame = render_swim_frame_3d(
            st,
            cfg,
            rig,
            width_px=500,
            height_px=500,
        )

        if cfg.render.save_frames_3d:
            cv2.imwrite(str(frames_dir / f"frame_{idx:06d}.png"), frame)

        if writer is None:
            writer_selection = open_mp4_writer(
                movie_path,
                fps=fps_3d,
                frame_size=(frame.shape[1], frame.shape[0]),
            )
            writer = writer_selection.writer
        writer.write(frame)
        frame_count += 1
        last_frame = frame

    if writer is not None:
        writer.release()

    if last_frame is not None:
        cv2.imwrite(str(out_dir / "swim3d_final.png"), last_frame)

    if writer_selection is None or last_frame is None:
        return None
    return VideoRenderResult(
        path=str(movie_path),
        selected_codec=writer_selection.selected_codec,
        attempted_codecs=writer_selection.attempted_codecs,
        fps=fps_3d,
        frame_size=(last_frame.shape[1], last_frame.shape[0]),
        frame_count=frame_count,
    )
