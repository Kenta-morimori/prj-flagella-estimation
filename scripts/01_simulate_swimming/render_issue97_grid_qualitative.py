#!/usr/bin/env python3
"""Render a 3x3 qualitative comparison movie for Issue #97 torque conditions."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
import logging
from pathlib import Path
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo

import cv2
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.pyplot as plt
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.sweeps import hook_overstretch
from sim_swim.render.render3d import (
    _flagella_colors,
    _hook_edges,
    _plot_flagella_helix_axis_3d,
    _plot_segments_3d,
    _resolve_view_range_um,
    _run_tumble_label,
    _select_frames,
)
from sim_swim.render.video_writer import VideoRenderResult, open_mp4_writer
from sim_swim.sim.core import SimulationState, Simulator
from sim_swim.sim.params import SimulationConfig


def _default_output_dir() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / "phase2_97"
        / "stage_f_grid_qual_3d_fp3_ft1p5_torque2p0_dur0p6"
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
    )


def _load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _setup_logger(log_path: Path) -> logging.Logger:
    logger = logging.getLogger(f"issue97.grid_qual.{log_path.parent.name}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger


def _git_info() -> dict[str, Any]:
    def run(cmd: list[str]) -> str:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True).strip()

    try:
        return {
            "commit": run(["git", "rev-parse", "HEAD"]),
            "commit_short": run(["git", "rev-parse", "--short", "HEAD"]),
            "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
            "is_clean": run(["git", "status", "--porcelain"]) == "",
        }
    except Exception:
        return {
            "commit": "unknown",
            "commit_short": "unknown",
            "branch": "unknown",
            "is_clean": False,
        }


def _condition_matrix() -> list[hook_overstretch.Condition]:
    args = argparse.Namespace(
        mode="torque-profile-grid",
        force_distributions=[
            "root_torque_segment_couples",
            "root_torque_axis_projection",
            "root_torque_hybrid_couples",
        ],
        fixed_attach_first_spring_scale=1.0,
        fixed_body_axis_angle_scale=1.0,
        fixed_first_second_spring_scale=1.0,
        fixed_attach_frame_position_scale=3.0,
        fixed_attach_frame_tangent_scale=1.5,
        torque_distribution_profiles=["diffusive", "uniform", "basal_unloading"],
    )
    return list(hook_overstretch.build_conditions(args))


def _short_distribution_label(value: str) -> str:
    mapping = {
        "root_torque_segment_couples": "segment_couples",
        "root_torque_axis_projection": "axis_projection",
        "root_torque_hybrid_couples": "hybrid_couples",
    }
    return mapping.get(value, value)


def _plot_cell(
    ax: plt.Axes,
    *,
    st: SimulationState,
    cfg: SimulationConfig,
    rig,
    title: str,
    fail_label: str,
) -> None:
    ax.set_facecolor("white")
    beads = st.bead_positions_um
    view_range = _resolve_view_range_um(cfg, rig)
    center = np.array(st.position_um, dtype=float)
    ax.set_xlim(center[0] - view_range, center[0] + view_range)
    ax.set_ylim(center[1] - view_range, center[1] + view_range)
    ax.set_zlim(center[2] - view_range, center[2] + view_range)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.grid(True)
    ax.set_title(title, fontsize=9, pad=6)

    _plot_segments_3d(
        ax,
        beads,
        rig.body_spring_edges,
        color=(0.35, 0.35, 0.35),
        linewidth=1.3,
    )
    body_pts = beads[np.concatenate(rig.body_layer_indices)]
    ax.scatter(
        body_pts[:, 0],
        body_pts[:, 1],
        body_pts[:, 2],
        color="k",
        s=5,
        depthshade=False,
    )

    colors = _flagella_colors(len(rig.flagella_indices))
    if cfg.render.render_flagella:
        for f_id, idxs in enumerate(rig.flagella_indices):
            color = colors[f_id % len(colors)] if colors else (0.1, 0.4, 0.7)
            pts = beads[idxs]
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=color, linewidth=1.8)
            ax.scatter(
                pts[:, 0],
                pts[:, 1],
                pts[:, 2],
                color=[color],
                s=4,
                depthshade=False,
            )
            if cfg.render.show_flagella_helix_axis_3d:
                _plot_flagella_helix_axis_3d(ax, beads, idxs, f_id, color)

    if rig.hook_triplets.size > 0:
        _plot_segments_3d(
            ax,
            beads,
            _hook_edges(rig.hook_triplets),
            color=(1.0, 0.85, 0.05),
            linewidth=2.2,
        )

    lines = [
        _run_tumble_label(st, cfg),
        f"t = {st.t:.3f} s",
        fail_label,
    ]
    ax.text2D(
        0.02,
        0.96,
        "\n".join(lines),
        transform=ax.transAxes,
        va="top",
        fontsize=8,
    )


def _render_grid_movie(
    *,
    states_by_condition: list[list[SimulationState]],
    cfg_by_condition: list[SimulationConfig],
    rig_by_condition: list,
    condition_rows: list[dict[str, str]],
    out_dir: Path,
) -> VideoRenderResult:
    render_states_by_condition = [
        _select_frames(states, out_all_steps_3d=False, fps_hint=25.0)
        for states in states_by_condition
    ]
    frame_count = min(len(states) for states in render_states_by_condition)
    if frame_count <= 0:
        raise RuntimeError("No frames selected for grid render.")

    movie_path = out_dir / "grid_swim3d.mp4"
    writer = None
    writer_selection = None
    last_frame = None

    titles = [
        (
            f"{_short_distribution_label(row['force_distribution'])}\n"
            f"{row['torque_distribution_profile']}"
        )
        for row in condition_rows
    ]
    fail_labels = []
    for row in condition_rows:
        if row.get("final_shape_pass_nonbody", "") == "True":
            fail_labels.append("PASS")
        else:
            fail_t = row.get("first_fail_t_s", "")
            fail_c = row.get("first_fail_category_nonbody", "")
            fail_labels.append(f"FAIL {fail_c}@{fail_t[:6]}")

    for frame_idx in range(frame_count):
        fig = plt.figure(figsize=(12, 12))
        axes = [
            fig.add_subplot(3, 3, plot_index + 1, projection="3d")
            for plot_index in range(9)
        ]
        for idx, ax in enumerate(axes):
            _plot_cell(
                ax,
                st=render_states_by_condition[idx][frame_idx],
                cfg=cfg_by_condition[idx],
                rig=rig_by_condition[idx],
                title=titles[idx],
                fail_label=fail_labels[idx],
            )

        fig.tight_layout()
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        buf = np.asarray(canvas.buffer_rgba())
        frame = cv2.cvtColor(buf, cv2.COLOR_RGBA2BGR)
        plt.close(fig)

        if writer is None:
            writer_selection = open_mp4_writer(
                movie_path,
                fps=25.0,
                frame_size=(frame.shape[1], frame.shape[0]),
            )
            writer = writer_selection.writer
        writer.write(frame)
        last_frame = frame

    writer.release()
    if last_frame is not None:
        cv2.imwrite(str(out_dir / "grid_swim3d_final.png"), last_frame)
    if writer_selection is None or last_frame is None:
        raise RuntimeError("Failed to initialize video writer.")
    return VideoRenderResult(
        path=str(movie_path),
        selected_codec=writer_selection.selected_codec,
        attempted_codecs=writer_selection.attempted_codecs,
        fps=25.0,
        frame_size=(last_frame.shape[1], last_frame.shape[0]),
        frame_count=frame_count,
    )


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("conf/sim_swim.yaml"))
    parser.add_argument("--output-dir", type=Path, default=_default_output_dir())
    parser.add_argument("--duration-s", type=float, default=0.6)
    parser.add_argument("--dt-star", type=float, default=1.0e-4)
    parser.add_argument("--torque-nm", type=float, default=2.0e-20)
    parser.add_argument("--attach-seed", type=int, default=0)
    parser.add_argument("--phase-seed", type=int, default=0)
    parser.add_argument("--n-flagella", type=int, default=3)
    parser.add_argument("--progress-interval", type=int, default=5000)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    conditions = _condition_matrix()
    if args.dry_run:
        for condition in conditions:
            print(condition.condition_id)
        return

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=args.overwrite)
    logger = _setup_logger(output_dir / "run.log")
    logger.info("Starting issue97 grid qualitative render")

    base_cfg = _load_yaml(args.config)
    states_by_condition: list[list[SimulationState]] = []
    cfg_by_condition: list[SimulationConfig] = []
    rig_by_condition: list[Any] = []
    condition_rows: list[dict[str, str]] = []

    for index, condition in enumerate(conditions, start=1):
        logger.info("[%d/%d] %s", index, len(conditions), condition.condition_id)
        overrides = hook_overstretch._overrides_for_condition(args, condition)
        cfg = (
            SimulationConfig.from_dict(base_cfg)
            .with_overrides(overrides)
            .with_overrides(
                {
                    "render": {
                        "render_flagella": True,
                        "render_flagella_2d": False,
                        "show_flagella_helix_axis_3d": True,
                        "label_flagella": False,
                        "follow_camera_3d": True,
                        "save_frames_3d": False,
                        "save_frames_2d": False,
                    },
                    "output_sampling": {
                        "out_all_steps_3d": False,
                        "fps_out_3d": 25.0,
                    },
                }
            )
        )
        condition_dir = output_dir / condition.condition_id / "sim"
        condition_dir.mkdir(parents=True, exist_ok=True)
        simulator = Simulator(cfg)
        states = simulator.run(
            cfg.time.duration_s,
            step_summary_dir=condition_dir,
            stop_on_shape_fail=False,
            progress_interval=args.progress_interval,
        )
        states_by_condition.append(states)
        cfg_by_condition.append(cfg)
        rig_by_condition.append(simulator.rig)
        rows = list(
            csv.DictReader((condition_dir / "step_summary.csv").open(encoding="utf-8"))
        )
        summary = hook_overstretch._summary_row(
            cfg,
            condition,
            condition_dir,
            rows[-1],
            rows,
            helix_summary={},
            axis_center_summary=hook_overstretch._axis_center_phase_summary(
                list(
                    csv.DictReader(
                        (condition_dir / "flag_helix_axis_diagnostics.csv").open(
                            encoding="utf-8"
                        )
                    )
                )
            ),
        )
        condition_rows.append({k: str(v) for k, v in summary.items()})

    summary_path = output_dir / "summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=hook_overstretch.SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(condition_rows)

    render_result = _render_grid_movie(
        states_by_condition=states_by_condition,
        cfg_by_condition=cfg_by_condition,
        rig_by_condition=rig_by_condition,
        condition_rows=condition_rows,
        out_dir=output_dir,
    )

    manifest = {
        "git": _git_info(),
        "input": {
            "config": str(args.config),
            "duration_s": args.duration_s,
            "dt_star": args.dt_star,
            "torque_nm": args.torque_nm,
            "attach_seed": args.attach_seed,
            "phase_seed": args.phase_seed,
            "n_flagella": args.n_flagella,
        },
        "conditions": [condition.condition_id for condition in conditions],
        "render_video": {"grid_swim3d": render_result.to_manifest()},
        "outputs": {
            "root": str(output_dir),
            "summary_csv": str(summary_path),
            "log": str(output_dir / "run.log"),
        },
        "temporary_script": True,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    logger.info("Completed issue97 grid qualitative render: %s", render_result.path)
    print(output_dir)


if __name__ == "__main__":
    main()
