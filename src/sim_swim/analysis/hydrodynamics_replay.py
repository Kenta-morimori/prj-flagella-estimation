"""Opt-in fixed-camera 3D replay with reconstructed free-space RPY flow."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
from typing import Any

import cv2
import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.backends.backend_agg import FigureCanvasAgg

from sim_swim.analysis.flagella_count_behavior import (
    load_state_archive,
    normalize_base_overrides,
)
from sim_swim.analysis.hydrodynamics import load_hydro_archive, rpy_flow_velocity
from sim_swim.render.render3d import _select_frames, plot_swim_frame_3d
from sim_swim.render.video_writer import open_mp4_writer
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _load_cfg(root: Path, record: dict[str, Any]) -> SimulationConfig:
    manifest = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
    raw = (
        yaml.safe_load(Path(manifest["base_config"]).read_text(encoding="utf-8")) or {}
    )
    cfg = SimulationConfig.from_dict(raw).with_overrides(
        normalize_base_overrides(record["config_overrides"])
    )
    return replace(
        cfg, render=replace(cfg.render, follow_camera_3d=False, label_flagella=False)
    )


def _condition_record(root: Path, condition_id: str) -> dict[str, Any]:
    manifest = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
    return next(
        record
        for record in manifest["conditions"]
        if record["condition_id"] == condition_id
    )


def _flow_grid(state_positions_um: np.ndarray, view_range_um: float) -> np.ndarray:
    center = np.mean(state_positions_um, axis=0)
    offsets = np.linspace(-0.75 * view_range_um, 0.75 * view_range_um, 3)
    return np.asarray(
        [
            [center[0] + x, center[1] + y, center[2] + z]
            for x in offsets
            for y in offsets
            for z in offsets
        ]
    )


def render_condition(
    root: Path, condition_id: str, output_dir: Path, *, fps: float = 25.0
) -> Path:
    """Render one archive. Fixed camera is intentionally forced for overlays."""
    record = _condition_record(root, condition_id)
    cfg = _load_cfg(root, record)
    condition_dir = Path(record.get("output_dir", root / condition_id))
    states = _select_frames(
        load_state_archive(condition_dir / "state_archive.npz"), False, fps
    )
    hydro = load_hydro_archive(condition_dir / "hydro_archive.npz")
    simulator = Simulator(cfg)
    output_dir.mkdir(parents=True, exist_ok=True)
    movie = output_dir / f"{condition_id}_flow_overlay.mp4"
    figure = plt.figure(figsize=(6, 6), dpi=100)
    canvas = FigureCanvasAgg(figure)
    writer = None
    try:
        for state in states:
            figure.clear()
            index = int(np.argmin(np.abs(hydro.t_s - state.t)))
            grid_um = _flow_grid(state.bead_positions_um, cfg.render.view_range_um)
            velocity_um_s = 1e6 * rpy_flow_velocity(
                grid_um * 1e-6,
                hydro.positions_m[index],
                hydro.total_forces_N[index],
                bead_radius_m=hydro.bead_radius_m,
                viscosity_Pa_s=hydro.viscosity_Pa_s,
            )
            ax = figure.add_subplot(111, projection="3d")
            plot_swim_frame_3d(
                ax,
                state,
                cfg,
                simulator.rig,
                hide_ticks=True,
                title=condition_id,
                flow_vectors=(grid_um, velocity_um_s),
            )
            canvas.draw()
            frame = cv2.cvtColor(np.asarray(canvas.buffer_rgba()), cv2.COLOR_RGBA2BGR)
            if writer is None:
                selection = open_mp4_writer(
                    movie, fps=fps, frame_size=(frame.shape[1], frame.shape[0])
                )
                writer = selection.writer
            writer.write(frame)
    finally:
        if writer is not None:
            writer.release()
        plt.close(figure)
    (output_dir / f"{condition_id}_flow_overlay.json").write_text(
        json.dumps(
            {
                "condition_id": condition_id,
                "follow_camera_3d": False,
                "fps": fps,
                "grid_shape": [3, 3, 3],
                "movie": str(movie),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return movie


def overlay_manifest(condition_ids: list[str], *, fps: float) -> dict[str, Any]:
    """Return the fixed visualization contract for a campaign comparison."""
    return {
        "condition_ids": condition_ids,
        "follow_camera_3d": False,
        "fps": fps,
        "grid_shape": [3, 3, 3],
        "panel_count": len(condition_ids),
    }


def render_campaign_overlay(
    root: Path, condition_ids: list[str], output_dir: Path, *, fps: float = 25.0
) -> Path:
    """Render a synchronized fixed-camera horizontal multi-condition overlay."""
    if len(condition_ids) != 3:
        raise ValueError("campaign overlay requires exactly three condition IDs")
    assets = []
    for condition_id in condition_ids:
        record = _condition_record(root, condition_id)
        cfg = _load_cfg(root, record)
        condition_dir = Path(record.get("output_dir", root / condition_id))
        states = _select_frames(
            load_state_archive(condition_dir / "state_archive.npz"), False, fps
        )
        assets.append(
            (
                condition_id,
                cfg,
                Simulator(cfg),
                states,
                load_hydro_archive(condition_dir / "hydro_archive.npz"),
            )
        )
    frame_count = min(len(item[3]) for item in assets)
    if frame_count == 0:
        raise ValueError("campaign overlay has no frames")
    output_dir.mkdir(parents=True, exist_ok=True)
    movie = output_dir / "nflagella_phase0_flow_overlay.mp4"
    figure = plt.figure(figsize=(15, 5), dpi=100)
    canvas = FigureCanvasAgg(figure)
    writer = None
    try:
        for frame_index in range(frame_count):
            figure.clear()
            for panel_index, (condition_id, cfg, simulator, states, hydro) in enumerate(
                assets, start=1
            ):
                state = states[frame_index]
                index = int(np.argmin(np.abs(hydro.t_s - state.t)))
                grid_um = _flow_grid(state.bead_positions_um, cfg.render.view_range_um)
                velocity_um_s = 1e6 * rpy_flow_velocity(
                    grid_um * 1e-6,
                    hydro.positions_m[index],
                    hydro.total_forces_N[index],
                    bead_radius_m=hydro.bead_radius_m,
                    viscosity_Pa_s=hydro.viscosity_Pa_s,
                )
                axis = figure.add_subplot(1, 3, panel_index, projection="3d")
                plot_swim_frame_3d(
                    axis,
                    state,
                    cfg,
                    simulator.rig,
                    hide_ticks=True,
                    title=condition_id,
                    flow_vectors=(grid_um, velocity_um_s),
                )
            canvas.draw()
            frame = cv2.cvtColor(np.asarray(canvas.buffer_rgba()), cv2.COLOR_RGBA2BGR)
            if writer is None:
                writer = open_mp4_writer(
                    movie, fps=fps, frame_size=(frame.shape[1], frame.shape[0])
                ).writer
            writer.write(frame)
    finally:
        if writer is not None:
            writer.release()
        plt.close(figure)
    payload = {
        **overlay_manifest(condition_ids, fps=fps),
        "movie": str(movie),
        "frame_count": frame_count,
    }
    (output_dir / "nflagella_phase0_flow_overlay.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    return movie


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--condition-id")
    parser.add_argument("--campaign-nflagella-phase0", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fps", type=float, default=25.0)
    args = parser.parse_args(argv)
    if args.campaign_nflagella_phase0:
        print(
            render_campaign_overlay(
                args.run_dir,
                ["n1__phase0", "n2__phase0", "n3__phase0"],
                args.output_dir,
                fps=args.fps,
            )
        )
    elif args.condition_id:
        print(
            render_condition(
                args.run_dir, args.condition_id, args.output_dir, fps=args.fps
            )
        )
    else:
        parser.error("--condition-id or --campaign-nflagella-phase0 is required")
