"""Fixed-camera RPY flow, force-balance, and source-contribution replay."""

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
from sim_swim.analysis.hydrodynamics import (
    load_hydro_archive,
    rpy_flow_velocity,
    stokes_fluid_resistance,
    velocity_contributions,
)
from sim_swim.analysis.hydrodynamics_campaign import (
    FLOW_VOLUME_GRID_SIZE,
    select_qc_passed_conditions,
)
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


def _flow_grid(
    state_positions_um: np.ndarray, view_range_um: float, *, grid_size: int
) -> np.ndarray:
    center = np.mean(state_positions_um, axis=0)
    offsets = np.linspace(-0.75 * view_range_um, 0.75 * view_range_um, grid_size)
    return np.asarray(
        [
            [center[0] + x, center[1] + y, center[2] + z]
            for x in offsets
            for y in offsets
            for z in offsets
        ]
    )


def _visual_vectors(
    vectors: np.ndarray, *, length_um: float
) -> tuple[np.ndarray, float]:
    """Map physical vectors to visible arrows and return the reference magnitude."""
    norms = np.linalg.norm(vectors, axis=1)
    reference = max(float(np.quantile(norms, 0.95)), 1e-30)
    return vectors / reference * length_um, reference


def _source_velocity_vectors(hydro: Any, index: int) -> list[tuple[str, np.ndarray]]:
    """Velocity at body beads induced by body and each individual flagellum."""
    positions = hydro.positions_m[index]
    forces = hydro.total_forces_N[index]
    result: list[tuple[str, np.ndarray]] = []
    masks: list[tuple[str, np.ndarray]] = [("body", hydro.bead_is_body)]
    for flagellum_id in sorted(set(hydro.bead_flagella_id[~hydro.bead_is_body])):
        masks.append(
            (f"flagellum {int(flagellum_id)}", hydro.bead_flagella_id == flagellum_id)
        )
    for label, mask in masks:
        contribution = velocity_contributions(
            positions,
            forces,
            bead_radius_m=hydro.bead_radius_m,
            viscosity_Pa_s=hydro.viscosity_Pa_s,
            source_mask=mask,
        )["source"]
        result.append((label, contribution[hydro.bead_is_body] * 1e6))
    return result


def _draw_force_panel(
    ax: Any,
    state: Any,
    cfg: SimulationConfig,
    simulator: Simulator,
    hydro: Any,
    index: int,
    condition_id: str,
) -> float:
    plot_swim_frame_3d(
        ax,
        state,
        cfg,
        simulator.rig,
        hide_ticks=True,
        title=f"{condition_id}: force balance",
        show_legend=False,
    )
    positions = hydro.positions_m[index] * 1e6
    displayed_force, reference = _visual_vectors(
        hydro.total_forces_N[index], length_um=0.65
    )
    ax.quiver(
        *positions.T, *displayed_force.T, color="tab:orange", linewidth=0.8, alpha=0.85
    )
    ax.quiver(
        *positions.T,
        *stokes_fluid_resistance(displayed_force).T,
        color="tab:purple",
        linewidth=0.6,
        alpha=0.72,
    )
    ax.text2D(
        0.02,
        0.02,
        "orange: mechanical F_total\npurple: fluid resistance -F_total",
        transform=ax.transAxes,
        fontsize=7,
    )
    return reference


def _draw_source_panel(
    ax: Any,
    state: Any,
    cfg: SimulationConfig,
    simulator: Simulator,
    hydro: Any,
    index: int,
    condition_id: str,
) -> dict[str, float]:
    plot_swim_frame_3d(
        ax,
        state,
        cfg,
        simulator.rig,
        hide_ticks=True,
        title=f"{condition_id}: induced velocity",
        show_legend=False,
    )
    body_positions = hydro.positions_m[index][hydro.bead_is_body] * 1e6
    colors = ("tab:red", "tab:blue", "tab:green", "tab:purple", "tab:brown")
    references: dict[str, float] = {}
    for source_index, (label, vectors) in enumerate(
        _source_velocity_vectors(hydro, index)
    ):
        displayed, reference = _visual_vectors(vectors, length_um=0.55)
        ax.quiver(
            *body_positions.T,
            *displayed.T,
            color=colors[source_index % len(colors)],
            linewidth=0.8,
            alpha=0.78,
        )
        references[label] = reference
    ax.text2D(
        0.02,
        0.02,
        "arrows at body beads: body / each flagellum source",
        transform=ax.transAxes,
        fontsize=7,
    )
    return references


def render_condition(
    root: Path, condition_id: str, output_dir: Path, *, fps: float = 25.0
) -> Path:
    """Render one condition with dense 7^3 flow and force/source relation panels."""
    record = _condition_record(root, condition_id)
    cfg = _load_cfg(root, record)
    condition_dir = Path(record.get("output_dir", root / condition_id))
    states = _select_frames(
        load_state_archive(condition_dir / "state_archive.npz"), False, fps
    )
    hydro = load_hydro_archive(condition_dir / "hydro_archive.npz")
    simulator = Simulator(cfg)
    output_dir.mkdir(parents=True, exist_ok=True)
    movie = output_dir / f"{condition_id}_flow_force_overlay.mp4"
    figure = plt.figure(figsize=(15, 5), dpi=100)
    canvas = FigureCanvasAgg(figure)
    writer = None
    force_references: list[float] = []
    source_references: dict[str, list[float]] = {}
    try:
        for state in states:
            figure.clear()
            index = int(np.argmin(np.abs(hydro.t_s - state.t)))
            grid_um = _flow_grid(
                state.bead_positions_um,
                cfg.render.view_range_um,
                grid_size=FLOW_VOLUME_GRID_SIZE,
            )
            velocity_um_s = (
                rpy_flow_velocity(
                    grid_um * 1e-6,
                    hydro.positions_m[index],
                    hydro.total_forces_N[index],
                    bead_radius_m=hydro.bead_radius_m,
                    viscosity_Pa_s=hydro.viscosity_Pa_s,
                )
                * 1e6
            )
            flow = figure.add_subplot(1, 3, 1, projection="3d")
            plot_swim_frame_3d(
                flow,
                state,
                cfg,
                simulator.rig,
                hide_ticks=True,
                title=f"{condition_id}: RPY flow",
                show_legend=False,
                flow_vectors=(grid_um, velocity_um_s),
            )
            force = figure.add_subplot(1, 3, 2, projection="3d")
            force_references.append(
                _draw_force_panel(
                    force, state, cfg, simulator, hydro, index, condition_id
                )
            )
            source = figure.add_subplot(1, 3, 3, projection="3d")
            for label, reference in _draw_source_panel(
                source, state, cfg, simulator, hydro, index, condition_id
            ).items():
                source_references.setdefault(label, []).append(reference)
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
    (output_dir / f"{condition_id}_flow_force_overlay.json").write_text(
        json.dumps(
            {
                "condition_id": condition_id,
                "follow_camera_3d": False,
                "fps": fps,
                "flow_grid_shape": [FLOW_VOLUME_GRID_SIZE] * 3,
                "flow_arrow_policy": "normalized_direction; velocity unit um/s",
                "force_arrow_policy": "per-frame 95th-percentile scaling; force unit N",
                "force_arrow_reference_N_range": [
                    min(force_references),
                    max(force_references),
                ],
                "source_velocity_arrow_reference_um_s": {
                    label: [min(values), max(values)]
                    for label, values in source_references.items()
                },
                "stokes_force_balance": "F_hydro = -F_total",
                "movie": str(movie),
                "frame_count": len(states),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return movie


def render_all_qc_passed(
    root: Path, output_dir: Path, *, fps: float = 25.0
) -> list[Path]:
    """Render an individual video for every condition which passed strict QC."""
    selected, _skipped = select_qc_passed_conditions(root)
    return [
        render_condition(
            root, record["condition_id"], output_dir / record["condition_id"], fps=fps
        )
        for record in selected
    ]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--condition-id")
    source.add_argument("--all-qc-passed", action="store_true")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fps", type=float, default=25.0)
    args = parser.parse_args(argv)
    if args.all_qc_passed:
        for movie in render_all_qc_passed(args.run_dir, args.output_dir, fps=args.fps):
            print(movie)
    else:
        print(
            render_condition(
                args.run_dir, args.condition_id, args.output_dir, fps=args.fps
            )
        )


if __name__ == "__main__":
    main()
