"""Mac-side 3D plus nominal local-segment torque-weight replay for #203."""

from __future__ import annotations
import argparse
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import yaml
from matplotlib.backends.backend_agg import FigureCanvasAgg

from sim_swim.analysis.flagella_count_behavior import (
    load_state_archive,
    validate_replay_fps,
)
from sim_swim.render.render3d import _select_frames, plot_swim_frame_3d
from sim_swim.render.video_writer import open_mp4_writer
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def nominal_segment_weights(profile: str, segment_count: int) -> np.ndarray:
    """Return normalized nominal weights on the engine's bead-to-bead segments.

    #203 renders the configured nominal profile, not realized force.  Uniform is
    time-invariant; dynamic profiles intentionally require their recorded local
    twist state and are rejected rather than silently misrepresented.
    """
    if segment_count <= 0:
        raise ValueError("segment_count must be positive")
    if profile != "uniform":
        raise ValueError("composite replay currently requires uniform nominal weights")
    return np.full(segment_count, 1.0 / segment_count, dtype=float)


def _condition(root: Path, condition_id: str) -> tuple[dict[str, Any], dict[str, Any]]:
    manifest = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
    for record in manifest["conditions"]:
        if record["condition_id"] == condition_id:
            return manifest, record
    raise KeyError(f"condition not found: {condition_id}")


def render(root: Path, condition_id: str, output_dir: Path, fps: float = 10.0) -> Path:
    manifest, record = _condition(root, condition_id)
    base = Path(str(manifest["base_config"]))
    cfg = SimulationConfig.from_dict(
        yaml.safe_load(base.read_text(encoding="utf-8"))
    ).with_overrides(record["config_overrides"])
    if cfg.motor.force_distribution != "root_torque_segment_couples":
        raise ValueError("composite replay requires root_torque_segment_couples")
    output_dir.mkdir(parents=True, exist_ok=True)
    recorded = Path(str(record["output_dir"]))
    candidates = (recorded, root / recorded.name, root / "conditions" / recorded.name)
    raw = next(
        (candidate for candidate in candidates if candidate.is_dir()), candidates[-1]
    )
    states = load_state_archive(raw / "state_archive.npz")
    validate_replay_fps(states, fps)
    simulator = Simulator(cfg)
    frames = _select_frames(states, False, fps)
    figure = plt.figure(
        figsize=(11, max(4, 2.35 * len(simulator.rig.flagella_indices)))
    )
    canvas = FigureCanvasAgg(figure)
    movie_path = output_dir / f"{condition_id}_composite.mp4"
    selection = open_mp4_writer(
        movie_path, fps=fps, frame_size=(1100, int(figure.get_figheight() * 100))
    )
    try:
        for state in frames:
            figure.clear()
            ax3d = figure.add_axes([0.02, 0.08, 0.57, 0.86], projection="3d")
            plot_swim_frame_3d(
                ax3d,
                state,
                cfg,
                simulator.rig,
                title=f"{condition_id} | t={state.t:.3f} s",
                hide_ticks=True,
            )
            count = len(simulator.rig.flagella_indices)
            for flag_id, indices in enumerate(simulator.rig.flagella_indices):
                top = 0.94 - (flag_id + 1) * (0.86 / count)
                axis = figure.add_axes([0.66, top, 0.31, 0.86 / count - 0.035])
                weights = nominal_segment_weights(
                    cfg.motor.torque_distribution_profile, len(indices) - 1
                )
                axis.bar(np.arange(len(weights)), weights, color=f"C{flag_id}")
                axis.set(
                    title=f"Flagellum {flag_id}: sum={weights.sum():.3f}",
                    xlabel="local segment number",
                    ylabel="normalized nominal weight",
                    ylim=(0, max(float(weights.max()) * 1.25, 0.1)),
                )
                axis.grid(alpha=0.25)
            canvas.draw()
            rgba = np.asarray(canvas.buffer_rgba())
            selection.writer.write(rgba[:, :, :3])
    finally:
        selection.writer.release()
        plt.close(figure)
    (output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "kind": "phase2_issue203_composite_replay",
                "condition_id": condition_id,
                "profile": cfg.motor.torque_distribution_profile,
                "weight_unit": "local bead-to-bead segment",
                "weight_sum": 1.0,
                "movie": str(movie_path),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return movie_path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--condition-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fps", type=float, default=10.0)
    args = parser.parse_args(argv)
    print(render(args.run_dir, args.condition_id, args.output_dir, args.fps))
