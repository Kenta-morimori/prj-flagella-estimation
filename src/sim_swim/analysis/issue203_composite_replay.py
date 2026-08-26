"""Mac-side 3D plus nominal local-segment torque-weight replay for #203."""

from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import yaml
import cv2
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
        raise ValueError("nominal_segment_weights is only defined for uniform")
    return np.full(segment_count, 1.0 / segment_count, dtype=float)


def _profile_from_local_twist(profile: str, orientation: np.ndarray) -> np.ndarray:
    if profile == "uniform":
        return np.ones_like(orientation, dtype=float)
    if profile != "diffusive":
        raise ValueError(f"unsupported composite torque profile: {profile}")
    activity = np.abs(np.asarray(orientation, dtype=float))
    maximum = float(np.max(activity)) if activity.size else 0.0
    return np.ones_like(activity) if maximum <= 1e-12 else activity / maximum


def reconstructed_segment_weights(
    profile: str,
    segment_count: int,
    *,
    times_s: np.ndarray,
    dt_s: float,
    torque_Nm: float,
) -> list[np.ndarray]:
    """Reconstruct the engine's nominal local-twist segment weights.

    ``diffusive`` is dynamic but its local twist state is deterministic for the
    fixed RUN campaign.  The state archive does not store that internal state,
    so replay reproduces the exact update rule used by ``DynamicsEngine``.
    Returned vectors are normalized exactly as the force-couple implementation.
    """
    if segment_count <= 0 or dt_s <= 0:
        raise ValueError("segment_count and dt_s must be positive")
    if profile == "uniform":
        weight = nominal_segment_weights(profile, segment_count)
        return [weight.copy() for _ in times_s]
    if profile != "diffusive":
        raise ValueError(f"unsupported composite torque profile: {profile}")

    orientation = np.zeros(segment_count, dtype=float)
    current_step = 0
    result: list[np.ndarray] = []
    drive_rate = 2.0 * math.pi * 2.2 * (float(torque_Nm) / 2.0e-20)
    for t_s in np.asarray(times_s, dtype=float):
        target_step = max(current_step, int(round(float(t_s) / dt_s)))
        while current_step < target_step:
            lap = np.zeros_like(orientation)
            if segment_count > 1:
                lap[0] = orientation[1] - orientation[0]
                lap[-1] = orientation[-2] - orientation[-1]
                if segment_count > 2:
                    lap[1:-1] = (
                        orientation[:-2] - 2.0 * orientation[1:-1] + orientation[2:]
                    )
            orientation += dt_s * (80.0 * lap - 0.05 * orientation)
            orientation[0] += drive_rate * dt_s
            current_step += 1
        weight = _profile_from_local_twist(profile, orientation)
        result.append(weight / max(float(np.sum(weight)), 1e-12))
    return result


def composite_manifest_path(output_dir: Path, condition_id: str) -> Path:
    """Return the per-condition manifest path for a composite replay."""
    return output_dir / f"{condition_id}_composite_manifest.json"


def _written_frame_count(path: Path) -> int:
    """Return decodable frame count after a writer is released."""
    ffprobe = shutil.which("ffprobe")
    if ffprobe:
        completed = subprocess.run(
            [
                ffprobe,
                "-v",
                "error",
                "-count_frames",
                "-select_streams",
                "v:0",
                "-show_entries",
                "stream=nb_read_frames",
                "-of",
                "default=nokey=1:noprint_wrappers=1",
                str(path),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        try:
            return int(completed.stdout.strip())
        except ValueError:
            pass
    capture = cv2.VideoCapture(str(path))
    try:
        if not capture.isOpened():
            return 0
        return int(round(capture.get(cv2.CAP_PROP_FRAME_COUNT)))
    finally:
        capture.release()


def _write_frames_h264(
    movie_path: Path,
    frames: list[np.ndarray],
    *,
    fps: float,
    frame_size: tuple[int, int],
) -> tuple[str, tuple[str, ...]]:
    """Write RGB matplotlib frames through the shared H.264 BGR writer."""
    expected_shape = (frame_size[1], frame_size[0], 3)
    if any(
        image.shape != expected_shape or image.dtype != np.uint8 for image in frames
    ):
        received = [(image.shape, str(image.dtype)) for image in frames[:1]]
        raise ValueError(
            f"rendered frame does not match rawvideo shape {expected_shape}: {received}"
        )
    movie_path.unlink(missing_ok=True)
    selection = open_mp4_writer(movie_path, fps=fps, frame_size=frame_size)
    try:
        for image in frames:
            # Matplotlib emits RGB whereas the shared rawvideo contract is BGR.
            selection.writer.write(cv2.cvtColor(image, cv2.COLOR_RGB2BGR))
    finally:
        selection.writer.release()
    if _written_frame_count(movie_path) != len(frames):
        raise RuntimeError("H.264 MP4 writer produced no decodable frames")
    return selection.selected_codec, selection.attempted_codecs


def _condition(root: Path, condition_id: str) -> tuple[dict[str, Any], dict[str, Any]]:
    manifest = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
    for record in manifest["conditions"]:
        if record["condition_id"] == condition_id:
            return manifest, record
    raise KeyError(f"condition not found: {condition_id}")


def _condition_assets(
    root: Path, manifest: dict[str, Any], record: dict[str, Any], fps: float
) -> tuple[SimulationConfig, Simulator, list[Any]]:
    base = Path(str(manifest["base_config"]))
    cfg = SimulationConfig.from_dict(
        yaml.safe_load(base.read_text(encoding="utf-8"))
    ).with_overrides(record["config_overrides"])
    if cfg.motor.force_distribution != "root_torque_segment_couples":
        raise ValueError("composite replay requires root_torque_segment_couples")
    recorded = Path(str(record["output_dir"]))
    candidates = (recorded, root / recorded.name, root / "conditions" / recorded.name)
    raw = next(
        (candidate for candidate in candidates if candidate.is_dir()), candidates[-1]
    )
    states = load_state_archive(raw / "state_archive.npz")
    validate_replay_fps(states, fps)
    simulator = Simulator(cfg)
    return cfg, simulator, _select_frames(states, False, fps)


def _panel_weights(
    cfg: SimulationConfig, simulator: Simulator, frames: list[Any]
) -> list[list[np.ndarray]]:
    return [
        reconstructed_segment_weights(
            cfg.motor.torque_distribution_profile,
            len(indices) - 1,
            times_s=np.asarray([state.t for state in frames], dtype=float),
            dt_s=cfg.dt_s,
            torque_Nm=cfg.motor.torque_Nm,
        )
        for indices in simulator.rig.flagella_indices
    ]


def render(root: Path, condition_id: str, output_dir: Path, fps: float = 10.0) -> Path:
    manifest, record = _condition(root, condition_id)
    output_dir.mkdir(parents=True, exist_ok=True)
    cfg, simulator, frames = _condition_assets(root, manifest, record, fps)
    figure = plt.figure(
        figsize=(11, max(4, 2.35 * len(simulator.rig.flagella_indices)))
    )
    canvas = FigureCanvasAgg(figure)
    movie_path = output_dir / f"{condition_id}_composite.mp4"
    panel_weights = _panel_weights(cfg, simulator, frames)
    rendered_frames: list[np.ndarray] = []
    try:
        for frame_index, state in enumerate(frames):
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
                weights = panel_weights[flag_id][frame_index]
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
            rendered_frames.append(rgba[:, :, :3].copy())
    finally:
        plt.close(figure)
    frame_size = (rendered_frames[0].shape[1], rendered_frames[0].shape[0])
    selected_codec, attempted_codecs = _write_frames_h264(
        movie_path, rendered_frames, fps=fps, frame_size=frame_size
    )
    composite_manifest_path(output_dir, condition_id).write_text(
        json.dumps(
            {
                "kind": "phase2_issue203_composite_replay",
                "condition_id": condition_id,
                "profile": cfg.motor.torque_distribution_profile,
                "weight_unit": "local bead-to-bead segment",
                "weight_sum": 1.0,
                "weight_time_behavior": (
                    "time_invariant"
                    if cfg.motor.torque_distribution_profile == "uniform"
                    else "reconstructed_local_twist"
                ),
                "weight_reconstruction": {
                    "dt_s": cfg.dt_s,
                    "torque_Nm": cfg.motor.torque_Nm,
                },
                "selected_codec": selected_codec,
                "attempted_codecs": list(attempted_codecs),
                "frame_count": len(rendered_frames),
                "movie": str(movie_path),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return movie_path


def render_n_flagella_grid(
    root: Path, n_flagella: int, output_dir: Path, fps: float = 10.0
) -> Path:
    """Render the nine attach/phase seed conditions as one composite grid."""
    manifest = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
    records = sorted(
        (
            record
            for record in manifest["conditions"]
            if int((record.get("axis_values") or {}).get("n_flagella", -1))
            == n_flagella
        ),
        key=lambda record: (
            int((record.get("axis_values") or {}).get("attach_seed", -1)),
            int((record.get("axis_values") or {}).get("phase_seed", -1)),
        ),
    )
    if len(records) != 9:
        raise ValueError(f"Expected nine conditions for n_flagella={n_flagella}")
    assets = [_condition_assets(root, manifest, record, fps) for record in records]
    profiles = {cfg.motor.torque_distribution_profile for cfg, _, _ in assets}
    if len(profiles) != 1:
        raise ValueError("Grid conditions must use one torque profile")
    frame_count = min(len(frames) for _, _, frames in assets)
    if frame_count <= 0:
        raise RuntimeError("No composite grid frames selected")
    output_dir.mkdir(parents=True, exist_ok=True)
    movie_path = output_dir / f"nf{n_flagella:02d}_composite_grid.mp4"
    manifest_path = output_dir / f"nf{n_flagella:02d}_composite_grid_manifest.json"
    # A 16:9 canvas gives every 3x3 cell enough horizontal room for both the
    # 3D replay and its local-segment torque panel without oversized labels.
    figure = plt.figure(figsize=(16, 9))
    canvas = FigureCanvasAgg(figure)
    weights = [
        _panel_weights(cfg, simulator, frames[:frame_count])
        for cfg, simulator, frames in assets
    ]
    rendered_frames: list[np.ndarray] = []
    try:
        for frame_index in range(frame_count):
            figure.clear()
            for cell, (record, (cfg, simulator, frames), panel_weights) in enumerate(
                zip(records, assets, weights)
            ):
                row, col = divmod(cell, 3)
                x, y, width, height = col / 3, 1.0 - (row + 1) / 3, 1 / 3, 1 / 3
                condition_id = str(record["condition_id"])
                state = frames[frame_index]
                ax3d = figure.add_axes(
                    [x + 0.006, y + 0.025, width * 0.62, height * 0.93],
                    projection="3d",
                )
                plot_swim_frame_3d(
                    ax3d,
                    state,
                    cfg,
                    simulator.rig,
                    title=f"{condition_id}\nt={state.t:.2f} s",
                    hide_ticks=True,
                    show_status=False,
                    show_legend=False,
                )
                ax3d.set_title(f"{condition_id}\nt={state.t:.2f} s", fontsize=6, pad=2)
                count = len(simulator.rig.flagella_indices)
                for flag_id in range(count):
                    panel_height = height * 0.92 / count
                    top = y + height * 0.955 - (flag_id + 1) * panel_height
                    axis = figure.add_axes(
                        [x + width * 0.65, top, width * 0.34, panel_height - 0.008]
                    )
                    value = panel_weights[flag_id][frame_index]
                    axis.bar(np.arange(len(value)), value, color=f"C{flag_id}")
                    axis.set_ylim(0, max(float(value.max()) * 1.25, 0.1))
                    axis.set_xticks([])
                    axis.set_yticks([])
                    axis.set_title(
                        f"F{flag_id}  Σ={value.sum():.2f}", fontsize=5, pad=1
                    )
            canvas.draw()
            rendered_frames.append(np.asarray(canvas.buffer_rgba())[:, :, :3].copy())
    finally:
        plt.close(figure)
    frame_size = (rendered_frames[0].shape[1], rendered_frames[0].shape[0])
    selected_codec, attempted_codecs = _write_frames_h264(
        movie_path, rendered_frames, fps=fps, frame_size=frame_size
    )
    cfg = assets[0][0]
    manifest_path.write_text(
        json.dumps(
            {
                "kind": "phase2_issue203_composite_grid",
                "n_flagella": n_flagella,
                "profile": cfg.motor.torque_distribution_profile,
                "grid_shape": [3, 3],
                "condition_ids": [str(record["condition_id"]) for record in records],
                "weight_unit": "local bead-to-bead segment",
                "weight_sum": 1.0,
                "weight_time_behavior": (
                    "time_invariant"
                    if cfg.motor.torque_distribution_profile == "uniform"
                    else "reconstructed_local_twist"
                ),
                "selected_codec": selected_codec,
                "attempted_codecs": list(attempted_codecs),
                "frame_count": frame_count,
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
    parser.add_argument("--condition-id")
    parser.add_argument("--n-flagella", type=int)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fps", type=float, default=10.0)
    args = parser.parse_args(argv)
    if (args.condition_id is None) == (args.n_flagella is None):
        parser.error("specify exactly one of --condition-id or --n-flagella")
    if args.n_flagella is not None:
        print(
            render_n_flagella_grid(
                args.run_dir, args.n_flagella, args.output_dir, args.fps
            )
        )
    else:
        print(render(args.run_dir, str(args.condition_id), args.output_dir, args.fps))
