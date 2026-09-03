"""Render fixed-world, phase-seed RPY hydrodynamics comparison videos."""

from __future__ import annotations
import argparse
import json
from dataclasses import replace
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
    velocity_contributions,
)
from sim_swim.analysis.hydrodynamics_campaign import (
    FLOW_SLICE_GRID_SIZE,
    FLOW_VOLUME_GRID_SIZE,
    body_fixed_flow_slice,
    qc_result,
)
from sim_swim.render.render3d import _select_frames, plot_swim_frame_3d
from sim_swim.render.video_writer import open_mp4_writer
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _manifest(root: Path) -> dict[str, Any]:
    return json.loads((root / "run_manifest.json").read_text())


def _condition_dir(root: Path, record: dict[str, Any]) -> Path:
    candidate = Path(str(record.get("output_dir", "")))
    for path in (
        candidate,
        root / candidate.name,
        root / "conditions" / candidate.name,
    ):
        if path.is_dir():
            return path
    return root / "conditions" / str(record["condition_id"])


def _load_cfg(root: Path, record: dict[str, Any]) -> SimulationConfig:
    raw = yaml.safe_load(Path(_manifest(root)["base_config"]).read_text()) or {}
    cfg = SimulationConfig.from_dict(raw).with_overrides(
        normalize_base_overrides(record["config_overrides"])
    )
    return replace(
        cfg, render=replace(cfg.render, follow_camera_3d=False, label_flagella=False)
    )


def phase_seed_groups(root: Path) -> dict[tuple[int, int], list[dict[str, Any]]]:
    groups: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for record in _manifest(root)["conditions"]:
        axis = record["axis_values"]
        groups.setdefault(
            (int(axis["attach_seed"]), int(axis["n_flagella"])), []
        ).append(record)
    for key, records in groups.items():
        records.sort(key=lambda item: int(item["axis_values"]["phase_seed"]))
        if [int(item["axis_values"]["phase_seed"]) for item in records] != [0, 1, 2]:
            raise ValueError(f"{key}: expected phase seeds 0, 1, 2")
    if len(groups) != 12:
        raise ValueError(f"expected 12 attach/count groups, found {len(groups)}")
    return dict(sorted(groups.items()))


def _flow_grid(positions_um: np.ndarray, extent: float) -> np.ndarray:
    c = np.mean(positions_um, axis=0)
    offsets = np.linspace(-0.75 * extent, 0.75 * extent, FLOW_VOLUME_GRID_SIZE)
    return np.asarray(
        [
            [c[0] + x, c[1] + y, c[2] + z]
            for x in offsets
            for y in offsets
            for z in offsets
        ]
    )


def _bounds(hydros: list[Any]) -> list[list[float]]:
    p = np.concatenate([h.positions_m.reshape(-1, 3) * 1e6 for h in hydros])
    lo, hi = p.min(0), p.max(0)
    m = np.maximum((hi - lo) * 0.12, 1.0)
    return [[float(a), float(b)] for a, b in zip(lo - m, hi + m, strict=True)]


def _world_axes(ax: Any, b: list[list[float]]) -> None:
    ax.set_xlim(*b[0])
    ax.set_ylim(*b[1])
    ax.set_zlim(*b[2])
    ax.set_xlabel("x [µm]", fontsize=7)
    ax.set_ylabel("y [µm]", fontsize=7)
    ax.set_zlabel("z [µm]", fontsize=7)
    ax.tick_params(labelsize=6)


def _visible(v: np.ndarray) -> np.ndarray:
    ref = max(float(np.quantile(np.linalg.norm(v, axis=1), 0.95)), 1e-30)
    return v / ref * 0.55


def _source(
    ax: Any,
    state: Any,
    cfg: SimulationConfig,
    sim: Simulator,
    hydro: Any,
    index: int,
    b: list[list[float]],
    phase: int,
) -> None:
    plot_swim_frame_3d(
        ax,
        state,
        cfg,
        sim.rig,
        hide_ticks=False,
        title=f"phase seed {phase}: source contribution",
        show_legend=False,
    )
    p, f, body = (
        hydro.positions_m[index],
        hydro.total_forces_N[index],
        hydro.bead_is_body,
    )
    masks = [body] + [
        hydro.bead_flagella_id == i for i in sorted(set(hydro.bead_flagella_id[~body]))
    ]
    for color, mask in zip(
        ("tab:red", "tab:blue", "tab:green", "tab:purple", "tab:brown"),
        masks,
        strict=False,
    ):
        v = (
            velocity_contributions(
                p,
                f,
                bead_radius_m=hydro.bead_radius_m,
                viscosity_Pa_s=hydro.viscosity_Pa_s,
                source_mask=mask,
            )["source"][body]
            * 1e6
        )
        ax.quiver(
            *(p[body].T * 1e6), *_visible(v).T, color=color, linewidth=0.7, alpha=0.8
        )
    ax.text2D(
        0.02,
        0.02,
        "arrows at body beads: body / each flagellum source",
        transform=ax.transAxes,
        fontsize=6,
    )
    _world_axes(ax, b)


def _slice(ax: Any, hydro: Any, index: int, phase: int) -> None:
    points, v, r = body_fixed_flow_slice(hydro, index)
    p = hydro.positions_m[index]
    c = np.mean(p[hydro.bead_is_body], 0)
    beads = (p - c) @ r * 1e6
    ax.quiver(
        points[:, 0],
        points[:, 1],
        v[:, 0],
        v[:, 1],
        color="tab:blue",
        alpha=0.6,
        width=0.002,
    )
    ax.scatter(
        beads[hydro.bead_is_body, 0], beads[hydro.bead_is_body, 1], s=5, color="black"
    )
    ax.scatter(
        beads[~hydro.bead_is_body, 0],
        beads[~hydro.bead_is_body, 1],
        s=4,
        color="tab:green",
    )
    ax.set(
        title=f"phase seed {phase}: body-fixed RPY slice",
        xlabel="long axis [µm]",
        ylabel="radial axis [µm]",
        aspect="equal",
    )
    ax.tick_params(labelsize=6)


def render_phase_seed_group(
    root: Path,
    attach_seed: int,
    n_flagella: int,
    output_dir: Path,
    *,
    fps: float = 25.0,
) -> tuple[Path, dict[str, Any]]:
    records = phase_seed_groups(root)[(attach_seed, n_flagella)]
    valid = {}
    omitted = []
    for rec in records:
        phase = int(rec["axis_values"]["phase_seed"])
        d = _condition_dir(root, rec)
        passed, qc = qc_result(d)
        if not passed:
            omitted.append(
                {
                    "condition_id": rec["condition_id"],
                    "phase_seed": phase,
                    "reason": "strict_shape_qc_not_passed",
                    "qc": qc,
                }
            )
            continue
        cfg = _load_cfg(root, rec)
        hydro = load_hydro_archive(d / "hydro_archive.npz")
        valid[phase] = (
            hydro,
            _select_frames(load_state_archive(d / "state_archive.npz"), False, fps),
            cfg,
            Simulator(cfg),
        )
    b = _bounds([v[0] for v in valid.values()])
    output_dir.mkdir(parents=True, exist_ok=True)
    movie = output_dir / f"as{attach_seed:03d}__nf{n_flagella:02d}_phase_seeds_flow.mp4"
    fig = plt.figure(figsize=(15, 15), dpi=100)
    canvas = FigureCanvasAgg(fig)
    writer = None
    count = max(len(v[1]) for v in valid.values())
    try:
        for fi in range(count):
            fig.clear()
            gs = fig.add_gridspec(3, 3)
            for row, phase in enumerate((0, 1, 2)):
                if phase not in valid:
                    ax = fig.add_subplot(gs[row, :])
                    ax.axis("off")
                    ax.text(
                        0.5,
                        0.5,
                        "QC failed; visualization omitted",
                        ha="center",
                        va="center",
                        fontsize=16,
                    )
                    continue
                hydro, states, cfg, sim = valid[phase]
                state = states[min(fi, len(states) - 1)]
                idx = int(np.argmin(np.abs(hydro.t_s - state.t)))
                grid = _flow_grid(state.bead_positions_um, cfg.render.view_range_um)
                vel = (
                    rpy_flow_velocity(
                        grid * 1e-6,
                        hydro.positions_m[idx],
                        hydro.total_forces_N[idx],
                        bead_radius_m=hydro.bead_radius_m,
                        viscosity_Pa_s=hydro.viscosity_Pa_s,
                    )
                    * 1e6
                )
                ax = fig.add_subplot(gs[row, 0], projection="3d")
                plot_swim_frame_3d(
                    ax,
                    state,
                    cfg,
                    sim.rig,
                    hide_ticks=False,
                    title=f"phase seed {phase}: world RPY flow",
                    show_legend=False,
                    flow_vectors=(grid, vel),
                )
                _world_axes(ax, b)
                _source(
                    fig.add_subplot(gs[row, 1], projection="3d"),
                    state,
                    cfg,
                    sim,
                    hydro,
                    idx,
                    b,
                    phase,
                )
                _slice(fig.add_subplot(gs[row, 2]), hydro, idx, phase)
            fig.suptitle(
                f"attach seed {attach_seed}; n_flagella {n_flagella}; t={fi / fps:.3f} s",
                fontsize=12,
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
        plt.close(fig)
    metadata = {
        "attach_seed": attach_seed,
        "n_flagella": n_flagella,
        "movie": movie.name,
        "world_bounds_um": b,
        "phase_conditions": [
            {
                "phase_seed": int(r["axis_values"]["phase_seed"]),
                "condition_id": r["condition_id"],
            }
            for r in records
        ],
        "omitted": omitted,
        "frame_count": count,
    }
    manifest_path = output_dir.parent / "analysis_manifest.json"
    payload = (
        json.loads(manifest_path.read_text())
        if manifest_path.is_file()
        else {
            "kind": "free_space_rpy_phase_seed_flow_visualization",
            "input_campaign": str(root),
            "input_provenance": "hydro_archive.npz positions_m + total_forces_N",
            "layout": {
                "rows": "phase seeds 0, 1, 2",
                "columns": [
                    "world RPY flow",
                    "world source contribution",
                    "body-fixed long-axis slice",
                ],
            },
            "units": {"position": "µm", "velocity": "µm/s", "force": "N"},
            "world_flow_grid_shape": [FLOW_VOLUME_GRID_SIZE] * 3,
            "body_fixed_slice_grid_shape": [FLOW_SLICE_GRID_SIZE] * 2,
            "stokes_force_balance_verified": "F_hydro = -F_total; F_total + F_hydro = 0",
            "fps": fps,
            "groups": [],
        }
    )
    payload["groups"] = [
        item
        for item in payload["groups"]
        if (item["attach_seed"], item["n_flagella"]) != (attach_seed, n_flagella)
    ] + [metadata]
    payload["groups"].sort(key=lambda item: (item["attach_seed"], item["n_flagella"]))
    manifest_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    return movie, metadata


def render_all_phase_seed_groups(
    root: Path, output_dir: Path, *, fps: float = 25.0
) -> list[Path]:
    results = [
        render_phase_seed_group(root, a, n, output_dir, fps=fps)
        for a, n in phase_seed_groups(root)
    ]
    payload = {
        "kind": "free_space_rpy_phase_seed_flow_visualization",
        "input_campaign": str(root),
        "input_provenance": "hydro_archive.npz positions_m + total_forces_N",
        "layout": {
            "rows": "phase seeds 0, 1, 2",
            "columns": [
                "world RPY flow",
                "world source contribution",
                "body-fixed long-axis slice",
            ],
        },
        "units": {"position": "µm", "velocity": "µm/s", "force": "N"},
        "world_flow_grid_shape": [FLOW_VOLUME_GRID_SIZE] * 3,
        "body_fixed_slice_grid_shape": [FLOW_SLICE_GRID_SIZE] * 2,
        "stokes_force_balance_verified": "F_hydro = -F_total; F_total + F_hydro = 0",
        "fps": fps,
        "groups": [r for _, r in results],
    }
    (output_dir.parent / "analysis_manifest.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n"
    )
    return [m for m, _ in results]


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, required=True)
    p.add_argument("--phase-seed-groups", action="store_true", required=True)
    p.add_argument("--attach-seed", type=int)
    p.add_argument("--n-flagella", type=int)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--fps", type=float, default=25.0)
    a = p.parse_args(argv)
    if (a.attach_seed is None) != (a.n_flagella is None):
        p.error("--attach-seed and --n-flagella must be supplied together")
    if a.attach_seed is None:
        for movie in render_all_phase_seed_groups(a.run_dir, a.output_dir, fps=a.fps):
            print(movie)
    else:
        movie, _metadata = render_phase_seed_group(
            a.run_dir, a.attach_seed, a.n_flagella, a.output_dir, fps=a.fps
        )
        print(movie)


if __name__ == "__main__":
    main()
