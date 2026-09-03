"""Render fixed-world, phase-seed RPY hydrodynamics comparison videos."""

from __future__ import annotations
import argparse
import gc
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


def _axis_order(root: Path, axis_name: str) -> list[Any]:
    values: list[Any] = []
    for record in _manifest(root)["conditions"]:
        value = record["axis_values"].get(axis_name)
        if value not in values:
            values.append(value)
    if not values:
        raise ValueError(f"campaign has no axis {axis_name!r}")
    return values


def grouped_records(
    root: Path, *, row_axis: str, group_axes: tuple[str, ...]
) -> dict[tuple[Any, ...], list[dict[str, Any]]]:
    """Group arbitrary multi-run records without Issue- or seed-specific rules."""
    if row_axis in group_axes or not group_axes:
        raise ValueError("--row-axis must differ from at least one --group-axis")
    _axis_order(root, row_axis)
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for record in _manifest(root)["conditions"]:
        values = record["axis_values"]
        if row_axis not in values or any(axis not in values for axis in group_axes):
            raise ValueError("selected axes are not present in every condition")
        groups.setdefault(tuple(values[axis] for axis in group_axes), []).append(record)
    row_order = _axis_order(root, row_axis)
    for key, records in groups.items():
        records.sort(key=lambda item: row_order.index(item["axis_values"][row_axis]))
        actual = [item["axis_values"][row_axis] for item in records]
        if actual != row_order:
            raise ValueError(f"{key}: incomplete or duplicate {row_axis} rows")
    return groups


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
    row_label: str,
) -> None:
    plot_swim_frame_3d(
        ax,
        state,
        cfg,
        sim.rig,
        hide_ticks=False,
        title=f"{row_label}: source contribution",
        show_legend=False,
        show_status=False,
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


def _slice(ax: Any, hydro: Any, index: int, row_label: str) -> None:
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
        title=f"{row_label}: body-fixed RPY slice",
        xlabel="long axis [µm]",
        ylabel="radial axis [µm]",
        aspect="equal",
    )
    ax.tick_params(labelsize=6)


def render_group(
    root: Path,
    records: list[dict[str, Any]],
    *,
    row_axis: str,
    group_axes: tuple[str, ...],
    output_dir: Path,
    fps: float = 25.0,
) -> tuple[Path | None, dict[str, Any]]:
    valid = {}
    omitted = []
    for rec in records:
        row_value = rec["axis_values"][row_axis]
        d = _condition_dir(root, rec)
        passed, qc = qc_result(d)
        if not passed:
            omitted.append(
                {
                    "condition_id": rec["condition_id"],
                    "row_axis": row_axis,
                    "row_value": row_value,
                    "reason": "strict_shape_qc_not_passed",
                    "qc": qc,
                }
            )
            continue
        cfg = _load_cfg(root, rec)
        hydro = load_hydro_archive(d / "hydro_archive.npz")
        valid[row_value] = (
            hydro,
            _select_frames(load_state_archive(d / "state_archive.npz"), False, fps),
            cfg,
            Simulator(cfg),
        )
    group_ids = {axis: records[0]["axis_ids"][axis] for axis in group_axes}
    canonical_axes = tuple(
        sorted(
            group_axes, key=lambda axis: (axis != "n_flagella", group_axes.index(axis))
        )
    )
    group_stem = "__".join(group_ids[axis] for axis in canonical_axes)
    metadata = {
        "group_axes": dict(
            zip(
                group_axes,
                (records[0]["axis_values"][axis] for axis in group_axes),
                strict=True,
            )
        ),
        "group_ids": group_ids,
        "rows": [
            {
                "row_value": rec["axis_values"][row_axis],
                "condition_id": rec["condition_id"],
            }
            for rec in records
        ],
        "omitted": omitted,
    }
    if not valid:
        metadata["render_status"] = "all_rows_qc_failed"
        return None, metadata
    b = _bounds([v[0] for v in valid.values()])
    output_dir.mkdir(parents=True, exist_ok=True)
    movie = output_dir / f"{group_stem}_{row_axis}_flow.mp4"
    fig = plt.figure(figsize=(12, 12), dpi=80)
    canvas = FigureCanvasAgg(fig)
    writer = None
    count = max(len(v[1]) for v in valid.values())
    try:
        for fi in range(count):
            fig.clear()
            gs = fig.add_gridspec(len(records), 3)
            for row, rec in enumerate(records):
                row_value = rec["axis_values"][row_axis]
                row_label = f"{row_axis}={row_value}"
                if row_value not in valid:
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
                hydro, states, cfg, sim = valid[row_value]
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
                    title=f"{row_label}: world RPY flow",
                    show_legend=False,
                    show_status=False,
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
                    row_label,
                )
                _slice(fig.add_subplot(gs[row, 2]), hydro, idx, row_label)
            fig.suptitle(
                f"{group_stem}; t={fi / fps:.3f} s",
                fontsize=12,
            )
            canvas.draw()
            frame = cv2.cvtColor(np.asarray(canvas.buffer_rgba()), cv2.COLOR_RGBA2BGR)
            if writer is None:
                writer = open_mp4_writer(
                    movie, fps=fps, frame_size=(frame.shape[1], frame.shape[0])
                ).writer
            writer.write(frame)
            fig.clear()
            gc.collect()
    finally:
        if writer is not None:
            writer.release()
        plt.close(fig)
    metadata.update(
        {
            "movie": movie.name,
            "world_bounds_um": b,
            "frame_count": count,
            "render_status": "rendered",
        }
    )
    return movie, metadata


def render_grouped_flow_videos(
    root: Path,
    output_dir: Path,
    *,
    row_axis: str,
    group_axes: tuple[str, ...],
    fps: float = 25.0,
) -> list[Path]:
    results = [
        render_group(
            root,
            records,
            row_axis=row_axis,
            group_axes=group_axes,
            output_dir=output_dir,
            fps=fps,
        )
        for records in grouped_records(
            root, row_axis=row_axis, group_axes=group_axes
        ).values()
    ]
    payload = {
        "kind": "free_space_rpy_grouped_flow_visualization",
        "input_campaign": str(root),
        "input_provenance": "hydro_archive.npz positions_m + total_forces_N",
        "layout": {
            "row_axis": row_axis,
            "group_axes": list(group_axes),
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
    return [m for m, _ in results if m is not None]


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, required=True)
    p.add_argument("--row-axis", required=True)
    p.add_argument("--group-axis", action="append", required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--fps", type=float, default=25.0)
    a = p.parse_args(argv)
    for movie in render_grouped_flow_videos(
        a.run_dir,
        a.output_dir,
        row_axis=a.row_axis,
        group_axes=tuple(a.group_axis),
        fps=a.fps,
    ):
        print(movie)


if __name__ == "__main__":
    main()
