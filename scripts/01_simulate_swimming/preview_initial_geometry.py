#!/usr/bin/env python3
"""Render fixed-camera t=0 geometry from multi-run configs without simulation."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    geometry_preflight,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.model.builder import ModelBuilder
from sim_swim.render.render3d import _flagella_colors
from sim_swim.sim.params import SimulationConfig

_OVERVIEW_SLOTS = {
    1: (0,),
    2: (0, 3),
    3: (0, 2, 4),
    4: (0, 1, 3, 4),
    5: (0, 1, 2, 3, 4),
    6: (0, 1, 2, 3, 4, 5),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _plot_axial(ax, model, positions_over_b: np.ndarray, limit: float) -> None:
    """Project the body and flagella along the x body axis onto y–z."""
    center_layer = model.body_layer_indices[len(model.body_layer_indices) // 2]
    ring = positions_over_b[center_layer, 1:3]
    ax.plot(*np.vstack((ring, ring[0])).T, color="#9ca3af", linewidth=1.1)
    ax.scatter(*ring.T, s=18, color="#596579", zorder=3)
    colors = _flagella_colors(len(model.flagella_indices))
    for flag_id, (ids, attach_idx) in enumerate(
        zip(model.flagella_indices, model.flagella_attach_body_indices)
    ):
        flag = positions_over_b[ids, 1:3]
        attach = positions_over_b[int(attach_idx), 1:3]
        ax.plot(*flag.T, color=colors[flag_id], linewidth=2)
        ax.scatter(*flag.T, s=12, color=[colors[flag_id]], zorder=3)
        ax.plot(*np.vstack((attach, flag[0])).T, color="#111827", linewidth=2.4)
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), xlabel="y / b", ylabel="z / b")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(color="#e5e7eb", linewidth=0.5)
    ax.tick_params(labelsize=7)


def _plot_spatial(ax, model, positions_over_b: np.ndarray, limits: np.ndarray) -> None:
    for layer_index, ids in enumerate(model.body_layer_indices):
        closed = np.append(ids, ids[0])
        ax.plot(*positions_over_b[closed].T, color="#9ca3af", linewidth=0.75)
        if layer_index:
            previous = model.body_layer_indices[layer_index - 1]
            for before, after in zip(previous, ids):
                ax.plot(
                    *positions_over_b[[before, after]].T,
                    color="#9ca3af",
                    linewidth=0.65,
                )
    body = positions_over_b[model.body_indices]
    ax.scatter(*body.T, s=8, color="#596579", depthshade=False)
    colors = _flagella_colors(len(model.flagella_indices))
    for flag_id, (ids, attach_idx) in enumerate(
        zip(model.flagella_indices, model.flagella_attach_body_indices)
    ):
        flag = positions_over_b[ids]
        attach = positions_over_b[int(attach_idx)]
        ax.plot(*flag.T, color=colors[flag_id], linewidth=2)
        ax.scatter(*flag.T, s=10, color=[colors[flag_id]], depthshade=False)
        ax.plot(*np.vstack((attach, flag[0])).T, color="#111827", linewidth=2.4)
    ax.set(xlim=limits[0], ylim=limits[1], zlim=limits[2])
    ax.set_box_aspect(tuple(float(high - low) for low, high in limits))
    ax.view_init(elev=18, azim=-67)
    ax.set_proj_type("ortho")
    ax.set_xlabel("x / b", fontsize=8)
    ax.set_ylabel("y / b", fontsize=8)
    ax.set_zlabel("z / b", fontsize=8)
    ax.tick_params(labelsize=6)
    ax.grid(False)


def _render_axial_and_overview(items, output_dir: Path, stem: str) -> dict[str, object]:
    all_points = np.concatenate([positions for _, _, positions in items])
    axial_limit = float(np.max(np.abs(all_points[:, 1:3])) + 0.35)
    spatial_limits = np.column_stack(
        (all_points.min(axis=0) - 0.35, all_points.max(axis=0) + 0.35)
    )

    columns = min(4, len(items))
    rows = math.ceil(len(items) / columns)
    fig, axes = plt.subplots(
        rows, columns, figsize=(3.6 * columns, 3.55 * rows), squeeze=False
    )
    for ax, (condition, model, points) in zip(axes.flat, items):
        _plot_axial(ax, model, points, axial_limit)
        ax.set_title(condition["condition_id"].removesuffix("__rxaxis"), fontsize=9)
    for ax in axes.flat[len(items) :]:
        ax.axis("off")
    fig.suptitle("Initial geometry · axial projection at t=0", fontsize=13)
    fig.tight_layout()
    axial_path = output_dir / f"{stem}_axial_projection.png"
    fig.savefig(axial_path, dpi=200, facecolor="white")
    plt.close(fig)

    by_count = {}
    for item in items:
        n_flagella = int(item[0]["axis_values"]["n_flagella"])
        by_count.setdefault(n_flagella, []).append(item)
    representatives = []
    for n_flagella, candidates in sorted(by_count.items()):
        representatives.append(
            next(
                (
                    item
                    for item in candidates
                    if tuple(item[0]["axis_values"].get("attachment_slots", []))
                    == _OVERVIEW_SLOTS.get(n_flagella)
                ),
                candidates[0],
            )
        )
    overview_rows = math.ceil(len(representatives) / 2)
    figure = plt.figure(figsize=(14.5, 3.35 * overview_rows))
    for index, (condition, model, points) in enumerate(representatives):
        row, column = divmod(index, 2)
        spatial = figure.add_subplot(
            overview_rows, 4, row * 4 + column * 2 + 1, projection="3d"
        )
        axial = figure.add_subplot(overview_rows, 4, row * 4 + column * 2 + 2)
        _plot_spatial(spatial, model, points, spatial_limits)
        _plot_axial(axial, model, points, axial_limit)
        spatial.set_title(
            condition["condition_id"].removesuffix("__rxaxis") + " · 3D",
            fontsize=9,
        )
        axial.set_title("axial projection", fontsize=9)
    figure.suptitle("Initial geometry · t=0 (model preflight)", fontsize=13)
    figure.text(
        0.5,
        0.01,
        "gray: body · colors: flagella · black: hook",
        ha="center",
        fontsize=8,
    )
    figure.subplots_adjust(
        left=0.05, right=0.98, bottom=0.06, top=0.9, wspace=0.18, hspace=0.38
    )
    overview_path = output_dir / f"{stem}_all_counts_overview.png"
    figure.savefig(overview_path, dpi=200, facecolor="white")
    plt.close(figure)
    return {
        "axial_image": str(axial_path),
        "axial_image_sha256": _sha256(axial_path),
        "overview_image": str(overview_path),
        "overview_image_sha256": _sha256(overview_path),
        "overview_condition_ids": [item[0]["condition_id"] for item in representatives],
        "axial_projection": "body +x direction onto y-z, coordinates in b",
    }


def render_config(config_path: Path, output_dir: Path) -> dict[str, object]:
    campaign = normalize_campaign_config(load_yaml(config_path))
    conditions = build_campaign_conditions(campaign)
    preflight = geometry_preflight(campaign, conditions)
    base = load_yaml(Path(campaign["base_config"]))
    selected = [
        condition
        for condition in conditions
        if not condition["axis_values"].get("body_reaction_full_vector", False)
    ]
    items = []
    for condition in selected:
        cfg = SimulationConfig.from_dict(base).with_overrides(
            condition["config_overrides"]
        )
        model = ModelBuilder(cfg).build()
        items.append((condition, model, model.positions_m / cfg.b_m))
    columns = min(4, len(selected))
    rows = math.ceil(len(selected) / columns)
    fig = plt.figure(figsize=(4.3 * columns, 3.5 * rows), constrained_layout=True)
    for index, (condition, model, positions_over_b) in enumerate(items, start=1):
        positions = positions_over_b * float(base["scale"]["b_um"])
        ax = fig.add_subplot(rows, columns, index, projection="3d")
        body = positions[model.body_indices]
        ax.scatter(body[:, 0], body[:, 1], body[:, 2], s=15, c="#666666", alpha=0.8)
        for flag_id, (ids, attach) in enumerate(
            zip(model.flagella_indices, model.flagella_attach_body_indices)
        ):
            points = positions[ids]
            color = plt.cm.tab10(flag_id % 10)
            ax.plot(
                points[:, 0], points[:, 1], points[:, 2], color=color, linewidth=1.5
            )
            hook = positions[[int(attach), int(ids[0])]]
            ax.plot(hook[:, 0], hook[:, 1], hook[:, 2], color=color, linewidth=2.3)
            ax.scatter(points[0, 0], points[0, 1], points[0, 2], color=color, s=12)
        ax.set(xlim=(-7, 2), ylim=(-3, 3), zlim=(-3, 3))
        ax.set_box_aspect((9, 6, 6))
        ax.view_init(elev=20, azim=-65)
        ax.set_title(condition["condition_id"].removesuffix("__rxaxis"), fontsize=9)
        ax.set_xlabel("x (µm)", fontsize=7)
        ax.set_ylabel("y (µm)", fontsize=7)
        ax.set_zlabel("z (µm)", fontsize=7)
        ax.tick_params(labelsize=6)
    image_path = output_dir / f"{config_path.stem}_initial_geometry.png"
    fig.savefig(image_path, dpi=170)
    plt.close(fig)
    additional = _render_axial_and_overview(items, output_dir, config_path.stem)
    return {
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "image": str(image_path),
        "image_sha256": _sha256(image_path),
        "condition_count": len(conditions),
        "unique_shape_count": len(selected),
        "preflight": preflight,
        **additional,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    records = [render_config(path, args.output_dir) for path in args.config]
    manifest = {
        "kind": "initial_geometry_preview",
        "created_at_jst": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "simulation_executed": False,
        "camera": {
            "elev_deg": 20,
            "azim_deg": -65,
            "xlim_um": [-7, 2],
            "ylim_um": [-3, 3],
            "zlim_um": [-3, 3],
        },
        "axial_projection": "body +x direction onto y-z, coordinates in b",
        "campaigns": records,
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (args.output_dir / "run.log").write_text(
        f"Initial geometry preview only; {sum(r['unique_shape_count'] for r in records)} shapes; no simulation.\n",
        encoding="utf-8",
    )
    print(args.output_dir / "manifest.json")


if __name__ == "__main__":
    main()
