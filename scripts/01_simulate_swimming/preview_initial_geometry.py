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

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    geometry_preflight,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    columns = min(4, len(selected))
    rows = math.ceil(len(selected) / columns)
    fig = plt.figure(figsize=(4.3 * columns, 3.5 * rows), constrained_layout=True)
    for index, condition in enumerate(selected, start=1):
        cfg = SimulationConfig.from_dict(base).with_overrides(
            condition["config_overrides"]
        )
        model = ModelBuilder(cfg).build()
        positions = model.positions_m * 1e6
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
    return {
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "image": str(image_path),
        "image_sha256": _sha256(image_path),
        "condition_count": len(conditions),
        "unique_shape_count": len(selected),
        "preflight": preflight,
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
