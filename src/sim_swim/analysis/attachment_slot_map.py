"""Render the archived initial attachment topology in the body yz plane."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")

from sim_swim.analysis.flagella_count_behavior import load_state_archive
from sim_swim.render.render3d import _flagella_colors


def render_attachment_slot_map(replay_input: Path, output_dir: Path) -> dict[str, Path]:
    """Plot all archived initial rings; never infer attachment from a label alone."""

    import matplotlib.pyplot as plt

    manifest = json.loads(
        (replay_input / "run_manifest.json").read_text(encoding="utf-8")
    )
    records = list(manifest["conditions"])
    if not records:
        raise ValueError("attachment slot map requires completed conditions")
    output_dir.mkdir(parents=True, exist_ok=True)
    n_columns = 4
    n_rows = math.ceil(len(records) / n_columns)
    figure, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(22, 5.5 * n_rows),
        constrained_layout=True,
        squeeze=False,
    )
    mappings: list[dict[str, Any]] = []
    for axis, record in zip(axes.flat, records, strict=False):
        condition_id = str(record["condition_id"])
        actual = record["geometry"]["actual"]
        n_slots = int(actual["body_slots_per_layer"])
        n_layers = int(actual["body_layers"])
        if n_slots != 6 or n_layers % 2 != 1:
            raise ValueError(f"{condition_id}: expected hexagonal middle body layer")
        topology = list(actual["attachment_topology"])
        archive = Path(record["output_dir"]) / "state_archive.npz"
        if not archive.is_file():
            raise FileNotFoundError(archive)
        states = load_state_archive(archive)
        if not states:
            raise ValueError(f"{condition_id}: empty completed archive")
        positions = np.asarray(states[0].bead_positions_um, dtype=float)
        center_layer = n_layers // 2
        bead_ids = np.arange(center_layer * n_slots, (center_layer + 1) * n_slots)
        if (
            int(actual["body_beads"]) != n_layers * n_slots
            or positions.shape[0] < bead_ids[-1] + 1
        ):
            raise ValueError(f"{condition_id}: body topology does not match archive")
        ring = positions[bead_ids][:, 1:3]
        ring_center = ring.mean(axis=0)
        ring = ring - ring_center
        if not np.isfinite(ring).all():
            raise ValueError(f"{condition_id}: nonfinite archived ring")
        closed = np.vstack([ring, ring[0]])
        axis.plot(closed[:, 0], closed[:, 1], color="#394150", lw=1.5)
        axis.scatter(
            ring[:, 0], ring[:, 1], s=80, color="#d5d9e0", edgecolor="#394150", zorder=2
        )
        radius = max(float(np.linalg.norm(ring, axis=1).max()), 1e-6)
        for slot, point in enumerate(ring):
            label_position = point * 1.16
            axis.text(
                label_position[0],
                label_position[1],
                str(slot),
                ha="center",
                va="center",
                fontsize=10,
            )
        seen: set[int] = set()
        condition_map: list[dict[str, Any]] = []
        colors = _flagella_colors(len(topology))
        for flag_id, attachment in enumerate(topology):
            slot = int(attachment["slot"])
            bead = int(attachment["body_bead_index"])
            if (
                int(attachment["flag_id"]) != flag_id
                or int(attachment["layer"]) != center_layer
                or not 0 <= slot < n_slots
                or bead != int(bead_ids[slot])
                or slot in seen
            ):
                raise ValueError(
                    f"{condition_id}: archived F/slot/body bead mapping is inconsistent"
                )
            seen.add(slot)
            point = ring[slot]
            axis.scatter(
                [point[0]],
                [point[1]],
                s=165,
                color=[colors[flag_id]],
                edgecolor="black",
                zorder=3,
            )
            axis.annotate(
                f"F{flag_id}",
                xy=point,
                xytext=point * 1.53,
                color=colors[flag_id],
                fontsize=10,
                fontweight="bold",
                ha="center",
                va="center",
                arrowprops={"arrowstyle": "-", "color": colors[flag_id], "lw": 1.2},
            )
            condition_map.append(
                {
                    "flag_id": flag_id,
                    "slot": slot,
                    "body_bead_index": bead,
                    "color_rgb": list(colors[flag_id]),
                }
            )
        expected_slots = [
            int(value) for value in record["axis_values"]["attachment_slots"]
        ]
        if [item["slot"] for item in condition_map] != expected_slots:
            raise ValueError(
                f"{condition_id}: slot axis differs from archived geometry"
            )
        axis.set_title(
            condition_id + ("  (full ring)" if len(topology) == 6 else ""), fontsize=12
        )
        axis.set_xlabel("+y [µm] →")
        axis.set_ylabel("+z [µm] ↑")
        axis.set_aspect("equal")
        axis.set_xlim(-2.05 * radius, 2.05 * radius)
        axis.set_ylim(-2.05 * radius, 2.05 * radius)
        axis.grid(alpha=0.15)
        mappings.append(
            {
                "condition_id": condition_id,
                "mapping": condition_map,
                "full_ring_rotation_equivalent": len(topology) == 6,
                "archive": str(archive.resolve()),
                "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
                "initial_state_t_s": float(states[0].t),
            }
        )
    for axis in axes.flat[len(records) :]:
        axis.axis("off")
    figure.suptitle(
        "Initial attachment slots | rear view along +x | middle hexagonal body layer\n"
        "F colors match the existing 3D replay; slot labels are 0-based. n=6 occupies all slots.",
        fontsize=15,
    )
    image_path = output_dir / "attachment_slots_all_conditions.png"
    figure.savefig(image_path, dpi=220)
    plt.close(figure)
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "kind": "initial_attachment_slot_map",
                "view": "rear_to_front_along_positive_x",
                "horizontal_axis": "+y",
                "vertical_axis": "+z",
                "body_layer": "middle",
                "flag_color_source": "sim_swim.render.render3d._flagella_colors",
                "replay_input_manifest": str(
                    (replay_input / "run_manifest.json").resolve()
                ),
                "replay_input_manifest_sha256": hashlib.sha256(
                    (replay_input / "run_manifest.json").read_bytes()
                ).hexdigest(),
                "condition_count": len(mappings),
                "conditions": mappings,
                "image": str(image_path.resolve()),
                "image_sha256": hashlib.sha256(image_path.read_bytes()).hexdigest(),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "attachment_slot_map": image_path,
        "attachment_slot_map_manifest": manifest_path,
    }
