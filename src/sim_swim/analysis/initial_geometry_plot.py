"""Plot archived t=0 hex attachment geometry for common model evaluation."""

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


_BALANCED_SLOTS = {
    1: (0,),
    2: (0, 3),
    3: (0, 2, 4),
    4: (0, 1, 3, 4),
    5: (0, 1, 2, 3, 4),
    6: (0, 1, 2, 3, 4, 5),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _initial_records(
    replay_input: Path, b_um: float
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    manifest = json.loads(
        (replay_input / "run_manifest.json").read_text(encoding="utf-8")
    )
    records = list(manifest["conditions"])
    if not records or not math.isfinite(b_um) or b_um <= 0.0:
        raise ValueError(
            "Initial geometry requires completed conditions and positive b_um"
        )
    unique: dict[tuple[int, tuple[int, ...], str], dict[str, Any]] = {}
    sources: list[dict[str, Any]] = []
    for record in records:
        condition_id = str(record["condition_id"])
        actual = dict(record["geometry"]["actual"])
        if (
            int(actual["body_slots_per_layer"]) != 6
            or int(actual["body_layers"]) % 2 != 1
        ):
            raise ValueError(
                f"{condition_id}: expected hexagonal body with a middle layer"
            )
        slots = tuple(int(value) for value in record["axis_values"]["attachment_slots"])
        topology = list(actual["attachment_topology"])
        if slots != tuple(int(item["slot"]) for item in topology):
            raise ValueError(
                f"{condition_id}: archived attachment slots differ from the axis"
            )
        middle = int(actual["body_layers"]) // 2
        for flag_id, item in enumerate(topology):
            slot = int(item["slot"])
            if (
                int(item["flag_id"]) != flag_id
                or int(item["layer"]) != middle
                or not 0 <= slot < 6
                or int(item["body_bead_index"]) != middle * 6 + slot
            ):
                raise ValueError(
                    f"{condition_id}: invalid archived attachment topology"
                )
        archive = Path(str(record["output_dir"])) / "state_archive.npz"
        if not archive.is_file():
            raise FileNotFoundError(archive)
        archive_sha256 = _sha256(archive)
        expected_sha256 = dict(record.get("artifact_sha256", {}) or {}).get(
            "state_archive.npz"
        )
        if expected_sha256 is not None and archive_sha256 != str(expected_sha256):
            raise ValueError(f"{condition_id}: state_archive.npz SHA-256 mismatch")
        states = load_state_archive(archive)
        if not states or abs(float(states[0].t)) > 1e-12:
            raise ValueError(f"{condition_id}: completed archive has no t=0 state")
        positions = np.asarray(states[0].bead_positions_um, dtype=float) / b_um
        body_count = int(actual["body_beads"])
        flag_counts = [int(value) for value in actual["flagellum_beads"]]
        if (
            not np.isfinite(positions).all()
            or positions.shape != (body_count + sum(flag_counts), 3)
            or len(slots) != len(flag_counts)
            or len(slots) != int(record["axis_values"]["n_flagella"])
        ):
            raise ValueError(f"{condition_id}: invalid archived t=0 bead geometry")
        shape_hash = hashlib.sha256(positions.astype("<f8").tobytes()).hexdigest()
        key = (len(slots), slots, shape_hash)
        source = {
            "condition_id": condition_id,
            "archive": str(archive.resolve()),
            "archive_sha256": archive_sha256,
            "t_s": float(states[0].t),
            "geometry_sha256": shape_hash,
        }
        sources.append(source)
        if key not in unique:
            unique[key] = {
                "n_flagella": len(slots),
                "attachment_slots": list(slots),
                "geometry_sha256": shape_hash,
                "condition_ids": [],
                "positions": positions,
                "body_count": body_count,
                "flag_counts": flag_counts,
                "topology": topology,
            }
        unique[key]["condition_ids"].append(condition_id)
    panels = sorted(
        unique.values(),
        key=lambda panel: (
            panel["n_flagella"],
            tuple(panel["attachment_slots"]),
            min(panel["condition_ids"]),
        ),
    )
    return panels, sources


def _plot_panel(
    axis_3d: Any, axis_axial: Any, panel: dict[str, Any], limits: np.ndarray
) -> None:
    points = panel["positions"]
    n_body = panel["body_count"]
    body = points[:n_body]
    n_layers = n_body // 6
    for layer in range(n_layers):
        ids = np.arange(layer * 6, (layer + 1) * 6)
        closed = np.append(ids, ids[0])
        axis_3d.plot(*points[closed].T, color="#9ca3af", lw=0.75, alpha=0.6)
        if layer:
            for slot in range(6):
                axis_3d.plot(
                    *points[[ids[slot] - 6, ids[slot]]].T,
                    color="#9ca3af",
                    lw=0.65,
                    alpha=0.5,
                )
    axis_3d.scatter(*body.T, s=7, color="#596579", depthshade=False)
    ring = points[(n_layers // 2) * 6 : (n_layers // 2 + 1) * 6, 1:3]
    axis_axial.plot(*np.vstack([ring, ring[0]]).T, color="#9ca3af", lw=1.1)
    axis_axial.scatter(*ring.T, s=16, color="#596579", zorder=3)
    colors = _flagella_colors(panel["n_flagella"])
    offset = n_body
    for flag_id, count in enumerate(panel["flag_counts"]):
        flag = points[offset : offset + count]
        attach = points[int(panel["topology"][flag_id]["body_bead_index"])]
        color = colors[flag_id]
        axis_3d.plot(*flag.T, color=color, lw=2)
        axis_3d.scatter(*flag.T, s=10, color=[color], depthshade=False)
        axis_3d.plot(*np.vstack([attach, flag[0]]).T, color="#111827", lw=2.4)
        axis_axial.plot(*flag[:, 1:3].T, color=color, lw=2)
        axis_axial.scatter(*flag[:, 1:3].T, s=12, color=[color], zorder=3)
        axis_axial.plot(
            *np.vstack([attach[1:3], flag[0, 1:3]]).T, color="#111827", lw=2.4
        )
        offset += count
    axis_3d.set(xlim=limits[0], ylim=limits[1], zlim=limits[2])
    axis_3d.set_box_aspect(tuple(float(high - low) for low, high in limits))
    axis_3d.view_init(elev=18, azim=-67)
    axis_3d.set_proj_type("ortho")
    axis_3d.set_xlabel("x / b", fontsize=8)
    axis_3d.set_ylabel("y / b", fontsize=8)
    axis_3d.set_zlabel("z / b", fontsize=8)
    axis_3d.tick_params(labelsize=6)
    axis_3d.grid(False)
    axis_axial.set(xlim=limits[1], ylim=limits[2], xlabel="y / b", ylabel="z / b")
    axis_axial.set_aspect("equal", adjustable="box")
    axis_axial.tick_params(labelsize=7)
    axis_axial.grid(color="#e5e7eb", lw=0.5)


def _save_grid(
    panels: list[dict[str, Any]], path: Path, limits: np.ndarray, *, overview: bool
) -> None:
    import matplotlib.pyplot as plt

    columns = 2 if overview else len(panels)
    rows = math.ceil(len(panels) / columns) if overview else 1
    figure = plt.figure(
        figsize=(14.5, 3.35 * rows) if overview else (max(5.8, 3.55 * columns), 5.8)
    )
    for index, panel in enumerate(panels):
        row = index // columns
        column = index % columns
        side = (
            figure.add_subplot(rows, 4, row * 4 + column * 2 + 1, projection="3d")
            if overview
            else figure.add_subplot(2, columns, column + 1, projection="3d")
        )
        axial = (
            figure.add_subplot(rows, 4, row * 4 + column * 2 + 2)
            if overview
            else figure.add_subplot(2, columns, columns + column + 1)
        )
        _plot_panel(side, axial, panel, limits)
        label = f"n={panel['n_flagella']} · slots={''.join(map(str, panel['attachment_slots']))}"
        side.set_title(label + " · 3D", fontsize=9)
        axial.set_title("axial projection", fontsize=9)
    figure.suptitle("Archived initial geometry · t=0", fontsize=13)
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
    figure.savefig(path, dpi=200, facecolor="white")
    plt.close(figure)


def render_initial_geometry(
    replay_input: Path, output_dir: Path, *, b_um: float
) -> dict[str, Path]:
    """Write standard, archive-derived overview and per-count initial PNGs."""

    panels, sources = _initial_records(replay_input, b_um)
    all_points = np.concatenate([panel["positions"] for panel in panels])
    limits = np.column_stack(
        [all_points.min(axis=0) - 0.35, all_points.max(axis=0) + 0.35]
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    by_count: dict[int, list[dict[str, Any]]] = {}
    for panel in panels:
        by_count.setdefault(panel["n_flagella"], []).append(panel)
    representatives = [
        next(
            (
                panel
                for panel in by_count[n]
                if tuple(panel["attachment_slots"]) == _BALANCED_SLOTS.get(n)
            ),
            by_count[n][0],
        )
        for n in sorted(by_count)
    ]
    overview_path = output_dir / "all_counts_overview.png"
    _save_grid(representatives, overview_path, limits, overview=True)
    outputs = {"initial_geometry_overview": overview_path}
    panel_files: list[dict[str, Any]] = []
    for n, candidates in sorted(by_count.items()):
        for page_index, start in enumerate(range(0, len(candidates), 4), start=1):
            page = candidates[start : start + 4]
            path = output_dir / f"nf{n:02d}_page{page_index:02d}.png"
            _save_grid(page, path, limits, overview=False)
            outputs[f"initial_geometry_nf{n:02d}_page{page_index:02d}"] = path
            for panel in page:
                panel_files.append(
                    {
                        "condition_ids": panel["condition_ids"],
                        "n_flagella": n,
                        "attachment_slots": panel["attachment_slots"],
                        "geometry_sha256": panel["geometry_sha256"],
                        "image": str(path.resolve()),
                    }
                )
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "kind": "archived_initial_geometry",
                "t_s": 0.0,
                "b_um": b_um,
                "camera_3d": {"elevation_deg": 18, "azimuth_deg": -67},
                "condition_count": len(sources),
                "panel_count": len(panels),
                "sources": sources,
                "panels": panel_files,
                "images": {path.name: _sha256(path) for path in outputs.values()},
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    outputs["initial_geometry_manifest"] = manifest_path
    return outputs
