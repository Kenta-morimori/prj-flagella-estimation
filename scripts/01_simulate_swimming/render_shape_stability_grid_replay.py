#!/usr/bin/env python3
"""Render and/or plot Phase 2 shape-stability grid comparisons from sweep outputs."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
import logging
from pathlib import Path
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo

from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.pyplot as plt
import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.cli_profiles import (
    key_value_args_to_cli_args,
    split_config_key,
)
from sim_swim.sim.params import SimulationConfig

CANONICAL_TORQUE_DISTRIBUTION_CONDITION_IDS = (
    "segment_couples_diffusive_fp3_ft1p5",
    "segment_couples_uniform_fp3_ft1p5",
    "axis_projection_diffusive_fp3_ft1p5",
    "axis_projection_uniform_fp3_ft1p5",
)

METRIC_FIELDS = (
    "condition_id",
    "final_shape_pass_nonbody",
    "first_fail_t_s",
    "first_fail_category_nonbody",
    "hook_len_rel_err_max",
    "max_flag_bond_rel_err",
    "body_roll_net_abs_revolutions",
    "body_roll_direction_consistency",
    "axis_center_net_abs_revolutions_mean",
    "axis_center_direction_consistency_mean",
    "axis_center_body_relative_net_abs_revolutions_mean",
    "axis_center_body_relative_direction_consistency_mean",
    "axis_center_to_body_roll_ratio_mean",
    "body_shape_pass",
    "body_fail_category",
    "body_spring_max_stretch_ratio",
    "body_bend_max_error_deg",
    "body_centerline_max_deviation_um",
    "body_triangle_area_ratio_min",
)


def _default_output_dir() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / "phase2_replay"
        / "shape_stability_grid"
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
    )


def _load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _setup_logger(log_path: Path) -> logging.Logger:
    logger = logging.getLogger(f"shape_stability_grid_replay.{log_path.parent.name}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_path, encoding="utf-8")
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    return logger


def _git_info() -> dict[str, Any]:
    def run(cmd: list[str]) -> str:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True).strip()

    try:
        return {
            "commit": run(["git", "rev-parse", "HEAD"]),
            "commit_short": run(["git", "rev-parse", "--short", "HEAD"]),
            "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
            "is_clean": run(["git", "status", "--porcelain"]) == "",
        }
    except Exception:
        return {
            "commit": "unknown",
            "commit_short": "unknown",
            "branch": "unknown",
            "is_clean": False,
        }


def _short_distribution_label(value: str) -> str:
    mapping = {
        "root_torque_segment_couples": "segment_couples",
        "root_torque_axis_projection": "axis_projection",
    }
    return mapping.get(value, value)


def _display_condition_label(condition_label: str) -> str:
    replacements = {
        "as": "attach_seed",
        "ps": "phase_seed",
        "nf": "n_flagella",
    }
    parts = []
    for part in condition_label.split(","):
        stripped = part.strip()
        name, separator, value = stripped.partition("=")
        display_name = replacements.get(name, name)
        parts.append(f"{display_name}{separator}{value}" if separator else stripped)
    return ", ".join(parts)


def _label_for_row(row: dict[str, str]) -> str:
    condition_label = str(row.get("condition_label", "")).strip()
    if condition_label:
        return _display_condition_label(condition_label)
    condition_id = row.get("condition_id", "")
    if condition_id in CANONICAL_TORQUE_DISTRIBUTION_CONDITION_IDS:
        return (
            f"{_short_distribution_label(row['force_distribution'])}\n"
            f"{row['torque_distribution_profile']}"
        )
    tangent_mode = row.get("local_attach_frame_tangent_mode", "")
    if tangent_mode and tangent_mode != "vector":
        return f"{condition_id}\n{tangent_mode}"
    return condition_id


def _has_first_fail(row: dict[str, str]) -> bool:
    fail_t = str(row.get("first_fail_t_s", "")).strip().lower()
    return fail_t not in {"", "nan", "none"}


def _row_passes_nonbody(row: dict[str, str]) -> bool:
    return row.get("final_shape_pass_nonbody", "") == "True" and not _has_first_fail(
        row
    )


def _fail_label(row: dict[str, str]) -> str:
    fail_t = str(row.get("first_fail_t_s", ""))
    fail_c = str(row.get("first_fail_category_nonbody", ""))
    if _has_first_fail(row):
        return f"FAIL {fail_c}@{fail_t[:6]}"
    if _row_passes_nonbody(row):
        return "PASS"
    return f"FAIL {fail_c}@{fail_t[:6]}"


def _replay_status_lines(st: Any, cfg: Any, fail_label: str) -> list[str]:
    from sim_swim.render.render3d import _run_tumble_label

    return [
        _run_tumble_label(st, cfg),
        f"t = {st.t:.3f} s",
        f"motor torque / flag = {cfg.motor_torque_Nm:.2e} N m",
        fail_label,
    ]


def _auto_grid_shape(n_conditions: int) -> tuple[int, int]:
    if n_conditions <= 0:
        raise ValueError("n_conditions must be positive")
    n_cols = int(np.ceil(np.sqrt(n_conditions)))
    n_rows = int(np.ceil(n_conditions / n_cols))
    return n_rows, n_cols


def _grid_layout_for_rows(
    condition_rows: list[dict[str, str]],
) -> tuple[int, int, list[tuple[int, int]]]:
    row_indexes = [
        int(float(row["grid_row_index"]))
        for row in condition_rows
        if str(row.get("grid_row_index", "")).strip() != ""
    ]
    col_indexes = [
        int(float(row["grid_col_index"]))
        for row in condition_rows
        if str(row.get("grid_col_index", "")).strip() != ""
    ]
    if row_indexes and col_indexes:
        n_rows = max(row_indexes) + 1
        n_cols = max(col_indexes) + 1
        return (
            n_rows,
            n_cols,
            [
                (
                    int(float(row["grid_row_index"])),
                    int(float(row["grid_col_index"])),
                )
                for row in condition_rows
            ],
        )

    n_rows, n_cols = _auto_grid_shape(len(condition_rows))
    return (
        n_rows,
        n_cols,
        [divmod(plot_index, n_cols) for plot_index in range(len(condition_rows))],
    )


def _page_index_groups(n_conditions: int, max_panels_per_grid: int) -> list[list[int]]:
    if n_conditions <= 0:
        raise ValueError("n_conditions must be positive")
    max_panels = max(1, int(max_panels_per_grid))
    return [
        list(range(start, min(start + max_panels, n_conditions)))
        for start in range(0, n_conditions, max_panels)
    ]


def _float_or_nan(value: str | None) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")


def _condition_records(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    records = manifest.get("conditions", [])
    return {str(record["condition_id"]): dict(record) for record in records}


def _ordered_rows(
    rows: list[dict[str, str]],
    manifest: dict[str, Any] | None = None,
) -> list[dict[str, str]]:
    if manifest is not None:
        condition_order = [
            str(condition_id)
            for condition_id in manifest.get("condition_order", []) or []
        ]
        if condition_order:
            row_by_id = {row.get("condition_id", ""): row for row in rows}
            ordered = [
                row_by_id[condition_id]
                for condition_id in condition_order
                if condition_id in row_by_id
            ]
            if ordered:
                return ordered
    condition_ids = [row.get("condition_id", "") for row in rows]
    if all(
        condition_id in condition_ids
        for condition_id in CANONICAL_TORQUE_DISTRIBUTION_CONDITION_IDS
    ):
        return [
            row
            for condition_id in CANONICAL_TORQUE_DISTRIBUTION_CONDITION_IDS
            for row in rows
            if row.get("condition_id") == condition_id
        ]
    if manifest is not None:
        manifest_conditions = manifest.get("conditions", []) or []
        if manifest_conditions:
            row_by_id = {row.get("condition_id", ""): row for row in rows}
            ordered = [
                row_by_id[str(record.get("condition_id", ""))]
                for record in manifest_conditions
                if str(record.get("condition_id", "")) in row_by_id
            ]
            if ordered:
                return ordered
    return [row for row in rows if row.get("condition_id")]


def _load_inputs(
    input_dir: Path,
) -> tuple[list[dict[str, str]], dict[str, dict[str, Any]], Path]:
    summary_path = input_dir / "summary.csv"
    manifest_path = input_dir / "run_manifest.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing summary.csv under {input_dir}")
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing run_manifest.json under {input_dir}")
    rows = _load_csv_rows(summary_path)
    manifest = _load_json(manifest_path)
    records = _condition_records(manifest)
    base_cfg_path = Path(str(manifest.get("base_config") or manifest["config"]))
    ordered_rows = _ordered_rows(rows, manifest)
    if not ordered_rows:
        raise RuntimeError(f"No condition rows found in summary.csv under {input_dir}")
    for row in ordered_rows:
        condition_id = row["condition_id"]
        if condition_id not in records:
            raise RuntimeError(f"Missing {condition_id} in run_manifest.json")
        archive_path = _archive_path(input_dir, records[condition_id])
        if not archive_path.exists():
            raise FileNotFoundError(
                f"Missing state archive for {condition_id}: {archive_path}"
            )
    return ordered_rows, records, base_cfg_path


def _archive_path(input_dir: Path, condition_record: dict[str, Any]) -> Path:
    output_dir = condition_record.get("output_dir")
    if output_dir:
        return Path(str(output_dir)) / "state_archive.npz"
    return input_dir / str(condition_record["condition_id"]) / "state_archive.npz"


def _build_cfg(
    *,
    base_cfg_path: Path,
    condition_record: dict[str, Any],
    fps_out_3d: float,
) -> SimulationConfig:
    raw_cfg = _load_yaml(base_cfg_path)
    cfg = SimulationConfig.from_dict(raw_cfg).with_overrides(
        condition_record["config_overrides"]
    )
    return cfg.with_overrides(
        {
            "render": {
                "render_flagella": True,
                "render_flagella_2d": False,
                "show_flagella_helix_axis_3d": True,
                "label_flagella": False,
                "follow_camera_3d": True,
                "save_frames_3d": False,
                "save_frames_2d": False,
            },
            "output_sampling": {
                "out_all_steps_3d": False,
                "fps_out_3d": fps_out_3d,
            },
        }
    )


def _plot_cell(
    ax: plt.Axes,
    *,
    st: Any,
    cfg: Any,
    rig: Any,
    title: str,
    fail_label: str,
) -> None:
    from sim_swim.render.render3d import (
        _flagella_colors,
        _hook_edges,
        _plot_flagella_helix_axis_3d,
        _plot_segments_3d,
        _resolve_view_range_um,
    )

    ax.set_facecolor("white")
    beads = st.bead_positions_um
    view_range = _resolve_view_range_um(cfg, rig)
    center = np.array(st.position_um, dtype=float)
    ax.set_xlim(center[0] - view_range, center[0] + view_range)
    ax.set_ylim(center[1] - view_range, center[1] + view_range)
    ax.set_zlim(center[2] - view_range, center[2] + view_range)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.grid(True)
    ax.set_title(title, fontsize=9, pad=6)

    _plot_segments_3d(
        ax,
        beads,
        rig.body_spring_edges,
        color=(0.35, 0.35, 0.35),
        linewidth=1.3,
    )
    body_pts = beads[np.concatenate(rig.body_layer_indices)]
    ax.scatter(
        body_pts[:, 0],
        body_pts[:, 1],
        body_pts[:, 2],
        color="k",
        s=5,
        depthshade=False,
    )

    colors = _flagella_colors(len(rig.flagella_indices))
    if cfg.render.render_flagella:
        for f_id, idxs in enumerate(rig.flagella_indices):
            color = colors[f_id % len(colors)] if colors else (0.1, 0.4, 0.7)
            pts = beads[idxs]
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=color, linewidth=1.8)
            ax.scatter(
                pts[:, 0],
                pts[:, 1],
                pts[:, 2],
                color=[color],
                s=4,
                depthshade=False,
            )
            if cfg.render.show_flagella_helix_axis_3d:
                _plot_flagella_helix_axis_3d(ax, beads, idxs, f_id, color)

    if rig.hook_triplets.size > 0:
        _plot_segments_3d(
            ax,
            beads,
            _hook_edges(rig.hook_triplets),
            color=(1.0, 0.85, 0.05),
            linewidth=2.2,
        )

    ax.text2D(
        0.02,
        0.96,
        "\n".join(_replay_status_lines(st, cfg, fail_label)),
        transform=ax.transAxes,
        va="top",
        fontsize=8,
    )


def _render_grid_movie(
    *,
    states_by_condition: list[list[Any]],
    cfg_by_condition: list[SimulationConfig],
    rig_by_condition: list[Any],
    condition_rows: list[dict[str, str]],
    out_dir: Path,
    fps_out_3d: float,
    max_panels_per_grid: int,
) -> Any:
    import cv2

    from sim_swim.render.render3d import _select_frames
    from sim_swim.render.video_writer import VideoRenderResult, open_mp4_writer

    render_states_by_condition = [
        _select_frames(states, out_all_steps_3d=False, fps_hint=fps_out_3d)
        for states in states_by_condition
    ]
    frame_count = min(len(states) for states in render_states_by_condition)
    if frame_count <= 0:
        raise RuntimeError("No frames selected for grid render.")

    titles = [_label_for_row(row) for row in condition_rows]
    fail_labels = [_fail_label(row) for row in condition_rows]
    page_index_groups = _page_index_groups(
        len(condition_rows),
        max_panels_per_grid=max_panels_per_grid,
    )
    is_single_page = len(page_index_groups) == 1
    cleanup_patterns = ["grid_swim3d_page*.mp4", "grid_swim3d_page*_final.png"]
    if not is_single_page:
        cleanup_patterns.extend(["grid_swim3d.mp4", "grid_swim3d_final.png"])
    for pattern in cleanup_patterns:
        for stale_path in out_dir.glob(pattern):
            stale_path.unlink()
    page_manifests: list[dict[str, Any]] = []

    for page_number, page_indexes in enumerate(page_index_groups, start=1):
        page_rows = [condition_rows[index] for index in page_indexes]
        page_titles = [titles[index] for index in page_indexes]
        page_fail_labels = [fail_labels[index] for index in page_indexes]
        n_rows, n_cols, subplot_positions = _grid_layout_for_rows(page_rows)
        stem = "grid_swim3d" if is_single_page else f"grid_swim3d_page{page_number:02d}"
        movie_path = out_dir / f"{stem}.mp4"
        final_png_path = out_dir / f"{stem}_final.png"
        writer = None
        writer_selection = None
        last_frame = None

        for frame_idx in range(frame_count):
            fig = plt.figure(figsize=(4.8 * n_cols, 4.8 * n_rows))
            axes = {}
            for row_index, col_index in subplot_positions:
                axes[(row_index, col_index)] = fig.add_subplot(
                    n_rows,
                    n_cols,
                    row_index * n_cols + col_index + 1,
                    projection="3d",
                )
            for page_idx, (row_index, col_index) in enumerate(subplot_positions):
                condition_idx = page_indexes[page_idx]
                ax = axes[(row_index, col_index)]
                _plot_cell(
                    ax,
                    st=render_states_by_condition[condition_idx][frame_idx],
                    cfg=cfg_by_condition[condition_idx],
                    rig=rig_by_condition[condition_idx],
                    title=page_titles[page_idx],
                    fail_label=page_fail_labels[page_idx],
                )
            fig.tight_layout()
            canvas = FigureCanvasAgg(fig)
            canvas.draw()
            buf = np.asarray(canvas.buffer_rgba())
            frame = cv2.cvtColor(buf, cv2.COLOR_RGBA2BGR)
            plt.close(fig)

            if writer is None:
                writer_selection = open_mp4_writer(
                    movie_path,
                    fps=fps_out_3d,
                    frame_size=(frame.shape[1], frame.shape[0]),
                )
                writer = writer_selection.writer
            writer.write(frame)
            last_frame = frame

        if writer is not None:
            writer.release()
        if last_frame is not None:
            cv2.imwrite(str(final_png_path), last_frame)
        if writer_selection is None or last_frame is None:
            raise RuntimeError("Failed to initialize video writer.")
        render_result = VideoRenderResult(
            path=str(movie_path),
            selected_codec=writer_selection.selected_codec,
            attempted_codecs=writer_selection.attempted_codecs,
            fps=fps_out_3d,
            frame_size=(last_frame.shape[1], last_frame.shape[0]),
            frame_count=frame_count,
        )
        page_manifests.append(
            {
                "page": page_number,
                "conditions": [
                    condition_rows[index]["condition_id"] for index in page_indexes
                ],
                "final_png": str(final_png_path),
                "grid_shape": [n_rows, n_cols],
                "video": render_result.to_manifest(),
            }
        )

    return {
        "max_panels_per_grid": max(1, int(max_panels_per_grid)),
        "page_count": len(page_manifests),
        "pages": page_manifests,
    }


def _write_metrics(
    *,
    rows: list[dict[str, str]],
    out_dir: Path,
) -> Path:
    metrics_path = out_dir / "shape_stability_metrics.csv"
    with metrics_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METRIC_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in METRIC_FIELDS})
    return metrics_path


def _metric_series(rows: list[dict[str, str]], field: str) -> list[float]:
    return [_float_or_nan(row.get(field)) for row in rows]


def _axis_n_flagella_value(row: dict[str, str]) -> float:
    return _float_or_nan(row.get("axis_n_flagella_value"))


def _is_n_flagella_axis_campaign(rows: list[dict[str, str]]) -> bool:
    return bool(rows) and all(np.isfinite(_axis_n_flagella_value(row)) for row in rows)


def _seed_value(row: dict[str, str], name: str) -> float:
    value = _float_or_nan(row.get(f"axis_{name}_value") or row.get(name))
    return value if np.isfinite(value) else 0.0


def _seed_offsets(rows: list[dict[str, str]]) -> dict[tuple[float, float], float]:
    seed_keys = sorted(
        {
            (_seed_value(row, "attach_seed"), _seed_value(row, "phase_seed"))
            for row in rows
        }
    )
    if len(seed_keys) <= 1:
        return {key: 0.0 for key in seed_keys}
    offsets = np.linspace(-0.24, 0.24, len(seed_keys))
    return {key: float(offset) for key, offset in zip(seed_keys, offsets)}


def _plot_metrics_by_n_flagella(
    *,
    rows: list[dict[str, str]],
    out_path: Path,
) -> None:
    from matplotlib.lines import Line2D

    duration_s = max(_float_or_nan(row.get("duration_s")) for row in rows)
    metric_values = {
        "qc_time": [
            duration_s
            if _row_passes_nonbody(row)
            else _float_or_nan(row.get("first_fail_t_s"))
            for row in rows
        ],
        "max_flag_bond": [
            _float_or_nan(row.get("max_flag_bond_rel_err")) for row in rows
        ],
        "axis_lab": _metric_series(rows, "axis_center_net_abs_revolutions_mean"),
        "axis_body": _metric_series(
            rows, "axis_center_body_relative_net_abs_revolutions_mean"
        ),
        "body_roll": _metric_series(rows, "body_roll_net_abs_revolutions"),
        "axis_body_ratio": _metric_series(rows, "axis_center_to_body_roll_ratio_mean"),
    }
    panels = [
        (
            "QC duration",
            "Time to first QC failure or run end [s]",
            metric_values["qc_time"],
        ),
        (
            "Flagellar bond deformation",
            "Maximum flagellar bond relative error [-]",
            metric_values["max_flag_bond"],
        ),
        (
            "Helix-axis rotation in lab frame",
            "Mean rotation [revolutions]",
            metric_values["axis_lab"],
        ),
        (
            "Helix-axis rotation relative to body",
            "Mean body-relative rotation [revolutions]",
            metric_values["axis_body"],
        ),
        (
            "Body roll",
            "Body roll [revolutions]",
            metric_values["body_roll"],
        ),
        (
            "Flagella-to-body rotation ratio",
            "Mean helix-axis rotation / body roll [-]",
            metric_values["axis_body_ratio"],
        ),
    ]
    seed_offsets = _seed_offsets(rows)
    x_values = [
        _axis_n_flagella_value(row)
        + seed_offsets[
            (_seed_value(row, "attach_seed"), _seed_value(row, "phase_seed"))
        ]
        for row in rows
    ]
    n_values = sorted({_axis_n_flagella_value(row) for row in rows})
    pass_color = "#2f855a"
    fail_color = "#c05621"

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for panel_index, (ax, (title, y_label, values)) in enumerate(
        zip(axes.flat, panels)
    ):
        for row, x_value, value in zip(rows, x_values, values):
            if _row_passes_nonbody(row):
                ax.scatter(
                    x_value,
                    value,
                    color=pass_color,
                    marker="o",
                    s=28,
                    zorder=3,
                )
            else:
                ax.scatter(
                    x_value,
                    value,
                    color=fail_color,
                    marker="x",
                    s=48,
                    linewidths=1.8,
                    zorder=4,
                )
                if panel_index == 0:
                    attach_seed = _seed_value(row, "attach_seed")
                    phase_seed = _seed_value(row, "phase_seed")
                    ax.annotate(
                        f"attach_seed={attach_seed:g}, phase_seed={phase_seed:g}",
                        (x_value, value),
                        xytext=(-4, 6),
                        textcoords="offset points",
                        ha="right",
                        va="bottom",
                        fontsize=7,
                    )
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Number of flagella (n_flagella)")
        ax.set_ylabel(y_label)
        ax.set_xticks(n_values)
        ax.grid(axis="y", alpha=0.25)

    bond_ax = axes.flat[1]
    bond_ax.axhline(1.0, color="#555555", linestyle="--", linewidth=1.0)
    bond_ax.text(
        0.98,
        1.0,
        "QC limit = 1.0",
        transform=bond_ax.get_yaxis_transform(),
        ha="right",
        va="bottom",
        fontsize=8,
        color="#444444",
    )
    fig.suptitle("Shape stability across flagella counts and seed conditions")
    fig.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                color="none",
                markerfacecolor=pass_color,
                markeredgecolor=pass_color,
                label="PASS: no QC failure during run",
            ),
            Line2D(
                [0],
                [0],
                marker="x",
                color=fail_color,
                linestyle="none",
                markeredgewidth=1.8,
                label="FAIL: QC threshold exceeded at least once",
            ),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=2,
        frameon=False,
    )
    fig.text(
        0.5,
        0.015,
        "Each marker is one attach_seed x phase_seed condition; horizontal offsets separate seed combinations.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0.0, 0.05, 1.0, 0.91))
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def _plot_metrics_as_bars(
    *,
    rows: list[dict[str, str]],
    out_path: Path,
) -> None:
    labels = [_label_for_row(row) for row in rows]
    colors = ["#2f855a" if _row_passes_nonbody(row) else "#c05621" for row in rows]
    duration_s = max(_float_or_nan(row.get("duration_s")) for row in rows)
    first_fail = [
        duration_s
        if _row_passes_nonbody(row)
        else _float_or_nan(row.get("first_fail_t_s"))
        for row in rows
    ]
    panels = [
        ("first_fail_t_s_or_duration", first_fail),
        (
            "max_flag_bond_rel_err",
            [_float_or_nan(row.get("max_flag_bond_rel_err")) for row in rows],
        ),
        (
            "axis_center_net_abs_revolutions_mean",
            _metric_series(rows, "axis_center_net_abs_revolutions_mean"),
        ),
        (
            "axis_center_body_relative_net_abs_revolutions_mean",
            _metric_series(rows, "axis_center_body_relative_net_abs_revolutions_mean"),
        ),
        (
            "body_roll_net_abs_revolutions",
            _metric_series(rows, "body_roll_net_abs_revolutions"),
        ),
        (
            "axis_center_to_body_roll_ratio_mean",
            _metric_series(rows, "axis_center_to_body_roll_ratio_mean"),
        ),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    for ax, (title, values) in zip(axes.flat, panels):
        ax.bar(labels, values, color=colors)
        ax.set_title(title, fontsize=10)
        ax.tick_params(axis="x", labelrotation=15)
    fig.tight_layout()
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def _plot_metrics(
    *,
    rows: list[dict[str, str]],
    out_dir: Path,
) -> Path:
    out_path = out_dir / "shape_stability_metrics.png"
    if _is_n_flagella_axis_campaign(rows):
        _plot_metrics_by_n_flagella(rows=rows, out_path=out_path)
    else:
        _plot_metrics_as_bars(rows=rows, out_path=out_path)
    return out_path


def _replay_defaults(config_path: Path | None) -> dict[str, Any]:
    if config_path is None:
        return {}
    return dict((_load_yaml(config_path).get("replay") or {}))


def _config_run_dir(config_path: Path | None) -> Path | None:
    if config_path is None:
        return None
    output_cfg = dict((_load_yaml(config_path).get("output") or {}))
    if bool(output_cfg.get("timestamp_subdir", True)):
        return None
    base_dir = output_cfg.get("base_dir")
    return Path(str(base_dir)) if base_dir else None


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    config_from_key, parser_argv = split_config_key(raw_argv)
    parser_argv = key_value_args_to_cli_args(parser_argv)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--run-dir", type=Path, default=None)
    parser.add_argument("--input-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--mode",
        choices=("both", "plot-only", "render-only"),
        default=None,
    )
    parser.add_argument("--fps-out-3d", type=float, default=None)
    parser.add_argument("--max-panels-per-grid", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(parser_argv)
    if config_from_key is not None and args.config is not None:
        parser.error("Use either config=PATH or --config PATH (not both)")
    args.config = config_from_key or args.config

    if args.input_dir is not None and args.run_dir is not None:
        parser.error("Use either run_dir=PATH or input_dir=PATH (not both)")
    if args.input_dir is None and args.run_dir is None:
        args.run_dir = _config_run_dir(args.config)
    if args.input_dir is None:
        args.input_dir = args.run_dir
    if args.input_dir is None:
        parser.error(
            "--input-dir or --run-dir is required when output.timestamp_subdir is true"
        )

    replay_cfg = _replay_defaults(args.config)
    if args.mode is None:
        args.mode = str(replay_cfg.get("mode") or "both")
    if args.mode not in {"both", "plot-only", "render-only"}:
        parser.error(f"Invalid replay mode: {args.mode}")
    if args.fps_out_3d is None:
        args.fps_out_3d = float(replay_cfg.get("fps_out_3d") or 25.0)
    if args.max_panels_per_grid is None:
        args.max_panels_per_grid = int(replay_cfg.get("max_panels_per_grid") or 9)
    args.max_panels_per_grid = max(1, int(args.max_panels_per_grid))
    if args.output_dir is None:
        if args.run_dir is not None:
            output_subdir = str(replay_cfg.get("output_subdir") or "replay")
            args.output_dir = args.run_dir / output_subdir
        else:
            args.output_dir = _default_output_dir()
    return args


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    rows, records, base_cfg_path = _load_inputs(args.input_dir)
    if args.dry_run:
        for row in rows:
            print(row["condition_id"])
        return

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=args.overwrite)
    logger = _setup_logger(output_dir / "run.log")
    logger.info("Starting shape-stability replay render/plot")
    logger.info("input_dir=%s", args.input_dir)
    logger.info("mode=%s", args.mode)

    metrics_path = None
    metrics_plot_path = None
    render_result = None

    if args.mode in {"both", "plot-only"}:
        metrics_path = _write_metrics(rows=rows, out_dir=output_dir)
        metrics_plot_path = _plot_metrics(rows=rows, out_dir=output_dir)
        logger.info("Wrote metrics outputs: %s %s", metrics_path, metrics_plot_path)

    if args.mode in {"both", "render-only"}:
        from sim_swim.analysis.flagella_count_behavior import load_state_archive
        from sim_swim.sim.core import Simulator

        states_by_condition: list[list[Any]] = []
        cfg_by_condition: list[SimulationConfig] = []
        rig_by_condition: list[Any] = []
        for row in rows:
            condition_id = row["condition_id"]
            record = records[condition_id]
            cfg = _build_cfg(
                base_cfg_path=base_cfg_path,
                condition_record=record,
                fps_out_3d=args.fps_out_3d,
            )
            states = load_state_archive(_archive_path(args.input_dir, record))
            simulator = Simulator(cfg)
            states_by_condition.append(states)
            cfg_by_condition.append(cfg)
            rig_by_condition.append(simulator.rig)
        render_result = _render_grid_movie(
            states_by_condition=states_by_condition,
            cfg_by_condition=cfg_by_condition,
            rig_by_condition=rig_by_condition,
            condition_rows=rows,
            out_dir=output_dir,
            fps_out_3d=args.fps_out_3d,
            max_panels_per_grid=args.max_panels_per_grid,
        )
        logger.info(
            "Rendered grid videos: pages=%d",
            int(render_result["page_count"]),
        )

    manifest = {
        "git": _git_info(),
        "tool": "render_shape_stability_grid_replay",
        "input": {
            "input_dir": str(args.input_dir),
            "base_config": str(base_cfg_path),
            "mode": args.mode,
            "fps_out_3d": args.fps_out_3d,
            "max_panels_per_grid": args.max_panels_per_grid,
        },
        "conditions": [row["condition_id"] for row in rows],
        "outputs": {
            "root": str(output_dir),
            "metrics_csv": str(metrics_path) if metrics_path is not None else "",
            "metrics_png": str(metrics_plot_path)
            if metrics_plot_path is not None
            else "",
            "render_log": str(output_dir / "run.log"),
        },
    }
    if render_result is not None:
        manifest["render_video"] = {"grid_swim3d": render_result}
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(output_dir)


if __name__ == "__main__":
    main()
