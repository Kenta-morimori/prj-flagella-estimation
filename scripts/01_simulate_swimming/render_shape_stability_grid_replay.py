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


def _label_for_row(row: dict[str, str]) -> str:
    condition_label = str(row.get("condition_label", "")).strip()
    if condition_label:
        return condition_label
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


def _fail_label(row: dict[str, str]) -> str:
    if row.get("final_shape_pass_nonbody", "") == "True":
        return "PASS"
    fail_t = row.get("first_fail_t_s", "")
    fail_c = row.get("first_fail_category_nonbody", "")
    return f"FAIL {fail_c}@{fail_t[:6]}"


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
        _run_tumble_label,
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
        "\n".join([_run_tumble_label(st, cfg), f"t = {st.t:.3f} s", fail_label]),
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
) -> Any:
    import cv2

    from sim_swim.render.render3d import _select_frames
    from sim_swim.render.video_writer import open_mp4_writer

    render_states_by_condition = [
        _select_frames(states, out_all_steps_3d=False, fps_hint=fps_out_3d)
        for states in states_by_condition
    ]
    frame_count = min(len(states) for states in render_states_by_condition)
    if frame_count <= 0:
        raise RuntimeError("No frames selected for grid render.")

    movie_path = out_dir / "grid_swim3d.mp4"
    writer = None
    writer_selection = None
    last_frame = None
    titles = [_label_for_row(row) for row in condition_rows]
    fail_labels = [_fail_label(row) for row in condition_rows]
    n_conditions = len(condition_rows)
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
        subplot_positions = [
            (
                int(float(row["grid_row_index"])),
                int(float(row["grid_col_index"])),
            )
            for row in condition_rows
        ]
    else:
        n_cols = min(3, max(1, int(np.ceil(np.sqrt(n_conditions)))))
        n_rows = int(np.ceil(n_conditions / n_cols))
        subplot_positions = [
            divmod(plot_index, n_cols) for plot_index in range(n_conditions)
        ]

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
        for idx, (row_index, col_index) in enumerate(subplot_positions):
            ax = axes[(row_index, col_index)]
            _plot_cell(
                ax,
                st=render_states_by_condition[idx][frame_idx],
                cfg=cfg_by_condition[idx],
                rig=rig_by_condition[idx],
                title=titles[idx],
                fail_label=fail_labels[idx],
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

    writer.release()
    if last_frame is not None:
        cv2.imwrite(str(out_dir / "grid_swim3d_final.png"), last_frame)
    if writer_selection is None or last_frame is None:
        raise RuntimeError("Failed to initialize video writer.")
    from sim_swim.render.video_writer import VideoRenderResult

    return VideoRenderResult(
        path=str(movie_path),
        selected_codec=writer_selection.selected_codec,
        attempted_codecs=writer_selection.attempted_codecs,
        fps=fps_out_3d,
        frame_size=(last_frame.shape[1], last_frame.shape[0]),
        frame_count=frame_count,
    )


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


def _plot_metrics(
    *,
    rows: list[dict[str, str]],
    out_dir: Path,
) -> Path:
    labels = [_label_for_row(row) for row in rows]
    colors = [
        "#2f855a" if row.get("final_shape_pass_nonbody", "") == "True" else "#c05621"
        for row in rows
    ]
    duration_s = max(_float_or_nan(row.get("duration_s")) for row in rows)
    first_fail = [
        duration_s
        if row.get("final_shape_pass_nonbody", "") == "True"
        else _float_or_nan(row.get("first_fail_t_s"))
        for row in rows
    ]
    max_flag_bond = [_float_or_nan(row.get("max_flag_bond_rel_err")) for row in rows]
    axis_rev = _metric_series(rows, "axis_center_net_abs_revolutions_mean")
    axis_relative_rev = _metric_series(
        rows, "axis_center_body_relative_net_abs_revolutions_mean"
    )
    body_roll_rev = _metric_series(rows, "body_roll_net_abs_revolutions")
    axis_body_ratio = _metric_series(rows, "axis_center_to_body_roll_ratio_mean")

    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    panels = [
        ("first_fail_t_s_or_duration", first_fail),
        ("max_flag_bond_rel_err", max_flag_bond),
        ("axis_center_net_abs_revolutions_mean", axis_rev),
        ("axis_center_body_relative_net_abs_revolutions_mean", axis_relative_rev),
        ("body_roll_net_abs_revolutions", body_roll_rev),
        ("axis_center_to_body_roll_ratio_mean", axis_body_ratio),
    ]
    for ax, (title, values) in zip(axes.flat, panels):
        ax.bar(labels, values, color=colors)
        ax.set_title(title, fontsize=10)
        ax.tick_params(axis="x", labelrotation=15)
    fig.tight_layout()
    out_path = out_dir / "shape_stability_metrics.png"
    fig.savefig(out_path, dpi=160)
    plt.close(fig)
    return out_path


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=_default_output_dir())
    parser.add_argument(
        "--mode",
        choices=("both", "plot-only", "render-only"),
        default="both",
    )
    parser.add_argument("--fps-out-3d", type=float, default=25.0)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args(argv)


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
        )
        logger.info("Rendered video: %s", render_result.path)

    manifest = {
        "git": _git_info(),
        "tool": "render_shape_stability_grid_replay",
        "input": {
            "input_dir": str(args.input_dir),
            "base_config": str(base_cfg_path),
            "mode": args.mode,
            "fps_out_3d": args.fps_out_3d,
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
        manifest["render_video"] = {"grid_swim3d": render_result.to_manifest()}
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    print(output_dir)


if __name__ == "__main__":
    main()
