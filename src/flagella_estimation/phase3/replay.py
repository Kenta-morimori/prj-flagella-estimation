"""Replay helpers for Phase 3 common clip datasets."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg
import cv2
import numpy as np

from sim_swim.analysis.flagella_count_behavior import load_state_archive
from sim_swim.render.body2d import BodyCapsuleRenderConfig, render_body_capsule_frame
from sim_swim.render.video_writer import open_mp4_writer


@dataclass(frozen=True)
class ReplayConfig:
    dataset_dir: Path
    output_dir: Path | None = None
    n_flagella: int | None = None
    run_id: str | None = None
    clip_index: int | None = None
    time_band: str | None = None
    qc_label: str | None = None
    training_candidate: bool | None = None
    max_clips: int = 12
    frames_per_clip: int = 4
    mp4_fps: float = 5.0
    panel_width_px: int = 640
    panel_height_px: int = 260


def default_replay_output_dir(dataset_dir: Path) -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return dataset_dir / "replay" / now.strftime("%Y%m%d_%H%M%S")


def render_contact_sheet(cfg: ReplayConfig) -> Path:
    records = _load_metadata_jsonl(cfg.dataset_dir / "clip_metadata.jsonl")
    selected = [record for record in records if _matches(record, cfg)]
    selected = selected[: cfg.max_clips]
    if not selected:
        raise ValueError("no clips match replay filters")

    output_dir = cfg.output_dir or default_replay_output_dir(cfg.dataset_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(
        len(selected),
        cfg.frames_per_clip,
        figsize=(2.4 * cfg.frames_per_clip, 2.2 * len(selected)),
        squeeze=False,
    )
    for row_index, record in enumerate(selected):
        clip = record["clip"]
        qc = record["qc"]
        provenance = record["provenance"]
        labels = record["labels"]
        arr = np.load(_resolve_path(cfg.dataset_dir, clip["output_path"]))
        frame_indices = np.linspace(
            0, arr.shape[0] - 1, num=min(cfg.frames_per_clip, arr.shape[0]), dtype=int
        )
        title = (
            f"nf={labels['n_flagella']} run={provenance['run_id']} "
            f"clip={clip['clip_index']} {clip.get('time_band', '')} "
            f"qc={qc.get('qc_label', qc.get('status'))}"
        )
        for col_index in range(cfg.frames_per_clip):
            ax = axes[row_index][col_index]
            ax.axis("off")
            if col_index < len(frame_indices):
                ax.imshow(arr[frame_indices[col_index]], cmap="gray", vmin=0, vmax=255)
            if col_index == 0:
                ax.set_title(title, fontsize=8, loc="left")
    fig.tight_layout()
    sheet_path = output_dir / "contact_sheet.png"
    fig.savefig(sheet_path, dpi=150)
    plt.close(fig)

    manifest = {
        "pipeline_name": "phase3_clip_replay_contact_sheet",
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "input_dataset": str(cfg.dataset_dir),
        "output_dir": str(output_dir),
        "filters": {
            "n_flagella": cfg.n_flagella,
            "run_id": cfg.run_id,
            "clip_index": cfg.clip_index,
            "time_band": cfg.time_band,
            "qc_label": cfg.qc_label,
            "training_candidate": cfg.training_candidate,
            "max_clips": cfg.max_clips,
            "frames_per_clip": cfg.frames_per_clip,
        },
        "clip_count": len(selected),
        "clip_ids": [record["clip"]["clip_id"] for record in selected],
        "outputs": {"contact_sheet_png": str(sheet_path)},
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output_dir / "run.log").write_text(
        "\n".join(
            [
                f"created_at={manifest['created_at']}",
                f"input_dataset={cfg.dataset_dir}",
                f"output_dir={output_dir}",
                f"clip_count={len(selected)}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return output_dir


def render_3d_2d_grid_mp4(cfg: ReplayConfig) -> Path:
    """Render a human-review MP4 grid with source 3D and body-only 2D panels."""

    records = _select_records(cfg)
    output_dir = cfg.output_dir or default_replay_output_dir(cfg.dataset_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    clips = [_load_replay_clip(record, cfg) for record in records]
    if not clips:
        raise ValueError("no clips match replay filters")

    cols = 2 if len(clips) > 1 else 1
    rows = int(np.ceil(len(clips) / cols))
    cell_w = int(cfg.panel_width_px)
    cell_h = int(cfg.panel_height_px)
    frame_size = (cell_w * cols, cell_h * rows)
    video_path = output_dir / "3d_2d_grid.mp4"
    writer_selection = open_mp4_writer(
        video_path,
        fps=max(float(cfg.mp4_fps), 1.0e-6),
        frame_size=frame_size,
    )
    frame_count = max(len(clip["states"]) for clip in clips)

    for frame_index in range(frame_count):
        canvas = np.full((frame_size[1], frame_size[0], 3), 255, dtype=np.uint8)
        for clip_index, clip in enumerate(clips):
            states = clip["states"]
            state = states[min(frame_index, len(states) - 1)]
            row = clip_index // cols
            col = clip_index % cols
            y0 = row * cell_h
            x0 = col * cell_w
            panel = _render_clip_pair_panel(
                state,
                clip["record"],
                cell_w=cell_w,
                cell_h=cell_h,
                image_size_px=clip["image_size_px"],
                pixel_size_um=clip["pixel_size_um"],
            )
            canvas[y0 : y0 + cell_h, x0 : x0 + cell_w] = panel
        writer_selection.writer.write(canvas)
    writer_selection.writer.release()

    manifest = {
        "pipeline_name": "phase3_clip_replay_3d_2d_grid",
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "input_dataset": str(cfg.dataset_dir),
        "output_dir": str(output_dir),
        "filters": _filters_dict(cfg),
        "clip_count": len(clips),
        "clip_ids": [clip["record"]["clip"]["clip_id"] for clip in clips],
        "outputs": {"mp4_grid": str(video_path)},
        "video": {
            "fps": cfg.mp4_fps,
            "frame_count": frame_count,
            "frame_size": list(frame_size),
            "selected_codec": writer_selection.selected_codec,
            "attempted_codecs": writer_selection.attempted_codecs,
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (output_dir / "run.log").write_text(
        "\n".join(
            [
                f"created_at={manifest['created_at']}",
                f"input_dataset={cfg.dataset_dir}",
                f"output_dir={output_dir}",
                f"clip_count={len(clips)}",
                f"mp4_grid={video_path}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return output_dir


def _matches(record: dict[str, Any], cfg: ReplayConfig) -> bool:
    clip = record["clip"]
    qc = record["qc"]
    provenance = record["provenance"]
    labels = record["labels"]
    if cfg.n_flagella is not None and int(labels["n_flagella"]) != cfg.n_flagella:
        return False
    if cfg.run_id is not None and provenance.get("run_id") != cfg.run_id:
        return False
    if cfg.clip_index is not None and int(clip["clip_index"]) != cfg.clip_index:
        return False
    if cfg.time_band is not None and clip.get("time_band") != cfg.time_band:
        return False
    if cfg.qc_label is not None and qc.get("qc_label") != cfg.qc_label:
        return False
    if (
        cfg.training_candidate is not None
        and bool(qc.get("training_candidate", True)) != cfg.training_candidate
    ):
        return False
    return True


def _select_records(cfg: ReplayConfig) -> list[dict[str, Any]]:
    records = _load_metadata_jsonl(cfg.dataset_dir / "clip_metadata.jsonl")
    selected = [record for record in records if _matches(record, cfg)]
    selected = selected[: cfg.max_clips]
    if not selected:
        raise ValueError("no clips match replay filters")
    return selected


def _filters_dict(cfg: ReplayConfig) -> dict[str, Any]:
    return {
        "n_flagella": cfg.n_flagella,
        "run_id": cfg.run_id,
        "clip_index": cfg.clip_index,
        "time_band": cfg.time_band,
        "qc_label": cfg.qc_label,
        "training_candidate": cfg.training_candidate,
        "max_clips": cfg.max_clips,
        "frames_per_clip": cfg.frames_per_clip,
    }


def _load_replay_clip(record: dict[str, Any], cfg: ReplayConfig) -> dict[str, Any]:
    clip = record["clip"]
    source_path = _resolve_path(cfg.dataset_dir, record["source_video"]["source_path"])
    states = load_state_archive(source_path)
    start = int(clip["source_frame_start"])
    end_inclusive = int(clip["source_frame_end"])
    clip_states = states[start : end_inclusive + 1]
    if not clip_states:
        raise ValueError(f"empty source clip window: {clip['clip_id']}")
    if cfg.frames_per_clip > 0 and len(clip_states) > cfg.frames_per_clip:
        indices = np.linspace(
            0,
            len(clip_states) - 1,
            num=cfg.frames_per_clip,
            dtype=int,
        )
        clip_states = [clip_states[int(index)] for index in indices]
    normalization = record.get("normalization", {})
    crop_size = normalization.get("crop_size_px", [96, 96])
    image_size_px = int(crop_size[0] if isinstance(crop_size, list) else crop_size)
    pixel_size_um = float(normalization.get("pixel_size_um") or 0.1)
    return {
        "record": record,
        "states": clip_states,
        "image_size_px": image_size_px,
        "pixel_size_um": pixel_size_um,
    }


def _render_clip_pair_panel(
    state: Any,
    record: dict[str, Any],
    *,
    cell_w: int,
    cell_h: int,
    image_size_px: int,
    pixel_size_um: float,
) -> np.ndarray:
    label_h = 34
    body = record["labels"]
    clip = record["clip"]
    qc = record["qc"]
    provenance = record["provenance"]
    title = (
        f"nf={body['n_flagella']} run={provenance['run_id']} "
        f"clip={clip['clip_index']} {clip.get('time_band', '')} "
        f"qc={qc.get('qc_label', qc.get('status'))}"
    )
    panel = np.full((cell_h, cell_w, 3), 255, dtype=np.uint8)
    cv2.putText(
        panel,
        title[:115],
        (8, 22),
        cv2.FONT_HERSHEY_SIMPLEX,
        0.48,
        (35, 35, 35),
        1,
        cv2.LINE_AA,
    )
    sub_h = cell_h - label_h
    sub_w = cell_w // 2
    frame3d = _render_3d_frame(state, width=sub_w, height=sub_h)
    frame2d_gray, _ = render_body_capsule_frame(
        state,
        BodyCapsuleRenderConfig(
            image_size_px=image_size_px,
            pixel_size_um=pixel_size_um,
            tracking_center=True,
        ),
    )
    frame2d = cv2.cvtColor(frame2d_gray, cv2.COLOR_GRAY2BGR)
    frame2d = cv2.resize(
        frame2d, (cell_w - sub_w, sub_h), interpolation=cv2.INTER_NEAREST
    )
    panel[label_h:, :sub_w] = frame3d
    panel[label_h:, sub_w:] = frame2d
    cv2.line(panel, (sub_w, label_h), (sub_w, cell_h), (220, 220, 220), 1)
    return panel


def _render_3d_frame(state: Any, *, width: int, height: int) -> np.ndarray:
    fig = plt.figure(figsize=(width / 100.0, height / 100.0), dpi=100)
    ax = fig.add_subplot(111, projection="3d")
    beads = np.asarray(state.bead_positions_um, dtype=float)
    center = np.asarray(state.position_um, dtype=float)
    view_range = 4.0
    ax.set_xlim(center[0] - view_range, center[0] + view_range)
    ax.set_ylim(center[1] - view_range, center[1] + view_range)
    ax.set_zlim(center[2] - view_range, center[2] + view_range)
    ax.set_box_aspect((1, 1, 1))
    if beads.ndim == 2 and beads.shape[1] >= 3 and beads.size:
        ax.scatter(
            beads[:, 0],
            beads[:, 1],
            beads[:, 2],
            c="black",
            s=5,
            depthshade=False,
        )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xlabel("3D")
    ax.grid(True, alpha=0.25)
    fig.tight_layout(pad=0.1)
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())
    frame = cv2.cvtColor(buf, cv2.COLOR_RGBA2BGR)
    plt.close(fig)
    if frame.shape[:2] != (height, width):
        frame = cv2.resize(frame, (width, height), interpolation=cv2.INTER_AREA)
    return frame


def _load_metadata_jsonl(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            records.append(json.loads(line))
    return records


def _resolve_path(dataset_dir: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    for candidate in (path, dataset_dir / path):
        if candidate.is_file():
            return candidate
    return dataset_dir / path
