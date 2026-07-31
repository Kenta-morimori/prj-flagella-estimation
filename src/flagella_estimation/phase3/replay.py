"""Replay helpers for Phase 3 common clip datasets."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

import matplotlib
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cv2
import numpy as np

from sim_swim.analysis.flagella_count_behavior import load_state_archive
from flagella_estimation.phase3.render import select_frames
from flagella_estimation.phase3.pipeline import _environment_info, _git_info
from sim_swim.render.body2d import BodyCapsuleRenderConfig, render_body_capsule_frame
from sim_swim.render.grid_movie import (
    auto_grid_layout,
    compose_grid_frame,
    write_mp4_grid,
)
from sim_swim.render.render3d import render_swim_frame_3d
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


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
    max_clips: int | None = None
    clips_per_video: int = 12
    frames_per_clip: int = 4
    panel_layout: str = "3d+2d"
    panel_width_px: int = 640
    panel_height_px: int = 260


def default_replay_output_dir(dataset_dir: Path) -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    base = dataset_dir / "replay" / now.strftime("%Y%m%d_%H%M%S")
    if not base.exists():
        return base
    for index in range(1, 1000):
        candidate = base.with_name(f"{base.name}_{index:03d}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(
        f"could not allocate unique replay output directory under {base}"
    )


def render_contact_sheet(cfg: ReplayConfig) -> Path:
    records = _load_metadata_jsonl(cfg.dataset_dir / "clip_metadata.jsonl")
    selected = [record for record in records if _matches(record, cfg)]
    if cfg.max_clips is not None:
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
        "git": _git_info(),
        "environment": _replay_environment_info(),
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
    if cfg.clips_per_video <= 0:
        raise ValueError("clips_per_video must be positive")
    output_dir = cfg.output_dir or default_replay_output_dir(cfg.dataset_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    dataset_render = _load_dataset_render(cfg.dataset_dir)

    record_batches = [
        records[start : start + cfg.clips_per_video]
        for start in range(0, len(records), cfg.clips_per_video)
    ]
    video_entries = []
    for batch_index, record_batch in enumerate(record_batches, start=1):
        clips = [
            _load_replay_clip(record, cfg, dataset_render) for record in record_batch
        ]
        layout = auto_grid_layout(
            len(clips),
            cell_width_px=cfg.panel_width_px,
            cell_height_px=cfg.panel_height_px,
            max_cols=2,
        )
        video_name = (
            "3d_2d_grid.mp4"
            if len(record_batches) == 1
            else f"3d_2d_grid_{batch_index:03d}.mp4"
        )
        video_path = output_dir / video_name
        frame_count = max(len(clip["states"]) for clip in clips)
        fps = _resolve_replay_fps(clips)

        def frames() -> Any:
            for frame_index in range(frame_count):
                panels = []
                for clip in clips:
                    states = clip["states"]
                    state = states[min(frame_index, len(states) - 1)]
                    panels.append(
                        _render_clip_pair_panel(
                            state,
                            clip["record"],
                            clip["sim_cfg"],
                            clip["rig"],
                            cell_w=layout.cell_width_px,
                            cell_h=layout.cell_height_px,
                            body_render_cfg=clip["body_render_cfg"],
                            panel_layout=cfg.panel_layout,
                        )
                    )
                yield compose_grid_frame(panels, layout)

        render_result = write_mp4_grid(video_path, frames=frames(), fps=fps)
        video_entries.append(
            {
                "batch_index": batch_index,
                "path": str(video_path),
                "clip_count": len(clips),
                "clip_ids": [clip["record"]["clip"]["clip_id"] for clip in clips],
                **render_result.to_manifest(),
            }
        )

    manifest = {
        "pipeline_name": "phase3_clip_replay_3d_2d_grid",
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "input_dataset": str(cfg.dataset_dir),
        "output_dir": str(output_dir),
        "filters": _filters_dict(cfg),
        "clip_count": len(records),
        "clip_ids": [record["clip"]["clip_id"] for record in records],
        "outputs": {
            "mp4_grid": video_entries[0]["path"],
            "mp4_grids": [entry["path"] for entry in video_entries],
        },
        "panel_layout": cfg.panel_layout,
        "video": {
            key: value
            for key, value in video_entries[0].items()
            if key not in {"batch_index", "path", "clip_count", "clip_ids"}
        },
        "video_count": len(video_entries),
        "videos": video_entries,
        "git": _git_info(),
        "environment": _replay_environment_info(),
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
                f"clip_count={len(records)}",
                f"video_count={len(video_entries)}",
                *[f"mp4_grid={entry['path']}" for entry in video_entries],
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
    if cfg.max_clips is not None:
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
        "clips_per_video": cfg.clips_per_video,
        "frames_per_clip": cfg.frames_per_clip,
        "panel_layout": cfg.panel_layout,
    }


def _load_replay_clip(
    record: dict[str, Any],
    cfg: ReplayConfig,
    dataset_render: dict[str, Any] | None = None,
) -> dict[str, Any]:
    clip = record["clip"]
    source_path = _resolve_path(cfg.dataset_dir, record["source_video"]["source_path"])
    frame_rate_hz = float(clip.get("frame_rate_hz") or 25.0)
    states = select_frames(load_state_archive(source_path), frame_rate_hz)
    start = int(clip["source_frame_start"])
    end_inclusive = int(clip["source_frame_end"])
    clip_states = states[start : end_inclusive + 1]
    if not clip_states:
        raise ValueError(f"empty source clip window: {clip['clip_id']}")
    normalization = record.get("normalization", {})
    crop_size = normalization.get("crop_size_px", [96, 96])
    image_size_px = int(crop_size[0] if isinstance(crop_size, list) else crop_size)
    pixel_size_um = float(normalization.get("pixel_size_um") or 0.1)
    render = dataset_render or _load_dataset_render(cfg.dataset_dir)
    body_render_cfg = BodyCapsuleRenderConfig(
        image_size_px=image_size_px,
        pixel_size_um=pixel_size_um,
        body_length_um=float(render.get("body_length_um", 2.0)),
        body_width_um=float(render.get("body_width_um", 1.0)),
        body_intensity=int(render.get("body_intensity", 60)),
        background_intensity=int(render.get("background_intensity", 255)),
        tracking_center=bool(render.get("tracking_center", True)),
    )
    sim_cfg = _load_simulation_config_for_record(record, cfg.dataset_dir)
    sim_cfg = sim_cfg.with_overrides(
        {
            "output_sampling": {
                "out_all_steps_3d": False,
                "fps_out_3d": frame_rate_hz,
                "fps_out_2d": frame_rate_hz,
            },
            "render": {
                "render_flagella": True,
                "label_flagella": True,
                "follow_camera_3d": True,
            },
        }
    )
    rig = Simulator(sim_cfg).rig
    return {
        "record": record,
        "states": clip_states,
        "body_render_cfg": body_render_cfg,
        "sim_cfg": sim_cfg,
        "rig": rig,
    }


def _render_clip_pair_panel(
    state: Any,
    record: dict[str, Any],
    sim_cfg: SimulationConfig,
    rig: Any,
    *,
    cell_w: int,
    cell_h: int,
    body_render_cfg: BodyCapsuleRenderConfig,
    panel_layout: str,
) -> np.ndarray:
    label_h = 34
    body = record["labels"]
    clip = record["clip"]
    qc = record["qc"]
    provenance = record["provenance"]
    title = (
        f"nf={body['n_flagella']} run={provenance['run_id']} "
        f"clip={clip['clip_index']} {clip.get('time_band', '')} "
        f"qc={_qc_display_label(qc)}"
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
    panels = _render_content_panels(
        state,
        record,
        sim_cfg,
        rig,
        panel_layout=panel_layout,
        total_width_px=cell_w,
        height_px=sub_h,
        body_render_cfg=body_render_cfg,
    )
    x0 = 0
    for content_index, content in enumerate(panels):
        panel[label_h:, x0 : x0 + content.shape[1]] = content
        x0 += content.shape[1]
        if content_index < len(panels) - 1:
            cv2.line(panel, (x0, label_h), (x0, cell_h), (220, 220, 220), 1)
    return panel


def _render_content_panels(
    state: Any,
    record: dict[str, Any],
    sim_cfg: SimulationConfig,
    rig: Any,
    *,
    panel_layout: str,
    total_width_px: int,
    height_px: int,
    body_render_cfg: BodyCapsuleRenderConfig,
) -> list[np.ndarray]:
    modes = panel_layout.split("+")
    if modes not in (["3d", "2d"], ["3d"], ["2d"]):
        raise ValueError("panel_layout must be one of: 3d+2d, 3d, 2d")
    widths = _split_width(total_width_px, len(modes))
    panels = []
    for mode, width in zip(modes, widths):
        if mode == "3d":
            panels.append(
                render_swim_frame_3d(
                    state,
                    sim_cfg,
                    rig,
                    width_px=width,
                    height_px=height_px,
                    extra_status_lines=_shape_status_lines(record, state),
                    hide_ticks=True,
                )
            )
        else:
            frame2d_gray, _ = render_body_capsule_frame(
                state,
                body_render_cfg,
            )
            frame2d = cv2.cvtColor(frame2d_gray, cv2.COLOR_GRAY2BGR)
            panels.append(
                cv2.resize(frame2d, (width, height_px), interpolation=cv2.INTER_NEAREST)
            )
    return panels


def _split_width(total_width: int, count: int) -> list[int]:
    base = total_width // count
    widths = [base for _ in range(count)]
    widths[-1] += total_width - sum(widths)
    return widths


def _shape_status_lines(record: dict[str, Any], state: Any) -> list[str]:
    qc = record["qc"]
    first_fail = qc.get("run_first_fail_t_s")
    if first_fail is None:
        return ["shape first-fail: none", f"shape status: {qc.get('qc_label')}"]
    first_fail_s = float(first_fail)
    relation = (
        "after first-fail" if float(state.t) >= first_fail_s else "before first-fail"
    )
    return [
        f"shape first-fail: t = {first_fail_s:.3f} s",
        f"shape status: {relation}",
    ]


def _qc_display_label(qc: dict[str, Any]) -> str:
    label = str(qc.get("qc_label", qc.get("status")))
    if label == "diagnostic":
        reason = str(qc.get("exclusion_reason") or "post_first_fail")
        return f"diagnostic {reason} training_excluded"
    return label


def _resolve_replay_fps(clips: list[dict[str, Any]]) -> float:
    return float(clips[0]["record"]["clip"].get("frame_rate_hz") or 25.0)


def _load_simulation_config_for_record(
    record: dict[str, Any], dataset_dir: Path
) -> SimulationConfig:
    manifest_path = dataset_dir / "manifest.json"
    input_dataset = None
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            input_dataset = Path(str(manifest.get("input_dataset", "")))
        except json.JSONDecodeError:
            input_dataset = None
    run_id = str(record["provenance"]["run_id"])
    candidates = []
    if input_dataset is not None:
        candidates.append(input_dataset / "configs" / f"{run_id}.yaml")
    raw_run_dir = record["provenance"].get("raw_run_dir")
    if raw_run_dir:
        raw_path = Path(str(raw_run_dir))
        candidates.extend(
            [
                raw_path / "config.yaml",
                raw_path.parent / "configs" / f"{run_id}.yaml",
                raw_path.parent.parent
                / "dataset"
                / "v1_r1_duration_3s"
                / "configs"
                / f"{run_id}.yaml",
            ]
        )
    for path in candidates:
        if path.is_file():
            raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
            return SimulationConfig.from_dict(raw)
    return SimulationConfig.from_dict(
        yaml.safe_load(Path("conf/sim_swim.yaml").read_text(encoding="utf-8")) or {}
    )


def _load_metadata_jsonl(path: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            records.append(json.loads(line))
    return records


def _replay_environment_info() -> dict[str, Any]:
    return {
        **_environment_info(),
        "opencv": cv2.__version__,
        "matplotlib": matplotlib.__version__,
    }


def _load_dataset_render(dataset_dir: Path) -> dict[str, Any]:
    manifest_path = dataset_dir / "manifest.json"
    if not manifest_path.is_file():
        return {}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    render = manifest.get("render", {})
    return render if isinstance(render, dict) else {}


def _resolve_path(dataset_dir: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path
    for candidate in (path, dataset_dir / path):
        if candidate.is_file():
            return candidate
    return dataset_dir / path
