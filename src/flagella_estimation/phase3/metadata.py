"""Metadata builders for Phase 3 common clip records."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from flagella_estimation.phase3.render import FrameGeometry
from flagella_estimation.phase3.windows import FrameWindow


SCHEMA_VERSION = "phase3_clip_metadata/v0"
PIPELINE_NAME = "phase3_gt_passthrough"
PIPELINE_VERSION = "0.1.0"


def build_gt_passthrough_metadata(
    *,
    dataset_id: str,
    source_video_id: str,
    source_path: Path,
    frame_rate_hz: float,
    source_frame_count: int,
    source_duration_s: float,
    run_id: str,
    raw_run_dir: Path,
    n_flagella: int,
    track_id: str,
    group_key: str,
    clip_id: str,
    clip_index: int,
    window: FrameWindow,
    window_policy: str,
    output_path: Path,
    crop_size_px: int,
    pixel_size_um: float,
    frame_geometries: list[FrameGeometry],
) -> dict[str, Any]:
    """Build a #127-compatible metadata object for one pseudo-GT clip."""

    frame_count = window.frame_count
    t_start_s = window.start / frame_rate_hz
    t_end_s = window.end / frame_rate_hz
    frames = []
    for local_index, geometry in enumerate(frame_geometries):
        source_frame_index = window.start + local_index
        frames.append(
            {
                "frame_id": f"{source_video_id}:f{source_frame_index:06d}",
                "clip_frame_index": local_index,
                "source_frame_index": source_frame_index,
                "t_s": source_frame_index / frame_rate_hz,
                "bbox_xywh_px": list(geometry.bbox_xywh_px),
                "crop_xywh_px": list(geometry.crop_xywh_px),
                "center_xy_px": list(geometry.center_xy_px),
                "body_axis_angle_rad": geometry.body_axis_angle_rad,
                "body_length_px": geometry.body_length_px,
                "body_width_px": geometry.body_width_px,
                "detection_confidence": 1.0,
                "qc_status": "pass",
            }
        )

    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": dataset_id,
        "source_video": {
            "source_video_id": source_video_id,
            "source_kind": "phase2_pseudo",
            "source_path": str(source_path),
            "frame_rate_hz": frame_rate_hz,
            "width_px": crop_size_px,
            "height_px": crop_size_px,
            "duration_s": source_duration_s,
            "frame_count": source_frame_count,
            "codec_fourcc": None,
            "file_size_bytes": source_path.stat().st_size
            if source_path.exists()
            else None,
        },
        "processing_mode": "gt_passthrough",
        "provenance": {
            "pipeline_name": PIPELINE_NAME,
            "pipeline_version": PIPELINE_VERSION,
            "model_id": "phase2_flagella_count_behavior_v1",
            "dataset_version": "v1",
            "run_id": run_id,
            "render_id": "state_archive_numpy_v1",
            "condition_id": run_id,
            "raw_run_dir": str(raw_run_dir),
        },
        "track": {
            "track_id": track_id,
            "group_key": group_key,
            "source_track_id": track_id,
            "source_frame_start": 0,
            "source_frame_end": source_frame_count,
            "t_start_s": 0.0,
            "t_end_s": source_duration_s,
        },
        "clip": {
            "clip_id": clip_id,
            "clip_index": clip_index,
            "window_policy": window_policy,
            "output_path": str(output_path),
            "frame_count": frame_count,
            "frame_rate_hz": frame_rate_hz,
            "duration_s": frame_count / frame_rate_hz,
            "source_frame_start": window.start,
            "source_frame_end": window.end,
            "t_start_s": t_start_s,
            "t_end_s": t_end_s,
        },
        "normalization": {
            "crop_size_px": [crop_size_px, crop_size_px],
            "centering_mode": "body_center",
            "scale_mode": "fixed_um_per_px",
            "rotation_mode": "none",
            "pixel_size_um": pixel_size_um,
        },
        "frames": frames,
        "labels": {
            "n_flagella": n_flagella,
            "label_source": "phase2_gt",
        },
        "qc": {
            "status": "pass",
            "exclusion_reason": None,
            "detection_confidence_min": 1.0,
            "tracking_gap_count": 0,
            "notes": "pseudo GT passthrough from Phase 2 state archive",
        },
    }
