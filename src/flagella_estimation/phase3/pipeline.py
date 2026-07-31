"""Phase 3 pseudo-GT passthrough pipeline."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo

import numpy as np
import yaml

from flagella_estimation.phase3.metadata import (
    SCHEMA_VERSION,
    build_gt_passthrough_metadata,
)
from flagella_estimation.phase3.render import (
    RENDER_MODE,
    render_clip_array,
    select_frames,
)
from flagella_estimation.phase3.splits import (
    assign_grouped_splits,
    assert_no_group_leakage,
)
from flagella_estimation.phase3.windows import generate_windows
from sim_swim.analysis.flagella_count_behavior import load_state_archive
from sim_swim.render.body2d import BodyCapsuleRenderConfig, body_capsule_render_id


@dataclass(frozen=True)
class Phase3Config:
    dataset_id: str
    input_dataset: Path
    output_dir: Path
    config_path: Path | None = None
    cli_overrides: tuple[str, ...] = ()
    duration_s: float = 0.5
    window_policy: str = "non_overlap"
    overlap_stride_fraction: float = 0.5
    fps_out: float = 25.0
    image_size_px: int = 96
    pixel_size_um: float = 0.1
    body_length_um: float = 2.0
    body_width_um: float = 1.0
    body_intensity: int = 60
    allowed_n_flagella: tuple[int, ...] = (1, 2, 3)
    max_per_class: int | None = None
    baseline_torque_Nm: float = 2.0e-20
    require_use_for_ml_candidate: bool = True
    source_require_use_for_ml_candidate: bool = True
    dataset_version: str = "v1"
    dataset_revision: str | None = None
    output_name: str = "phase3_gt_passthrough_v1"


def _now_jst() -> datetime:
    return datetime.now(ZoneInfo("Asia/Tokyo"))


def default_output_dir() -> Path:
    return default_named_output_dir("phase3_gt_passthrough_v1")


def default_named_output_dir(name: str) -> Path:
    now = _now_jst()
    return Path("outputs") / now.strftime("%Y-%m-%d") / now.strftime("%H%M%S") / name


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _to_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _to_float(value: Any) -> float:
    return float(str(value).strip())


def _optional_float(value: Any) -> float | None:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    if np.isnan(parsed):
        return None
    return parsed


def _window_qc(
    *,
    window_end_s: float,
    run_shape_pass: bool,
    run_first_fail_t_s: float | None,
) -> tuple[bool, bool, bool, bool, str | None, str]:
    if run_first_fail_t_s is None:
        if run_shape_pass:
            return True, True, True, False, None, "pre_first_fail"
        return False, False, False, True, "run_shape_fail", "run_fail"
    is_pre_first_fail = window_end_s <= run_first_fail_t_s + 1.0e-12
    if is_pre_first_fail:
        return True, True, True, False, None, "pre_first_fail"
    return False, False, False, True, "post_or_contains_first_fail", "diagnostic"


def _load_source_manifest(input_dataset: Path) -> dict[str, Any]:
    path = input_dataset / "dataset_manifest.json"
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def _source_metadata(source_manifest: dict[str, Any]) -> dict[str, Any]:
    effective_campaign = source_manifest.get("effective_campaign_config")
    if not isinstance(effective_campaign, dict):
        return {}
    metadata = effective_campaign.get("metadata")
    return metadata if isinstance(metadata, dict) else {}


def _dataset_summary(
    metadata_records: list[dict[str, Any]], split_rows: list[dict[str, Any]]
) -> dict[str, Any]:
    groups = {record["track"]["group_key"] for record in metadata_records}
    classes = sorted(
        {int(record["labels"]["n_flagella"]) for record in metadata_records}
    )
    qc_labels = sorted({str(record["qc"]["qc_label"]) for record in metadata_records})
    split_groups: dict[str, set[str]] = {}
    for row in split_rows:
        split_groups.setdefault(str(row["split"]), set()).add(str(row["group_key"]))
    return {
        "clip_count": len(metadata_records),
        "group_count": len(groups),
        "class_count": len(classes),
        "classes": json.dumps(classes),
        "qc_labels": json.dumps(qc_labels),
        "training_candidate_clip_count": sum(
            1 for record in metadata_records if record["qc"]["training_candidate"]
        ),
        "diagnostic_clip_count": sum(
            1 for record in metadata_records if record["qc"]["diagnostic_only"]
        ),
        "train_group_count": len(split_groups.get("train", set())),
        "val_group_count": len(split_groups.get("val", set())),
        "test_group_count": len(split_groups.get("test", set())),
    }


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


def _environment_info() -> dict[str, Any]:
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
    }


def load_config(path: Path | None, overrides: list[str] | None = None) -> Phase3Config:
    raw: dict[str, Any] = {}
    if path is not None:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    data = _apply_overrides(raw, overrides or [])
    clip = dict(data.get("clip", {}) or {})
    output_sampling = dict(data.get("output_sampling", {}) or {})
    render = dict(data.get("render", {}) or {})
    filters = dict(data.get("filters", {}) or {})
    freeze = dict(data.get("freeze", {}) or {})
    metadata = dict(data.get("metadata", {}) or {})
    output_name = str(data.get("output_name", "phase3_gt_passthrough_v1"))
    return Phase3Config(
        dataset_id=str(data.get("dataset_id", "phase3_gt_passthrough_v1")),
        input_dataset=Path(str(data.get("input_dataset", ""))),
        output_dir=Path(
            str(data.get("output_dir") or default_named_output_dir(output_name))
        ),
        config_path=path,
        cli_overrides=tuple(overrides or ()),
        duration_s=float(clip.get("duration_s", 0.5)),
        window_policy=str(clip.get("window_policy", "non_overlap")),
        overlap_stride_fraction=float(clip.get("overlap_stride_fraction", 0.5)),
        fps_out=float(output_sampling.get("fps_out", 25.0)),
        image_size_px=int(render.get("image_size_px", 96)),
        pixel_size_um=float(render.get("pixel_size_um", 0.1)),
        body_length_um=float(render.get("body_length_um", 2.0)),
        body_width_um=float(render.get("body_width_um", 1.0)),
        body_intensity=int(render.get("body_intensity", 60)),
        allowed_n_flagella=tuple(
            int(v) for v in filters.get("allowed_n_flagella", [1, 2, 3])
        ),
        max_per_class=(
            None
            if filters.get("max_per_class") in (None, "")
            else int(filters.get("max_per_class"))
        ),
        baseline_torque_Nm=float(freeze.get("baseline_torque_Nm", 2.0e-20)),
        require_use_for_ml_candidate=bool(
            filters.get("require_use_for_ml_candidate", True)
        ),
        source_require_use_for_ml_candidate=bool(
            filters.get(
                "source_require_use_for_ml_candidate",
                filters.get("require_use_for_ml_candidate", True),
            )
        ),
        dataset_version=str(metadata.get("dataset_version", "v1")),
        dataset_revision=(
            None
            if metadata.get("dataset_revision") in (None, "")
            else str(metadata.get("dataset_revision"))
        ),
        output_name=output_name,
    )


def _apply_overrides(raw: dict[str, Any], overrides: list[str]) -> dict[str, Any]:
    data = dict(raw)
    for item in overrides:
        if "=" not in item:
            raise ValueError(f"Invalid override; expected KEY=VALUE: {item}")
        key, value = item.split("=", 1)
        node = data
        parts = key.split(".")
        for part in parts[:-1]:
            child = node.setdefault(part, {})
            if not isinstance(child, dict):
                raise ValueError(f"Override path conflict: {key}")
            node = child
        node[parts[-1]] = yaml.safe_load(value)
    return data


def validate_training_candidate(
    row: dict[str, str], cfg: Phase3Config
) -> tuple[bool, str | None]:
    n_flagella = int(_to_float(row.get("n_flagella", "0")))
    if n_flagella not in cfg.allowed_n_flagella:
        return False, "n_flagella_not_in_mvp_scope"
    if cfg.source_require_use_for_ml_candidate and not _to_bool(
        row.get("use_for_ml_candidate")
    ):
        return False, "not_use_for_ml_candidate"
    torque = _to_float(row.get("torque_Nm", "nan"))
    if not np.isclose(torque, cfg.baseline_torque_Nm, rtol=1.0e-9, atol=0.0):
        return False, "torque_variation_diagnostic_only"
    return True, None


def select_samples(
    rows: list[dict[str, str]], cfg: Phase3Config
) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    per_class: dict[int, int] = {}
    for row in rows:
        ok, _reason = validate_training_candidate(row, cfg)
        if not ok:
            continue
        n_flagella = int(_to_float(row["n_flagella"]))
        current_count = per_class.get(n_flagella, 0)
        if cfg.max_per_class is not None and current_count >= cfg.max_per_class:
            continue
        selected.append(row)
        per_class[n_flagella] = current_count + 1
    return selected


def build_clip_dataset(cfg: Phase3Config) -> Path:
    summary_path = cfg.input_dataset / "summary.csv"
    if not summary_path.is_file():
        raise FileNotFoundError(f"summary.csv not found: {summary_path}")
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    clips_dir = cfg.output_dir / "clips"
    clips_dir.mkdir(parents=True, exist_ok=True)

    all_rows = _read_csv(summary_path)
    source_manifest = _load_source_manifest(cfg.input_dataset)
    source_metadata = _source_metadata(source_manifest)
    dataset_version = str(
        cfg.dataset_version or source_metadata.get("dataset_version") or "v1"
    )
    dataset_revision = cfg.dataset_revision
    if dataset_revision is None and source_metadata.get("dataset_revision") is not None:
        dataset_revision = str(source_metadata["dataset_revision"])
    render_config = BodyCapsuleRenderConfig(
        image_size_px=cfg.image_size_px,
        pixel_size_um=cfg.pixel_size_um,
        body_length_um=cfg.body_length_um,
        body_width_um=cfg.body_width_um,
        body_intensity=cfg.body_intensity,
        tracking_center=True,
    )
    render_id = body_capsule_render_id(render_config)
    selected_rows = select_samples(all_rows, cfg)
    metadata_records: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    split_rows: list[dict[str, Any]] = []

    group_labels = {
        f"phase2:{dataset_version}:{row['sample_id']}": int(
            _to_float(row["n_flagella"])
        )
        for row in selected_rows
    }
    split_by_group = assign_grouped_splits(
        group_labels.keys(),
        group_labels=group_labels,
    )

    for row in selected_rows:
        sample_id = str(row["sample_id"])
        n_flagella = int(_to_float(row["n_flagella"]))
        raw_run_dir = Path(str(row["raw_dir"]))
        archive_path = raw_run_dir / "state_archive.npz"
        if not archive_path.is_file():
            qc_rows.append(
                {
                    "sample_id": sample_id,
                    "n_flagella": n_flagella,
                    "status": "fail",
                    "exclusion_reason": "missing_state_archive",
                    "clip_count": 0,
                }
            )
            continue

        states = select_frames(load_state_archive(archive_path), cfg.fps_out)
        windows = generate_windows(
            source_frame_count=len(states),
            frame_rate_hz=cfg.fps_out,
            duration_s=cfg.duration_s,
            policy=cfg.window_policy,
            overlap_stride_fraction=cfg.overlap_stride_fraction,
        )
        group_key = f"phase2:{dataset_version}:{sample_id}"
        source_duration_s = len(states) / cfg.fps_out
        run_shape_pass = _to_bool(row.get("shape_pass", "true"))
        run_first_fail_t_s = _optional_float(row.get("first_fail_t_s"))
        run_training_count = 0
        run_diagnostic_count = 0
        for clip_index, window in enumerate(windows):
            clip_id = f"{sample_id}_c{clip_index:04d}"
            output_path = clips_dir / f"{clip_id}.npy"
            clip_array, geometries = render_clip_array(
                states[window.start : window.end],
                image_size_px=cfg.image_size_px,
                pixel_size_um=cfg.pixel_size_um,
                body_length_um=cfg.body_length_um,
                body_width_um=cfg.body_width_um,
                body_intensity=cfg.body_intensity,
            )
            (
                window_shape_pass,
                is_pre_first_fail,
                training_candidate,
                diagnostic_only,
                exclusion_reason,
                qc_label,
            ) = _window_qc(
                window_end_s=window.end / cfg.fps_out,
                run_shape_pass=run_shape_pass,
                run_first_fail_t_s=run_first_fail_t_s,
            )
            run_training_count += int(training_candidate)
            run_diagnostic_count += int(diagnostic_only)
            np.save(output_path, clip_array)
            metadata = build_gt_passthrough_metadata(
                dataset_id=cfg.dataset_id,
                source_video_id=sample_id,
                source_path=archive_path,
                frame_rate_hz=cfg.fps_out,
                source_frame_count=len(states),
                source_duration_s=source_duration_s,
                run_id=sample_id,
                raw_run_dir=raw_run_dir,
                n_flagella=n_flagella,
                track_id=f"{sample_id}:gt_track_0000",
                group_key=group_key,
                clip_id=clip_id,
                clip_index=clip_index,
                window=window,
                window_policy=cfg.window_policy,
                output_path=output_path,
                crop_size_px=cfg.image_size_px,
                pixel_size_um=cfg.pixel_size_um,
                frame_geometries=geometries,
                render_id=render_id,
                dataset_version=dataset_version,
                dataset_revision=dataset_revision,
                run_shape_pass=run_shape_pass,
                run_first_fail_t_s=run_first_fail_t_s,
                window_shape_pass=window_shape_pass,
                is_pre_first_fail=is_pre_first_fail,
                training_candidate=training_candidate,
                diagnostic_only=diagnostic_only,
                exclusion_reason=exclusion_reason,
                qc_label=qc_label,
            )
            metadata_records.append(metadata)
            split_rows.append(
                {
                    "clip_id": clip_id,
                    "sample_id": sample_id,
                    "group_key": group_key,
                    "split": split_by_group[group_key],
                    "n_flagella": n_flagella,
                    "clip_index": clip_index,
                    "t_start_s": window.start / cfg.fps_out,
                    "t_end_s": window.end / cfg.fps_out,
                    "qc_label": qc_label,
                    "training_candidate": training_candidate,
                    "exclusion_reason": exclusion_reason,
                }
            )
        qc_rows.append(
            {
                "sample_id": sample_id,
                "n_flagella": n_flagella,
                "status": "pass" if windows else "fail",
                "exclusion_reason": None if windows else "no_complete_window",
                "clip_count": len(windows),
                "run_shape_pass": run_shape_pass,
                "run_first_fail_t_s": run_first_fail_t_s,
                "training_candidate_clip_count": run_training_count,
                "diagnostic_clip_count": run_diagnostic_count,
            }
        )

    assert_no_group_leakage(
        {"group_key": str(row["group_key"]), "split": str(row["split"])}
        for row in split_rows
    )
    metadata_path = cfg.output_dir / "clip_metadata.jsonl"
    with metadata_path.open("w", encoding="utf-8") as handle:
        for record in metadata_records:
            handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")
    _write_csv(
        cfg.output_dir / "split_summary.csv",
        split_rows,
        [
            "clip_id",
            "sample_id",
            "group_key",
            "split",
            "n_flagella",
            "clip_index",
            "t_start_s",
            "t_end_s",
            "qc_label",
            "training_candidate",
            "exclusion_reason",
        ],
    )
    _write_csv(
        cfg.output_dir / "qc_summary.csv",
        qc_rows,
        [
            "sample_id",
            "n_flagella",
            "status",
            "exclusion_reason",
            "clip_count",
            "run_shape_pass",
            "run_first_fail_t_s",
            "training_candidate_clip_count",
            "diagnostic_clip_count",
        ],
    )
    dataset_summary = _dataset_summary(metadata_records, split_rows)
    _write_csv(
        cfg.output_dir / "dataset_summary.csv",
        [dataset_summary],
        list(dataset_summary),
    )
    manifest = {
        "pipeline_name": "phase3_gt_passthrough",
        "schema_version": SCHEMA_VERSION,
        "created_at": _now_jst().isoformat(),
        "dataset_id": cfg.dataset_id,
        "dataset_version": dataset_version,
        "dataset_revision": dataset_revision,
        "input_dataset": str(cfg.input_dataset),
        "source_dataset_manifest": str(cfg.input_dataset / "dataset_manifest.json"),
        "output_dir": str(cfg.output_dir),
        "clip": {
            "duration_s": cfg.duration_s,
            "window_policy": cfg.window_policy,
            "overlap_stride_fraction": cfg.overlap_stride_fraction,
        },
        "output_sampling": {
            "fps_out": cfg.fps_out,
        },
        "render": {
            "render_mode": RENDER_MODE,
            "render_id": render_id,
            "projection": "orthographic",
            "body_deformation_rendered": False,
            "rendered_objects": ["body"],
            "excluded_objects": ["flagella"],
            "body_shape": "capsule",
            "body_length_definition": "end_to_end",
            "image_size_px": cfg.image_size_px,
            "pixel_size_um": cfg.pixel_size_um,
            "body_length_um": cfg.body_length_um,
            "body_width_um": cfg.body_width_um,
            "body_intensity": cfg.body_intensity,
            "background_intensity": 255,
            "tracking_center": True,
        },
        "filters": {
            "allowed_n_flagella": list(cfg.allowed_n_flagella),
            "require_use_for_ml_candidate": cfg.require_use_for_ml_candidate,
            "source_require_use_for_ml_candidate": (
                cfg.source_require_use_for_ml_candidate
            ),
            "baseline_torque_Nm": cfg.baseline_torque_Nm,
            "max_per_class": cfg.max_per_class,
        },
        "window_qc_policy": {
            "first_fail_policy": "strict_pre_first_fail",
            "training_candidate_rule": (
                "window_end_s <= run_first_fail_t_s or no run first-fail"
            ),
            "diagnostic_rule": (
                "first-fail-containing and later windows are diagnostic-only"
            ),
            "warmup_s": [0.0, 0.5, 1.0],
            "early_clips_kept": True,
        },
        "invocation": {
            "config_path": str(cfg.config_path)
            if cfg.config_path is not None
            else None,
            "cli_overrides": list(cfg.cli_overrides),
        },
        "outputs": {
            "clip_metadata_jsonl": str(metadata_path),
            "clips_dir": str(clips_dir),
            "split_summary_csv": str(cfg.output_dir / "split_summary.csv"),
            "qc_summary_csv": str(cfg.output_dir / "qc_summary.csv"),
            "dataset_summary_csv": str(cfg.output_dir / "dataset_summary.csv"),
        },
        "sample_count": len(selected_rows),
        "clip_count": len(metadata_records),
        "dataset_summary": dataset_summary,
        "git": _git_info(),
        "environment": _environment_info(),
    }
    (cfg.output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    (cfg.output_dir / "run.log").write_text(
        "\n".join(
            [
                f"created_at={manifest['created_at']}",
                f"input_dataset={cfg.input_dataset}",
                f"output_dir={cfg.output_dir}",
                f"sample_count={len(selected_rows)}",
                f"clip_count={len(metadata_records)}",
                "training_candidate_clip_count="
                f"{dataset_summary['training_candidate_clip_count']}",
                f"diagnostic_clip_count={dataset_summary['diagnostic_clip_count']}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return cfg.output_dir
