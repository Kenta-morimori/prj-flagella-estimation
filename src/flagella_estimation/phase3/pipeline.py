"""Phase 3 pseudo-GT passthrough pipeline."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
import subprocess
from typing import Any
from zoneinfo import ZoneInfo

import numpy as np
import yaml

from flagella_estimation.phase3.metadata import (
    SCHEMA_VERSION,
    build_gt_passthrough_metadata,
)
from flagella_estimation.phase3.render import render_clip_array, select_frames
from flagella_estimation.phase3.splits import (
    assign_grouped_splits,
    assert_no_group_leakage,
)
from flagella_estimation.phase3.windows import generate_windows
from sim_swim.analysis.flagella_count_behavior import load_state_archive


@dataclass(frozen=True)
class Phase3Config:
    dataset_id: str
    input_dataset: Path
    output_dir: Path
    duration_s: float = 0.5
    window_policy: str = "non_overlap"
    overlap_stride_fraction: float = 0.5
    frame_rate_hz: float = 25.0
    crop_size_px: int = 96
    pixel_size_um: float = 0.1
    allowed_n_flagella: tuple[int, ...] = (1, 2, 3)
    max_per_class: int | None = None
    baseline_torque_Nm: float = 2.0e-20
    require_use_for_ml_candidate: bool = True


def _now_jst() -> datetime:
    return datetime.now(ZoneInfo("Asia/Tokyo"))


def default_output_dir() -> Path:
    now = _now_jst()
    return (
        Path("outputs")
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / "phase3_gt_passthrough_v1"
    )


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


def load_config(path: Path | None, overrides: list[str] | None = None) -> Phase3Config:
    raw: dict[str, Any] = {}
    if path is not None:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    data = _apply_overrides(raw, overrides or [])
    clip = dict(data.get("clip", {}) or {})
    filters = dict(data.get("filters", {}) or {})
    freeze = dict(data.get("freeze", {}) or {})
    return Phase3Config(
        dataset_id=str(data.get("dataset_id", "phase3_gt_passthrough_v1")),
        input_dataset=Path(str(data.get("input_dataset", ""))),
        output_dir=Path(str(data.get("output_dir") or default_output_dir())),
        duration_s=float(clip.get("duration_s", 0.5)),
        window_policy=str(clip.get("window_policy", "non_overlap")),
        overlap_stride_fraction=float(clip.get("overlap_stride_fraction", 0.5)),
        frame_rate_hz=float(clip.get("frame_rate_hz", 25.0)),
        crop_size_px=int(clip.get("crop_size_px", 96)),
        pixel_size_um=float(clip.get("pixel_size_um", 0.1)),
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
    if cfg.require_use_for_ml_candidate and not _to_bool(
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
    selected_rows = select_samples(all_rows, cfg)
    metadata_records: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    split_rows: list[dict[str, Any]] = []

    group_labels = {
        f"phase2:v1:{row['sample_id']}": int(_to_float(row["n_flagella"]))
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

        states = select_frames(load_state_archive(archive_path), cfg.frame_rate_hz)
        windows = generate_windows(
            source_frame_count=len(states),
            frame_rate_hz=cfg.frame_rate_hz,
            duration_s=cfg.duration_s,
            policy=cfg.window_policy,
            overlap_stride_fraction=cfg.overlap_stride_fraction,
        )
        group_key = f"phase2:v1:{sample_id}"
        source_duration_s = len(states) / cfg.frame_rate_hz
        for clip_index, window in enumerate(windows):
            clip_id = f"{sample_id}_c{clip_index:04d}"
            output_path = clips_dir / f"{clip_id}.npy"
            clip_array, geometries = render_clip_array(
                states[window.start : window.end],
                image_size_px=cfg.crop_size_px,
                pixel_size_um=cfg.pixel_size_um,
            )
            np.save(output_path, clip_array)
            metadata = build_gt_passthrough_metadata(
                dataset_id=cfg.dataset_id,
                source_video_id=sample_id,
                source_path=archive_path,
                frame_rate_hz=cfg.frame_rate_hz,
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
                crop_size_px=cfg.crop_size_px,
                pixel_size_um=cfg.pixel_size_um,
                frame_geometries=geometries,
            )
            metadata_records.append(metadata)
            split_rows.append(
                {
                    "clip_id": clip_id,
                    "sample_id": sample_id,
                    "group_key": group_key,
                    "split": split_by_group[group_key],
                    "n_flagella": n_flagella,
                }
            )
        qc_rows.append(
            {
                "sample_id": sample_id,
                "n_flagella": n_flagella,
                "status": "pass" if windows else "fail",
                "exclusion_reason": None if windows else "no_complete_window",
                "clip_count": len(windows),
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
        ["clip_id", "sample_id", "group_key", "split", "n_flagella"],
    )
    _write_csv(
        cfg.output_dir / "qc_summary.csv",
        qc_rows,
        ["sample_id", "n_flagella", "status", "exclusion_reason", "clip_count"],
    )
    manifest = {
        "pipeline_name": "phase3_gt_passthrough",
        "schema_version": SCHEMA_VERSION,
        "created_at": _now_jst().isoformat(),
        "dataset_id": cfg.dataset_id,
        "input_dataset": str(cfg.input_dataset),
        "output_dir": str(cfg.output_dir),
        "clip": {
            "duration_s": cfg.duration_s,
            "window_policy": cfg.window_policy,
            "overlap_stride_fraction": cfg.overlap_stride_fraction,
            "frame_rate_hz": cfg.frame_rate_hz,
            "crop_size_px": cfg.crop_size_px,
            "pixel_size_um": cfg.pixel_size_um,
        },
        "filters": {
            "allowed_n_flagella": list(cfg.allowed_n_flagella),
            "require_use_for_ml_candidate": cfg.require_use_for_ml_candidate,
            "baseline_torque_Nm": cfg.baseline_torque_Nm,
        },
        "outputs": {
            "clip_metadata_jsonl": str(metadata_path),
            "clips_dir": str(clips_dir),
            "split_summary_csv": str(cfg.output_dir / "split_summary.csv"),
            "qc_summary_csv": str(cfg.output_dir / "qc_summary.csv"),
        },
        "sample_count": len(selected_rows),
        "clip_count": len(metadata_records),
        "git": _git_info(),
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
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return cfg.output_dir
