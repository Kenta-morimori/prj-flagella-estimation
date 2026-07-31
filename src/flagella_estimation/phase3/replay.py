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
import numpy as np


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
