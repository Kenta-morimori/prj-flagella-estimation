"""Phase 4 loader for Phase 3 common clip datasets."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np


VALID_SPLITS = {"train", "val", "test"}


@dataclass(frozen=True)
class Phase4ClipSample:
    clip_id: str
    clip_path: Path
    split: str
    n_flagella: int
    group_key: str
    frame_count: int
    frame_shape: tuple[int, int]
    metadata: dict[str, Any]


@dataclass(frozen=True)
class Phase4DatasetAudit:
    dataset_dir: Path
    sample_count: int
    split_counts: dict[str, int]
    class_counts: dict[int, int]
    split_class_counts: dict[tuple[str, int], int]
    group_count: int


def load_phase3_common_clip_dataset(dataset_dir: Path) -> list[Phase4ClipSample]:
    """Load Phase 3 common clip records without starting model training."""

    dataset_dir = Path(dataset_dir)
    manifest = _load_json(dataset_dir / "manifest.json")
    records = _load_metadata_jsonl(dataset_dir / "clip_metadata.jsonl")
    split_rows = _load_split_summary(dataset_dir / "split_summary.csv")
    split_by_clip_id = {row["clip_id"]: row for row in split_rows}

    expected_clip_count = int(manifest.get("clip_count", len(records)))
    if expected_clip_count != len(records):
        raise ValueError(
            f"manifest clip_count={expected_clip_count} but metadata has {len(records)}"
        )
    if len(split_by_clip_id) != len(split_rows):
        raise ValueError("split_summary.csv contains duplicate clip_id values")

    metadata_clip_ids = {
        str(_required_mapping(record, "clip")["clip_id"]) for record in records
    }
    split_clip_ids = set(split_by_clip_id)
    if metadata_clip_ids != split_clip_ids:
        missing_from_split = sorted(metadata_clip_ids - split_clip_ids)
        missing_from_metadata = sorted(split_clip_ids - metadata_clip_ids)
        raise ValueError(
            "metadata and split_summary clip_id sets differ: "
            f"missing_from_split={missing_from_split}, "
            f"missing_from_metadata={missing_from_metadata}"
        )

    samples: list[Phase4ClipSample] = []
    seen_clip_ids: set[str] = set()
    for record in records:
        clip = _required_mapping(record, "clip")
        labels = _required_mapping(record, "labels")
        track = _required_mapping(record, "track")
        normalization = _required_mapping(record, "normalization")

        clip_id = str(clip["clip_id"])
        if clip_id in seen_clip_ids:
            raise ValueError(f"duplicate clip_id in metadata: {clip_id}")
        seen_clip_ids.add(clip_id)
        if clip_id not in split_by_clip_id:
            raise ValueError(f"clip_id missing from split_summary.csv: {clip_id}")

        split_row = split_by_clip_id[clip_id]
        split = str(split_row["split"])
        if split not in VALID_SPLITS:
            raise ValueError(f"unsupported split for {clip_id}: {split}")

        group_key = str(track["group_key"])
        if str(split_row["group_key"]) != group_key:
            raise ValueError(f"group_key mismatch for {clip_id}")

        n_flagella = int(labels["n_flagella"])
        if int(split_row["n_flagella"]) != n_flagella:
            raise ValueError(f"n_flagella mismatch for {clip_id}")

        clip_path = _resolve_clip_path(dataset_dir, str(clip["output_path"]))
        arr = np.load(clip_path)
        if arr.dtype != np.uint8:
            raise ValueError(f"{clip_id} dtype must be uint8, got {arr.dtype}")
        if arr.ndim != 3:
            raise ValueError(f"{clip_id} must have shape (T, H, W), got {arr.shape}")

        frame_count = int(clip["frame_count"])
        if arr.shape[0] != frame_count:
            raise ValueError(
                f"{clip_id} frame_count mismatch: metadata={frame_count}, array={arr.shape[0]}"
            )
        crop_size = normalization.get("crop_size_px")
        if crop_size is not None and list(arr.shape[1:]) != list(crop_size):
            raise ValueError(
                f"{clip_id} crop_size_px mismatch: {crop_size} vs {arr.shape[1:]}"
            )
        if len(record.get("frames", [])) != frame_count:
            raise ValueError(f"{clip_id} frames metadata length mismatch")

        samples.append(
            Phase4ClipSample(
                clip_id=clip_id,
                clip_path=clip_path,
                split=split,
                n_flagella=n_flagella,
                group_key=group_key,
                frame_count=frame_count,
                frame_shape=(int(arr.shape[1]), int(arr.shape[2])),
                metadata=record,
            )
        )

    return samples


def audit_phase4_clip_dataset(samples: list[Phase4ClipSample]) -> Phase4DatasetAudit:
    """Return counts and raise if group_key leakage is detected."""

    split_by_group: dict[str, str] = {}
    split_counts: dict[str, int] = {}
    class_counts: dict[int, int] = {}
    split_class_counts: dict[tuple[str, int], int] = {}
    dataset_dir = samples[0].clip_path.parents[1] if samples else Path(".")

    for sample in samples:
        previous_split = split_by_group.setdefault(sample.group_key, sample.split)
        if previous_split != sample.split:
            raise ValueError(
                f"group_key leakage: {sample.group_key} appears in "
                f"{previous_split} and {sample.split}"
            )
        split_counts[sample.split] = split_counts.get(sample.split, 0) + 1
        class_counts[sample.n_flagella] = class_counts.get(sample.n_flagella, 0) + 1
        key = (sample.split, sample.n_flagella)
        split_class_counts[key] = split_class_counts.get(key, 0) + 1

    return Phase4DatasetAudit(
        dataset_dir=dataset_dir,
        sample_count=len(samples),
        split_counts=split_counts,
        class_counts=class_counts,
        split_class_counts=split_class_counts,
        group_count=len(split_by_group),
    )


def _load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    return json.loads(path.read_text(encoding="utf-8"))


def _load_metadata_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    records: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            records.append(json.loads(line))
    return records


def _load_split_summary(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _required_mapping(record: dict[str, Any], key: str) -> dict[str, Any]:
    value = record.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"metadata field must be an object: {key}")
    return value


def _resolve_clip_path(dataset_dir: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    candidates = [path] if path.is_absolute() else [path, dataset_dir / path]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(raw_path)
