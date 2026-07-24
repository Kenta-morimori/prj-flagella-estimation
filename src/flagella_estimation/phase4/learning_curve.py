"""Grouped learning curves for the Phase 4 diagnostic baseline."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
import itertools
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo

import numpy as np
import yaml

from flagella_estimation.phase4.baseline import (
    classification_metrics,
    confusion_matrix,
    extract_clip_features,
    fit_nearest_centroid,
    predict_nearest_centroid,
)
from flagella_estimation.phase4.dataset import (
    Phase4ClipSample,
    audit_phase4_clip_dataset,
    load_phase3_common_clip_dataset,
)
from flagella_estimation.phase4.training import (
    Phase4BaselineConfig,
    validate_frozen_dataset,
)


@dataclass(frozen=True)
class Phase4LearningCurveConfig:
    dataset_dir: Path
    output_dir: Path
    config_path: Path | None = None
    cli_overrides: tuple[str, ...] = ()
    allowed_n_flagella: tuple[int, ...] = (1, 2, 3)
    required_dataset_version: str = "v1"
    required_clip_duration_s: float = 0.5
    required_window_policy: str = "non_overlap"
    development_splits: tuple[str, ...] = ("train", "val")
    protected_split: str = "test"
    holdout_groups_per_class: int = 1
    train_groups_per_class: tuple[int, ...] = ()
    repeats: int = 20
    seed: int = 0


@dataclass(frozen=True)
class GroupFeature:
    group_key: str
    split: str
    n_flagella: int
    features: np.ndarray
    clip_count: int


def default_output_dir() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / "phase4_grouped_learning_curve"
    )


def load_learning_curve_config(
    path: Path | None, overrides: list[str] | None = None
) -> Phase4LearningCurveConfig:
    raw: dict[str, Any] = {}
    if path is not None:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    data = _apply_overrides(raw, overrides or [])
    curve = dict(data.get("learning_curve", {}) or {})
    freeze = dict(data.get("freeze", {}) or {})
    requested_sizes = curve.get("train_groups_per_class", []) or []
    return Phase4LearningCurveConfig(
        dataset_dir=Path(str(data.get("dataset_dir", ""))),
        output_dir=Path(str(data.get("output_dir") or default_output_dir())),
        config_path=path,
        cli_overrides=tuple(overrides or ()),
        allowed_n_flagella=tuple(
            int(value) for value in freeze.get("allowed_n_flagella", [1, 2, 3])
        ),
        required_dataset_version=str(freeze.get("dataset_version", "v1")),
        required_clip_duration_s=float(freeze.get("clip_duration_s", 0.5)),
        required_window_policy=str(freeze.get("window_policy", "non_overlap")),
        development_splits=tuple(
            str(value) for value in curve.get("development_splits", ["train", "val"])
        ),
        protected_split=str(curve.get("protected_split", "test")),
        holdout_groups_per_class=int(curve.get("holdout_groups_per_class", 1)),
        train_groups_per_class=tuple(int(value) for value in requested_sizes),
        repeats=int(curve.get("repeats", 20)),
        seed=int(data.get("seed", 0)),
    )


def aggregate_group_features(
    samples: list[Phase4ClipSample],
) -> list[GroupFeature]:
    """Average clip features within each independent group."""

    by_group: dict[str, list[tuple[Phase4ClipSample, np.ndarray]]] = {}
    for sample in samples:
        feature_row = extract_clip_features(np.load(sample.clip_path))
        by_group.setdefault(sample.group_key, []).append((sample, feature_row))

    groups: list[GroupFeature] = []
    for group_key in sorted(by_group):
        rows = by_group[group_key]
        splits = {sample.split for sample, _ in rows}
        labels = {sample.n_flagella for sample, _ in rows}
        if len(splits) != 1:
            raise ValueError(f"group_key leakage while aggregating: {group_key}")
        if len(labels) != 1:
            raise ValueError(f"group_key has multiple labels: {group_key}")
        groups.append(
            GroupFeature(
                group_key=group_key,
                split=next(iter(splits)),
                n_flagella=next(iter(labels)),
                features=np.mean([features for _, features in rows], axis=0),
                clip_count=len(rows),
            )
        )
    return groups


def evaluate_grouped_learning_curve(cfg: Phase4LearningCurveConfig) -> Path:
    """Evaluate disjoint grouped train/holdout partitions without using test."""

    if cfg.repeats < 1:
        raise ValueError("repeats must be at least 1")
    source_manifest = json.loads(
        (cfg.dataset_dir / "manifest.json").read_text(encoding="utf-8")
    )
    samples = load_phase3_common_clip_dataset(cfg.dataset_dir)
    audit_phase4_clip_dataset(samples)
    validate_frozen_dataset(
        samples,
        source_manifest,
        Phase4BaselineConfig(
            dataset_dir=cfg.dataset_dir,
            output_dir=cfg.output_dir,
            allowed_n_flagella=cfg.allowed_n_flagella,
            required_dataset_version=cfg.required_dataset_version,
            required_clip_duration_s=cfg.required_clip_duration_s,
            required_window_policy=cfg.required_window_policy,
            seed=cfg.seed,
        ),
    )
    groups = aggregate_group_features(samples)
    classes = np.asarray(sorted(cfg.allowed_n_flagella), dtype=np.int64)
    development = [group for group in groups if group.split in cfg.development_splits]
    protected = [group for group in groups if group.split == cfg.protected_split]
    validate_group_sets(development, protected, classes)
    if cfg.holdout_groups_per_class < 1:
        raise ValueError("holdout_groups_per_class must be at least 1")

    development_by_class = {
        int(class_id): sorted(
            (group for group in development if group.n_flagella == int(class_id)),
            key=lambda group: group.group_key,
        )
        for class_id in classes
    }
    max_balanced_size = (
        min(len(rows) for rows in development_by_class.values())
        - cfg.holdout_groups_per_class
    )
    if max_balanced_size < 1:
        raise ValueError(
            "development pool needs at least one train group plus "
            "holdout_groups_per_class for every class"
        )
    requested_sizes = (
        cfg.train_groups_per_class
        if cfg.train_groups_per_class
        else tuple(range(1, max_balanced_size + 1))
    )
    if any(size < 1 or size > max_balanced_size for size in requested_sizes):
        raise ValueError(
            f"train_groups_per_class must be within 1..{max_balanced_size}"
        )

    curve_rows: list[dict[str, Any]] = []
    confusion_rows: list[dict[str, Any]] = []
    for train_size in requested_sizes:
        partitions = _select_balanced_partitions(
            development_by_class,
            train_size=train_size,
            holdout_size=cfg.holdout_groups_per_class,
            repeats=cfg.repeats,
            seed=cfg.seed,
        )
        for repeat_index, (selected, holdout) in enumerate(partitions):
            train_features = np.stack([group.features for group in selected])
            train_labels = np.asarray(
                [group.n_flagella for group in selected], dtype=np.int64
            )
            evaluation_features = np.stack([group.features for group in holdout])
            evaluation_labels = np.asarray(
                [group.n_flagella for group in holdout], dtype=np.int64
            )
            model = fit_nearest_centroid(train_features, train_labels)
            predicted = predict_nearest_centroid(model, evaluation_features)
            metrics = classification_metrics(
                evaluation_labels, predicted, model.classes
            )
            matrix = confusion_matrix(evaluation_labels, predicted, model.classes)
            curve_rows.append(
                {
                    "train_groups_per_class": train_size,
                    "repeat_index": repeat_index,
                    "train_group_count": len(selected),
                    "train_clip_count": sum(group.clip_count for group in selected),
                    "evaluation_group_count": len(holdout),
                    "accuracy": metrics["accuracy"],
                    "balanced_accuracy": metrics["balanced_accuracy"],
                    "macro_f1": metrics["macro_f1"],
                    "selected_group_keys": json.dumps(
                        [group.group_key for group in selected],
                        ensure_ascii=False,
                        separators=(",", ":"),
                    ),
                    "holdout_group_keys": json.dumps(
                        [group.group_key for group in holdout],
                        ensure_ascii=False,
                        separators=(",", ":"),
                    ),
                }
            )
            for row_index, actual in enumerate(model.classes):
                actual_count = int(matrix[row_index, :].sum())
                recall = (
                    float(matrix[row_index, row_index] / actual_count)
                    if actual_count
                    else 0.0
                )
                for column_index, predicted_class in enumerate(model.classes):
                    confusion_rows.append(
                        {
                            "train_groups_per_class": train_size,
                            "repeat_index": repeat_index,
                            "actual_n_flagella": int(actual),
                            "predicted_n_flagella": int(predicted_class),
                            "count": int(matrix[row_index, column_index]),
                            "actual_class_recall": recall,
                        }
                    )

    summary_rows = _summarize_curve(curve_rows, confusion_rows, classes)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(cfg.output_dir / "learning_curve.csv", curve_rows)
    _write_csv(cfg.output_dir / "learning_curve_summary.csv", summary_rows)
    _write_csv(cfg.output_dir / "confusion.csv", confusion_rows)
    _write_manifest(
        cfg=cfg,
        groups=groups,
        classes=classes,
        requested_sizes=requested_sizes,
        curve_rows=curve_rows,
    )
    return cfg.output_dir


def validate_group_sets(
    development: list[GroupFeature],
    evaluation: list[GroupFeature],
    classes: np.ndarray,
) -> None:
    development_keys = {group.group_key for group in development}
    evaluation_keys = {group.group_key for group in evaluation}
    overlap = sorted(development_keys & evaluation_keys)
    if overlap:
        raise ValueError(f"development/evaluation group leakage: {overlap}")
    required = set(int(value) for value in classes)
    development_classes = {group.n_flagella for group in development}
    evaluation_classes = {group.n_flagella for group in evaluation}
    if development_classes != required:
        raise ValueError(
            f"development classes={sorted(development_classes)} "
            f"but required={sorted(required)}"
        )
    if evaluation_classes != required:
        raise ValueError(
            f"evaluation classes={sorted(evaluation_classes)} "
            f"but required={sorted(required)}"
        )


def _select_balanced_partitions(
    development_by_class: dict[int, list[GroupFeature]],
    *,
    train_size: int,
    holdout_size: int,
    repeats: int,
    seed: int,
) -> list[tuple[list[GroupFeature], list[GroupFeature]]]:
    combination_count = math.prod(
        math.comb(len(groups), train_size)
        * math.comb(len(groups) - train_size, holdout_size)
        for groups in development_by_class.values()
    )
    target_count = min(repeats, combination_count)
    if combination_count <= 10_000:
        per_class = []
        for class_id in sorted(development_by_class):
            groups = development_by_class[class_id]
            class_partitions = []
            for train_groups in itertools.combinations(groups, train_size):
                train_keys = {group.group_key for group in train_groups}
                remaining = [
                    group for group in groups if group.group_key not in train_keys
                ]
                for holdout_groups in itertools.combinations(remaining, holdout_size):
                    class_partitions.append((train_groups, holdout_groups))
            per_class.append(class_partitions)
        all_partitions = []
        for class_parts in itertools.product(*per_class):
            selected = list(
                itertools.chain.from_iterable(parts[0] for parts in class_parts)
            )
            holdout = list(
                itertools.chain.from_iterable(parts[1] for parts in class_parts)
            )
            all_partitions.append((selected, holdout))
        rng = np.random.default_rng(seed + train_size)
        order = rng.permutation(len(all_partitions))[:target_count]
        return [all_partitions[int(index)] for index in order]

    rng = np.random.default_rng(seed + train_size)
    partitions: list[tuple[list[GroupFeature], list[GroupFeature]]] = []
    seen: set[tuple[tuple[str, ...], tuple[str, ...]]] = set()
    max_attempts = target_count * 100
    for _ in range(max_attempts):
        selected = []
        holdout = []
        for class_id in sorted(development_by_class):
            groups = development_by_class[class_id]
            indices = rng.choice(
                len(groups),
                size=train_size + holdout_size,
                replace=False,
            )
            selected.extend(
                groups[int(index)] for index in sorted(indices[:train_size])
            )
            holdout.extend(groups[int(index)] for index in sorted(indices[train_size:]))
        key = (
            tuple(sorted(group.group_key for group in selected)),
            tuple(sorted(group.group_key for group in holdout)),
        )
        if key in seen:
            continue
        seen.add(key)
        partitions.append((selected, holdout))
        if len(partitions) == target_count:
            return partitions
    raise RuntimeError("could not sample the requested number of unique partitions")


def _summarize_curve(
    curve_rows: list[dict[str, Any]],
    confusion_rows: list[dict[str, Any]],
    classes: np.ndarray,
) -> list[dict[str, Any]]:
    summary_rows = []
    train_sizes = sorted({int(row["train_groups_per_class"]) for row in curve_rows})
    for train_size in train_sizes:
        matching = [
            row
            for row in curve_rows
            if int(row["train_groups_per_class"]) == train_size
        ]
        summary: dict[str, Any] = {
            "train_groups_per_class": train_size,
            "repeat_count": len(matching),
        }
        for metric in ("accuracy", "balanced_accuracy", "macro_f1"):
            values = np.asarray([float(row[metric]) for row in matching])
            summary[f"{metric}_mean"] = float(values.mean())
            summary[f"{metric}_std"] = float(values.std())
            summary[f"{metric}_p02_5"] = float(np.quantile(values, 0.025))
            summary[f"{metric}_p97_5"] = float(np.quantile(values, 0.975))
        for class_id in classes:
            recalls = [
                float(row["actual_class_recall"])
                for row in confusion_rows
                if int(row["train_groups_per_class"]) == train_size
                and int(row["actual_n_flagella"]) == int(class_id)
                and int(row["predicted_n_flagella"]) == int(class_id)
            ]
            values = np.asarray(recalls)
            prefix = f"n_flagella_{int(class_id)}_recall"
            summary[f"{prefix}_mean"] = float(values.mean())
            summary[f"{prefix}_p02_5"] = float(np.quantile(values, 0.025))
            summary[f"{prefix}_p97_5"] = float(np.quantile(values, 0.975))
        summary_rows.append(summary)
    return summary_rows


def _write_manifest(
    *,
    cfg: Phase4LearningCurveConfig,
    groups: list[GroupFeature],
    classes: np.ndarray,
    requested_sizes: tuple[int, ...],
    curve_rows: list[dict[str, Any]],
) -> None:
    created_at = datetime.now(ZoneInfo("Asia/Tokyo")).isoformat()
    split_counts = {
        split: sum(group.split == split for group in groups)
        for split in sorted({group.split for group in groups})
    }
    manifest = {
        "pipeline_name": "phase4_grouped_learning_curve",
        "pipeline_version": "0.1.0",
        "created_at": created_at,
        "input_dataset": str(cfg.dataset_dir),
        "output_dir": str(cfg.output_dir),
        "classes": [int(value) for value in classes],
        "freeze": {
            "allowed_n_flagella": list(cfg.allowed_n_flagella),
            "dataset_version": cfg.required_dataset_version,
            "clip_duration_s": cfg.required_clip_duration_s,
            "window_policy": cfg.required_window_policy,
        },
        "group_count": len(groups),
        "group_counts_by_split": split_counts,
        "learning_curve": {
            "development_splits": list(cfg.development_splits),
            "protected_split": cfg.protected_split,
            "holdout_groups_per_class": cfg.holdout_groups_per_class,
            "train_groups_per_class": list(requested_sizes),
            "requested_repeats": cfg.repeats,
            "evaluated_subset_count": len(curve_rows),
            "seed": cfg.seed,
            "unit": "unique track.group_key per n_flagella class",
        },
        "invocation": {
            "config_path": str(cfg.config_path) if cfg.config_path else None,
            "cli_overrides": list(cfg.cli_overrides),
        },
        "outputs": {
            "learning_curve": str(cfg.output_dir / "learning_curve.csv"),
            "learning_curve_summary": str(
                cfg.output_dir / "learning_curve_summary.csv"
            ),
            "confusion": str(cfg.output_dir / "confusion.csv"),
        },
        "git": _git_info(),
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
    }
    (cfg.output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    (cfg.output_dir / "run.log").write_text(
        "\n".join(
            [
                f"created_at={created_at}",
                f"input_dataset={cfg.dataset_dir}",
                f"output_dir={cfg.output_dir}",
                f"group_count={len(groups)}",
                f"train_groups_per_class={list(requested_sizes)}",
                f"evaluated_subset_count={len(curve_rows)}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty CSV: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _git_info() -> dict[str, Any]:
    def run(command: list[str]) -> str:
        return subprocess.check_output(
            command, stderr=subprocess.STDOUT, text=True
        ).strip()

    try:
        return {
            "commit": run(["git", "rev-parse", "HEAD"]),
            "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
            "is_clean": run(["git", "status", "--porcelain"]) == "",
        }
    except Exception:
        return {"commit": "unknown", "branch": "unknown", "is_clean": False}


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
