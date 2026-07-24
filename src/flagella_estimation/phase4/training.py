"""Training workflow for the Phase 4 deterministic baseline."""

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

from flagella_estimation.phase4.baseline import (
    FEATURE_NAMES,
    NearestCentroidModel,
    classification_metrics,
    confusion_matrix,
    extract_clip_features,
    fit_nearest_centroid,
    predict_nearest_centroid,
)
from flagella_estimation.phase4.dataset import (
    Phase4ClipSample,
    Phase4DatasetAudit,
    audit_phase4_clip_dataset,
    load_phase3_common_clip_dataset,
)


@dataclass(frozen=True)
class Phase4BaselineConfig:
    dataset_dir: Path
    output_dir: Path
    config_path: Path | None = None
    cli_overrides: tuple[str, ...] = ()
    allowed_n_flagella: tuple[int, ...] = (1, 2, 3)
    required_dataset_version: str = "v1"
    required_clip_duration_s: float = 0.5
    required_window_policy: str = "non_overlap"
    seed: int = 0


def default_output_dir() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / "phase4_baseline_v1"
    )


def load_baseline_config(
    path: Path | None, overrides: list[str] | None = None
) -> Phase4BaselineConfig:
    raw: dict[str, Any] = {}
    if path is not None:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    data = _apply_overrides(raw, overrides or [])
    freeze = dict(data.get("freeze", {}) or {})
    return Phase4BaselineConfig(
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
        seed=int(data.get("seed", 0)),
    )


def train_baseline_classifier(cfg: Phase4BaselineConfig) -> Path:
    """Train and evaluate the baseline using the dataset's fixed grouped splits."""

    source_manifest = json.loads(
        (cfg.dataset_dir / "manifest.json").read_text(encoding="utf-8")
    )
    samples = load_phase3_common_clip_dataset(cfg.dataset_dir)
    audit = audit_phase4_clip_dataset(samples)
    _validate_frozen_dataset(samples, source_manifest, cfg)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    features = np.stack(
        [extract_clip_features(np.load(sample.clip_path)) for sample in samples]
    )
    labels = np.asarray([sample.n_flagella for sample in samples], dtype=np.int64)
    train_mask = np.asarray([sample.split == "train" for sample in samples])
    train_classes = set(int(value) for value in labels[train_mask])
    required_classes = set(cfg.allowed_n_flagella)
    if train_classes != required_classes:
        raise ValueError(
            f"train split classes={sorted(train_classes)} "
            f"but required={sorted(required_classes)}"
        )

    model = fit_nearest_centroid(features[train_mask], labels[train_mask])
    predictions = predict_nearest_centroid(model, features)
    metrics_by_split: dict[str, dict[str, float | int]] = {}
    matrices: dict[str, np.ndarray] = {}
    for split in ("train", "val", "test"):
        split_mask = np.asarray([sample.split == split for sample in samples])
        if not split_mask.any():
            raise ValueError(f"dataset has no samples in required split: {split}")
        metrics_by_split[split] = classification_metrics(
            labels[split_mask], predictions[split_mask], model.classes
        )
        matrices[split] = confusion_matrix(
            labels[split_mask], predictions[split_mask], model.classes
        )

    _write_outputs(
        cfg=cfg,
        samples=samples,
        features=features,
        labels=labels,
        predictions=predictions,
        model=model,
        metrics_by_split=metrics_by_split,
        matrices=matrices,
        audit=audit,
        source_manifest=source_manifest,
    )
    return cfg.output_dir


def _validate_frozen_dataset(
    samples: list[Phase4ClipSample],
    source_manifest: dict[str, Any],
    cfg: Phase4BaselineConfig,
) -> None:
    if not samples:
        raise ValueError("dataset has no clips")
    actual_classes = {sample.n_flagella for sample in samples}
    required_classes = set(cfg.allowed_n_flagella)
    if actual_classes != required_classes:
        raise ValueError(
            f"dataset classes={sorted(actual_classes)} "
            f"but required={sorted(required_classes)}"
        )

    manifest_clip = source_manifest.get("clip", {})
    if not np.isclose(
        float(manifest_clip.get("duration_s")),
        cfg.required_clip_duration_s,
        rtol=0.0,
        atol=1.0e-9,
    ):
        raise ValueError("manifest clip duration is outside freeze")
    if manifest_clip.get("window_policy") != cfg.required_window_policy:
        raise ValueError("manifest window policy is outside freeze")

    for sample in samples:
        provenance = sample.metadata.get("provenance", {})
        clip = sample.metadata.get("clip", {})
        qc = sample.metadata.get("qc", {})
        if provenance.get("dataset_version") != cfg.required_dataset_version:
            raise ValueError(f"{sample.clip_id} dataset_version is outside freeze")
        if clip.get("window_policy") != cfg.required_window_policy:
            raise ValueError(f"{sample.clip_id} window policy is outside freeze")
        if qc.get("status") != "pass":
            raise ValueError(f"{sample.clip_id} qc.status must be pass")


def _write_outputs(
    *,
    cfg: Phase4BaselineConfig,
    samples: list[Phase4ClipSample],
    features: np.ndarray,
    labels: np.ndarray,
    predictions: np.ndarray,
    model: NearestCentroidModel,
    metrics_by_split: dict[str, dict[str, float | int]],
    matrices: dict[str, np.ndarray],
    audit: Phase4DatasetAudit,
    source_manifest: dict[str, Any],
) -> None:
    np.savez(
        cfg.output_dir / "model.npz",
        classes=model.classes,
        feature_names=np.asarray(FEATURE_NAMES),
        feature_mean=model.feature_mean,
        feature_scale=model.feature_scale,
        centroids=model.centroids,
    )
    prediction_rows = []
    for sample, actual, predicted, feature_row in zip(
        samples, labels, predictions, features, strict=True
    ):
        row: dict[str, Any] = {
            "clip_id": sample.clip_id,
            "group_key": sample.group_key,
            "split": sample.split,
            "n_flagella": int(actual),
            "predicted_n_flagella": int(predicted),
            "correct": bool(actual == predicted),
        }
        row.update(
            {
                f"feature_{name}": float(value)
                for name, value in zip(FEATURE_NAMES, feature_row, strict=True)
            }
        )
        prediction_rows.append(row)
    _write_csv(cfg.output_dir / "predictions.csv", prediction_rows)

    confusion_rows = []
    for split, matrix in matrices.items():
        for row_index, actual in enumerate(model.classes):
            for column_index, predicted in enumerate(model.classes):
                confusion_rows.append(
                    {
                        "split": split,
                        "actual_n_flagella": int(actual),
                        "predicted_n_flagella": int(predicted),
                        "count": int(matrix[row_index, column_index]),
                    }
                )
    _write_csv(cfg.output_dir / "confusion_matrix.csv", confusion_rows)

    metrics = {
        "baseline_name": "standardized_nearest_centroid_v1",
        "classes": [int(value) for value in model.classes],
        "feature_names": list(FEATURE_NAMES),
        "metrics_by_split": metrics_by_split,
    }
    _write_json(cfg.output_dir / "metrics.json", metrics)

    created_at = datetime.now(ZoneInfo("Asia/Tokyo")).isoformat()
    manifest = {
        "pipeline_name": "phase4_baseline_classifier",
        "pipeline_version": "0.1.0",
        "created_at": created_at,
        "input_dataset": str(cfg.dataset_dir),
        "output_dir": str(cfg.output_dir),
        "seed": cfg.seed,
        "freeze": {
            "allowed_n_flagella": list(cfg.allowed_n_flagella),
            "dataset_version": cfg.required_dataset_version,
            "clip_duration_s": cfg.required_clip_duration_s,
            "window_policy": cfg.required_window_policy,
        },
        "dataset_audit": {
            "sample_count": audit.sample_count,
            "group_count": audit.group_count,
            "split_counts": audit.split_counts,
            "class_counts": {
                str(key): value for key, value in audit.class_counts.items()
            },
        },
        "source_manifest": {
            "dataset_id": source_manifest.get("dataset_id"),
            "pipeline_name": source_manifest.get("pipeline_name"),
            "schema_version": source_manifest.get("schema_version"),
            "git": source_manifest.get("git"),
        },
        "invocation": {
            "config_path": str(cfg.config_path) if cfg.config_path else None,
            "cli_overrides": list(cfg.cli_overrides),
        },
        "outputs": {
            "model": str(cfg.output_dir / "model.npz"),
            "metrics": str(cfg.output_dir / "metrics.json"),
            "predictions": str(cfg.output_dir / "predictions.csv"),
            "confusion_matrix": str(cfg.output_dir / "confusion_matrix.csv"),
        },
        "git": _git_info(),
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
    }
    _write_json(cfg.output_dir / "manifest.json", manifest)
    (cfg.output_dir / "run.log").write_text(
        "\n".join(
            [
                f"created_at={created_at}",
                f"input_dataset={cfg.dataset_dir}",
                f"output_dir={cfg.output_dir}",
                f"sample_count={audit.sample_count}",
                f"group_count={audit.group_count}",
                f"test_macro_f1={metrics_by_split['test']['macro_f1']}",
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


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


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
