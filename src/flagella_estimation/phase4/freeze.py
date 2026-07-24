"""Machine-readable Phase 4 dataset freeze audit."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from flagella_estimation.phase4.dataset import (
    Phase4ClipSample,
    audit_phase4_clip_dataset,
    load_phase3_common_clip_dataset,
)


@dataclass(frozen=True)
class DatasetFreezePolicy:
    allowed_n_flagella: tuple[int, ...] = (1, 2, 3)
    dataset_version: str = "v1"
    model_ids: tuple[str, ...] = ("phase2_flagella_count_behavior_v1",)
    source_model_ids: tuple[str, ...] = ("flag_spring2p25_body2p5_candidate",)
    render_ids: tuple[str, ...] = ("state_archive_numpy_v1",)
    pipeline_name: str = "phase3_gt_passthrough"
    schema_version: str = "phase3_clip_metadata/v0"
    source_kind: str = "phase2_pseudo"
    processing_mode: str = "gt_passthrough"
    label_source: str = "phase2_gt"
    clip_duration_s: float = 0.5
    window_policy: str = "non_overlap"
    baseline_torque_Nm: float = 2.0e-20
    require_use_for_ml_candidate: bool = True
    group_key_prefix: str = "phase2:v1:"


@dataclass(frozen=True)
class DatasetFreezeAudit:
    status: str
    dataset_dir: Path
    errors: tuple[str, ...]
    warnings: tuple[str, ...]
    observed: dict[str, Any]
    policy: DatasetFreezePolicy

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "dataset_dir": str(self.dataset_dir),
            "errors": list(self.errors),
            "warnings": list(self.warnings),
            "observed": self.observed,
            "policy": asdict(self.policy),
        }


def audit_phase4_dataset_freeze(
    dataset_dir: Path, policy: DatasetFreezePolicy
) -> DatasetFreezeAudit:
    """Load a Phase 3 dataset and return all freeze violations."""

    dataset_dir = Path(dataset_dir)
    try:
        manifest = json.loads(
            (dataset_dir / "manifest.json").read_text(encoding="utf-8")
        )
        samples = load_phase3_common_clip_dataset(dataset_dir)
        audit = audit_phase4_clip_dataset(samples)
    except (
        FileNotFoundError,
        KeyError,
        TypeError,
        ValueError,
        json.JSONDecodeError,
    ) as exc:
        return DatasetFreezeAudit(
            status="FAIL",
            dataset_dir=dataset_dir,
            errors=(f"dataset contract failed: {exc}",),
            warnings=(),
            observed={},
            policy=policy,
        )
    return audit_loaded_dataset_freeze(
        dataset_dir=dataset_dir,
        manifest=manifest,
        samples=samples,
        policy=policy,
        group_count=audit.group_count,
    )


def audit_loaded_dataset_freeze(
    *,
    dataset_dir: Path,
    manifest: dict[str, Any],
    samples: list[Phase4ClipSample],
    policy: DatasetFreezePolicy,
    group_count: int | None = None,
) -> DatasetFreezeAudit:
    """Audit already-loaded samples against one accepted freeze policy."""

    errors: list[str] = []
    warnings: list[str] = []
    manifest_clip = manifest.get("clip", {})
    manifest_filters = manifest.get("filters", {})
    _expect(
        errors,
        "manifest.pipeline_name",
        manifest.get("pipeline_name"),
        policy.pipeline_name,
    )
    _expect(
        errors,
        "manifest.schema_version",
        manifest.get("schema_version"),
        policy.schema_version,
    )
    _expect(
        errors,
        "manifest.clip.window_policy",
        manifest_clip.get("window_policy"),
        policy.window_policy,
    )
    _expect_float(
        errors,
        "manifest.clip.duration_s",
        manifest_clip.get("duration_s"),
        policy.clip_duration_s,
    )
    _expect(
        errors,
        "manifest.filters.allowed_n_flagella",
        sorted(manifest_filters.get("allowed_n_flagella", [])),
        sorted(policy.allowed_n_flagella),
    )
    _expect(
        errors,
        "manifest.filters.require_use_for_ml_candidate",
        manifest_filters.get("require_use_for_ml_candidate"),
        policy.require_use_for_ml_candidate,
    )
    _expect_float(
        errors,
        "manifest.filters.baseline_torque_Nm",
        manifest_filters.get("baseline_torque_Nm"),
        policy.baseline_torque_Nm,
    )

    observed = _observed_values(samples, group_count=group_count)
    _expect(
        errors, "classes", observed["n_flagella"], sorted(policy.allowed_n_flagella)
    )
    _expect(
        errors,
        "dataset_versions",
        observed["dataset_versions"],
        [policy.dataset_version],
    )
    _expect(errors, "model_ids", observed["model_ids"], sorted(policy.model_ids))
    _expect(errors, "render_ids", observed["render_ids"], sorted(policy.render_ids))
    _expect(errors, "source_kinds", observed["source_kinds"], [policy.source_kind])
    _expect(
        errors,
        "processing_modes",
        observed["processing_modes"],
        [policy.processing_mode],
    )
    _expect(errors, "label_sources", observed["label_sources"], [policy.label_source])
    _expect(
        errors, "window_policies", observed["window_policies"], [policy.window_policy]
    )
    _expect(errors, "qc_statuses", observed["qc_statuses"], ["pass"])
    invalid_group_keys = [
        sample.group_key
        for sample in samples
        if not sample.group_key.startswith(policy.group_key_prefix)
    ]
    if invalid_group_keys:
        errors.append(
            f"group_key prefix mismatch: expected {policy.group_key_prefix!r}, "
            f"found={sorted(set(invalid_group_keys))}"
        )

    observed["source_provenance"] = _audit_source_provenance(
        dataset_dir=Path(dataset_dir),
        phase3_manifest=manifest,
        samples=samples,
        policy=policy,
        errors=errors,
    )

    return DatasetFreezeAudit(
        status="PASS" if not errors else "FAIL",
        dataset_dir=Path(dataset_dir),
        errors=tuple(errors),
        warnings=tuple(warnings),
        observed=observed,
        policy=policy,
    )


def assert_phase4_dataset_freeze(
    *,
    dataset_dir: Path,
    manifest: dict[str, Any],
    samples: list[Phase4ClipSample],
    policy: DatasetFreezePolicy,
) -> None:
    audit = audit_loaded_dataset_freeze(
        dataset_dir=dataset_dir,
        manifest=manifest,
        samples=samples,
        policy=policy,
    )
    if audit.status != "PASS":
        raise ValueError("dataset freeze failed: " + "; ".join(audit.errors))


def _observed_values(
    samples: list[Phase4ClipSample], *, group_count: int | None
) -> dict[str, Any]:
    def values(path: tuple[str, ...]) -> list[Any]:
        observed = set()
        for sample in samples:
            value: Any = sample.metadata
            for key in path:
                value = value.get(key) if isinstance(value, dict) else None
            observed.add(value)
        return sorted(observed, key=lambda item: str(item))

    return {
        "sample_count": len(samples),
        "group_count": group_count
        if group_count is not None
        else len({sample.group_key for sample in samples}),
        "n_flagella": sorted({sample.n_flagella for sample in samples}),
        "dataset_versions": values(("provenance", "dataset_version")),
        "model_ids": values(("provenance", "model_id")),
        "render_ids": values(("provenance", "render_id")),
        "source_kinds": values(("source_video", "source_kind")),
        "processing_modes": values(("processing_mode",)),
        "label_sources": values(("labels", "label_source")),
        "window_policies": values(("clip", "window_policy")),
        "qc_statuses": values(("qc", "status")),
        "splits": sorted({sample.split for sample in samples}),
    }


def _audit_source_provenance(
    *,
    dataset_dir: Path,
    phase3_manifest: dict[str, Any],
    samples: list[Phase4ClipSample],
    policy: DatasetFreezePolicy,
    errors: list[str],
) -> dict[str, Any]:
    observed: dict[str, Any] = {}
    source_dir = _resolve_existing_path(
        phase3_manifest.get("input_dataset"), relative_to=dataset_dir
    )
    if source_dir is None:
        errors.append(
            "source provenance unavailable: manifest.input_dataset does not exist"
        )
        return observed

    source_manifest_path = source_dir / "dataset_manifest.json"
    if not source_manifest_path.is_file():
        errors.append(
            f"source provenance unavailable: {source_manifest_path} does not exist"
        )
        return observed
    try:
        source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"source dataset manifest is invalid: {exc}")
        return observed

    effective_campaign = source_manifest.get("effective_campaign_config")
    if not isinstance(effective_campaign, dict):
        effective_campaign = {}
    effective_metadata = effective_campaign.get("metadata")
    if not isinstance(effective_metadata, dict):
        effective_metadata = {}
    observed.update(
        {
            "dataset_dir": str(source_dir),
            "dataset_manifest_path": str(source_manifest_path),
            "dataset_manifest_sha256": _sha256(source_manifest_path),
            "dataset_id": source_manifest.get("dataset_id"),
            "dataset_version": effective_metadata.get("dataset_version"),
            "model_id": effective_metadata.get("model_id"),
        }
    )
    _expect(
        errors,
        "source.dataset_id",
        source_manifest.get("dataset_id"),
        policy.dataset_version,
    )
    _expect(
        errors,
        "source.dataset_version",
        effective_metadata.get("dataset_version"),
        policy.dataset_version,
    )
    source_model_id = effective_metadata.get("model_id")
    if source_model_id not in policy.source_model_ids:
        errors.append(
            "source.model_id: "
            f"expected one of {sorted(policy.source_model_ids)!r}, "
            f"observed={source_model_id!r}"
        )

    raw_source_samples = source_manifest.get("samples")
    if not isinstance(raw_source_samples, list):
        errors.append("source.samples: expected a list in dataset_manifest.json")
        raw_source_samples = []
    source_samples = {
        record.get("sample_id"): record
        for record in raw_source_samples
        if isinstance(record, dict) and record.get("sample_id")
    }
    run_to_label: dict[str, int] = {}
    for sample in samples:
        run_id = sample.metadata.get("provenance", {}).get("run_id")
        if not isinstance(run_id, str) or not run_id:
            errors.append(f"clip {sample.clip_id}: provenance.run_id is required")
            continue
        previous = run_to_label.setdefault(run_id, sample.n_flagella)
        if previous != sample.n_flagella:
            errors.append(
                f"source run {run_id}: conflicting Phase 3 labels "
                f"{previous!r} and {sample.n_flagella!r}"
            )

    config_hashes: dict[str, str] = {}
    torques: set[float] = set()
    switching_values: set[bool] = set()
    brownian_values: set[bool] = set()
    source_classes: set[int] = set()
    for run_id, phase3_label in sorted(run_to_label.items()):
        record = source_samples.get(run_id)
        if record is None:
            errors.append(f"source run {run_id}: not found in dataset_manifest.json")
            continue
        _expect(
            errors, f"source run {run_id}.status", record.get("status"), "completed"
        )
        config_path = _resolve_existing_path(
            record.get("config_path"), relative_to=source_dir
        )
        if config_path is None or not config_path.is_file():
            errors.append(f"source run {run_id}: resolved config is unavailable")
            continue
        try:
            config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
        except (OSError, yaml.YAMLError) as exc:
            errors.append(f"source run {run_id}: resolved config is invalid: {exc}")
            continue
        if not isinstance(config, dict):
            errors.append(
                f"source run {run_id}: resolved config must contain a mapping"
            )
            continue

        config_hashes[run_id] = _sha256(config_path)
        switching = config.get("motor", {}).get("enable_switching")
        brownian = config.get("brownian", {}).get("enabled")
        torque = config.get("motor", {}).get("torque_Nm")
        source_n_flagella = config.get("flagella", {}).get("n_flagella")
        if isinstance(switching, bool):
            switching_values.add(switching)
        if isinstance(brownian, bool):
            brownian_values.add(brownian)
        try:
            torques.add(float(torque))
        except (TypeError, ValueError):
            pass
        try:
            source_classes.add(int(source_n_flagella))
        except (TypeError, ValueError):
            pass

        _expect(errors, f"source run {run_id}.motor.enable_switching", switching, False)
        _expect(errors, f"source run {run_id}.brownian.enabled", brownian, False)
        _expect_float(
            errors,
            f"source run {run_id}.motor.torque_Nm",
            torque,
            policy.baseline_torque_Nm,
        )
        _expect(
            errors,
            f"source run {run_id}.flagella.n_flagella",
            source_n_flagella,
            phase3_label,
        )
        _expect(
            errors,
            f"source run {run_id}.manifest.n_flagella",
            record.get("n_flagella"),
            phase3_label,
        )

    observed.update(
        {
            "selected_run_count": len(run_to_label),
            "resolved_config_count": len(config_hashes),
            "resolved_config_sha256": config_hashes,
            "torque_Nm": sorted(torques),
            "motor_enable_switching": sorted(switching_values),
            "brownian_enabled": sorted(brownian_values),
            "n_flagella": sorted(source_classes),
        }
    )
    return observed


def _resolve_existing_path(value: Any, *, relative_to: Path) -> Path | None:
    if not isinstance(value, str) or not value:
        return None
    path = Path(value).expanduser()
    candidates = (
        [path] if path.is_absolute() else [Path.cwd() / path, relative_to / path]
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _expect(errors: list[str], field: str, actual: Any, expected: Any) -> None:
    if actual != expected:
        errors.append(f"{field}: expected={expected!r}, observed={actual!r}")


def _expect_float(errors: list[str], field: str, actual: Any, expected: float) -> None:
    try:
        matches = np.isclose(float(actual), expected, rtol=1.0e-9, atol=0.0)
    except (TypeError, ValueError):
        matches = False
    if not matches:
        errors.append(f"{field}: expected={expected!r}, observed={actual!r}")
