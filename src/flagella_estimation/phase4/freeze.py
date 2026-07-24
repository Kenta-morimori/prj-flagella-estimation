"""Machine-readable Phase 4 dataset freeze audit."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np

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
    behavior_regime: str = "RUN_fixed"
    brownian_policy: str = "excluded"
    torque_variation_policy: str = "excluded"
    n_flagella_4_policy: str = "diagnostic_only"


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
    warnings = [
        "Phase 3 manifest does not embed the raw Phase 2 config hash; "
        "physical-regime exclusions rely on dataset_version/model_id registry assertions."
    ]
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

    required_assertions = {
        "behavior_regime": (policy.behavior_regime, "RUN_fixed"),
        "brownian_policy": (policy.brownian_policy, "excluded"),
        "torque_variation_policy": (policy.torque_variation_policy, "excluded"),
        "n_flagella_4_policy": (policy.n_flagella_4_policy, "diagnostic_only"),
    }
    for name, (actual, expected) in required_assertions.items():
        _expect(errors, f"registry_assertion.{name}", actual, expected)

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
