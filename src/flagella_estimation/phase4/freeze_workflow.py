"""CLI-oriented workflow for Phase 4 dataset freeze audits."""

from __future__ import annotations

from dataclasses import asdict, dataclass
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

from flagella_estimation.phase4.freeze import (
    DatasetFreezeAudit,
    DatasetFreezePolicy,
    audit_phase4_dataset_freeze,
)


@dataclass(frozen=True)
class Phase4FreezeAuditConfig:
    dataset_dir: Path
    output_dir: Path
    policy: DatasetFreezePolicy
    config_path: Path | None = None
    cli_overrides: tuple[str, ...] = ()


def default_output_dir() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / "phase4_dataset_freeze_audit"
    )


def load_freeze_audit_config(
    path: Path | None, overrides: list[str] | None = None
) -> Phase4FreezeAuditConfig:
    raw: dict[str, Any] = {}
    if path is not None:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    data = _apply_overrides(raw, overrides or [])
    freeze = dict(data.get("freeze", {}) or {})
    policy = DatasetFreezePolicy(
        allowed_n_flagella=tuple(
            int(value) for value in freeze.get("allowed_n_flagella", [1, 2, 3])
        ),
        dataset_version=str(freeze.get("dataset_version", "v1")),
        model_ids=tuple(
            str(value)
            for value in freeze.get("model_ids", ["phase2_flagella_count_behavior_v1"])
        ),
        source_model_ids=tuple(
            str(value)
            for value in freeze.get(
                "source_model_ids", ["flag_spring2p25_body2p5_candidate"]
            )
        ),
        render_ids=tuple(
            str(value) for value in freeze.get("render_ids", ["state_archive_numpy_v1"])
        ),
        pipeline_name=str(freeze.get("pipeline_name", "phase3_gt_passthrough")),
        schema_version=str(freeze.get("schema_version", "phase3_clip_metadata/v0")),
        source_kind=str(freeze.get("source_kind", "phase2_pseudo")),
        processing_mode=str(freeze.get("processing_mode", "gt_passthrough")),
        label_source=str(freeze.get("label_source", "phase2_gt")),
        clip_duration_s=float(freeze.get("clip_duration_s", 0.5)),
        window_policy=str(freeze.get("window_policy", "non_overlap")),
        baseline_torque_Nm=float(freeze.get("baseline_torque_Nm", 2.0e-20)),
        require_use_for_ml_candidate=bool(
            freeze.get("require_use_for_ml_candidate", True)
        ),
        group_key_prefix=str(freeze.get("group_key_prefix", "phase2:v1:")),
    )
    return Phase4FreezeAuditConfig(
        dataset_dir=Path(str(data.get("dataset_dir", ""))),
        output_dir=Path(str(data.get("output_dir") or default_output_dir())),
        policy=policy,
        config_path=path,
        cli_overrides=tuple(overrides or ()),
    )


def run_freeze_audit(
    cfg: Phase4FreezeAuditConfig,
) -> tuple[Path, DatasetFreezeAudit]:
    """Run the audit and write PASS or FAIL artifacts."""

    if cfg.output_dir.resolve() == cfg.dataset_dir.resolve():
        raise ValueError("output_dir must differ from dataset_dir")
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    audit = audit_phase4_dataset_freeze(cfg.dataset_dir, cfg.policy)
    created_at = datetime.now(ZoneInfo("Asia/Tokyo")).isoformat()
    audit_path = cfg.output_dir / "freeze_audit.json"
    audit_path.write_text(
        json.dumps(audit.to_dict(), ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    manifest = {
        "pipeline_name": "phase4_dataset_freeze_audit",
        "pipeline_version": "0.1.0",
        "created_at": created_at,
        "status": audit.status,
        "input_dataset": str(cfg.dataset_dir),
        "output_dir": str(cfg.output_dir),
        "policy": asdict(cfg.policy),
        "invocation": {
            "config_path": str(cfg.config_path) if cfg.config_path else None,
            "cli_overrides": list(cfg.cli_overrides),
        },
        "outputs": {"freeze_audit": str(audit_path)},
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
                f"status={audit.status}",
                f"error_count={len(audit.errors)}",
                f"warning_count={len(audit.warnings)}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return cfg.output_dir, audit


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
