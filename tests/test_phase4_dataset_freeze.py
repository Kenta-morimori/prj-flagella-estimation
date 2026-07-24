from __future__ import annotations

import json
from pathlib import Path

import pytest

from phase4_test_utils import write_phase4_fixture_dataset
from flagella_estimation.phase4.freeze import DatasetFreezePolicy
from flagella_estimation.phase4.freeze_workflow import (
    Phase4FreezeAuditConfig,
    load_freeze_audit_config,
    run_freeze_audit,
)


@pytest.mark.light
def test_dataset_freeze_audit_writes_pass_artifacts(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "phase3_dataset"
    output_dir = tmp_path / "freeze_audit"
    write_phase4_fixture_dataset(dataset_dir)

    result_dir, audit = run_freeze_audit(
        Phase4FreezeAuditConfig(
            dataset_dir=dataset_dir,
            output_dir=output_dir,
            policy=DatasetFreezePolicy(),
        )
    )

    assert result_dir == output_dir
    assert audit.status == "PASS"
    assert audit.errors == ()
    assert audit.observed["n_flagella"] == [1, 2, 3]
    assert audit.observed["group_count"] == 9
    assert len(audit.warnings) == 1
    assert {path.name for path in output_dir.iterdir()} == {
        "freeze_audit.json",
        "manifest.json",
        "run.log",
    }
    saved = json.loads((output_dir / "freeze_audit.json").read_text(encoding="utf-8"))
    assert saved["status"] == "PASS"
    assert saved["policy"]["behavior_regime"] == "RUN_fixed"
    assert saved["policy"]["brownian_policy"] == "excluded"


@pytest.mark.light
def test_dataset_freeze_audit_records_all_policy_failures(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "phase3_dataset"
    output_dir = tmp_path / "freeze_audit"
    write_phase4_fixture_dataset(dataset_dir, dataset_version="v2")
    manifest_path = dataset_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["filters"]["baseline_torque_Nm"] = 3.0e-20
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    _, audit = run_freeze_audit(
        Phase4FreezeAuditConfig(
            dataset_dir=dataset_dir,
            output_dir=output_dir,
            policy=DatasetFreezePolicy(),
        )
    )

    assert audit.status == "FAIL"
    assert any("baseline_torque_Nm" in error for error in audit.errors)
    assert any("dataset_versions" in error for error in audit.errors)
    assert any("group_key prefix mismatch" in error for error in audit.errors)
    saved = json.loads((output_dir / "freeze_audit.json").read_text(encoding="utf-8"))
    assert saved["status"] == "FAIL"
    assert len(saved["errors"]) >= 3


@pytest.mark.light
def test_load_freeze_audit_config_supports_registry_assertions(
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "freeze.yaml"
    config_path.write_text(
        "dataset_dir: input\noutput_dir: output\nfreeze: {}\n",
        encoding="utf-8",
    )
    cfg = load_freeze_audit_config(
        config_path,
        [
            "dataset_dir=override_input",
            "output_dir=override_output",
            "freeze.registry_assertions.behavior_regime=RUN_fixed",
            "freeze.registry_assertions.brownian=excluded",
        ],
    )
    assert cfg.dataset_dir == Path("override_input")
    assert cfg.output_dir == Path("override_output")
    assert cfg.policy.behavior_regime == "RUN_fixed"
    assert cfg.policy.brownian_policy == "excluded"
