from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

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
    assert audit.warnings == ()
    provenance = audit.observed["source_provenance"]
    assert provenance["selected_run_count"] == 9
    assert provenance["resolved_config_count"] == 9
    assert provenance["motor_enable_switching"] == [False]
    assert provenance["brownian_enabled"] == [False]
    assert provenance["torque_Nm"] == [2.0e-20]
    assert len(provenance["dataset_manifest_sha256"]) == 64
    assert all(
        len(digest) == 64 for digest in provenance["resolved_config_sha256"].values()
    )
    assert {path.name for path in output_dir.iterdir()} == {
        "freeze_audit.json",
        "manifest.json",
        "run.log",
    }
    saved = json.loads((output_dir / "freeze_audit.json").read_text(encoding="utf-8"))
    assert saved["status"] == "PASS"
    assert saved["policy"]["source_model_ids"] == ["flag_spring2p25_body2p5_candidate"]


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
def test_dataset_freeze_audit_rejects_brownian_source_config(
    tmp_path: Path,
) -> None:
    dataset_dir = tmp_path / "phase3_dataset"
    write_phase4_fixture_dataset(dataset_dir)
    source_dir = dataset_dir.parent / f"{dataset_dir.name}_phase2_source"
    config_path = source_dir / "configs" / "nf01_train.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    config["brownian"]["enabled"] = True
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    _, audit = run_freeze_audit(
        Phase4FreezeAuditConfig(
            dataset_dir=dataset_dir,
            output_dir=tmp_path / "freeze_audit",
            policy=DatasetFreezePolicy(),
        )
    )

    assert audit.status == "FAIL"
    assert any("brownian.enabled" in error for error in audit.errors)


@pytest.mark.light
def test_dataset_freeze_audit_rejects_missing_source_manifest(
    tmp_path: Path,
) -> None:
    dataset_dir = tmp_path / "phase3_dataset"
    write_phase4_fixture_dataset(dataset_dir)
    source_dir = dataset_dir.parent / f"{dataset_dir.name}_phase2_source"
    (source_dir / "dataset_manifest.json").unlink()

    _, audit = run_freeze_audit(
        Phase4FreezeAuditConfig(
            dataset_dir=dataset_dir,
            output_dir=tmp_path / "freeze_audit",
            policy=DatasetFreezePolicy(),
        )
    )

    assert audit.status == "FAIL"
    assert any("source provenance unavailable" in error for error in audit.errors)


@pytest.mark.light
def test_load_freeze_audit_config_supports_source_model_ids(
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
            "freeze.source_model_ids=[model_a, model_b]",
        ],
    )
    assert cfg.dataset_dir == Path("override_input")
    assert cfg.output_dir == Path("override_output")
    assert cfg.policy.source_model_ids == ("model_a", "model_b")
