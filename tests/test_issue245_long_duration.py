from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest

from sim_swim.analysis.model_development_evaluation import (
    _render_replays,
    build_evaluation,
    collect_rows,
)
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.analysis.parallel_job import load_parallel_job, resolve_execution


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "conf/phase2_multi_run/2010_hex_project_long_duration_2s_issue245.yaml"
SCREEN_CONFIG = (
    ROOT
    / "conf/phase2_multi_run/2010_hex_project_attachment_patterns_1tau_issue245.yaml"
)
SCREEN_JOB = (
    ROOT
    / "conf/phase2_parallel/issue245_2010_hex_long_duration/attachment_screen_job.yaml"
)
MAIN_JOB = ROOT / "conf/phase2_parallel/issue245_2010_hex_long_duration/job.yaml"


def _load_execution_target_module():
    path = ROOT / "tools/codex/issue_execution_target.py"
    spec = importlib.util.spec_from_file_location("issue_execution_target_245", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _summary() -> dict:
    metrics = {
        "hook_angle_err_max_deg": {"max": 61.0},
        "local_attach_first_rel_err": {"max": 0.01},
        "hook_len_rel_err_max": {"max": 0.01},
        "flag_bond_rel_err_max": {"max": 0.03},
        "flag_bend_err_max_deg": {"max": 1.0},
        "flag_torsion_err_max_deg": {"max": 1.0},
        "body_spring_max_stretch_ratio": {"max": 0.04},
        "motor_force_balance_residual_ratio": {"max": 1.0e-12},
        "motor_torque_balance_residual_ratio": {"max": 0.01},
    }
    return {
        "execution": {"status": "completed"},
        "gates": {
            "finite": {"any_fail": False},
            "shape_body": {"any_fail": False},
            "shape_nonbody": {
                "any_fail": True,
                "first_failure_category": "hook",
                "first_failure_t_s": 0.00004,
            },
        },
        "all_step_metrics": metrics,
    }


def _write_long_run(root: Path) -> Path:
    config = normalize_campaign_config(load_yaml(CONFIG))
    run_dir = root / "synchronized"
    run_dir.mkdir(parents=True)
    records, rows = [], []
    for condition in build_campaign_conditions(config):
        directory = run_dir / "conditions" / condition["condition_id"]
        directory.mkdir(parents=True)
        (directory / "run_summary.json").write_text(json.dumps(_summary()))
        (directory / "performance.json").write_text(
            json.dumps({"wall_time_s": 10.0, "steps_per_s": 100.0})
        )
        (directory / "state_archive.npz").write_bytes(b"portable archive")
        record = {
            **condition,
            "output_dir": f"/cs10/original/{condition['condition_id']}",
            "time": {"dt_star": 1.0e-4, "dt_internal_s": 4.0e-6, "duration_s": 2.0},
        }
        records.append(record)
        rows.append({"condition_id": condition["condition_id"], "completed": "True"})
    with (run_dir / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["condition_id", "completed"])
        writer.writeheader()
        writer.writerows(rows)
    (run_dir / "run_manifest.json").write_text(
        json.dumps(
            {
                "base_config": "conf/sim_swim_2010_hex.yaml",
                "campaign_config": str(CONFIG.relative_to(ROOT)),
                "model_profile": load_yaml(ROOT / "conf/sim_swim_2010_hex.yaml")[
                    "model_profile"
                ],
                "git": {"commit": "test", "is_clean": True},
                "conditions": records,
            }
        )
    )
    return run_dir


def test_issue245_static_condition_contract() -> None:
    config = normalize_campaign_config(load_yaml(CONFIG))
    conditions = build_campaign_conditions(config)
    assert len(conditions) == 13
    assert [c["condition_id"] for c in conditions] == [
        "nf01__slots0",
        "nf02__slots01",
        "nf02__slots02",
        "nf02__slots03",
        "nf03__slots012",
        "nf03__slots013",
        "nf03__slots014",
        "nf03__slots024",
        "nf04__slots0123",
        "nf04__slots0124",
        "nf04__slots0134",
        "nf05__slots01234",
        "nf06__slots012345",
    ]
    assert [c["axis_values"]["attachment_slots"] for c in conditions] == [
        [0],
        [0, 1],
        [0, 2],
        [0, 3],
        [0, 1, 2],
        [0, 1, 3],
        [0, 1, 4],
        [0, 2, 4],
        [0, 1, 2, 3],
        [0, 1, 2, 4],
        [0, 1, 3, 4],
        [0, 1, 2, 3, 4],
        [0, 1, 2, 3, 4, 5],
    ]
    overrides = config["base_overrides"]
    assert overrides["time"]["duration"] == {"value": 2.0, "unit": "s"}
    assert overrides["time"]["integration"]["dt_star"] == 1.0e-4
    assert overrides["motor"]["torque_Nm"] == 2.5e-20
    assert overrides["output"]["archive_interval_s"] == 0.001
    assert config["output"]["save_state_archive"] is True
    assert config["base_overrides"]["seed"]["phase_seed"] == 0
    assert config["development_evaluation"]["long_duration"]["expected_steps"] == 500000


def test_issue245_cs10_issue_form_metadata_maps_to_execution_label() -> None:
    module = _load_execution_target_module()
    body = "\n".join(
        (
            "### Heavy/runtime execution target",
            "",
            "cs10_user_run",
            "",
            "### cs10 execution mode",
            "",
            "parallel",
            "",
            "### Parallel job config path",
            "",
            "conf/phase2_parallel/issue245_2010_hex_long_duration/job.yaml",
            "",
            "### Parallel worker plan",
            "",
            "max_workers: auto (effective 8 workers), cs10_qualified",
        )
    )
    assert module.execution_label_from_issue_body(body) == "execution:cs10"
    assert "max_workers: auto (effective 8 workers), cs10_qualified" in body


def test_issue245_parallel_jobs_use_cs10_qualified_auto_workers() -> None:
    qualification = load_parallel_job(SCREEN_JOB)
    main = load_parallel_job(MAIN_JOB)
    assert qualification.condition_ids == main.condition_ids
    for job in (qualification, main):
        execution = resolve_execution(job, None)
        assert job.max_workers == "auto"
        assert execution.worker_policy == "cs10_qualified"
        assert execution.max_workers == 8


def test_long_duration_collects_portable_archives_and_window_qc(tmp_path: Path) -> None:
    run_dir = _write_long_run(tmp_path)
    rows, _ = collect_rows(config=load_yaml(CONFIG), run_dirs=[run_dir])
    assert len(rows) == 13
    assert {row["screen_status"] for row in rows} == {"pass"}
    assert sum(bool(row["full_ring_rotation_equivalent"]) for row in rows) == 1
    assert all(len(row["state_archive_sha256"]) == 64 for row in rows)
    outputs = build_evaluation(
        config_path=CONFIG, run_dirs=[run_dir], output_dir=tmp_path / "evaluation"
    )
    assert outputs["window_qc_csv"].is_file()
    assert outputs["attachment_pattern_heatmap"].is_file()
    manifest = json.loads(outputs["manifest"].read_text())
    assert manifest["full_ring_rotation_equivalent_condition_ids"] == [
        "nf06__slots012345",
    ]


def test_attachment_pattern_replay_pages_each_count(
    monkeypatch, tmp_path: Path
) -> None:
    run_dir = _write_long_run(tmp_path)
    rows, _ = collect_rows(config=load_yaml(CONFIG), run_dirs=[run_dir])
    calls: list[list[str]] = []

    def fake_replay_main(args: list[str]) -> None:
        calls.append(args)

    import sim_swim.analysis.phase2_replay as replay_module

    monkeypatch.setattr(replay_module, "main", fake_replay_main)
    _render_replays(
        rows,
        replay_input=tmp_path / "replay_input",
        output_dir=tmp_path / "evaluation",
        stage="long_duration",
    )
    assert len(calls) == 6
    selected = [
        arg
        for call in calls
        for arg in call
        if arg.startswith("nf") and "__slots" in arg
    ]
    assert selected == [row["condition_id"] for row in rows]
    assert all("--camera-3d" in call and "--camera-2d" in call for call in calls)


def test_long_duration_rejects_missing_or_mismatched_archive(tmp_path: Path) -> None:
    run_dir = _write_long_run(tmp_path)
    archive = run_dir / "conditions/nf01__slots0/state_archive.npz"
    archive.unlink()
    with pytest.raises(FileNotFoundError, match="state_archive.npz"):
        collect_rows(config=load_yaml(CONFIG), run_dirs=[run_dir])

    run_dir = _write_long_run(tmp_path / "again")
    manifest_path = run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    record = manifest["conditions"][0]
    archive = run_dir / "conditions" / record["condition_id"] / "state_archive.npz"
    record["artifact_sha256"] = {
        "state_archive.npz": hashlib.sha256(b"wrong").hexdigest()
    }
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        collect_rows(config=load_yaml(CONFIG), run_dirs=[run_dir])


def test_long_duration_rejects_partial_only_checkpoint(tmp_path: Path) -> None:
    run_dir = _write_long_run(tmp_path)
    condition = run_dir / "conditions/nf01__slots0"
    (condition / "state_archive.npz").unlink()
    (condition / "state_archive.partial.npz").write_bytes(b"checkpoint only")
    (condition / "progress.json").write_text(
        json.dumps({"status": "partial"}), encoding="utf-8"
    )
    (condition / "diagnostic_samples.csv").write_text("step\n25\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Partial checkpoint"):
        collect_rows(config=load_yaml(CONFIG), run_dirs=[run_dir])


def test_long_duration_manifest_marks_partial_artifacts_diagnostic_only(
    tmp_path: Path,
) -> None:
    run_dir = _write_long_run(tmp_path)
    outputs = build_evaluation(
        config_path=CONFIG, run_dirs=[run_dir], output_dir=tmp_path / "evaluation"
    )
    manifest = json.loads(outputs["manifest"].read_text())
    policy = manifest["long_duration_artifact_policy"]
    assert policy["required_completed_artifacts"] == [
        "run_summary.json",
        "performance.json",
        "state_archive.npz",
    ]
    assert policy["partial_checkpoint_policy"].startswith("diagnostic-only")
    assert "run.log" in policy["excluded_operational_logs"]
