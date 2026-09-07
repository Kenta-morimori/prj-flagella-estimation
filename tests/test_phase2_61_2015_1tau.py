from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
from sim_swim.analysis.cli_profiles import args_from_profile, load_profile

from sim_swim.analysis.issue61_2015_1tau import analyze
from sim_swim.analysis.issue61_2015_1tau import _first_threshold_crossing
from sim_swim.analysis.parallel_job import (
    _aggregate_stage_a_campaign,
    build_plan,
    load_parallel_job,
    resolve_execution,
)
from sim_swim.analysis.sweeps import stage_a_2015


TORQUES = (1.0e-21, 2.5e-20, 1.0e-19)
ROOT = Path(__file__).parents[1]
PARALLEL_JOB = ROOT / "conf/phase2_parallel/issue61_2015_1tau/job.yaml"


def test_issue61_profile_fixes_three_tracking_1tau_conditions() -> None:
    profile = load_profile(
        ROOT / "conf/phase2_sweeps/2015_issue61_1tau_tracking_stability.yaml"
    )
    args = stage_a_2015._parse_args(args_from_profile(profile))

    assert args.campaign_issue == 61
    assert args.profiles == ["project"]
    assert args.motor_torques_nm == list(TORQUES)
    assert args.link_reference_torque is True
    assert args.dt_star == pytest.approx(1.0e-5)
    assert args.duration_tau == pytest.approx(1.0)


def test_issue61_parallel_job_has_three_isolated_tracking_shards(
    tmp_path: Path,
) -> None:
    job = load_parallel_job(PARALLEL_JOB)
    execution = resolve_execution(job, None)
    plan = build_plan(job, execution, tmp_path / "parallel")

    assert job.is_stage_a_campaign_job
    assert execution.worker_policy == "cs10_qualified"
    assert execution.max_workers == 3
    assert [item["task_id"] for item in plan["configs"]] == [
        "project_torque_1em21",
        "project_torque_2p5em20",
        "project_torque_1em19",
    ]
    assert len({item["output_dir"] for item in plan["configs"]}) == 3
    assert [item["overrides"] for item in plan["configs"]] == [
        ["motor_torques_nm=1.0e-21"],
        ["motor_torques_nm=2.5e-20"],
        ["motor_torques_nm=1.0e-19"],
    ]


def test_stage_a_torque_condition_ids_are_lossless() -> None:
    assert stage_a_2015._physical_torque_condition_id("project", 1.0e-21, 3) == (
        "project_torque_1em21"
    )
    assert stage_a_2015._physical_torque_condition_id("project", 2.5e-20, 3) == (
        "project_torque_2p5em20"
    )


def test_reference_evidence_hashes_the_source_manifest(tmp_path: Path) -> None:
    manifest = tmp_path / "source-manifest.json"
    manifest.write_text('{"source": "fixture"}\n', encoding="utf-8")
    evidence = tmp_path / "evidence.json"
    evidence.write_text(
        json.dumps(
            [
                {
                    "label": "fixture",
                    "source_run_root": "/source/run",
                    "manifest_path": str(manifest),
                }
            ]
        ),
        encoding="utf-8",
    )

    records = stage_a_2015._reference_evidence(evidence)
    assert records[0]["source_run_root"] == "/source/run"
    assert len(records[0]["manifest_sha256"]) == 64


def _campaign(tmp_path: Path, *, fail_metric: str | None = None) -> Path:
    root = tmp_path / "campaign"
    root.mkdir()
    conditions = []
    rows = []
    for index, torque in enumerate(TORQUES):
        condition_id = f"project_torque_{index}"
        condition_dir = root / condition_id
        condition_dir.mkdir()
        (condition_dir / "trajectory.csv").write_text("t_s\n0\n", encoding="utf-8")
        (condition_dir / "state_archive.npz").write_bytes(b"fixture")
        (condition_dir / "run_summary.json").write_text(
            json.dumps(
                {
                    "execution": {"status": "completed"},
                    "gates": {
                        name: {"status": "available", "any_fail": False}
                        for name in ("finite", "shape_nonbody", "shape_body")
                    },
                }
            ),
            encoding="utf-8",
        )
        row = {
            "condition_id": condition_id,
            "status": "completed",
            "completion_pass": "True",
            "finite_pass_all": "True",
            "wall_time_s": "10.0",
            "steps_per_s": "100.0",
        }
        for metric in (
            "body_spring_max_stretch_ratio",
            "body_length_rel_drift_max",
            "body_width_rel_drift_max",
            "body_cross_section_area_rel_drift_max",
            "max_flag_bond_rel_err",
            "max_hook_len_rel_err",
            "max_hook_angle_err_deg",
            "max_flag_bend_err_deg",
            "max_flag_torsion_err_deg",
            "max_flag_helix_radius_abs_err_over_b",
            "max_flag_helix_pitch_rel_err",
            "max_motor_force_balance_residual_ratio",
            "max_motor_torque_balance_residual_ratio",
        ):
            row[metric] = "0.0"
        if index == 0 and fail_metric:
            row[fail_metric] = "999.0"
        rows.append(row)
        conditions.append(
            {
                "condition_id": condition_id,
                "motor_torque_Nm": torque,
                "output_dir": str(condition_dir),
                "time": {"tau_s": 1.0, "dt_internal_s": 1e-5, "total_steps": 1_000_000},
                "config_overrides": {
                    "motor.torque_Nm": torque,
                    "motor.reference_torque_Nm": torque,
                    "time.scale_policy": "reference_torque",
                },
            }
        )
    with (root / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (root / "run_manifest.json").write_text(
        json.dumps(
            {
                "kind": "stage_a_2015",
                "issue": 61,
                "stage": "motor_on",
                "profiles": ["project"],
                "duration_tau": 1.0,
                "dt_star": 1e-5,
                "link_reference_torque": True,
                "reference_evidence": [{"manifest_sha256": "a" * 64}],
                "conditions": conditions,
            }
        ),
        encoding="utf-8",
    )
    return root


def test_stage_a_parallel_aggregate_builds_issue61_analysis_input(
    tmp_path: Path,
) -> None:
    job = load_parallel_job(PARALLEL_JOB)
    root = tmp_path / "parallel"
    manifest = build_plan(job, resolve_execution(job, None), root)
    root.mkdir()
    manifest["output_root"] = str(root)
    for record, torque in zip(manifest["configs"], TORQUES, strict=True):
        child_root = Path(record["output_dir"]) / "2026-09-06" / record["task_id"]
        # A single-torque shard intentionally keeps the historical child ID.
        condition_id = "project"
        condition_dir = child_root / condition_id
        condition_dir.mkdir(parents=True)
        row = {
            "condition_id": condition_id,
            "status": "completed",
            "completion_pass": "True",
            "finite_pass_all": "True",
            "wall_time_s": "10.0",
            "steps_per_s": "100.0",
        }
        for metric in (
            "body_spring_max_stretch_ratio",
            "body_length_rel_drift_max",
            "body_width_rel_drift_max",
            "body_cross_section_area_rel_drift_max",
            "max_flag_bond_rel_err",
            "max_hook_len_rel_err",
            "max_hook_angle_err_deg",
            "max_flag_bend_err_deg",
            "max_flag_torsion_err_deg",
            "max_flag_helix_radius_abs_err_over_b",
            "max_flag_helix_pitch_rel_err",
            "max_motor_force_balance_residual_ratio",
            "max_motor_torque_balance_residual_ratio",
        ):
            row[metric] = "0.0"
        with (child_root / "summary.csv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        (condition_dir / "run_summary.json").write_text(
            json.dumps(
                {
                    "execution": {"status": "completed"},
                    "gates": {
                        "finite": {
                            "status": "available",
                            "any_fail": False,
                            "final_pass": True,
                        },
                        "shape_nonbody": {
                            "status": "available",
                            "any_fail": False,
                            "final_pass": True,
                        },
                        "shape_body": {
                            "status": "available",
                            "any_fail": False,
                            "final_pass": True,
                        },
                    },
                }
            ),
            encoding="utf-8",
        )
        condition = {
            "condition_id": condition_id,
            "profile": "project",
            "motor_torque_Nm": torque,
            "output_dir": str(condition_dir),
            "config_overrides": {
                "motor.torque_Nm": torque,
                "motor.reference_torque_Nm": torque,
                "time.scale_policy": "reference_torque",
            },
            "time": {"tau_s": 1.0, "dt_internal_s": 1.0e-5, "total_steps": 100000},
        }
        child = {
            "kind": "stage_a_2015",
            "issue": 61,
            "stage": "motor_on",
            "duration_tau": 1.0,
            "dt_star": 1.0e-5,
            "comparison_role": "issue61_2015_project_1tau_tracking_stability",
            "motor_enabled": True,
            "diagonal_braces_enabled": False,
            "link_reference_torque": True,
            "base_config": "conf/sim_swim_2015.yaml",
            "reference_evidence": [],
            "motor_torques_Nm": [torque],
            "conditions": [condition],
            "performance_json": str(child_root / "performance.json"),
        }
        (child_root / "run_manifest.json").write_text(
            json.dumps(child), encoding="utf-8"
        )
        (child_root / "performance.json").write_text(
            json.dumps({"kind": "stage_a_2015_performance", "conditions": [{}]}),
            encoding="utf-8",
        )
        record["status"] = "succeeded"

    campaign = _aggregate_stage_a_campaign(job, manifest)
    assert (campaign / "conditions/project_torque_2p5em20").is_symlink()
    aggregate_manifest = json.loads((campaign / "run_manifest.json").read_text())
    assert [item["condition_id"] for item in aggregate_manifest["conditions"]] == [
        "project_torque_1em21",
        "project_torque_2p5em20",
        "project_torque_1em19",
    ]
    assert all(
        item["source_condition_id"] == "project"
        for item in aggregate_manifest["conditions"]
    )
    result = analyze(
        run_root=campaign,
        threshold_contract=ROOT / "conf/phase2_validation/2015_stage_a_thresholds.yaml",
        output_dir=tmp_path / "analysis",
    )
    assert (
        json.loads((result / "issue61_decision.json").read_text())["strict_pass_count"]
        == 3
    )


def test_stage_a_parallel_aggregate_rejects_task_torque_mismatch(
    tmp_path: Path,
) -> None:
    job = load_parallel_job(PARALLEL_JOB)
    root = tmp_path / "parallel"
    manifest = build_plan(job, resolve_execution(job, None), root)
    root.mkdir()
    manifest["output_root"] = str(root)
    manifest["configs"][0]["task_id"] = "project_torque_1em19"
    for record, torque in zip(manifest["configs"], TORQUES, strict=True):
        child_root = Path(record["output_dir"]) / "2026-09-06" / "child"
        condition_dir = child_root / "project"
        condition_dir.mkdir(parents=True)
        (condition_dir / "run_summary.json").write_text(
            json.dumps({"execution": {"status": "completed"}}), encoding="utf-8"
        )
        (child_root / "summary.csv").write_text(
            "condition_id,status\nproject,completed\n", encoding="utf-8"
        )
        (child_root / "performance.json").write_text(
            json.dumps({"conditions": [{}]}), encoding="utf-8"
        )
        (child_root / "run_manifest.json").write_text(
            json.dumps(
                {
                    "kind": "stage_a_2015",
                    "issue": 61,
                    "stage": "motor_on",
                    "duration_tau": 1.0,
                    "dt_star": 1e-5,
                    "comparison_role": "issue61_2015_project_1tau_tracking_stability",
                    "motor_enabled": True,
                    "diagonal_braces_enabled": False,
                    "link_reference_torque": True,
                    "base_config": "conf/sim_swim_2015.yaml",
                    "reference_evidence": [],
                    "performance_json": str(child_root / "performance.json"),
                    "conditions": [
                        {
                            "condition_id": "project",
                            "profile": "project",
                            "motor_torque_Nm": torque,
                            "output_dir": str(condition_dir),
                            "config_overrides": {
                                "motor.torque_Nm": torque,
                                "motor.reference_torque_Nm": torque,
                                "time.scale_policy": "reference_torque",
                            },
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        record["status"] = "succeeded"
    with pytest.raises(RuntimeError, match="task/torque mismatch"):
        _aggregate_stage_a_campaign(job, manifest)


def test_issue61_analysis_records_all_pass_and_blocks_promotion(tmp_path: Path) -> None:
    root = _campaign(tmp_path)
    output = analyze(
        run_root=root,
        threshold_contract=Path("conf/phase2_validation/2015_stage_a_thresholds.yaml"),
        output_dir=tmp_path / "analysis",
    )
    decision = json.loads(
        (output / "issue61_decision.json").read_text(encoding="utf-8")
    )
    assert decision["status"] == "pass"
    assert decision["strict_pass_count"] == 3
    assert "no supported-profile promotion" in decision["handoff"]


def test_issue61_analysis_records_first_threshold_failure(tmp_path: Path) -> None:
    root = _campaign(tmp_path, fail_metric="max_flag_bond_rel_err")
    output = analyze(
        run_root=root,
        threshold_contract=Path("conf/phase2_validation/2015_stage_a_thresholds.yaml"),
        output_dir=tmp_path / "analysis",
    )
    rows = list(csv.DictReader((output / "issue61_summary.csv").open(encoding="utf-8")))
    assert rows[0]["first_failing_criterion"] == "max_flag_bond_rel_err"
    decision = json.loads(
        (output / "issue61_decision.json").read_text(encoding="utf-8")
    )
    assert decision["status"] == "fail"
    assert "do not promote" in decision["handoff"]


def test_issue61_streams_first_body_drift_crossing(tmp_path: Path) -> None:
    diagnostics = tmp_path / "body_constraint_diagnostics.csv"
    diagnostics.write_text(
        "step,t_s,body_length_um\n0,0.0,2.0\n10,0.1,2.01\n20,0.2,2.03\n",
        encoding="utf-8",
    )
    assert _first_threshold_crossing(
        tmp_path, criterion="body_length_rel_drift_max", limit=0.01
    ) == {"criterion": "body_length_rel_drift_max", "t_s": 0.2, "step": 20}


def test_issue61_analysis_rejects_non_tracking_manifest(tmp_path: Path) -> None:
    root = _campaign(tmp_path)
    manifest_path = root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["link_reference_torque"] = False
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="link_reference_torque"):
        analyze(
            run_root=root,
            threshold_contract=Path(
                "conf/phase2_validation/2015_stage_a_thresholds.yaml"
            ),
            output_dir=tmp_path / "analysis",
        )
