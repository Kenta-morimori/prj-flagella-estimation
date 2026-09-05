from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
from sim_swim.analysis.cli_profiles import args_from_profile, load_profile

from sim_swim.analysis.issue61_2015_1tau import analyze
from sim_swim.analysis.sweeps import stage_a_2015


TORQUES = (1.0e-21, 2.5e-20, 1.0e-19)
ROOT = Path(__file__).parents[1]


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
