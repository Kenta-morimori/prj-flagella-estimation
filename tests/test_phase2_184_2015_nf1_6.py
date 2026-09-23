from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from sim_swim.analysis.issue184_2015_nf1_6 import analyze


ROOT = Path(__file__).parents[1]
IDS = tuple(f"nf{index:02d}" for index in range(1, 7))


def _campaign(tmp_path: Path, *, fail: bool = False) -> Path:
    root = tmp_path / "campaign"
    conditions_root = root / "conditions"
    conditions_root.mkdir(parents=True)
    conditions = []
    rows = []
    for index, condition_id in enumerate(IDS, start=1):
        condition_dir = conditions_root / condition_id
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
        (condition_dir / "body_constraint_diagnostics.csv").write_text(
            "step,t_s,body_spring_max_stretch_ratio,body_length_um,body_width_mean_um,body_width_min_um,body_width_max_um,body_cross_section_area_min_um2,body_cross_section_area_max_um2\n"
            "0,0,0,2,1,1,1,1,1\n",
            encoding="utf-8",
        )
        (condition_dir / "step_summary.csv").write_text(
            "step,t_s,flag_bond_rel_err_max,hook_len_rel_err_max,hook_angle_err_max_deg,flag_bend_err_max_deg,flag_torsion_err_max_deg,flag_helix_radius_abs_err_over_b_max,flag_helix_pitch_rel_err_max,motor_force_balance_residual_ratio,motor_torque_balance_residual_ratio\n"
            f"0,0,0,0,0,0,0,0,{1 if fail and index == 1 else 0},0,0\n",
            encoding="utf-8",
        )
        rows.append(
            {
                "condition_id": condition_id,
                "wall_time_s": "10",
                "steps_per_s": "100",
                "body_spring_max_stretch_ratio": "0",
                "flag_bond_rel_err_max": "0",
                "hook_len_rel_err_max": "0",
                "hook_angle_err_max_deg": "0",
                "flag_bend_err_max_deg": "0",
                "flag_torsion_err_max_deg": "0",
                "flag_helix_radius_abs_err_over_b_max": "0",
                "flag_helix_pitch_rel_err_max": "1" if fail and index == 1 else "0",
                "motor_force_balance_residual_ratio": "0",
                "motor_torque_balance_residual_ratio": "0",
            }
        )
        attachments = [
            {
                "flag_id": flag_id,
                "body_bead_index": 12 + flag_id,
                "layer": 2,
                "slot": flag_id,
            }
            for flag_id in range(index)
        ]
        conditions.append(
            {
                "condition_id": condition_id,
                "output_dir": str(condition_dir),
                "axis_values": {"n_flagella": index},
                "config_overrides": {
                    "flagella": {"placement_mode": "seeded_surface"},
                    "seed": {"attach_seed": 0, "phase_seed": 0},
                    "time": {"scale_policy": "reference_torque"},
                    "motor": {"torque_Nm": 2.5e-20, "reference_torque_Nm": 2.5e-20},
                },
                "geometry_preflight": {
                    "placement_mode": "seeded_surface",
                    "attach_seed": 0,
                    "phase_seed": 0,
                    "attachments": attachments,
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
                "kind": "generic_multi_run",
                "parallel_aggregate": True,
                "condition_order": list(IDS),
                "conditions": conditions,
            }
        ),
        encoding="utf-8",
    )
    (root / "campaign_completion.json").write_text(
        '{"status":"completed","exit_code":0}\n', encoding="utf-8"
    )
    return root


def test_issue184_analysis_accepts_complete_seeded_surface_campaign(
    tmp_path: Path,
) -> None:
    root = _campaign(tmp_path)
    output = analyze(
        run_root=root,
        threshold_contract=ROOT / "conf/phase2_validation/2015_stage_a_thresholds.yaml",
        output_dir=tmp_path / "analysis",
    )
    decision = json.loads((output / "issue184_decision.json").read_text())
    assert decision["status"] == "pass"
    assert decision["strict_pass_count"] == 6


def test_issue184_analysis_records_strict_failure_and_rejects_bad_topology(
    tmp_path: Path,
) -> None:
    root = _campaign(tmp_path, fail=True)
    output = analyze(
        run_root=root,
        threshold_contract=ROOT / "conf/phase2_validation/2015_stage_a_thresholds.yaml",
        output_dir=tmp_path / "analysis",
    )
    rows = list(csv.DictReader((output / "issue184_summary.csv").open()))
    assert rows[0]["first_failing_criterion"] == "max_flag_helix_pitch_rel_err"
    manifest_path = root / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["conditions"][3]["config_overrides"]["flagella"]["placement_mode"] = (
        "seeded_center_layer"
    )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="not seeded_surface"):
        analyze(
            run_root=root,
            threshold_contract=ROOT
            / "conf/phase2_validation/2015_stage_a_thresholds.yaml",
            output_dir=tmp_path / "rejected",
        )
