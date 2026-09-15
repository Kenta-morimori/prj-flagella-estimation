from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from sim_swim.analysis.issue184_2015_nf1_6 import EXPECTED_IDS
from sim_swim.analysis.issue184_2015_nf1_6_provisional import analyze


ROOT = Path(__file__).parents[1]


def _source(root: Path, condition_id: str, *, fail: bool = False) -> Path:
    index = int(condition_id[-2:])
    placement = (
        "seeded_surface" if condition_id in {"nf04", "nf05"} else "seeded_center_layer"
    )
    condition_dir = root / condition_id
    condition_dir.mkdir(parents=True)
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
        f"0,0,0,0,0,0,0,0,{1 if fail else 0},0,0\n",
        encoding="utf-8",
    )
    topology = [
        {"flag_id": flag, "body_bead_index": 12 + flag, "layer": 2, "slot": flag}
        for flag in range(index)
    ]
    condition = {
        "condition_id": condition_id,
        "output_dir": str(condition_dir),
        "axis_values": {"n_flagella": index},
        "config_overrides": {
            "flagella": {"placement_mode": placement},
            "seed": {"attach_seed": 0, "phase_seed": 0},
            "time": {
                "scale_policy": "reference_torque",
                "duration": {"value": 10.0, "unit": "tau"},
                "integration": {"dt_star": 1.0e-5},
            },
            "motor": {"torque_Nm": 2.5e-20, "reference_torque_Nm": 2.5e-20},
        },
        "geometry": {"actual": {"attachment_topology": topology}},
    }
    (root / "run_manifest.json").write_text(
        json.dumps(
            {
                "kind": "generic_multi_run",
                "git": {"commit": "fixture"},
                "condition_order": [condition_id],
                "conditions": [condition],
            }
        ),
        encoding="utf-8",
    )
    row = {
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
        "flag_helix_pitch_rel_err_max": "1" if fail else "0",
        "motor_force_balance_residual_ratio": "0",
        "motor_torque_balance_residual_ratio": "0",
    }
    with (root / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    return root


def _sources(tmp_path: Path, *, fail: bool = False) -> dict[str, Path]:
    return {
        condition_id: _source(
            tmp_path / condition_id, condition_id, fail=fail and condition_id == "nf04"
        )
        for condition_id in EXPECTED_IDS
    }


def test_provisional_analysis_records_mixed_topology_and_strict_qc(
    tmp_path: Path,
) -> None:
    output = analyze(
        sources=_sources(tmp_path, fail=True),
        threshold_contract=ROOT / "conf/phase2_validation/2015_stage_a_thresholds.yaml",
        output_dir=tmp_path / "analysis",
    )
    decision = json.loads((output / "issue184_provisional_decision.json").read_text())
    rows = list(csv.DictReader((output / "issue184_provisional_summary.csv").open()))
    assert decision["status"] == "provisional"
    assert decision["strict_status"] == "fail"
    assert rows[0]["placement_mode"] == "seeded_center_layer"
    assert rows[3]["placement_mode"] == "seeded_surface"
    assert rows[3]["first_failing_criterion"] == "max_flag_helix_pitch_rel_err"


def test_provisional_analysis_rejects_missing_or_wrong_provenance(
    tmp_path: Path,
) -> None:
    sources = _sources(tmp_path)
    with pytest.raises(ValueError, match="exactly nf01"):
        analyze(
            sources={key: value for key, value in sources.items() if key != "nf06"},
            threshold_contract=ROOT
            / "conf/phase2_validation/2015_stage_a_thresholds.yaml",
            output_dir=tmp_path / "missing",
        )
    manifest_path = sources["nf04"] / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["conditions"][0]["config_overrides"]["flagella"]["placement_mode"] = (
        "seeded_center_layer"
    )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="placement mode mismatch"):
        analyze(
            sources=sources,
            threshold_contract=ROOT
            / "conf/phase2_validation/2015_stage_a_thresholds.yaml",
            output_dir=tmp_path / "wrong",
        )
