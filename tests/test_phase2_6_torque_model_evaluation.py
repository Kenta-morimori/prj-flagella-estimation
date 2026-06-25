from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import pytest

from sim_swim.sim.core import SimulationState


def _load_eval_script():
    script_path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "01_simulate_swimming"
        / "run_phase2_6_torque_model_evaluation.py"
    )
    spec = importlib.util.spec_from_file_location(
        "phase2_6_torque_model_evaluation", script_path
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _state(t: float, position_um: tuple[float, float, float]) -> SimulationState:
    return SimulationState(
        t=t,
        position_um=position_um,
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(0.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=[],
    )


def test_body_motion_summary_reports_displacement_and_speed() -> None:
    script = _load_eval_script()
    summary = script._summarize_body_motion(
        [_state(0.0, (0.0, 0.0, 0.0)), _state(0.5, (1.5, 0.0, 0.0))]
    )

    assert summary["body_displacement_um"] == pytest.approx(1.5)
    assert summary["body_path_length_um"] == pytest.approx(1.5)
    assert summary["body_mean_speed_um_s"] == pytest.approx(3.0)
    assert summary["body_path_mean_speed_um_s"] == pytest.approx(3.0)
    assert summary["body_axis_angle_change_deg"] == pytest.approx(0.0)


def test_auto_scale_torque_selection_includes_upper_and_first_fail() -> None:
    script = _load_eval_script()
    torques = [0.5e-20, 2.0e-20, 4.0e-20, 1.0e-19]
    rows = [
        {
            "torque_Nm": 0.5e-20,
            "helix_retention_pass": False,
            "first_fail_category": "motor_no_rotation",
        },
        {
            "torque_Nm": 2.0e-20,
            "helix_retention_pass": True,
            "first_fail_category": "none",
        },
        {
            "torque_Nm": 4.0e-20,
            "helix_retention_pass": False,
            "first_fail_category": "flag",
        },
        {
            "torque_Nm": 1.0e-19,
            "helix_retention_pass": False,
            "first_fail_category": "flag",
        },
    ]

    selected = script._select_auto_scale_torques(rows, torques=torques)

    assert selected == [2.0e-20, 4.0e-20, 1.0e-19]


def test_torque_upper_bound_is_rejected(tmp_path: Path, monkeypatch) -> None:
    script = _load_eval_script()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_phase2_6_torque_model_evaluation",
            "--torques",
            "2.0e-19",
            "--output-dir",
            str(tmp_path / "upper_bound"),
        ],
    )

    with pytest.raises(SystemExit, match="1.0e-19"):
        script.main()


def test_phase2_6_eval_smoke_writes_summary(tmp_path: Path, monkeypatch) -> None:
    script = _load_eval_script()
    out_dir = tmp_path / "phase2_6_eval_smoke"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_phase2_6_torque_model_evaluation",
            "--torques",
            "1.0e-20",
            "--scale-torques",
            "none",
            "--duration",
            "0.001",
            "--dt-star",
            "1.0e-3",
            "--max-conditions",
            "1",
            "--output-dir",
            str(out_dir),
        ],
    )

    script.main()

    summary_csv = out_dir / "phase2_6_torque_model_evaluation_summary.csv"
    assert summary_csv.is_file()
    with summary_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 1
    row = rows[0]
    assert row["force_distribution"] == "root_torque_segment_couples"
    assert row["n_flagella"] == "1"
    assert row["stub_mode"] == "full_flagella"
    assert row["local_hook_scale"] == "1.0"
    assert row["local_spring_scale"] == "1.0"
    assert row["local_bend_scale"] == "1.0"
    assert row["local_torsion_scale"] == "1.0"

    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "completed"
    assert manifest["condition_count"] == 1
