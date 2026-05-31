from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest


def _load_sweep_script():
    script_path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "01_simulate_swimming"
        / "run_motor_scale_sweep.py"
    )
    spec = importlib.util.spec_from_file_location(
        "motor_scale_sweep_script", script_path
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_phase24_local_hook_scale_sweep_reproduces_hook_gate(
    tmp_path: Path, monkeypatch
) -> None:
    """P2-4-003: standard local_hook_scale sweep keeps hook first-fail traceable."""
    sweep_script = _load_sweep_script()
    out_dir = tmp_path / "phase24_hook_sweep"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_motor_scale_sweep",
            "--target",
            "local_hook_scale",
            "--values",
            "1,8",
            "--torques",
            "1.2e-21,4e-21",
            "--duration",
            "0.02",
            "--body-stiffness-scale",
            "50",
            "--output-dir",
            str(out_dir),
        ],
    )

    sweep_script.main()

    summary_csv = out_dir / "local_hook_scale_sweep_summary.csv"
    assert summary_csv.is_file()
    with summary_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 4
    by_condition = {(float(row["torque_Nm"]), float(row["value"])): row for row in rows}

    for scale in (1.0, 8.0):
        row = by_condition[(1.2e-21, scale)]
        assert row["body_stiffness_scale"] == "50.0"
        assert row["shape_pass_nonbody"] == "True"
        assert row["body_shape_pass"] == "True"
        assert row["shape_pass"] == "True"
        assert row["first_fail_category"] == "none"
        assert row["first_fail_category_nonbody"] == "none"
        assert row["body_fail_category"] == "none"
        assert float(row["hook_len_rel_err_max"]) < 1.0
        assert float(row["local_attach_first_rel_err"]) < 1.0

    for scale in (1.0, 8.0):
        row = by_condition[(4.0e-21, scale)]
        assert row["body_stiffness_scale"] == "50.0"
        assert row["finite_pass"] == "True"
        assert row["shape_pass_nonbody"] == "False"
        assert row["shape_pass"] == "False"
        assert row["first_fail_category"] == "hook"
        assert row["first_fail_category_nonbody"] == "hook"
        assert row["body_fail_category"] in {
            "none",
            "body_nonfinite",
            "body_spring",
            "body_bend",
            "body_centerline",
            "body_area",
        }
        assert float(row["hook_len_rel_err_max"]) > 1.0
        assert float(row["local_attach_first_rel_err"]) > 1.0
        assert row["body_shape_pass"] in {"True", "False"}
        assert row["body_spring_max_stretch_ratio"] != ""


def test_body_stiffness_scale_must_be_positive(tmp_path: Path, monkeypatch) -> None:
    """body stiffness override は無効値を argparse error として早期拒否する。"""
    sweep_script = _load_sweep_script()

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_motor_scale_sweep",
            "--body-stiffness-scale",
            "0",
            "--output-dir",
            str(tmp_path / "invalid_body_stiffness"),
        ],
    )

    with pytest.raises(SystemExit) as excinfo:
        sweep_script.main()

    assert excinfo.value.code == 2
