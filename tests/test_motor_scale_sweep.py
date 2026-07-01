from __future__ import annotations

import csv
from pathlib import Path

import pytest

from sim_swim.analysis.sweeps import motor_scale


def test_phase24_local_hook_scale_sweep_reproduces_hook_gate(
    tmp_path: Path,
) -> None:
    """P2-4-003: standard local_hook_scale sweep keeps hook first-fail traceable."""
    out_dir = tmp_path / "phase24_hook_sweep"

    motor_scale.main(
        [
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

    summary_csv = out_dir / "summary.csv"
    assert summary_csv.is_file()
    with summary_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 4
    by_condition = {(float(row["torque_Nm"]), float(row["value"])): row for row in rows}

    for scale in (1.0, 8.0):
        row = by_condition[(1.2e-21, scale)]
        assert row["body_stiffness_scale"] == "50.0"
        assert row["dt_star"] == "0.001"
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


def test_body_stiffness_scale_must_be_positive(tmp_path: Path) -> None:
    """body stiffness override は無効値を argparse error として早期拒否する。"""
    with pytest.raises(SystemExit) as excinfo:
        motor_scale.main(
            [
                "--body-stiffness-scale",
                "0",
                "--output-dir",
                str(tmp_path / "invalid_body_stiffness"),
            ]
        )

    assert excinfo.value.code == 2


def test_stub_mode_full_flagella_summary_columns(tmp_path: Path) -> None:
    """P2-5-004: sweep helper can run the full_flagella single-flagellum baseline."""
    out_dir = tmp_path / "phase25_full_flagella_sweep"

    motor_scale.main(
        [
            "--target",
            "local_hook_scale",
            "--values",
            "8",
            "--torques",
            "1.2e-21,4e-21",
            "--duration",
            "0.02",
            "--stub-mode",
            "full_flagella",
            "--output-dir",
            str(out_dir),
        ],
    )

    summary_csv = out_dir / "summary.csv"
    with summary_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 2
    by_torque = {float(row["torque_Nm"]): row for row in rows}

    safe = by_torque[1.2e-21]
    assert safe["stub_mode"] == "full_flagella"
    assert safe["n_flagella"] == "1"
    assert safe["dt_star"] == "0.001"
    assert safe["shape_pass"] == "True"
    assert safe["first_fail_category"] == "none"

    fail = by_torque[4.0e-21]
    assert fail["stub_mode"] == "full_flagella"
    assert fail["n_flagella"] == "1"
    assert fail["shape_pass"] == "False"
    assert fail["first_fail_category"] == "flag"
    assert fail["first_fail_category_nonbody"] == "flag"
