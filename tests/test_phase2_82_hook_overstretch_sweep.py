from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace


def _load_sweep_script():
    script_path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "01_simulate_swimming"
        / "run_phase2_82_hook_overstretch_sweep.py"
    )
    spec = importlib.util.spec_from_file_location("phase2_82_sweep", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_phase2_82_body_first_grid_conditions() -> None:
    script = _load_sweep_script()
    args = SimpleNamespace(
        mode="body-first-grid",
        attach_first_spring_scales=[1.0, 2.0],
        body_axis_angle_scales=[1.0, 3.0],
        fixed_first_second_spring_scale=1.0,
    )

    conditions = script.build_conditions(args)

    assert [condition.condition_id for condition in conditions] == [
        "af1_axis1_fs1",
        "af1_axis3_fs1",
        "af2_axis1_fs1",
        "af2_axis3_fs1",
    ]
    assert all(condition.mode == "body-first-grid" for condition in conditions)
    assert conditions[-1].scales == {
        "local_attach_first_spring_scale": 2.0,
        "local_attach_first_body_axis_angle_scale": 3.0,
        "local_first_second_spring_scale": 1.0,
    }


def test_phase2_82_first_second_grid_conditions() -> None:
    script = _load_sweep_script()
    args = SimpleNamespace(
        mode="first-second-grid",
        fixed_attach_first_spring_scale=2.0,
        fixed_body_axis_angle_scale=1.5,
        first_second_spring_scales=[1.0, 2.0, 3.0],
    )

    conditions = script.build_conditions(args)

    assert [condition.condition_id for condition in conditions] == [
        "af2_axis1p5_fs1",
        "af2_axis1p5_fs2",
        "af2_axis1p5_fs3",
    ]
    assert all(condition.mode == "first-second-grid" for condition in conditions)
    assert conditions[1].scales == {
        "local_attach_first_spring_scale": 2.0,
        "local_attach_first_body_axis_angle_scale": 1.5,
        "local_first_second_spring_scale": 2.0,
    }


def test_phase2_82_summary_row_records_fail_and_max_hook_events(
    tmp_path: Path,
) -> None:
    script = _load_sweep_script()
    cfg = SimpleNamespace(
        time=SimpleNamespace(duration_s=0.5),
        dt_star=1.0e-4,
        motor_torque_Nm=2.5e-20,
        flagella=SimpleNamespace(n_flagella=3),
        motor=SimpleNamespace(
            local_attach_first_spring_scale=3.0,
            local_attach_first_body_axis_angle_scale=1.25,
            local_first_second_spring_scale=1.25,
        ),
    )
    condition = script.Condition(
        condition_id="af3_axis1p25_fs1p25",
        mode="first-second-grid",
        description="test",
        scales={},
    )
    rows = [
        {
            "t_s": "0.0",
            "shape_pass_nonbody": "True",
            "hook_len_rel_err_max": "0.1",
            "hook_len_rel_err_max_flag_id": "0",
        },
        {
            "t_s": "0.2",
            "shape_pass_nonbody": "False",
            "hook_len_rel_err_max": "1.01",
            "hook_len_rel_err_max_flag_id": "1",
            "local_attach_first_rel_err": "0.9",
            "local_first_second_rel_err": "0.2",
        },
        {
            "t_s": "0.5",
            "shape_pass_nonbody": "False",
            "hook_len_rel_err_max": "1.4",
            "hook_len_rel_err_max_flag_id": "2",
            "hook_len_rel_err_max_attach_body_bead_index": "8",
            "hook_len_rel_err_max_flag_first_bead_index": "25",
            "hook_len_rel_err_max_len_over_b": "0.6",
        },
    ]

    row = script._summary_row(
        cfg,
        condition,
        tmp_path,
        rows[-1],
        rows,
        helix_summary={},
    )

    assert row["first_fail_t_s"] == "0.2"
    assert row["first_fail_hook_len_rel_err_max"] == "1.01"
    assert row["first_fail_hook_len_rel_err_max_flag_id"] == "1"
    assert row["max_hook_len_rel_err_t_s"] == "0.5"
    assert row["max_hook_len_rel_err"] == "1.4"
    assert row["max_hook_len_rel_err_flag_id"] == "2"
    assert row["max_hook_len_rel_err_attach_body_bead_index"] == "8"
    assert row["max_hook_len_rel_err_flag_first_bead_index"] == "25"
