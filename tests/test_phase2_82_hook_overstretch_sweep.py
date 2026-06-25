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
