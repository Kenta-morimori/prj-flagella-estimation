from __future__ import annotations

from pathlib import Path

import pytest

from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.analysis.torque_dt_stability import build_plan
from sim_swim.analysis.torque_dt_stability_campaign import (
    _assert_condition,
    _validate_campaign_contract,
)
from sim_swim.analysis.sweeps.generic_multi_run import run_campaign
from sim_swim.sim.params import SimulationConfig


CONFIG = Path("conf/phase2_multi_run/2010_project_tau_policy_torque_dt_0p05s.yaml")


def test_issue199_plan_expands_fixed_real_time_policy_grid() -> None:
    raw = load_yaml(CONFIG)
    _validate_campaign_contract(raw)
    plan = build_plan(CONFIG)

    assert plan["duration"] == {"value": 0.05, "unit": "s"}
    assert len(plan["conditions"]) == 12
    by_id = {row["condition_id"]: row for row in plan["conditions"]}
    fixed = by_id["tau_fixed_control_T1e-21_dt1e-3"]
    linked = by_id["torque_linked_tau_T1e-21_dt1e-3"]
    assert fixed["time"]["tau_s"] == pytest.approx(1.0)
    assert linked["time"]["tau_s"] == pytest.approx(1.0)
    assert fixed["time"]["total_steps"] == linked["time"]["total_steps"] == 50
    assert by_id["torque_linked_tau_T1e-19_dt1e-4"]["time"]["total_steps"] == 50001


def test_issue199_conditions_enforce_the_declared_time_policy() -> None:
    raw = load_yaml(CONFIG)
    base = load_yaml(Path(raw["base_config"]))
    plan = build_plan(CONFIG)
    for condition in plan["conditions"]:
        cfg = SimulationConfig.from_dict(base).with_overrides(
            condition["config_overrides"]
        )
        _assert_condition(cfg, condition)
        assert cfg.torque_for_forces_Nm == pytest.approx(abs(cfg.motor_torque_Nm))


def test_issue199_dry_run_has_no_output_side_effects(
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert run_campaign(["--campaign-config", str(CONFIG), "--dry-run"]) == Path()
    lines = capsys.readouterr().out.splitlines()
    assert len(lines) == 12
    assert lines[0].startswith("tau_fixed_control_T")


def test_issue199_contract_rejects_extra_torque() -> None:
    raw = load_yaml(CONFIG)
    raw["torques_Nm"].append(1.2e-18)
    with pytest.raises(ValueError, match="fixed 3x2 torque-dt grid"):
        _validate_campaign_contract(raw)
