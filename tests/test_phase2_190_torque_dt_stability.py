from pathlib import Path

import pytest

from sim_swim.analysis.torque_dt_stability import build_plan
from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.sim.params import SimulationConfig


CONFIG = Path("conf/phase2_sweeps/2010_project_torque_linked_dt_stability.yaml")


def test_2010_project_keeps_legacy_time_scale_by_default() -> None:
    cfg = SimulationConfig.from_dict(load_yaml(Path("conf/sim_swim_2010.yaml")))
    assert cfg.tau_s == pytest.approx(1.0)
    assert cfg.time_scale_policy == "legacy_fixed_tau_s_1"


def test_issue_190_plan_links_all_torques_and_derives_tau() -> None:
    plan = build_plan(CONFIG)
    assert len(plan["conditions"]) == 12
    high = next(
        row for row in plan["conditions"] if row["condition_id"] == "T1p2e-18_dt1e-4"
    )
    assert high["time"]["tau_s"] == pytest.approx(1.0e-3 * (1.0e-6**3) / 1.2e-18)
    assert high["time"]["total_steps"] == 10000
    assert high["time"]["motor_torque_Nm"] == pytest.approx(1.2e-18)
    assert high["time"]["reference_torque_Nm"] == pytest.approx(1.2e-18)
    assert high["time"]["torque_for_forces_Nm"] == pytest.approx(1.2e-18)
    assert high["comparison_archive_interval_s"] == pytest.approx(
        high["time"]["duration_s"] / 200
    )
    assert "output.archive_interval_s=" in high["execution_command"]


def test_issue_190_condition_ids_preserve_non_power_of_ten_torques() -> None:
    ids = {row["condition_id"] for row in build_plan(CONFIG)["conditions"]}

    assert "T2p5e-20_dt1e-4" in ids
    assert "T1p2e-18_dt1e-5" in ids
