from __future__ import annotations

from pathlib import Path

import pytest

from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.analysis.reference_torque_comparison import build_plan
from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig


def _raw(base_config: str) -> dict[str, object]:
    return {
        "base_config": base_config,
        "reference_torque_Nm": 1.2e-18,
        "policies": ["fixed-reference", "tracking-reference"],
        "time_bases": ["same-real-time", "same-dimensionless-time"],
        "motor_torque_scales": [0.5, 1.0],
        "durations": {"same_real_time_s": 0.01, "same_dimensionless_tau": 1.0},
    }


def test_plan_separates_fixed_and_tracking_scales_and_manifest() -> None:
    plan = build_plan(_raw("conf/sim_swim_2015.yaml"))
    conditions = {row["condition_id"]: row for row in plan["conditions"]}
    fixed = conditions["fixed_reference_same_real_time_motor_scale_0.5"]
    tracking = conditions["tracking_reference_same_real_time_motor_scale_0.5"]

    assert fixed["equivalence"]["physical_properties_fixed"] is True
    assert fixed["equivalence"]["reference_torque_Nm"] == pytest.approx(1.2e-18)
    assert fixed["equivalence"]["torque_for_forces_Nm"] == pytest.approx(1.2e-18)
    assert tracking["equivalence"]["reference_linked_properties"] is True
    assert tracking["equivalence"]["reference_torque_Nm"] == pytest.approx(0.6e-18)
    assert tracking["equivalence"]["torque_for_forces_Nm"] == pytest.approx(0.6e-18)
    assert tracking["equivalence"]["tau_tracks_reference_torque"] is True
    assert fixed["time"]["torque_for_forces_source"] == (
        "motor.torque_for_forces_override_Nm"
    )
    assert tracking["time"]["torque_for_forces_source"] == ("motor.reference_torque_Nm")
    assert tracking["time"]["material_coefficients"]["bend_k_Nm"] == pytest.approx(
        fixed["time"]["material_coefficients"]["bend_k_Nm"] * 0.5
    )
    assert tracking["time"]["paper_notation"]["tau"]["s"] == pytest.approx(
        tracking["time"]["tau_s"]
    )
    assert tracking["time"]["paper_notation"]["delta_t"]["s"] == pytest.approx(
        tracking["time"]["dt_internal_s"]
    )


def test_2010_project_default_tau_tracks_reference_torque() -> None:
    raw = _raw("conf/sim_swim_2010.yaml")
    raw["reference_torque_Nm"] = 2.5e-20
    plan = build_plan(raw)
    tracking = next(
        row
        for row in plan["conditions"]
        if row["comparison_policy"] == "tracking-reference"
    )
    assert tracking["equivalence"]["tau_tracks_reference_torque"] is True
    assert tracking["time"]["tau_s"] != pytest.approx(1.0)


def test_2010_project_keeps_legacy_tau_when_explicit() -> None:
    raw = _raw("conf/sim_swim_2010.yaml")
    raw["reference_torque_Nm"] = 2.5e-20
    raw["time_scale_policy"] = "legacy_fixed_tau_s_1"
    plan = build_plan(raw)
    tracking = next(
        row
        for row in plan["conditions"]
        if row["comparison_policy"] == "tracking-reference"
    )
    assert plan["time_scale_policy"] == "legacy_fixed_tau_s_1"
    assert tracking["equivalence"]["tau_tracks_reference_torque"] is False
    assert tracking["time"]["tau_s"] == pytest.approx(1.0)


def test_engine_material_coefficients_follow_force_scale_not_motor_drive() -> None:
    raw = load_yaml(Path("conf/sim_swim_2015.yaml"))
    raw["model_profile"]["implementation_status"] = "supported"
    base = SimulationConfig.from_dict(raw)
    fixed = base.with_overrides(
        {
            "motor": {
                "torque_Nm": 0.6e-18,
                "reference_torque_Nm": 1.2e-18,
                "torque_for_forces_override_Nm": 1.2e-18,
            }
        }
    )
    tracking = base.with_overrides(
        {
            "motor": {
                "torque_Nm": 0.6e-18,
                "reference_torque_Nm": 0.6e-18,
                "torque_for_forces_override_Nm": 0.0,
            }
        }
    )
    fixed_engine = DynamicsEngine(ModelBuilder(fixed).build(), fixed)
    tracking_engine = DynamicsEngine(ModelBuilder(tracking).build(), tracking)

    assert fixed.motor_torque_Nm == pytest.approx(tracking.motor_torque_Nm)
    assert tracking.tau_s == pytest.approx(fixed.tau_s * 2.0)
    assert tracking_engine.spring_h == pytest.approx(fixed_engine.spring_h * 0.5)
    assert tracking_engine.k_bend == pytest.approx(fixed_engine.k_bend * 0.5)
    assert tracking_engine.k_torsion == pytest.approx(fixed_engine.k_torsion * 0.5)
    assert tracking_engine.k_hook == pytest.approx(fixed_engine.k_hook * 0.5)
    assert tracking_engine.repulsion_A == pytest.approx(fixed_engine.repulsion_A * 0.5)
