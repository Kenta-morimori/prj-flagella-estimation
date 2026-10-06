from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    geometry_preflight,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.analysis.parallel_job import (
    build_plan,
    load_parallel_job,
    resolve_execution,
)
from sim_swim.dynamics.forces import (
    _attach_body_support,
    compute_root_torque_segment_couples_forces,
)
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("model_name,expected", [("hex", 13), ("project", 3)])
def test_local_reaction_contract_and_force_balance(
    model_name: str, expected: int, tmp_path: Path
) -> None:
    config_name = "2010_hex_project" if model_name == "hex" else "2010_project"
    campaign = normalize_campaign_config(
        load_yaml(
            ROOT
            / f"conf/phase2_multi_run/{config_name}_local_reaction_1tau_issue255.yaml"
        )
    )
    conditions = build_campaign_conditions(campaign)
    preflight = geometry_preflight(campaign, conditions)
    job = load_parallel_job(
        ROOT
        / f"conf/phase2_parallel/issue255_motor_reaction/{model_name}_local_1tau_job.yaml"
    )
    execution = resolve_execution(job, None)
    plan = build_plan(job, execution, tmp_path / "plan")
    assert len(conditions) == len(preflight) == len(plan["configs"]) == expected
    assert execution.max_workers == min(8, expected)
    assert len({record["output_dir"] for record in plan["configs"]}) == expected

    base = load_yaml(ROOT / campaign["base_config"])
    old_campaign = normalize_campaign_config(
        load_yaml(
            ROOT
            / f"conf/phase2_multi_run/{config_name}_motor_reaction_1tau_issue255.yaml"
        )
    )
    old_full = {
        record["condition_id"].removesuffix("__rxfull"): record
        for record in build_campaign_conditions(old_campaign)
        if record["condition_id"].endswith("__rxfull")
    }
    for condition in conditions:
        cid = condition["condition_id"]
        assert cid in old_full
        local = SimulationConfig.from_dict(base).with_overrides(
            condition["config_overrides"]
        )
        global_full = SimulationConfig.from_dict(base).with_overrides(
            old_full[cid]["config_overrides"]
        )
        assert local.motor.body_reaction_full_vector
        assert local.motor.body_reaction_support == "attach_one_ring"
        assert global_full.motor.body_reaction_support == "all_body"
        assert local.motor.force_distribution == global_full.motor.force_distribution
        assert local.motor.torque_Nm == global_full.motor.torque_Nm
        assert local.time == global_full.time
        model = ModelBuilder(local).build()
        np.testing.assert_array_equal(
            model.positions_m, ModelBuilder(global_full).build().positions_m
        )
        assert preflight[cid]["initial_hook_body_axis_error_deg"] < 1e-6
        assert preflight[cid]["initial_helix_rear_angle_max_deg"] < 1e-6
        weights = [np.ones(len(ids) - 1) for ids in model.flagella_indices]
        common = dict(
            positions_m=model.positions_m,
            flagella_indices=model.flagella_indices,
            body_indices=model.body_indices,
            torque_per_flag=local.motor_torque_Nm * model.torque_signs,
            segment_weights=weights,
            full_vector_body_reaction=True,
        )
        global_forces, _ = compute_root_torque_segment_couples_forces(**common)
        local_forces, diag = compute_root_torque_segment_couples_forces(
            **common,
            body_reaction_support="attach_one_ring",
            flagella_attach_body_indices=model.flagella_attach_body_indices,
            body_ring_edges=model.body_ring_edges,
            body_vertical_edges=model.body_vertical_edges,
        )
        flag = np.concatenate(model.flagella_indices)
        np.testing.assert_array_equal(local_forces[flag], global_forces[flag])
        assert diag.reaction_support_bead_counts == (5,) * len(model.flagella_indices)
        assert not diag.reaction_fallback_used
        support = set()
        for attach in model.flagella_attach_body_indices:
            support.update(
                _attach_body_support(
                    attach_index=int(attach),
                    body_ring_edges=model.body_ring_edges,
                    body_vertical_edges=model.body_vertical_edges,
                )
            )
        assert np.all(local_forces[list(set(model.body_indices) - support)] == 0)
        scale = abs(local.motor_torque_Nm) * len(model.flagella_indices)
        assert np.linalg.norm(local_forces.sum(axis=0)) < 1e-8 * scale / local.b_m
        total_torque = np.sum(np.cross(model.positions_m, local_forces), axis=0)
        assert np.linalg.norm(total_torque) < 1e-8 * scale


def test_local_reaction_fails_without_body_support() -> None:
    base = load_yaml(ROOT / "conf/sim_swim_2010.yaml")
    cfg = SimulationConfig.from_dict(base)
    model = ModelBuilder(cfg).build()
    with pytest.raises(RuntimeError, match="Local body reaction solver failed"):
        compute_root_torque_segment_couples_forces(
            positions_m=model.positions_m,
            flagella_indices=model.flagella_indices,
            body_indices=model.body_indices,
            torque_per_flag=cfg.motor_torque_Nm * model.torque_signs,
            segment_weights=[np.ones(len(ids) - 1) for ids in model.flagella_indices],
            full_vector_body_reaction=True,
            body_reaction_support="attach_one_ring",
            flagella_attach_body_indices=model.flagella_attach_body_indices,
            body_ring_edges=np.zeros((0, 2), dtype=int),
            body_vertical_edges=np.zeros((0, 2), dtype=int),
        )


@pytest.mark.parametrize("model_name,expected", [("hex", 13), ("project", 3)])
def test_local_reaction_2s_preserves_1tau_conditions(
    model_name: str, expected: int, tmp_path: Path
) -> None:
    config_name = "2010_hex_project" if model_name == "hex" else "2010_project"
    campaign_dir = ROOT / "conf/phase2_multi_run"
    screen = normalize_campaign_config(
        load_yaml(campaign_dir / f"{config_name}_local_reaction_1tau_issue255.yaml")
    )
    main = normalize_campaign_config(
        load_yaml(campaign_dir / f"{config_name}_local_reaction_2s_issue255.yaml")
    )
    assert main["development_evaluation"]["stage"] == "long_duration"
    assert (
        main["development_evaluation"]["qc"] == screen["development_evaluation"]["qc"]
    )
    assert main["sweep"] == screen["sweep"]
    assert main["base_config"] == screen["base_config"]
    main_overrides = dict(main["base_overrides"])
    main_overrides["time"] = dict(main_overrides["time"])
    main_overrides["time"]["duration"] = screen["base_overrides"]["time"]["duration"]
    assert main_overrides == screen["base_overrides"]
    conditions = build_campaign_conditions(main)
    preflight = geometry_preflight(main, conditions)
    job = load_parallel_job(
        ROOT
        / f"conf/phase2_parallel/issue255_motor_reaction/{model_name}_local_2s_job.yaml"
    )
    execution = resolve_execution(job, None)
    plan = build_plan(job, execution, tmp_path / "plan")
    assert len(conditions) == len(preflight) == len(plan["configs"]) == expected
    assert execution.max_workers == min(8, expected)
    assert len({record["output_dir"] for record in plan["configs"]}) == expected
    base = load_yaml(ROOT / main["base_config"])
    for condition in conditions:
        cid = condition["condition_id"]
        cfg = SimulationConfig.from_dict(base).with_overrides(
            condition["config_overrides"]
        )
        assert cfg.total_steps == 500_000
        assert cfg.motor.body_reaction_support == "attach_one_ring"
        assert preflight[cid]["initial_helix_rear_angle_max_deg"] <= 1e-6
        assert preflight[cid]["initial_hook_angle_max_abs_error_deg"] <= 1e-6
        assert preflight[cid]["initial_hook_body_axis_error_deg"] <= 1e-6
        assert preflight[cid]["initial_hook_force_norm_N"] <= 1e-18
        assert preflight[cid]["initial_min_outward_projection_m2"] > 0
        assert (
            preflight[cid]["initial_min_nonattached_bead_distance_m"]
            >= 2 * cfg.scale.bead_radius_a_over_b * cfg.b_m
        )


@pytest.mark.parametrize(
    "overrides",
    [
        {"body_reaction_support": "unknown"},
        {"body_reaction_support": "attach_one_ring"},
        {
            "body_reaction_support": "attach_one_ring",
            "body_reaction_full_vector": True,
            "force_distribution": "triplet",
        },
    ],
)
def test_local_reaction_rejects_invalid_settings(overrides: dict[str, object]) -> None:
    base = load_yaml(ROOT / "conf/sim_swim_2010.yaml")
    base["motor"].update(overrides)
    with pytest.raises(ValueError, match="motor.body_reaction_support"):
        SimulationConfig.from_dict(base)
