from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.dynamics.forces import compute_root_torque_segment_couples_forces
from sim_swim.dynamics.frame_potential import attachment_frame_energy_forces
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]


def geometry():
    campaign = normalize_campaign_config(
        load_yaml(
            ROOT
            / "conf/phase2_multi_run/2010_hex_project_local_reaction_1tau_issue255.yaml"
        )
    )
    condition = next(
        c
        for c in build_campaign_conditions(campaign)
        if c["condition_id"] == "nf03__slots024"
    )
    cfg = SimulationConfig.from_dict(
        load_yaml(ROOT / campaign["base_config"])
    ).with_overrides(condition["config_overrides"])
    model = ModelBuilder(cfg).build()
    return cfg, model, DynamicsEngine(model, cfg)


def frame_args(model, engine):
    return (
        model,
        engine.hook_attach_layer_indices,
        (
            engine.hook_attach_first_frame_local_m,
            engine.hook_first_second_frame_local_m,
        ),
        (
            engine.hook_attach_first_rest_lengths_m,
            engine.hook_first_second_rest_lengths_m,
        ),
        (engine.k_hook * 0.25, engine.k_hook * 0.1),
    )


def test_frame_gradient_conservation_and_covariance():
    cfg, model, engine = geometry()
    args = frame_args(model, engine)
    pos = model.positions_m.copy()
    pos[model.hook_triplets[:, 1]] += np.array([0.08, -0.04, 0.03]) * cfg.b_m
    pos[model.body_layer_indices[-1][0]] += np.array([0.03, 0.02, -0.01]) * cfg.b_m
    energy, force = attachment_frame_energy_forces(pos, *args)
    epsilon = cfg.b_m * 1e-6
    numerical = np.zeros_like(pos)
    relevant = np.unique(
        np.concatenate([*model.body_layer_indices, model.hook_triplets.ravel()])
    )
    for i in relevant:
        for axis in range(3):
            high, low = pos.copy(), pos.copy()
            high[i, axis] += epsilon
            low[i, axis] -= epsilon
            numerical[i, axis] = -(
                attachment_frame_energy_forces(high, *args)[0]
                - attachment_frame_energy_forces(low, *args)[0]
            ) / (2 * epsilon)
    np.testing.assert_allclose(
        force, numerical, rtol=1e-5, atol=np.max(abs(force)) * 1e-7
    )
    assert np.linalg.norm(force.sum(0)) < np.linalg.norm(force) * 1e-12
    assert (
        np.linalg.norm(np.cross(pos, force).sum(0))
        < np.linalg.norm(force) * cfg.b_m * 1e-11
    )
    rotation, _ = np.linalg.qr(np.random.default_rng(1).normal(size=(3, 3)))
    transformed = pos @ rotation.T + np.array([2, -3, 4]) * cfg.b_m
    rotated_energy, rotated_force = attachment_frame_energy_forces(transformed, *args)
    assert rotated_energy == pytest.approx(energy, rel=1e-12)
    np.testing.assert_allclose(
        rotated_force, force @ rotation.T, atol=np.max(abs(force)) * 1e-12
    )
    initial_energy, initial_force = attachment_frame_energy_forces(
        model.positions_m, *args
    )
    assert initial_energy < 1e-45
    assert np.max(abs(initial_force)) < 1e-24


def test_frame_degenerate_fails():
    _, model, engine = geometry()
    with pytest.raises(RuntimeError, match="degenerate"):
        attachment_frame_energy_forces(
            np.zeros_like(model.positions_m), *frame_args(model, engine)
        )


@pytest.mark.parametrize("deform", [False, True])
def test_drive_correction_full_vector_balance_and_minimum_norm(deform):
    cfg, model, _ = geometry()
    pos = model.positions_m.copy()
    if deform:
        ids = model.flagella_indices[0]
        # One nearly axial segment and a distorted helix.
        from sim_swim.dynamics.forces import _principal_axis_or_none

        axis = _principal_axis_or_none(pos[ids])
        perpendicular = np.cross(axis, np.array([0.0, 1.0, 0.0]))
        pos[ids[1]] = pos[ids[0]] + (0.1 * axis + 0.001 * perpendicular) * cfg.b_m
    weights = [np.geomspace(1, 0.001, len(ids) - 1) for ids in model.flagella_indices]
    common = dict(
        positions_m=pos,
        flagella_indices=model.flagella_indices,
        body_indices=model.body_indices,
        torque_per_flag=np.full(3, cfg.motor_torque_Nm),
        segment_weights=weights,
        full_vector_body_reaction=True,
        body_reaction_support="attach_one_ring",
        flagella_attach_body_indices=model.flagella_attach_body_indices,
        body_ring_edges=model.body_ring_edges,
        body_vertical_edges=model.body_vertical_edges,
    )
    legacy, _ = compute_root_torque_segment_couples_forces(**common)
    corrected, diag = compute_root_torque_segment_couples_forces(
        **common, segment_torque_correction="minimum_norm"
    )
    np.testing.assert_array_equal(
        legacy,
        compute_root_torque_segment_couples_forces(
            **common, segment_torque_correction="none"
        )[0],
    )
    metrics = dict(diag.drive_metrics)
    assert metrics["motor_drive_transverse_ratio_max"] < 1e-9
    assert metrics["motor_drive_axial_error_ratio_max"] < 1e-9
    assert np.linalg.norm(corrected.sum(0)) < np.linalg.norm(corrected) * 1e-10
    assert (
        np.linalg.norm(np.cross(pos - pos.mean(0), corrected).sum(0))
        < cfg.motor_torque_Nm * 1e-8
    )
    assert diag.reaction_support_bead_counts == (5, 5, 5)
    # Minimum-norm correction is orthogonal to every force/torque-free perturbation.
    from sim_swim.dynamics.forces import _cross_matrix

    ids = model.flagella_indices[0]
    arms = (pos[ids] - pos[ids[0]]) / cfg.b_m
    operator = np.vstack(
        (
            np.tile(np.eye(3), (1, len(ids))),
            np.hstack([_cross_matrix(arm) for arm in arms]),
        )
    )
    perturbation = np.random.default_rng(2).normal(size=3 * len(ids))
    null = perturbation - np.linalg.pinv(operator) @ (operator @ perturbation)
    delta = (corrected[ids] - legacy[ids]).ravel()
    assert abs(delta @ null) < np.linalg.norm(delta) * np.linalg.norm(null) * 1e-9


def test_engine_correction_telemetry():
    cfg, _, _ = geometry()
    cfg = replace(
        cfg,
        motor=replace(
            cfg.motor,
            attach_frame_reaction="energy_gradient",
            segment_torque_correction="minimum_norm",
        ),
    )
    model = ModelBuilder(cfg).build()
    engine = DynamicsEngine(model, cfg)
    diag = engine.step(cfg.dt_star)
    metrics = diag.force_balance_diagnostics
    assert metrics["force_evaluation_t_s"] == pytest.approx(0, abs=1e-12)
    assert metrics["motor_drive_transverse_ratio_max"] < 1e-9
    assert metrics["motor_reaction_support_count_min"] == 5
    assert metrics["motor_reaction_solver_success"] == 1
    assert metrics["motor_reaction_fallback_used"] == 0


@pytest.mark.parametrize(
    "override",
    [
        {"motor.attach_frame_reaction": "invalid"},
        {"motor.segment_torque_correction": "invalid"},
        {
            "motor.attach_frame_reaction": "energy_gradient",
            "motor.local_attach_frame_tangent_mode": "basal_bearing",
        },
    ],
)
def test_config_rejects_invalid_candidates(override):
    cfg, _, _ = geometry()
    with pytest.raises(ValueError):
        cfg.with_overrides(
            {"motor": {k.removeprefix("motor."): v for k, v in override.items()}}
        )


def test_corrected_drive_degenerate_solver_fails():
    cfg, model, _ = geometry()
    pos = model.positions_m.copy()
    ids = model.flagella_indices[0]
    pos[ids] = pos[ids[0]]
    with pytest.raises(RuntimeError, match="degenerate"):
        compute_root_torque_segment_couples_forces(
            pos,
            [ids],
            model.body_indices,
            np.array([cfg.motor_torque_Nm]),
            [np.ones(len(ids) - 1)],
            full_vector_body_reaction=True,
            segment_torque_correction="minimum_norm",
        )


@pytest.mark.parametrize("profile", ["2010_hex_project", "2010_project"])
def test_all_initial_shapes_conservative_frame(profile):
    campaign = normalize_campaign_config(
        load_yaml(
            ROOT / f"conf/phase2_multi_run/{profile}_local_reaction_1tau_issue255.yaml"
        )
    )
    base = load_yaml(ROOT / campaign["base_config"])
    for condition in build_campaign_conditions(campaign):
        cfg = SimulationConfig.from_dict(base).with_overrides(
            condition["config_overrides"]
        )
        model = ModelBuilder(cfg).build()
        engine = DynamicsEngine(model, cfg)
        _, forces = attachment_frame_energy_forces(
            model.positions_m, *frame_args(model, engine)
        )
        assert np.max(abs(forces)) < 1e-24


def test_frame_drive_campaigns_preflight_and_stage_contract(tmp_path):
    from sim_swim.analysis.multi_run_campaign import geometry_preflight
    from sim_swim.analysis.parallel_job import (
        build_plan,
        load_parallel_job,
        resolve_execution,
    )

    campaigns = []
    expected = {"nf03__slots024", "nf04__slots0134", "nf05__slots01234"}
    for stage in ("1tau", "1s"):
        path = (
            ROOT
            / f"conf/phase2_multi_run/2010_hex_project_frame_drive_{stage}_issue255.yaml"
        )
        campaign = normalize_campaign_config(load_yaml(path))
        conditions = build_campaign_conditions(campaign)
        assert {c["condition_id"] for c in conditions} == expected
        assert len(geometry_preflight(campaign, conditions)) == 3
        job = load_parallel_job(
            ROOT
            / f"conf/phase2_parallel/issue255_motor_reaction/hex_frame_drive_{stage}_job.yaml"
        )
        execution = resolve_execution(job, None)
        plan = build_plan(job, execution, tmp_path / stage)
        assert execution.max_workers == 3
        assert len({r["output_dir"] for r in plan["configs"]}) == 3
        base = SimulationConfig.from_dict(load_yaml(ROOT / campaign["base_config"]))
        configs = {
            c["condition_id"]: base.with_overrides(c["config_overrides"])
            for c in conditions
        }
        for cfg in configs.values():
            assert cfg.motor.attach_frame_reaction == "energy_gradient"
            assert cfg.motor.segment_torque_correction == "minimum_norm"
            assert cfg.motor.body_reaction_support == "attach_one_ring"
            assert not cfg.brownian.enabled
            assert cfg.total_steps == (10000 if stage == "1tau" else 250000)
        campaigns.append(configs)
    for cid in expected:
        short, long = campaigns[0][cid], campaigns[1][cid]
        assert short.motor == long.motor
        assert short.flagella == long.flagella
        assert short.potentials == long.potentials
        assert short.seed == long.seed
        assert short.dt_star == long.dt_star
        np.testing.assert_array_equal(
            ModelBuilder(short).build().positions_m,
            ModelBuilder(long).build().positions_m,
        )


def test_diagnostic_writer_and_extremum_time(tmp_path):
    import csv

    from sim_swim.analysis.online_run_summary import OnlineRunSummary
    from sim_swim.sim.debug_summary import StepSummaryRecorder

    cfg, model, engine = geometry()
    recorder = StepSummaryRecorder(model, cfg, tmp_path)
    diag = engine.step(cfg.dt_star)
    recorder.record(0, engine.t_star, diag)
    path = recorder.write_csv()
    row = next(csv.DictReader(path.open()))
    assert float(row["force_evaluation_t_s"]) == pytest.approx(0.0, abs=1e-15)
    assert float(row["motor_reaction_support_count_min"]) == 5
    summary = OnlineRunSummary(expected_steps=1)
    summary.record(recorder.last_row)
    assert summary.maximum_times["component_frame_torque_norm_Nm"] == pytest.approx(
        0.0, abs=1e-15
    )
