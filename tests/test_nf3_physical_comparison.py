"""Independent physical hypotheses: preserve controls and fail closed."""

from copy import deepcopy
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    geometry_preflight,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.analysis.parallel_job import load_parallel_job, resolve_execution
from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.dynamics.forces import (
    compute_root_torque_segment_couples_forces,
    compute_segment_repulsion_forces,
    continuous_segment_repulsion_potential,
)
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]
SHORT = ROOT / "conf/phase2_multi_run/2010_hex_nf3_physics_1tau_issue263.yaml"
LONG = ROOT / "conf/phase2_multi_run/2010_hex_nf3_physics_1p5s_issue263.yaml"


def configs(path=SHORT):
    campaign = normalize_campaign_config(load_yaml(path))
    rows = build_campaign_conditions(campaign)
    base = SimulationConfig.from_dict(load_yaml(ROOT / campaign["base_config"]))
    return campaign, rows, [base.with_overrides(r["config_overrides"]) for r in rows]


def test_five_independent_controls_and_duration_only_extension():
    short, rows, cfgs = configs()
    long, longrows, longcfgs = configs(LONG)
    assert [r["condition_id"] for r in rows] == [r["condition_id"] for r in longrows]
    assert len(rows) == 5
    baseline = asdict(cfgs[0])
    for cfg, (section, key, value) in zip(
        cfgs[1:],
        [
            ("motor", "segment_force_placement", "root_three"),
            ("potentials", "spring_spring_repulsion", None),
            ("flagella", "initial_root_tangent_angle_deg", 39.276),
            ("potentials", "spring_spring_repulsion", None),
        ],
    ):
        expected = deepcopy(baseline)
        if cfg is cfgs[2]:
            expected["potentials"][key]["force_curve"] = "continuous_cutoff"
        elif cfg is cfgs[4]:
            expected["potentials"][key]["body_flagella_enabled"] = True
        else:
            expected[section][key] = value
        assert asdict(cfg) == expected
    for left, right in zip(cfgs, longcfgs):
        a, b = asdict(left), asdict(right)
        a.pop("time")
        b.pop("time")
        assert a == b
        assert left.time.duration_s == pytest.approx(0.04)
        assert right.time.duration_s == pytest.approx(1.5)
        assert right.tau_s * right.dt_star == pytest.approx(4e-6)
    assert short["development_evaluation"]["stage"] == "short_screen"
    assert long["development_evaluation"]["stage"] == "long_duration"
    for suffix in ["1tau", "1p5s"]:
        job = load_parallel_job(
            ROOT
            / f"conf/phase2_parallel/issue263_initial_bundle/nf3_physics_{suffix}_job.yaml"
        )
        assert resolve_execution(job, None).max_workers == 5
        assert job.condition_ids == tuple(r["condition_id"] for r in rows)
    assert len(geometry_preflight(short, rows)) == 5


@pytest.mark.parametrize("deformed", [False, True])
def test_local_motor_preserves_full_torque_body_reaction_and_zero_net_force(deformed):
    _, _, cfgs = configs()
    cfg = cfgs[0]
    m = ModelBuilder(cfg).build()
    p = m.positions_m.copy()
    if deformed:
        p[m.flagella_indices[0][1:4]] += (
            np.array([[0.01, 0.02, -0.01], [-0.02, 0.01, 0.02], [0.01, -0.02, 0.01]])
            * cfg.b_m
        )
    kwargs = dict(
        positions_m=p,
        flagella_indices=m.flagella_indices,
        body_indices=m.body_indices,
        torque_per_flag=np.full(3, cfg.motor_torque_Nm),
        segment_weights=[np.exp(-np.arange(10) / 3) for _ in range(3)],
        full_vector_body_reaction=True,
        body_reaction_support="attach_one_ring",
        flagella_attach_body_indices=m.flagella_attach_body_indices,
        body_ring_edges=m.body_ring_edges,
        body_vertical_edges=m.body_vertical_edges,
    )
    baseline, _ = compute_root_torque_segment_couples_forces(**kwargs)
    local, diag = compute_root_torque_segment_couples_forces(
        **kwargs, segment_force_placement="root_three"
    )
    np.testing.assert_array_equal(local[m.body_indices], baseline[m.body_indices])
    for ids in m.flagella_indices:
        np.testing.assert_array_equal(local[ids[3:]], 0.0)
        np.testing.assert_allclose(
            np.cross(p[ids] - p[ids[0]], local[ids]).sum(0),
            np.cross(p[ids] - p[ids[0]], baseline[ids]).sum(0),
            rtol=1e-12,
            atol=1e-32,
        )
        assert np.linalg.norm(local[ids].sum(0)) < 1e-25
    assert np.linalg.norm(np.cross(p, local).sum(0)) < 1e-31
    assert diag.degenerate_axis_count == 0
    assert (
        dict(diag.drive_metrics)["motor_root_localization_torque_error_ratio_max"]
        < 1e-12
    )


def test_degenerate_local_root_fails_instead_of_changing_torque_or_support():
    _, _, cfgs = configs()
    cfg = cfgs[0]
    m = ModelBuilder(cfg).build()
    p = m.positions_m.copy()
    ids = m.flagella_indices[0]
    p[ids[:3]] = p[ids[0]] + np.array([[0, 0, 0], [0.2, 0, 0], [0.4, 0, 0]]) * cfg.b_m
    with pytest.raises(RuntimeError, match="localization failed"):
        compute_root_torque_segment_couples_forces(
            positions_m=p,
            flagella_indices=m.flagella_indices[:1],
            body_indices=m.body_indices,
            torque_per_flag=np.array([cfg.motor_torque_Nm]),
            segment_weights=[np.ones(10)],
            full_vector_body_reaction=True,
            segment_force_placement="root_three",
        )


def repulsion(distance, curve):
    p = (
        np.array([[-0.3, 0, 0], [0.3, 0, 0], [0, distance, -0.3], [0, distance, 0.3]])
        * 1e-6
    )
    f = compute_segment_repulsion_forces(
        p,
        np.array([[0, 1], [2, 3]]),
        np.array([[0, 1]]),
        2.5e-20,
        0.2e-6,
        0.2e-6,
        force_curve=curve,
    )
    return p, f


@pytest.mark.parametrize("distance", [0.01, 0.08, 0.15, 0.199])
def test_continuous_repulsion_is_potential_gradient_and_conserves_momentum(distance):
    p, f = repulsion(distance, "continuous_cutoff")
    epsilon = 1e-12

    def u(d):
        return continuous_segment_repulsion_potential(d, 2.5e-20, 0.2e-6, 0.2e-6)

    gradient = -(u(distance * 1e-6 + epsilon) - u(distance * 1e-6 - epsilon)) / (
        2 * epsilon
    )
    assert f[2:, 1].sum() == pytest.approx(gradient, rel=1e-7, abs=1e-23)
    assert np.linalg.norm(f.sum(0)) < 1e-27
    assert np.linalg.norm(np.cross(p, f).sum(0)) < 1e-32


def test_cutoff_continuity_near_contact_scale_and_legacy_default():
    _, near = repulsion(1e-7, "continuous_cutoff")
    assert near[2:, 1].sum() == pytest.approx(0.125e-12, rel=1e-6)
    _, inside = repulsion(0.2 - 1e-9, "continuous_cutoff")
    assert inside[2:, 1].sum() < 1e-20
    for d in [0.2, 0.21]:
        _, f = repulsion(d, "continuous_cutoff")
        np.testing.assert_array_equal(f, 0)
    p, legacy = repulsion(0.1, "exponential_cutoff")
    f = compute_segment_repulsion_forces(
        p, np.array([[0, 1], [2, 3]]), np.array([[0, 1]]), 2.5e-20, 0.2e-6, 0.2e-6
    )
    np.testing.assert_array_equal(f, legacy)
    assert f[2:, 1].sum() == pytest.approx(0.125e-12 * np.exp(-0.5))


def test_root_translation_preserves_shape_phase_hook_and_symmetry():
    _, _, cfgs = configs()
    b = ModelBuilder(cfgs[0]).build()
    g = ModelBuilder(cfgs[3]).build()
    np.testing.assert_array_equal(
        b.positions_m[b.body_indices], g.positions_m[g.body_indices]
    )
    np.testing.assert_allclose(
        b.spring_rest_lengths_m, g.spring_rest_lengths_m, rtol=0, atol=1e-21
    )
    np.testing.assert_array_equal(
        b.flagella_initial_phases_rad, g.flagella_initial_phases_rad
    )
    for ids, attach in zip(g.flagella_indices, g.flagella_attach_body_indices):
        delta = g.positions_m[ids] - b.positions_m[ids]
        np.testing.assert_allclose(
            delta, np.broadcast_to(delta[0], delta.shape), rtol=0, atol=1e-21
        )
        assert np.linalg.norm(
            g.positions_m[ids[0]] - g.positions_m[attach]
        ) == pytest.approx(0.25e-6)
    for row in g.initial_geometry_metadata["filaments"]:
        assert row["attachment_motor_axis_distance_after_um"] == pytest.approx(
            0.1568060056
        )
    evaluation = DynamicsEngine(g, cfgs[3]).step(0.0)
    assert (
        np.linalg.norm(
            evaluation.spring_forces
            + evaluation.bend_forces
            + evaluation.torsion_forces
            + evaluation.hook_forces
        )
        < 1e-24
    )
    assert g.initial_geometry_metadata["minimum_bead_distance_um"] == pytest.approx(
        0.25
    )
    angle = 2 * np.pi / 3
    rotation = np.array(
        [
            [1, 0, 0],
            [0, np.cos(angle), -np.sin(angle)],
            [0, np.sin(angle), np.cos(angle)],
        ]
    )
    np.testing.assert_allclose(
        g.positions_m[g.flagella_indices[0]] @ rotation.T,
        g.positions_m[g.flagella_indices[1]],
        atol=1e-20,
    )


@pytest.mark.parametrize("angle", [float("nan"), 90.0, -90.0])
def test_invalid_root_angle_is_rejected(angle):
    _, _, cfgs = configs()
    with pytest.raises(ValueError, match="initial_root_tangent_angle_deg"):
        cfgs[0].with_overrides({"flagella": {"initial_root_tangent_angle_deg": angle}})


@pytest.mark.parametrize("index", [0, 1, 2, 3, 4])
def test_engine_bounded_steps_and_provenance(index):
    _, _, cfgs = configs()
    cfg = cfgs[index]
    model = ModelBuilder(cfg).build()
    engine = DynamicsEngine(model, cfg)
    for _ in range(4):
        d = engine.step(cfg.dt_star)
        assert np.isfinite(model.positions_m).all()
        assert np.isfinite(d.motor_forces).all()
    manifest = cfg.implementation_manifest()
    assert manifest["dynamics"]["implementation_status"] == "implemented"
