from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import yaml

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    normalize_campaign_config,
)
from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.dynamics.forces import compute_effective_vector_torque_forces
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.debug_summary import StepSummaryRecorder
from sim_swim.sim.params import SimulationConfig


def configurations(stage="1tau"):
    campaign = normalize_campaign_config(
        yaml.safe_load(
            Path(
                f"conf/phase2_multi_run/2010_hex_nf3_dj_{stage}_issue263.yaml"
            ).read_text()
        )
    )
    raw = yaml.safe_load(Path(campaign["base_config"]).read_text())
    return [
        SimulationConfig.from_dict(raw).with_overrides(c["config_overrides"])
        for c in build_campaign_conditions(campaign)
    ]


def drive(model, q, torque):
    return compute_effective_vector_torque_forces(
        positions_m=q,
        flagella_indices=model.flagella_indices,
        flagella_attach_body_indices=model.flagella_attach_body_indices,
        body_ring_edges=model.body_ring_edges,
        body_vertical_edges=model.body_vertical_edges,
        torque_per_flag=np.full(3, torque),
    )


def test_effective_drive_has_full_vector_and_dual_work():
    cfg = configurations()[1]
    model = ModelBuilder(cfg).build()
    q = model.positions_m
    F, diag = drive(model, q, cfg.motor_torque_Nm)
    assert not diag.reaction_fallback_used
    assert len(diag.reaction_support_bead_counts) == 3
    assert np.linalg.norm(F.sum(0)) < 1e-24
    assert np.linalg.norm(np.cross(q, F).sum(0)) < 1e-30
    for ids in model.flagella_indices:
        r = q[ids] - q[ids].mean(0)
        axis = np.linalg.svd(r)[2][0]
        if axis @ (q[ids[-1]] - q[ids[0]]) < 0:
            axis = -axis
        target = cfg.motor_torque_Nm * axis
        np.testing.assert_allclose(np.cross(r, F[ids]).sum(0), target, atol=1e-30)
        v = np.random.default_rng(12).normal(size=r.shape) * 1e-6
        inertia = sum((x @ x) * np.eye(3) - np.outer(x, x) for x in r)
        omega = np.linalg.solve(inertia, np.cross(r, v - v.mean(0)).sum(0))
        np.testing.assert_allclose(np.sum(F[ids] * v), target @ omega, atol=1e-30)
        residual = v - v.mean(0) - np.cross(omega, r)
        assert abs(np.sum(F[ids] * residual)) < 1e-30


def test_effective_drive_is_objective_and_unit_covariant():
    cfg = configurations()[1]
    m = ModelBuilder(cfg).build()
    q = m.positions_m
    f, _ = drive(m, q, cfg.motor_torque_Nm)
    Q = np.array([[0.0, -1, 0], [0, 0, 1], [-1, 0, 0]])
    assert np.linalg.det(Q) > 0
    g, _ = drive(m, q @ Q.T + np.array([3.0, -2.0, 5.0]) * 1e-6, cfg.motor_torque_Nm)
    np.testing.assert_allclose(g, f @ Q.T, atol=1e-24, rtol=1e-9)
    g, _ = drive(m, q * 1e6, cfg.motor_torque_Nm)
    np.testing.assert_allclose(g, f / 1e6, atol=1e-30, rtol=1e-9)


def test_parallel_segment_does_not_make_effective_drive_singular():
    cfg = configurations()[1]
    m = ModelBuilder(cfg).build()
    q = m.positions_m.copy()
    ids = m.flagella_indices[0]
    q[ids[1]] = q[ids[0]] + np.array([-0.58e-6, 0, 0])
    f, _ = drive(m, q, cfg.motor_torque_Nm)
    assert np.isfinite(f).all() and np.linalg.norm(f) < 1e-12
    q[ids] = np.c_[np.arange(len(ids)) * 1e-6, np.zeros((len(ids), 2))]
    with pytest.raises(RuntimeError, match="degenerate"):
        drive(m, q, cfg.motor_torque_Nm)


def test_matrix_geometry_references_and_duration_are_separate():
    cs = configurations()
    ms = [ModelBuilder(c).build() for c in cs]
    np.testing.assert_array_equal(ms[0].positions_m, ms[1].positions_m)
    np.testing.assert_array_equal(ms[2].positions_m, ms[3].positions_m)
    for m, c in zip(ms, cs):
        assert c.total_steps == 10000 and c.dt_s == pytest.approx(4e-6)
        assert c.potentials.spring_spring_repulsion.body_flagella_enabled
        lengths = np.linalg.norm(
            m.positions_m[m.spring_pairs[:, 0]] - m.positions_m[m.spring_pairs[:, 1]],
            axis=1,
        )
        np.testing.assert_allclose(lengths, m.spring_rest_lengths_m, atol=1e-20)
    for m, cfg in zip(ms[2:], cs[2:]):
        assert len(m.root_bending_rest_rad) == len(m.root_torsion_rest_rad) == 9
        e = DynamicsEngine(m, cfg)
        theta, phi = e._state_angles_rad()
        for f, ids in enumerate(m.flagella_indices):
            for row, tri in enumerate(m.bending_triplets):
                if m.bending_flag_ids[row] == f and tri[0] >= ids[3]:
                    assert theta[row] == pytest.approx(np.radians(142))
        before = dict(m.root_bending_rest_rad)
        m.positions_m[30, 0] += 1e-10
        assert m.root_bending_rest_rad == before
    for c in configurations("3p5s"):
        assert c.total_steps == 875000


@pytest.mark.parametrize("index", range(4))
def test_reference_aware_recorder_and_bounded_integration(index, tmp_path):
    c = configurations()[index]
    m = ModelBuilder(c).build()
    e = DynamicsEngine(m, c)
    recorder = StepSummaryRecorder(m, c, tmp_path)
    diag = e.step(0)
    recorder.record(0, 0, diag)
    assert recorder.last_row["flag_bend_err_max_deg"] < 1e-8
    assert recorder.last_row["flag_torsion_err_max_deg"] < 1e-8
    for k in range(32):
        diag = e.step(c.dt_star)
    recorder.record(32, 32 * c.dt_star, diag)
    assert np.isfinite(m.positions_m).all()
    assert recorder.last_row["flag_bend_err_max_deg"] < 5
    if index in (1, 3):
        assert all(np.all(a == 0) for a in e.local_twist_orientations)
    recorder.write_csv()


def test_transition_phase_perturbation_keeps_hook_and_bond_lengths():
    c = configurations()[3]
    base = ModelBuilder(c).build()
    changed = ModelBuilder(
        c.with_overrides({"flagella": {"initial_phase_offsets_deg": [5, 0, 0]}})
    ).build()
    np.testing.assert_allclose(
        base.positions_m[base.hook_triplets[:, 1]],
        changed.positions_m[changed.hook_triplets[:, 1]],
        atol=1e-20,
    )
    np.testing.assert_allclose(
        base.spring_rest_lengths_m, changed.spring_rest_lengths_m, atol=1e-20
    )
    assert np.max(abs(base.positions_m[30:41] - changed.positions_m[30:41])) > 1e-9


def test_transition_rejects_state_switching_and_preprocessing():
    c = configurations()[3]
    for changes in [
        {"motor": {"enable_switching": True}},
        {"flagella": {"initial_hook_force_neutral": True}},
    ]:
        with pytest.raises(ValueError, match="axis_transition_rear"):
            c.with_overrides(changes)
