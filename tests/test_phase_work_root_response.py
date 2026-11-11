"""Virtual-work and factor-isolation tests, not a swimming acceptance gate."""

from pathlib import Path

import numpy as np
import pytest
import yaml
from scipy.spatial.transform import Rotation

from sim_swim.analysis.model_development_evaluation import screen_status
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    normalize_campaign_config,
)
from sim_swim.analysis.online_run_summary import OnlineRunSummary
from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.dynamics.forces import _zero_net_force_vector_torque_drive
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.debug_summary import (
    BodyConstraintDiagnosticsRecorder,
    StepSummaryRecorder,
)
from sim_swim.sim.params import SimulationConfig
from sim_swim.sim.rotation_diagnostics import rigid_fit


def configurations(stage="1tau"):
    raw = yaml.safe_load(
        Path(f"conf/phase2_multi_run/2010_hex_nf3_mk_{stage}_issue263.yaml").read_text()
    )
    campaign = normalize_campaign_config(raw)
    base = SimulationConfig.from_dict(
        yaml.safe_load(Path(campaign["base_config"]).read_text())
    )
    return [
        base.with_overrides(c["config_overrides"])
        for c in build_campaign_conditions(campaign)
    ]


@pytest.fixture
def setup():
    cfg = configurations()[1]
    model = ModelBuilder(cfg).build()
    return cfg, model, DynamicsEngine(model, cfg).phase_drive


def oracle_phases(drive, positions):
    """Independent rigid-fit/SVD/rotation implementation of chi, no gradients."""
    phases = []
    for ref in drive.references:
        rb, _ = rigid_fit(ref.body_fit.reference, positions[ref.body_ids])
        rf, _ = rigid_fit(ref.flag_fit.reference, positions[ref.flag_ids])
        pts = positions[ref.flag_ids]
        e = np.linalg.svd(pts - pts.mean(0))[2][0]
        if e @ (pts[-1] - pts[0]) < 0:
            e = -e
        a = rf @ ref.output_frame[:, 2]
        cross = np.cross(a, e)
        sine = np.linalg.norm(cross)
        align = (
            np.eye(3)
            if sine < 1e-14
            else Rotation.from_rotvec(
                cross / sine * np.arctan2(sine, a @ e)
            ).as_matrix()
        )
        quat = Rotation.from_matrix(
            (rb @ ref.input_frame).T @ align @ rf @ ref.output_frame
        ).as_quat()
        z, w = quat[2:]
        phases.append((w + 1j * z) ** 2 / (w * w + z * z))
    return np.asarray(phases)


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_analytic_force_is_virtual_work_for_nonrigid_displacements(setup, seed):
    cfg, model, drive = setup
    rng = np.random.default_rng(seed)
    q = model.positions_m.copy() + rng.normal(size=model.positions_m.shape) * 8e-9
    torque = np.array([1.0, -1.0, 0.7]) * cfg.motor_torque_Nm
    force, diag = drive.forces(q, torque)
    assert np.linalg.norm(force.sum(0)) < 1e-26
    assert np.linalg.norm(np.cross(q, force).sum(0)) < 1e-31
    assert dict(diag.drive_metrics)["motor_drive_target_vector_error_ratio_max"] < 1e-10
    np.testing.assert_allclose(drive.phases(q), oracle_phases(drive, q), atol=1e-13)
    for _ in range(8):
        d = rng.normal(size=q.shape)
        d /= np.linalg.norm(d)
        eps = 5e-13
        virtual_work = (
            torque
            @ np.angle(
                oracle_phases(drive, q + eps * d) / oracle_phases(drive, q - eps * d)
            )
            / (2 * eps)
        )
        np.testing.assert_allclose(
            np.sum(force * d), virtual_work, atol=2e-22, rtol=2e-6
        )
    # A whole-coordinate difference checks individual force components as well.
    fd = np.zeros_like(q)
    if seed == 1:
        for i in range(len(q)):
            for j in range(3):
                plus, minus = q.copy(), q.copy()
                plus[i, j] += eps
                minus[i, j] -= eps
                fd[i, j] = (
                    torque
                    @ np.angle(oracle_phases(drive, plus) / oracle_phases(drive, minus))
                    / (2 * eps)
                )
        np.testing.assert_allclose(force, fd, atol=2e-22, rtol=2e-6)


def test_projection_torque_objectivity_unit_covariance_and_zero_input(setup):
    cfg, model, drive = setup
    q = model.positions_m.copy()
    rng = np.random.default_rng(8)
    q += rng.normal(size=q.shape) * 3e-9
    torques = np.array([cfg.motor_torque_Nm, 0.0, 0.0])
    f, _ = drive.forces(q, torques)
    ref = drive.references[0]
    values = drive._evaluate(q, ref)
    e, n = values[4:6]
    tf = np.cross(q[ref.flag_ids] - q[ref.flag_ids].mean(0), f[ref.flag_ids]).sum(0)
    np.testing.assert_allclose(tf, torques[0] * (n + e) / (1 + n @ e), atol=1e-31)
    assert tf @ e == pytest.approx(torques[0], rel=1e-10, abs=1e-31)
    assert tf @ n == pytest.approx(torques[0], rel=1e-10, abs=1e-31)
    assert abs(tf @ np.cross(n, e)) < 1e-31
    rotation = Rotation.random(random_state=rng).as_matrix()
    transformed = q @ rotation.T + np.array([2.0, -4.0, 3.0]) * 1e-6
    fg, _ = drive.forces(transformed, torques)
    np.testing.assert_allclose(fg, f @ rotation.T, atol=1e-25, rtol=1e-9)
    fg, _ = drive.forces(q * 1e6, torques)
    np.testing.assert_allclose(fg, f / 1e6, atol=1e-30, rtol=1e-9)
    zero, _ = drive.forces(q, np.zeros(3))
    np.testing.assert_array_equal(zero, np.zeros_like(q))


def test_transverse_phase_is_a_gauge_not_restoring_energy(setup):
    cfg, model, drive = setup
    q = model.positions_m.copy()
    q[model.flagella_indices[0][2], 2] += 2e-8
    torques = np.full(3, cfg.motor_torque_Nm)
    f, _ = drive.forces(q, torques)
    before = drive.phases(q)
    for ref in drive.references:
        ref.input_frame = (
            ref.input_frame @ Rotation.from_rotvec([0, 0, 0.37]).as_matrix()
        )
        ref.output_frame = (
            ref.output_frame @ Rotation.from_rotvec([0, 0, -0.53]).as_matrix()
        )
    fg, _ = drive.forces(q, torques)
    np.testing.assert_allclose(fg, f, atol=1e-26, rtol=1e-10)
    np.testing.assert_allclose(np.angle(drive.phases(q) / before), -0.9, atol=1e-12)


def test_full_spin_has_correct_work_without_wrapped_phase_jumps(setup):
    cfg, model, drive = setup
    q0 = model.positions_m.copy()
    ref = drive.references[0]
    axis = drive._evaluate(q0, ref)[4]
    center = q0[ref.flag_ids[0]]
    phases = []
    before = q0.copy()
    torques = np.array([cfg.motor_torque_Nm, 0, 0])
    for angle in np.linspace(0, 2 * np.pi, 73):
        q = q0.copy()
        rotation = Rotation.from_rotvec(axis * angle).as_matrix()
        q[ref.flag_ids] = (q0[ref.flag_ids] - center) @ rotation.T + center
        phases.append(drive.phases(q)[0])
        if len(phases) > 1:
            force, _ = drive.forces(before, torques)
            drive.record_displacement(before, q, torques, force)
        before = q
    assert np.angle(np.asarray(phases[1:]) / phases[:-1]).sum() == pytest.approx(
        2 * np.pi, abs=1e-11
    )
    assert drive.work_total_J == pytest.approx(
        2 * np.pi * cfg.motor_torque_Nm, rel=1e-10, abs=1e-31
    )


def test_not_just_a_minimum_norm_wrench_on_a_deformable_chain(setup):
    cfg, model, drive = setup
    q = model.positions_m.copy()
    ref = drive.references[0]
    q[ref.flag_ids[3], 1] += 8e-8
    torque = cfg.motor_torque_Nm
    f, _ = drive.forces(q, np.array([torque, 0, 0]))
    e, n = drive._evaluate(q, ref)[4:6]
    minimum, fail, _, _ = _zero_net_force_vector_torque_drive(
        positions_m=q,
        indices=ref.flag_ids,
        origin=q[ref.flag_ids].mean(0),
        target_torque_Nm=torque * (n + e) / (1 + n @ e),
    )
    assert fail == 0
    assert np.linalg.norm(
        f[ref.flag_ids] - minimum[ref.flag_ids]
    ) > 1e-4 * np.linalg.norm(f[ref.flag_ids])


@pytest.mark.parametrize("kind", ["line", "pca", "body", "antiparallel"])
def test_singular_fits_and_reverse_axes_reject_without_fallback(setup, kind):
    cfg, model, drive = setup
    q = model.positions_m.copy()
    ref = drive.references[0]
    if kind == "line":
        q[ref.flag_ids] = np.c_[
            np.arange(len(ref.flag_ids)) * 1e-6, np.zeros((len(ref.flag_ids), 2))
        ]
    elif kind == "pca":
        theta = np.linspace(0, 2 * np.pi, len(ref.flag_ids), endpoint=False)
        q[ref.flag_ids] = (
            np.c_[np.cos(theta), np.sin(theta), np.zeros_like(theta)] * 1e-6
        )
    elif kind == "body":
        q[ref.body_ids] = np.c_[
            np.arange(len(ref.body_ids)) * 1e-6, np.zeros((len(ref.body_ids), 2))
        ]
    else:
        e, n = drive._evaluate(q, ref)[4:6]
        cross = np.cross(e, -n)
        rotation = Rotation.from_rotvec(
            cross / np.linalg.norm(cross) * np.arccos(np.clip(e @ (-n), -1, 1))
        ).as_matrix()
        q[ref.flag_ids] = (q[ref.flag_ids] - q[ref.flag_ids].mean(0)) @ rotation.T + q[
            ref.flag_ids
        ].mean(0)
    with pytest.raises(RuntimeError, match="Phase drive"):
        drive.forces(q, np.full(3, cfg.motor_torque_Nm))


def test_matrix_keeps_geometry_material_and_duration_factors_separate():
    configs = configurations()
    models = [ModelBuilder(c).build() for c in configs]
    for cfg, model in zip(configs, models, strict=True):
        np.testing.assert_array_equal(model.positions_m, models[0].positions_m)
        assert model.root_bending_rest_rad == models[0].root_bending_rest_rad
        assert model.root_torsion_rest_rad == models[0].root_torsion_rest_rad
        assert cfg.dt_s == pytest.approx(4e-6) and cfg.total_steps == 10000
        assert cfg.motor.torque_Nm == cfg.motor.reference_torque_Nm == 2.5e-20
        assert cfg.fluid.viscosity_Pa_s == 0.001
        assert cfg.potentials.spring_spring_repulsion.body_flagella_enabled
    assert [c.potentials.bend.root_response for c in configs] == [
        "cosine",
        "cosine",
        "normalized_tangent",
        "normalized_tangent",
    ]
    for cfg in configurations("1s"):
        assert cfg.time.duration_s == pytest.approx(1.0)
        assert cfg.total_steps in (250000, 250001)
    assert (
        configurations()[0].implementation_manifest()["dynamics"]["force_distribution"]
        == "effective_vector_torque"
    )
    manifest = configs[1].implementation_manifest()["dynamics"]
    assert (
        manifest["drive_contract_version"] == 2
        and not manifest["diffusive_activity_used"]
    )


def test_K_is_passive_energy_gradient_with_fixed_natural_angles():
    cfg = configurations()[2]
    model = ModelBuilder(cfg).build()
    eng = DynamicsEngine(model, cfg)
    baseline_cfg = configurations()[0]
    baseline_model = ModelBuilder(baseline_cfg).build()
    control = DynamicsEngine(baseline_model, baseline_cfg)
    rng = np.random.default_rng(20)
    q = model.positions_m.copy() + rng.normal(size=model.positions_m.shape) * 1e-9
    model.positions_m = q.copy()
    baseline_model.positions_m = q.copy()
    force = eng.step(0).bend_forces - control.step(0).bend_forces
    roots = dict(model.root_bending_rest_rad)
    for row, theta in roots.items():
        assert eng.root_bend_factors[row] * np.sin(theta) ** 2 == pytest.approx(
            np.sin(np.radians(142)) ** 2
        )

    def energy(q):
        out = 0.0
        for row, theta0 in roots.items():
            a, b, c = q[model.bending_triplets[row]]
            u = a - b
            v = c - b
            cosine = u @ v / np.linalg.norm(u) / np.linalg.norm(v)
            out += (
                0.5
                * eng.k_bend
                * (eng.root_bend_factors[row] - 1)
                * (cosine - np.cos(theta0)) ** 2
            )
        return out

    for _ in range(8):
        d = rng.normal(size=q.shape)
        d /= np.linalg.norm(d)
        eps = 5e-13
        fd = -(energy(q + eps * d) - energy(q - eps * d)) / (2 * eps)
        np.testing.assert_allclose(np.sum(force * d), fd, rtol=2e-6, atol=1e-22)
    assert model.root_bending_rest_rad == roots
    nonroot = np.ones(len(model.bending_triplets), bool)
    nonroot[list(roots)] = False
    np.testing.assert_array_equal(eng.root_bend_factors[nonroot], 1.0)


@pytest.mark.parametrize("index", range(4))
def test_bounded_engine_telemetry_and_common_QC(index, tmp_path):
    cfg = configurations()[index]
    model = ModelBuilder(cfg).build()
    eng = DynamicsEngine(model, cfg)
    rec = StepSummaryRecorder(model, cfg, tmp_path)
    body_rec = BodyConstraintDiagnosticsRecorder(model, cfg, tmp_path, emit_csv=False)
    summary = OnlineRunSummary(expected_steps=32)
    for step in range(33):
        diag = eng.step(0 if step == 0 else cfg.dt_star)
        rec.record(step, step * cfg.dt_star, diag)
        summary.record(rec.last_row)
        body_rec.record(step, step * cfg.dt_s, diag)
        summary.record_body(body_rec.last_row)
    assert np.isfinite(model.positions_m).all()
    if index in (1, 3):
        metrics = diag.force_balance_diagnostics
        assert metrics["motor_drive_contract_version"] == 2
        assert np.isnan(metrics["motor_drive_axial_error_ratio_max"])
        assert metrics["motor_drive_output_projection_error_ratio_max"] < 1e-10
        assert metrics["motor_drive_target_vector_error_ratio_max"] < 1e-10
        assert metrics["motor_drive_phase_work_total_J"] > 0
        assert summary.maximum_times["motor_drive_phase_work_total_J"] == pytest.approx(
            32 * cfg.dt_s
        )
        assert all(np.all(a == 0) for a in eng.local_twist_orientations)
    qc = yaml.safe_load(
        Path("conf/phase2_multi_run/2010_hex_nf3_mk_1tau_issue263.yaml").read_text()
    )["development_evaluation"]["qc"]
    doc = summary.document(
        run_dir=tmp_path, completed=True, reason="", time_manifest={}, policy="compact"
    )
    assert screen_status(doc, qc=qc) == "pass"
    rec.write_csv()


def test_invalid_root_response_and_phase_reaction_modes_reject():
    cfg = configurations()[1]
    with pytest.raises(ValueError, match="without segment modifiers"):
        cfg.with_overrides({"motor": {"body_reaction_support": "all_body"}})
    with pytest.raises(ValueError, match="root_response"):
        cfg.with_overrides({"potentials": {"bend": {"root_response": "unknown"}}})
    with pytest.raises(ValueError, match="fixed axis_transition_rear"):
        cfg.with_overrides(
            {
                "potentials": {"bend": {"root_response": "normalized_tangent"}},
                "flagella": {
                    "initial_bundle_mode": "compact_rear",
                    "initial_hook_force_neutral": True,
                },
            }
        )
