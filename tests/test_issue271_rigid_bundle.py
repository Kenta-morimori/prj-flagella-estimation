import hashlib
import json
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
from sim_swim.model.builder import ModelBuilder
from sim_swim.model.rigid_common_bundle import PATTERNS
from sim_swim.sim.debug_summary import StepSummaryRecorder, _triplet_angles_rad
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]
STATE = ROOT / "conf/phase2_initial/2010_hex_rigid_common_bundle_v1.json"
CAMPAIGN = ROOT / "conf/phase2_multi_run/2010_hex_rigid_bundle_1tau_issue271.yaml"
OLD = (
    ROOT
    / "conf/phase2_multi_run/2010_hex_dense_count13_body_contact_1tau_issue263.yaml"
)
JOB = ROOT / "conf/phase2_parallel/issue271_rigid_bundle/short_1tau_job.yaml"


def _conditions(path):
    campaign = normalize_campaign_config(load_yaml(path))
    base = SimulationConfig.from_dict(load_yaml(ROOT / campaign["base_config"]))
    return campaign, [
        (row["condition_id"], base.with_overrides(row["config_overrides"]))
        for row in build_campaign_conditions(campaign)
    ]


def _shape_after(cfg, count, tmp_path):
    model = ModelBuilder(cfg).build()
    engine = DynamicsEngine(model, cfg)
    recorder = StepSummaryRecorder(model, cfg, tmp_path, emit_csv=False)
    rows = []
    for step in range(count):
        diag = engine.step(cfg.dt_star)
        recorder.record(step=step, t_star=engine.t_star, diag=diag)
        rows.append((recorder.last_row.copy(), diag))
    return rows


def test_all_13_rigid_geometries_preserve_full_shape_and_clearance():
    campaign, cases = _conditions(CAMPAIGN)
    base = SimulationConfig.from_dict(load_yaml(ROOT / campaign["base_config"]))
    digest = hashlib.sha256(STATE.read_bytes()).hexdigest()
    assert all(cfg.flagella.initial_bundle_state_sha256 == digest for _, cfg in cases)
    for slots in PATTERNS:
        cfg = base.with_overrides(
            {
                "flagella": {
                    "n_flagella": len(slots),
                    "placement_mode": "explicit_center_layer_slots",
                    "attachment_slots": list(slots),
                    "initial_bundle_mode": "rigid_common_axis_rear",
                    "initial_bundle_state_path": str(STATE),
                    "initial_bundle_state_sha256": digest,
                },
                "seed": {"attach_seed": 0, "phase_seed": 0},
            }
        )
        model = ModelBuilder(cfg).build()
        nominal = ModelBuilder(
            cfg.with_overrides(
                {
                    "flagella": {
                        "initial_bundle_mode": "none",
                        "initial_bundle_state_path": None,
                        "initial_bundle_state_sha256": None,
                    }
                }
            )
        ).build()
        np.testing.assert_array_equal(
            model.spring_rest_lengths_m, nominal.spring_rest_lengths_m
        )
        for ids in model.flagella_indices:
            assert len(ids) == 11
        state = model.initial_geometry_metadata
        assert state["whole_chain_rigid_shape_error_over_b"] < 1e-12
        assert state["hook_length_error_over_b"] < 1e-12
        assert state["min_self_segment_distance_um"] >= 0.2
        assert state["min_distal_body_segment_distance_um"] >= 0.2
        if len(slots) > 1:
            assert state["min_flag_segment_distance_um"] >= 0.2
            assert state["min_flag_bead_distance_um"] >= 0.2
        assert state["body_penetration_allowed_through_shaft_segment"] == 2


def test_short_job_has_five_separate_qualified_conditions_and_preflight():
    campaign, cases = _conditions(CAMPAIGN)
    assert [name for name, _ in cases] == [
        "nf01__slots0__base",
        "nf02__slots02__base",
        "nf03__slots024__base",
        "nf06__slots012345__base",
        "nf03__slots024__tilt1deg",
    ]
    assert campaign["development_evaluation"]["stage"] == "short_screen"
    assert campaign["development_evaluation"]["expected_condition_count"] == 5
    preflight = geometry_preflight(campaign, build_campaign_conditions(campaign))
    assert len(preflight) == 5
    assert all(row["initial_repulsion_force_norm_N"] == 0 for row in preflight.values())
    job = load_parallel_job(JOB)
    assert list(job.condition_ids) == [name for name, _ in cases]
    assert resolve_execution(job, None).max_workers == 5
    assert resolve_execution(job, None).worker_policy == "cs10_qualified"


def test_historical_failure_is_near_straight_torsion_instability(tmp_path):
    _, old = _conditions(OLD)
    for enabled in (True, False):
        cfg = next(
            cfg
            for name, cfg in old
            if name == ("nf02__slots02__bfon" if enabled else "nf02__slots02__bfoff")
        )
        model = ModelBuilder(cfg).build()
        rows = model.bending_triplets[model.bending_flag_ids >= 0]
        assert np.degrees(
            _triplet_angles_rad(model.positions_m, rows)
        ).max() == pytest.approx(179.2249531946, abs=1e-6)
        result = _shape_after(cfg, 2, tmp_path / str(enabled))
        assert result[0][0]["shape_pass_nonbody"]
        assert not result[1][0]["shape_pass_nonbody"]
        assert result[1][0]["flag_bond_rel_err_max"] == pytest.approx(
            3.0272269147, rel=1e-6
        )
        force = result[1][1]
        assert np.linalg.norm(force.repulsion_forces) == 0
        assert np.linalg.norm(force.torsion_forces) > 1e3 * np.linalg.norm(
            force.motor_forces
        )
    off = next(cfg for name, cfg in old if name == "nf02__slots02__bfoff")
    controlled_motor_off = off.with_overrides(
        {
            "motor": {
                "enabled": False,
                "allow_reference_torque_mismatch": True,
                "torque_for_forces_override_Nm": off.torque_for_forces_Nm,
            }
        }
    )
    controlled_motor_off = controlled_motor_off.with_overrides(
        {"time": {"integration": {"dt_star": off.dt_s / controlled_motor_off.tau_s}}}
    )
    assert controlled_motor_off.dt_s == pytest.approx(off.dt_s)
    assert controlled_motor_off.torque_for_forces_Nm == off.torque_for_forces_Nm
    assert not _shape_after(
        controlled_motor_off,
        2,
        tmp_path / "motor_off",
    )[-1][0]["shape_pass_nonbody"]
    assert _shape_after(
        off.with_overrides({"time": {"integration": {"dt_star": 5e-5}}}),
        4,
        tmp_path / "half_dt",
    )[-1][0]["shape_pass_nonbody"]
    assert _shape_after(
        off.with_overrides({"potentials": {"torsion": {"kt_over_T": 0.0}}}),
        2,
        tmp_path / "torsion_off_diagnostic_only",
    )[-1][0]["shape_pass_nonbody"]


def test_replacement_nf02_and_screen_cases_pass_first_ten_steps(tmp_path):
    _, cases = _conditions(CAMPAIGN)
    for name, cfg in cases:
        rows = _shape_after(cfg, 10, tmp_path / name)
        assert all(row["shape_pass_nonbody"] for row, _ in rows)
        assert max(row["flag_bond_rel_err_max"] for row, _ in rows) < 0.01


def test_asset_hash_mismatch_is_rejected():
    _, cases = _conditions(CAMPAIGN)
    cfg = cases[1][1].with_overrides(
        {"flagella": {"initial_bundle_state_sha256": "0" * 64}}
    )
    with pytest.raises(ValueError, match="SHA-256"):
        ModelBuilder(cfg).build()


def test_rehashed_internal_deformation_still_fails_recipe_check(tmp_path):
    _, cases = _conditions(CAMPAIGN)
    data = json.loads(STATE.read_text())
    record = next(
        row
        for row in data["conditions"]
        if row["attachment_slots"] == [0, 2] and row["variant"] == "base"
    )
    record["flagella_positions_over_b"][0][4][0] += 0.01
    altered = tmp_path / "changed.json"
    altered.write_text(json.dumps(data), encoding="utf-8")
    cfg = cases[1][1].with_overrides(
        {
            "flagella": {
                "initial_bundle_state_path": str(altered),
                "initial_bundle_state_sha256": hashlib.sha256(
                    altered.read_bytes()
                ).hexdigest(),
            }
        }
    )
    with pytest.raises(ValueError, match="pinned recipe"):
        ModelBuilder(cfg).build()
