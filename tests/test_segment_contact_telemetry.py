from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from sim_swim.analysis.online_run_summary import OnlineRunSummary
from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.dynamics.forces import compute_segment_repulsion_forces
from sim_swim.dynamics.segment_contact_diagnostics import SegmentContactTracker
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig


def crossing():
    # Endpoints are all outside cutoff, but segment interiors are only0.1um apart.
    return np.array([[-1, 0, 0], [1, 0, 0], [0, 0.1, -1], [0, 0.1, 1]]) * 1e-6


@pytest.mark.parametrize(
    "owners,body,group",
    [
        ([-1, -1, 0, 0], [True, True, False, False], "bf"),
        ([0, 0, 1, 1], [False] * 4, "inter_ff"),
        ([0, 0, 0, 0], [False] * 4, "self_ff"),
    ],
)
def test_actual_segment_contact_detects_interior_contact_missed_by_beads(
    owners, body, group
):
    q = crossing()
    springs = np.array([[0, 1], [2, 3]])
    pairs = np.array([[0, 1]])
    tracker = SegmentContactTracker(springs, pairs, np.array(body), np.array(owners))
    step = tracker.new_step()
    kw = dict(
        positions_m=q,
        spring_pairs=springs,
        segment_pair_indices=pairs,
        a_ss=2.5e-20,
        cutoff=0.2e-6,
        a_length=0.1e-6,
    )
    original = compute_segment_repulsion_forces(**kw)
    actual = compute_segment_repulsion_forces(**kw, contact_diagnostics=step)
    np.testing.assert_array_equal(actual, original)
    np.testing.assert_allclose(step.forces.sum(0), original, rtol=0, atol=0)
    m = tracker.metrics(step, evaluation_t_s=0.25, dt_s=4e-6)
    prefix = f"segment_contact_{group}_"
    assert m[prefix + "min_distance_um"] == pytest.approx(0.1)
    assert m[prefix + "active_pairs"] == 1
    assert m[prefix + "witness_s"] == pytest.approx(0.5)
    assert m[prefix + "witness_t"] == pytest.approx(0.5)
    assert m[prefix + "max_pair_force_N"] == pytest.approx(2.5e-13 * np.exp(-1))
    assert m[prefix + "max_bead_force_N"] == pytest.approx(
        m[prefix + "max_pair_force_N"] / 2
    )
    assert m[prefix + "first_active_t_s"] == 0.25
    assert m[prefix + "active_time_s"] == 4e-6
    assert np.linalg.norm(q[0] - q[2]) > 0.2e-6
    # Frozen queries do not mutate temporal history.
    tracker.metrics(step, evaluation_t_s=0.1, dt_s=0)
    later = tracker.metrics(tracker.new_step(), evaluation_t_s=0.26, dt_s=4e-6)
    assert later[prefix + "first_active_t_s"] == 0.25
    assert later[prefix + "active_time_s"] == 4e-6
    assert later["segment_contact_recorded_steps"] == 2


def test_empty_active_pair_list_and_cutoff_are_not_contact():
    q = crossing()
    springs = np.array([[0, 1], [2, 3]])
    for pairs in [np.empty((0, 2), int), np.array([[0, 1]])]:
        t = SegmentContactTracker(
            springs, pairs, np.zeros(4, bool), np.array([0, 0, 1, 1])
        )
        s = t.new_step()
        f = compute_segment_repulsion_forces(
            q, springs, pairs, 1.0, 0.1e-6, 0.1e-6, contact_diagnostics=s
        )
        m = t.metrics(s, evaluation_t_s=0, dt_s=4e-6)
        assert m["segment_contact_inter_ff_active_pairs"] == 0
        assert m["segment_contact_inter_ff_first_active_t_s"] == -1
        assert not np.any(f)
        if not len(pairs):
            assert m["segment_contact_inter_ff_min_distance_um"] == -1


def raw_config():
    raw = yaml.safe_load(Path("conf/sim_swim_2010_hex.yaml").read_text())
    raw.setdefault("flagella", {}).update(
        n_flagella=3,
        placement_mode="explicit_center_layer_slots",
        attachment_slots=[0, 2, 4],
        initial_bundle_mode="axis_transition_rear",
        initial_phase_offsets_deg=[0, 0, 0],
    )
    raw.setdefault("motor", {}).update(
        force_distribution="effective_vector_torque",
        enable_switching=False,
        body_reaction_full_vector=True,
        body_reaction_support="attach_one_ring",
        segment_torque_correction="none",
        segment_torque_distribution="fixed_weights",
        segment_force_placement="distributed",
    )
    raw.setdefault("potentials", {}).setdefault("torsion", {})["force_method"] = (
        "analytic"
    )
    raw.setdefault("brownian", {})["enabled"] = False
    return raw


def test_opt_in_engine_preserves_forces_and_trajectory_bitwise():
    raw = raw_config()
    on = deepcopy(raw)
    on.setdefault("output", {})["segment_contact_diagnostics"] = True
    c0, c1 = SimulationConfig.from_dict(raw), SimulationConfig.from_dict(on)
    assert not c0.output.segment_contact_diagnostics
    a, b = (
        DynamicsEngine(ModelBuilder(c0).build(), c0),
        DynamicsEngine(ModelBuilder(c1).build(), c1),
    )
    for _ in range(16):
        d0, d1 = a.step(1e-4), b.step(1e-4)
        np.testing.assert_array_equal(d0.total_forces, d1.total_forces)
        np.testing.assert_array_equal(d0.repulsion_forces, d1.repulsion_forces)
        np.testing.assert_array_equal(d0.positions_after_m, d1.positions_after_m)
    assert d1.force_balance_diagnostics["segment_contact_recorded_steps"] == 16
    b.step(0)
    assert b.segment_contact_tracker.steps == 16
    raw.setdefault("output", {})["segment_contact_diagnostics"] = "true"
    with pytest.raises(ValueError, match="segment_contact_diagnostics must be boolean"):
        SimulationConfig.from_dict(raw)


def test_contact_maximum_is_timestamped_at_force_evaluation():
    s = OnlineRunSummary(expected_steps=1)
    s.record(
        {
            "t_s": 0.1,
            "force_evaluation_t_s": 0.099996,
            "segment_contact_bf_max_pair_force_N": 1.0,
            "segment_contact_bf_active_time_s": 4e-6,
        }
    )
    assert s.maximum_times["segment_contact_bf_max_pair_force_N"] == 0.099996
    assert s.maximum_times["segment_contact_bf_active_time_s"] == 0.1


def test_contact_campaign_is_an_exact_one_factor_comparison():
    from sim_swim.analysis.multi_run_campaign import (
        build_campaign_conditions,
        normalize_campaign_config,
    )

    for stage in ("1tau", "1s"):
        path = Path(
            f"conf/phase2_multi_run/2010_hex_nf3_connector_contact_{stage}_issue263.yaml"
        )
        raw = yaml.safe_load(path.read_text())
        campaign = normalize_campaign_config(raw)
        conditions = build_campaign_conditions(campaign)
        assert [c["condition_id"] for c in conditions] == [
            "nf03__slots024__J_bfon_control",
            "nf03__slots024__J_bfoff_control",
        ]
        assert len(conditions) == 2
        a, b = [deepcopy(c["config_overrides"]) for c in conditions]
        assert (
            a["potentials"]["spring_spring_repulsion"].pop("body_flagella_enabled")
            is True
        )
        assert (
            b["potentials"]["spring_spring_repulsion"].pop("body_flagella_enabled")
            is False
        )
        assert a == b
        assert a["output"]["segment_contact_diagnostics"] is True
        assert a["motor"]["force_distribution"] == "effective_vector_torque"
        assert a["potentials"]["bend"]["root_response"] == "cosine"
        dev = raw["development_evaluation"]
        assert dev["geometry"]["require_initial_segment_clearance"]
        assert dev["expected_condition_count"] == 2
        assert dev["qc"]["max_motor_torque_balance_residual_ratio"] == 0.02


def test_opt_in_contact_metrics_reach_compact_run_summary(tmp_path):
    import json

    from sim_swim.sim.core import Simulator

    raw = raw_config()
    raw["time"] = {
        "duration_s": 0.000016,
        "dt_s": 0.001,
        "integration": {"dt_star": 0.0001},
        "scale_policy": "reference_torque",
    }
    raw.setdefault("output", {}).update(
        segment_contact_diagnostics=True, policy="compact"
    )
    cfg = SimulationConfig.from_dict(raw)
    Simulator(cfg).run(cfg.time.duration_s, step_summary_dir=tmp_path)
    summary = json.loads((tmp_path / "run_summary.json").read_text())
    metrics = summary["all_step_metrics"]
    assert metrics["segment_contact_available"]["final"] == 1
    assert metrics["segment_contact_recorded_steps"]["final"] == cfg.total_steps
    assert metrics["segment_contact_bf_candidate_pairs"]["final"] > 0
    assert metrics["segment_contact_bf_first_active_t_s"]["final"] == -1
    assert metrics["segment_contact_bf_min_distance_um"]["min"] > 0.2
