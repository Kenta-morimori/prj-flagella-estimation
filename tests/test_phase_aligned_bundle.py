import copy
import hashlib
import json
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
from sim_swim.model.builder import ModelBuilder
from sim_swim.model.initial_bundle import segment_distance_squared
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]
SHORT = (
    ROOT
    / "conf/phase2_multi_run/2010_hex_dense_count13_body_contact_1tau_issue263.yaml"
)
LONG = (
    ROOT / "conf/phase2_multi_run/2010_hex_dense_count13_body_contact_2s_issue263.yaml"
)
STATE = ROOT / "conf/phase2_initial/2010_hex_phase_aligned_bundle_v1.json"


def configurations(path=SHORT):
    campaign = normalize_campaign_config(load_yaml(path))
    base = SimulationConfig.from_dict(load_yaml(ROOT / campaign["base_config"]))
    conditions = build_campaign_conditions(campaign)
    return (
        campaign,
        conditions,
        [base.with_overrides(r["config_overrides"]) for r in conditions],
    )


def test_26_conditions_are_13_exact_geometry_and_mechanics_pairs():
    campaign, conditions, configs = configurations()
    assert len(configs) == 26
    assert sorted({c.flagella.n_flagella for c in configs}) == list(range(1, 7))
    preflight = geometry_preflight(campaign, conditions)
    for left, right, a, b in zip(
        configs[::2], configs[1::2], conditions[::2], conditions[1::2]
    ):
        da, db = asdict(left), asdict(right)
        assert (
            da["potentials"]["spring_spring_repulsion"].pop("body_flagella_enabled")
            is True
        )
        assert (
            db["potentials"]["spring_spring_repulsion"].pop("body_flagella_enabled")
            is False
        )
        assert da == db
        m = ModelBuilder(left).build()
        np.testing.assert_array_equal(
            m.positions_m, ModelBuilder(right).build().positions_m
        )
        nominal = ModelBuilder(
            left.with_overrides(
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
            m.spring_rest_lengths_m, nominal.spring_rest_lengths_m
        )
        for row in (a, b):
            p = preflight[row["condition_id"]]
            assert p["initial_geometry_state"]["common_phase_deg"] == 60
            assert p["initial_geometry_state"]["min_body_segment_distance_um"] > 0.21
            if left.flagella.n_flagella > 1:
                assert (
                    p["initial_geometry_state"]["min_flag_segment_distance_um"]
                    >= 0.21 - 1e-8
                )
            assert p["initial_hook_len_max_abs_error_m"] < 1e-15
            assert p["initial_repulsion_force_norm_N"] == 0
            assert p["initial_flag_bend_error_max_deg"] < 60
            assert p["initial_flag_torsion_error_max_deg"] < 120
            assert p["initial_elastic_force_norm_N"] > 0
    for label in ("1tau", "2s"):
        job = load_parallel_job(
            ROOT
            / f"conf/phase2_parallel/issue263_initial_bundle/dense_count13_contact_{label}_job.yaml"
        )
        assert len(job.condition_ids) == 26
        assert resolve_execution(job, None).max_workers == 8
        assert resolve_execution(job, None).worker_policy == "cs10_qualified"


def test_long_contract_changes_duration_only_and_retains_review_boundary():
    short, a, _ = configurations()
    long, b, _ = configurations(LONG)
    for x, y in zip(a, b):
        assert x["condition_id"] == y["condition_id"]
        left, right = (
            copy.deepcopy(x["config_overrides"]),
            copy.deepcopy(y["config_overrides"]),
        )
        assert left["time"].pop("duration") == {"value": 1.0, "unit": "tau"}
        assert right["time"].pop("duration") == {"value": 2.0, "unit": "s"}
        assert left == right
    assert (
        short["development_evaluation"]["long_duration"]["status"]
        == "requires_short_screen_review"
    )
    assert long["development_evaluation"]["long_duration"]["expected_steps"] == 500000


@pytest.mark.parametrize(
    "change,message",
    [
        ("bond", "bond lengths"),
        ("hook", "bond lengths"),
        ("phase", "common shape and phase"),
        ("body", "body geometry"),
    ],
)
def test_imported_geometry_cannot_silently_change_physics_or_phase(
    tmp_path, change, message
):
    _, _, configs = configurations()
    cfg = next(
        c
        for c in configs
        if c.flagella.n_flagella == 3 and c.flagella.attachment_slots == (0, 2, 4)
    )
    data = json.loads(STATE.read_text())
    record = next(r for r in data["conditions"] if r["attachment_slots"] == [0, 2, 4])
    points = np.array(record["flagella_positions_over_b"])
    if change == "bond":
        points[0, 3, 0] += 0.02
    elif change == "hook":
        points[0, :, 0] -= 0.02
    elif change == "phase":
        angle = 0.05
        rot = np.array(
            [
                [1, 0, 0],
                [0, np.cos(angle), -np.sin(angle)],
                [0, np.sin(angle), np.cos(angle)],
            ]
        )
        points[0] = points[0, 0] + (points[0] - points[0, 0]) @ rot.T
    else:
        data["body_positions_over_b"][0][0] += 0.1
    record["flagella_positions_over_b"] = points.tolist()
    path = tmp_path / "state.json"
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match=message):
        ModelBuilder(
            cfg.with_overrides(
                {
                    "flagella": {
                        "initial_bundle_state_path": str(path),
                        "initial_bundle_state_sha256": hashlib.sha256(
                            path.read_bytes()
                        ).hexdigest(),
                    }
                }
            )
        ).build()


@pytest.mark.parametrize("scale", [1e-6, 1.0, 1e6])
def test_true_segment_clearance_detects_interior_crossing_at_any_scale(scale):
    p = np.array([[-1.0, 0, 0]]) * scale
    q = -p
    r = np.array([[0.0, -1, 0]]) * scale
    s = -r
    assert segment_distance_squared(p, q, r, s)[0] == pytest.approx(0, abs=1e-20)


@pytest.mark.parametrize("enabled", [True, False])
def test_dense_nf6_initial_transient_can_take_two_internal_steps(enabled):
    _, _, configs = configurations()
    cfg = next(
        c
        for c in configs
        if c.flagella.n_flagella == 6
        and c.potentials.spring_spring_repulsion.body_flagella_enabled == enabled
    )
    model = ModelBuilder(cfg).build()
    engine = DynamicsEngine(model, cfg)
    for _ in range(2):
        diagnostics = engine.step(cfg.dt_star)
        assert np.isfinite(diagnostics.total_forces).all()
        assert np.isfinite(model.positions_m).all()


def test_precomputed_initial_stress_is_not_labelled_hook_neutral():
    _, _, configs = configurations()
    cfg = configs[0].with_overrides({"flagella": {"initial_hook_force_neutral": True}})
    with pytest.raises(ValueError, match="preprocessing"):
        ModelBuilder(cfg).build()


def test_imported_phase_and_asset_hash_are_retained_in_runtime_provenance():
    from sim_swim.sim.core import Simulator

    _, _, configs = configurations()
    simulator = Simulator(configs[-1])
    np.testing.assert_allclose(simulator.model.flagella_initial_phases_rad, np.pi / 3)
    summary = simulator.initial_geometry_summary
    assert (
        summary["initial_geometry_state"]["state_sha256"]
        == configs[-1].flagella.initial_bundle_state_sha256
    )
    assert summary["flagella"]["initial_phase_source"] == "imported common distal phase"
    assert all(
        "torsion_angle" in r["initial_geometry_failures"]
        or "bend_angle" in r["initial_geometry_failures"]
        for r in summary["per_flagellum"]
    )
    bad = configs[-1].with_overrides(
        {"flagella": {"initial_bundle_state_sha256": "0" * 64}}
    )
    with pytest.raises(ValueError, match="SHA-256"):
        ModelBuilder(bad).build()
