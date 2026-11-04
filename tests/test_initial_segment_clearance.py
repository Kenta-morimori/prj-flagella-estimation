from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest
import yaml

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    geometry_preflight,
    normalize_campaign_config,
)
from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.model.builder import ModelBuilder
from sim_swim.model.segment_contacts import (
    initial_body_cylinder_clearance,
    segment_clearance_metrics,
    select_repulsion_segment_pairs,
)
from sim_swim.sim.params import SimulationConfig


def matrix():
    campaign = normalize_campaign_config(
        yaml.safe_load(
            Path("conf/phase2_multi_run/2010_hex_nf3_dj_1tau_issue263.yaml").read_text()
        )
    )
    conditions = build_campaign_conditions(campaign)
    base = SimulationConfig.from_dict(
        yaml.safe_load(Path(campaign["base_config"]).read_text())
    )
    return campaign, conditions, base


def test_root_on_nonadjacent_body_edge_is_rejected_despite_bead_clearance():
    campaign, conditions, base = matrix()
    bad = deepcopy(conditions[:1])
    bad[0]["config_overrides"]["flagella"]["initial_bundle_mode"] = "compact_rear"
    cfg = base.with_overrides(bad[0]["config_overrides"])
    m = ModelBuilder(cfg).build()
    root = m.positions_m[30]
    np.testing.assert_allclose(
        root, (m.positions_m[6] + m.positions_m[12]) / 2, atol=1e-21
    )
    assert np.linalg.norm(root - m.positions_m[6]) > 2 * m.bead_radius_m
    assert np.linalg.norm(root - m.positions_m[12]) > 2 * m.bead_radius_m
    with pytest.raises(ValueError, match="Initial active segment overlap"):
        geometry_preflight(campaign, bad)
    legacy = deepcopy(campaign)
    legacy["development_evaluation"].pop("geometry")
    record = geometry_preflight(legacy, bad)[bad[0]["condition_id"]]
    assert record["initial_segment_clearance"]["status"] == "FAIL"
    assert not record["initial_segment_clearance"]["enforced"]


def test_four_conditions_pass_with_identical_coordinates_within_geometry():
    campaign, conditions, base = matrix()
    preflight = geometry_preflight(campaign, conditions)
    models = [
        ModelBuilder(base.with_overrides(c["config_overrides"])).build()
        for c in conditions
    ]
    for pair in ((0, 1), (2, 3)):
        np.testing.assert_array_equal(
            models[pair[0]].positions_m, models[pair[1]].positions_m
        )
    for cid, record in preflight.items():
        contact = record["initial_segment_clearance"]
        assert contact["enforced"] and contact["status"] == "PASS", cid
        assert record["initial_hook_force_norm_N"] < 1e-18
        assert record["initial_repulsion_force_norm_N"] == 0
        assert record["initial_elastic_force_norm_N"] < 1e-18
    m = models[0]
    assert m.initial_geometry_metadata["hook_from_outward_deg"] == 35
    assert (
        m.initial_geometry_metadata["initial_body_cylinder_clearance_lower_bound_um"]
        > 0.2
    )
    assert not m.root_bending_rest_rad and not m.root_torsion_rest_rad
    legacy = ModelBuilder(
        base.with_overrides(conditions[0]["config_overrides"]).with_overrides(
            {"flagella": {"initial_bundle_mode": "compact_rear"}}
        )
    ).build()
    for ids in m.flagella_indices:
        np.testing.assert_allclose(
            np.diff(m.positions_m[ids], axis=0),
            np.diff(legacy.positions_m[ids], axis=0),
            atol=2e-21,
        )


@pytest.mark.parametrize("body_enabled", [False, True])
def test_runtime_and_preflight_share_exact_pair_selection(body_enabled):
    campaign, conditions, base = matrix()
    cfg = base.with_overrides(conditions[0]["config_overrides"]).with_overrides(
        {
            "potentials": {
                "spring_spring_repulsion": {"body_flagella_enabled": body_enabled}
            }
        }
    )
    m = ModelBuilder(cfg).build()
    e = DynamicsEngine(m, cfg)
    mask = np.all(m.bead_is_body[m.spring_pairs], axis=1)
    selected, counts = select_repulsion_segment_pairs(
        m.segment_pair_indices, mask, body_enabled
    )
    # Independent formula protects the runtime refactor from changing interactions.
    pairs = m.segment_pair_indices
    i, j = mask[pairs[:, 0]], mask[pairs[:, 1]]
    expected = pairs[~(i & j) & (body_enabled | ~(i ^ j))]
    np.testing.assert_array_equal(selected, expected)
    np.testing.assert_array_equal(selected, e.segment_pair_indices_for_repulsion)
    assert counts == e.segment_repulsion_pair_counts
    assert (
        counts["body_flagella_active"] > 0
        if body_enabled
        else counts["body_flagella_active"] == 0
    )


@pytest.mark.parametrize("scale", [1, 1e-6, 1e6])
def test_segment_interior_overlap_is_detected_in_any_length_unit(scale):
    p = np.array([[-1, 0, 0], [1, 0, 0], [0, 0, 0], [0, 1, 0]], float) * scale
    springs = np.array([[0, 1], [2, 3]])
    result = segment_clearance_metrics(
        p, springs, np.array([[0, 1]]), np.array([True, False]), 0.2 * scale
    )
    g = result["groups"]["body_flagella"]
    assert result["status"] == "FAIL" and g["overlap_count"] == 1
    assert g["min_distance_m"] == 0
    assert g["witness"]["bead_indices"] == [[0, 1], [2, 3]]


def test_cylinder_bound_inspects_segment_interiors_and_is_covariant():
    body = np.array([[0, 0.5, 0], [1, 0.5, 0], [0, -0.5, 0], [1, -0.5, 0]])
    # Endpoints are outside, but the filament crosses the entire body interior.
    chain = np.array([[0.5, 1, 0], [0.5, -1, 0]])
    assert initial_body_cylinder_clearance(body, [chain], np.array([1, 0, 0])) < 0
    clear = np.array([[0.5, 0.75, 0], [-0.5, 0.75, 0]])
    a = initial_body_cylinder_clearance(body, [clear], np.array([1, 0, 0]))
    Q = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    b = initial_body_cylinder_clearance(
        body @ Q.T + 3, [clear @ Q.T + 3], Q @ np.array([1, 0, 0])
    )
    assert a > 0.2 and b == pytest.approx(a)


def test_strict_preflight_also_rejects_interflag_crossing(monkeypatch):
    campaign, conditions, base = matrix()
    m = ModelBuilder(base.with_overrides(conditions[2]["config_overrides"])).build()
    i, j = m.flagella_indices[:2]
    m.positions_m[j[4]] = m.positions_m[i[4]]
    monkeypatch.setattr(ModelBuilder, "build", lambda self: m)
    with pytest.raises(ValueError, match="Initial active segment overlap"):
        geometry_preflight(campaign, conditions[2:3])
