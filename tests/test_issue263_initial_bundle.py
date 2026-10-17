from dataclasses import asdict
from pathlib import Path

import numpy as np

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    geometry_preflight,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = (
    ROOT / "conf/phase2_multi_run/2010_hex_nf3_bundle_initial_screen_issue263.yaml"
)


def test_initial_bundle_screen_preserves_c_and_compacts_with_clear_hooks():
    campaign = normalize_campaign_config(load_yaml(CAMPAIGN))
    conditions = build_campaign_conditions(campaign)
    assert [row["condition_id"] for row in conditions] == [
        "nf03__slots024__C",
        "nf03__slots024__bundle",
        "nf03__slots024__bundle_phase5",
    ]
    base = SimulationConfig.from_dict(load_yaml(ROOT / campaign["base_config"]))
    configs = [base.with_overrides(row["config_overrides"]) for row in conditions]
    original = normalize_campaign_config(
        load_yaml(
            ROOT
            / "conf/phase2_multi_run/2010_hex_project_nf3_bundle_1tau_issue255.yaml"
        )
    )
    c_row = next(
        row
        for row in build_campaign_conditions(original)
        if row["condition_id"] == "nf03__slots024__C"
    )
    old_c = base.with_overrides(c_row["config_overrides"])
    old_c_data = asdict(old_c)
    new_c_data = asdict(configs[0])
    for data in (old_c_data, new_c_data):
        data["flagella"].pop("initial_bundle_mode")
        data["flagella"].pop("initial_phase_offsets_deg")
    assert new_c_data == old_c_data
    np.testing.assert_array_equal(
        ModelBuilder(configs[0]).build().positions_m,
        ModelBuilder(old_c).build().positions_m,
    )
    comparable = []
    for cfg in configs:
        data = asdict(cfg)
        data["flagella"].pop("initial_bundle_mode")
        data["flagella"].pop("initial_phase_offsets_deg")
        comparable.append(data)
    assert comparable[0] == comparable[1] == comparable[2]
    assert configs[2].flagella.initial_phase_offsets_deg == (0.0, 5.0, 0.0)

    preflight = geometry_preflight(campaign, conditions)
    control = preflight[conditions[0]["condition_id"]]
    for row in conditions[1:]:
        record = preflight[row["condition_id"]]
        assert (
            record["initial_bundle_radius_mean_um"]
            < control["initial_bundle_radius_mean_um"]
        )
        assert record["initial_min_nonattached_bead_distance_m"] >= 2e-7 - 1e-15
        assert record["initial_hook_len_max_abs_error_m"] < 1e-15
        assert record["initial_hook_force_norm_N"] < 1e-18
        assert record["initial_min_outward_projection_m2"] >= -1e-20
        np.testing.assert_allclose(record["initial_hook_from_outward_deg"], 90.0)

    nominal = ModelBuilder(configs[1]).build()
    perturbed = ModelBuilder(configs[2]).build()
    for flag_id in (0, 2):
        ids = nominal.flagella_indices[flag_id]
        np.testing.assert_allclose(
            nominal.positions_m[ids], perturbed.positions_m[ids], atol=1e-18
        )
    assert not np.allclose(
        nominal.positions_m[nominal.flagella_indices[1]],
        perturbed.positions_m[perturbed.flagella_indices[1]],
        atol=1e-18,
    )


def test_compact_initial_geometry_can_integrate_two_steps():
    campaign = normalize_campaign_config(load_yaml(CAMPAIGN))
    condition = build_campaign_conditions(campaign)[1]
    base = SimulationConfig.from_dict(load_yaml(ROOT / campaign["base_config"]))
    simulator = Simulator(base.with_overrides(condition["config_overrides"]))
    states = simulator.run(8e-6, output_policy="compact", archive_interval_s=4e-6)
    assert len(states) == 3
    assert np.isfinite(states[-1].bead_positions_um).all()
