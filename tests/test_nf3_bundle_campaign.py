import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from sim_swim.analysis.bundle_rotation import analyze
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]


def configurations():
    path = ROOT / "conf/phase2_multi_run/2010_hex_project_nf3_bundle_1tau_issue255.yaml"
    campaign = normalize_campaign_config(load_yaml(path))
    base = SimulationConfig.from_dict(load_yaml(ROOT / campaign["base_config"]))
    records = build_campaign_conditions(campaign)
    return (
        campaign,
        records,
        [base.with_overrides(r["config_overrides"]) for r in records],
    )


def test_four_conditions_change_only_two_modes_and_preserve_initial_geometry():
    _, _, configs = configurations()
    assert len(configs) == 4
    comparable = []
    positions = []
    for cfg in configs:
        assert cfg.flagella.n_flagella == 3
        assert cfg.flagella.attachment_slots == (0, 2, 4)
        assert cfg.potentials.torsion.force_method == "analytic"
        assert cfg.output.rotation_diagnostics
        sim = Simulator(cfg)
        positions.append(sim.model.positions_m.copy())
        d = asdict(cfg)
        d["motor"].pop("local_attach_frame_position_scale")
        d["motor"].pop("segment_torque_correction")
        comparable.append(d)
    assert all(d == comparable[0] for d in comparable)
    for p in positions[1:]:
        np.testing.assert_array_equal(p, positions[0])
    assert {
        (c.motor.local_attach_frame_position_scale, c.motor.segment_torque_correction)
        for c in configs
    } == {(f, t) for f in [1.0, 1.25] for t in ["none", "minimum_norm"]}


def test_internal_step_recorder_and_reusable_analysis(tmp_path):
    campaign, records, configs = configurations()
    # Bounded integration is an output/diagnostic test, not a physical screen.
    record = records[0]
    cfg = configs[0].with_overrides({"output": {"archive_interval_s": 0.00002}})
    sim = Simulator(cfg)
    run = tmp_path / "run"
    condition = run / "conditions" / record["condition_id"]
    condition.mkdir(parents=True)
    states = sim.run(
        0.00008,
        step_summary_dir=condition,
        output_policy="compact",
        archive_interval_s=0.00002,
    )
    assert len(states) >= 5
    rows = list(csv.DictReader((condition / "rotation_diagnostics.csv").open()))
    assert len(rows) >= 5 and float(rows[-1]["t_s"]) >= 0.00008 - 1e-12
    assert all(r["flag0_valid"] == "1" for r in rows)
    np.savez(
        condition / "state_archive.npz",
        t=np.array([s.t for s in states]),
        bead_positions_um=np.array([s.bead_positions_um for s in states]),
    )
    record["config_overrides"]["output"]["archive_interval_s"] = 0.00002
    (run / "run_manifest.json").write_text(
        json.dumps({"base_config": campaign["base_config"], "conditions": [record]})
    )
    evaluation = tmp_path / "evaluation"
    evaluation.mkdir()
    (evaluation / "summary.csv").write_text(
        "condition_id,screen_status\n" + record["condition_id"] + ",pass\n"
    )
    analyze(run, evaluation, tmp_path / "analysis")
    result = json.loads((tmp_path / "analysis/manifest.json").read_text())
    assert result["conditions"][0]["rotation_valid_fraction"] == 1.0
    assert result["conditions"][0]["state_count"] == len(states)
    assert result["acceptance_gate"] is False
    assert (tmp_path / "analysis/bundle_rotation.png").is_file()


def test_four_conditions_internal_torque_balance_on_deformed_geometry():
    from sim_swim.dynamics.engine import DynamicsEngine
    from sim_swim.model.builder import ModelBuilder

    _, _, configs = configurations()
    for cfg in configs:
        model = ModelBuilder(cfg).build()
        model.positions_m[model.flagella_indices[0][3]] += (
            np.array([0.004, -0.003, 0.002]) * cfg.b_m
        )
        engine = DynamicsEngine(model, cfg)
        diagnostics = engine.step(cfg.dt_star)
        metrics = diagnostics.force_balance_diagnostics
        assert metrics["component_torsion_torque_norm_Nm"] < 1e-30
        assert metrics["component_total_torque_norm_Nm"] < 1e-30
        assert metrics["motor_reaction_solver_success"] == 1
        assert metrics["motor_reaction_fallback_used"] == 0


def test_candidate_C_extension_changes_duration_only():
    from dataclasses import asdict

    _, _, configs = configurations()
    campaign = normalize_campaign_config(
        load_yaml(
            ROOT
            / "conf/phase2_multi_run/2010_hex_project_nf3_bundle_C_2s_issue255.yaml"
        )
    )
    conditions = build_campaign_conditions(campaign)
    assert len(conditions) == 1
    assert conditions[0]["condition_id"] == "nf03__slots024__C"
    cfg = SimulationConfig.from_dict(
        load_yaml(ROOT / campaign["base_config"])
    ).with_overrides(conditions[0]["config_overrides"])
    assert np.isclose(cfg.time.duration_s, 2.0, atol=1e-14, rtol=0)
    assert cfg.total_steps == 500000
    a, b = asdict(configs[2]), asdict(cfg)
    a.pop("time")
    b.pop("time")
    assert a == b
