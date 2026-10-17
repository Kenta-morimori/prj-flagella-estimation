import csv
import json
from copy import deepcopy
from dataclasses import asdict
from pathlib import Path

import numpy as np

from sim_swim.analysis.model_development_evaluation import (
    _render_replays,
    _row,
    build_evaluation,
)
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    geometry_preflight,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.analysis.parallel_job import load_parallel_job, resolve_execution
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = (
    ROOT / "conf/phase2_multi_run/2010_hex_nf3_bundle_initial_screen_issue263.yaml"
)
COUNT_SCREEN = ROOT / "conf/phase2_multi_run/2010_hex_nf3to6_bundle_1tau_issue263.yaml"
LONG_CAMPAIGN = ROOT / "conf/phase2_multi_run/2010_hex_nf3to6_bundle_2s_issue263.yaml"
COUNT_JOB = ROOT / "conf/phase2_parallel/issue263_initial_bundle/count_screen_job.yaml"
LONG_JOB = ROOT / "conf/phase2_parallel/issue263_initial_bundle/long_2s_job.yaml"


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


def test_common_evaluator_accepts_fixed_count_and_torque(tmp_path: Path):
    campaign = normalize_campaign_config(load_yaml(CAMPAIGN))
    record = build_campaign_conditions(campaign)[0]
    (tmp_path / "run_summary.json").write_text(
        json.dumps({"gates": {}, "all_step_metrics": {}}), encoding="utf-8"
    )
    (tmp_path / "performance.json").write_text(
        json.dumps({"wall_time_s": 1.0, "steps_per_s": 1000.0}),
        encoding="utf-8",
    )
    record["output_dir"] = str(tmp_path)
    record["time"] = {"dt_star": 1e-4, "dt_internal_s": 4e-6}
    row = _row(record, {"completed": "True"}, qc={}, stage="short_screen")
    assert row["n_flagella"] == 3
    assert json.loads(row["attachment_slots"]) == [0, 2, 4]
    assert row["torque_Nm"] == 2.5e-20
    (tmp_path / "state_archive.npz").write_bytes(b"test archive")
    record["time"]["duration_s"] = 2.0
    long_row = _row(record, {"completed": "True"}, qc={}, stage="long_duration")
    assert long_row["duration_s"] == 2.0
    assert long_row["full_ring_rotation_equivalent"] is False


def test_two_second_checkpoint_changes_only_duration():
    short_config = normalize_campaign_config(load_yaml(COUNT_SCREEN))
    short = build_campaign_conditions(short_config)
    long_config = normalize_campaign_config(load_yaml(LONG_CAMPAIGN))
    longer = build_campaign_conditions(long_config)
    assert [row["condition_id"] for row in short] == [
        "nf03__slots024__bundle",
        "nf04__slots0134__bundle",
        "nf05__slots01234__bundle",
        "nf06__slots012345__bundle",
    ]
    assert [row["condition_id"] for row in short] == [
        row["condition_id"] for row in longer
    ]
    for before, after in zip(short, longer, strict=True):
        old = deepcopy(before["config_overrides"])
        new = deepcopy(after["config_overrides"])
        assert old["time"]["duration"] == {"value": 1.0, "unit": "tau"}
        assert new["time"]["duration"] == {"value": 2.0, "unit": "s"}
        old["time"].pop("duration")
        new["time"].pop("duration")
        assert old == new
    assert long_config["development_evaluation"]["stage"] == "long_duration"
    assert (
        long_config["development_evaluation"]["long_duration"]["expected_steps"]
        == 500000
    )
    for path in (COUNT_JOB, LONG_JOB):
        job = load_parallel_job(path)
        execution = resolve_execution(job, None)
        assert job.condition_ids == tuple(row["condition_id"] for row in longer)
        assert execution.worker_policy == "cs10_qualified"
        assert execution.max_workers == 4
    preflight = geometry_preflight(short_config, short)
    for row in short:
        record = preflight[row["condition_id"]]
        assert record["initial_min_nonattached_bead_distance_m"] >= 2e-7 - 1e-15
        assert record["initial_hook_len_max_abs_error_m"] < 1e-15
        assert record["initial_hook_force_norm_N"] < 1e-18
        assert (
            row["axis_values"]["phase_offsets_deg"]
            == [0] * row["axis_values"]["n_flagella"]
        )


def test_long_evaluation_groups_fixed_geometry_conditions(monkeypatch, tmp_path: Path):
    import sim_swim.analysis.phase2_replay as replay

    calls = []
    monkeypatch.setattr(replay, "main", lambda args: calls.append(args))
    config = normalize_campaign_config(load_yaml(LONG_CAMPAIGN))
    rows = [
        {
            "condition_id": row["condition_id"],
            "n_flagella": row["axis_values"]["n_flagella"],
        }
        for row in build_campaign_conditions(config)
    ]
    _render_replays(
        rows,
        replay_input=tmp_path / "input",
        output_dir=tmp_path / "output",
        stage="long_duration",
    )
    assert len(calls) == 4
    for n_flagella, args in enumerate(calls, 3):
        assert args.count("--condition-id") == 1
        assert args[args.index("--camera-3d") + 1] == "fixed"
        assert args[args.index("--camera-2d") + 1] == "fixed"
        assert str(tmp_path / "output" / "replay" / f"nf{n_flagella:02d}") in args


def test_two_second_common_evaluation_accepts_four_portable_sources(tmp_path: Path):
    campaign = normalize_campaign_config(load_yaml(LONG_CAMPAIGN))
    source = tmp_path / "source"
    source.mkdir()
    records = []
    summary_rows = []
    metrics = {
        name: {"max": value}
        for name, value in {
            "local_attach_first_rel_err": 0.01,
            "hook_len_rel_err_max": 0.01,
            "flag_bond_rel_err_max": 0.01,
            "flag_bend_err_max_deg": 1.0,
            "flag_torsion_err_max_deg": 1.0,
            "motor_force_balance_residual_ratio": 1e-12,
            "motor_torque_balance_residual_ratio": 1e-12,
        }.items()
    }
    for condition in build_campaign_conditions(campaign):
        cid = condition["condition_id"]
        directory = source / "conditions" / cid
        directory.mkdir(parents=True)
        (directory / "run_summary.json").write_text(
            json.dumps(
                {
                    "execution": {"status": "completed"},
                    "gates": {
                        "finite": {"any_fail": False},
                        "shape_body": {"any_fail": False},
                        "shape_nonbody": {"any_fail": False},
                    },
                    "all_step_metrics": metrics,
                }
            )
        )
        (directory / "performance.json").write_text(
            json.dumps({"wall_time_s": 1.0, "steps_per_s": 500.0})
        )
        (directory / "state_archive.npz").write_bytes(b"portable archive")
        records.append(
            {
                **condition,
                "output_dir": f"/cs10/original/{cid}",
                "time": {"duration_s": 2.0},
            }
        )
        summary_rows.append({"condition_id": cid, "completed": "True"})
    with (source / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["condition_id", "completed"])
        writer.writeheader()
        writer.writerows(summary_rows)
    (source / "run_manifest.json").write_text(
        json.dumps(
            {
                "base_config": campaign["base_config"],
                "campaign_config": str(LONG_CAMPAIGN.relative_to(ROOT)),
                "model_profile": load_yaml(ROOT / campaign["base_config"])[
                    "model_profile"
                ],
                "git": {"commit": "test", "is_clean": True},
                "conditions": records,
            }
        )
    )
    output = build_evaluation(
        config_path=LONG_CAMPAIGN,
        run_dirs=[source],
        output_dir=tmp_path / "evaluation",
    )
    result = json.loads(output["manifest"].read_text())
    assert result["status_counts"] == {"pass": 4, "fail": 0}
    assert output["window_qc_csv"].is_file()
