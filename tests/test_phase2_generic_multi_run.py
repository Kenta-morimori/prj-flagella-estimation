from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.cli_profiles import list_profile_entries, load_profile_entry
from sim_swim.analysis.multi_run_campaign import (
    apply_campaign_cli_overrides,
    build_campaign_conditions,
    load_yaml,
)
from sim_swim.analysis.online_run_summary import OnlineRunSummary
from sim_swim.analysis.sweeps.generic_multi_run import (
    _condition_row,
    _summary_fieldnames,
    run_campaign,
)


def _load_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_generic_profile(
    path: Path, base_dir: Path, *, timestamp_subdir: bool
) -> None:
    path.write_text(
        "\n".join(
            [
                "kind: generic_multi_run",
                "metadata:",
                "  role: sweep",
                "  canonical: true",
                "base_config: conf/sim_swim_2010.yaml",
                "base_overrides: {}",
                "sweep:",
                "  axes:",
                "    torque:",
                "      key: motor.torque_Nm",
                "      short_name: torque",
                "      values: [1.5e-20, 2.0e-20]",
                "replay:",
                "  mode: both",
                "  fps_out_3d: 10.0",
                "  output_subdir: replay",
                "plot:",
                "  default_x_axis: torque",
                "  default_y_axis: null",
                "  metrics:",
                "    - first_fail_t_s",
                "output:",
                f"  base_dir: {base_dir}",
                f"  timestamp_subdir: {str(timestamp_subdir).lower()}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_generic_multi_run_profile_is_listed() -> None:
    entries = list_profile_entries(role="sweep", canonical_only=True)
    paths = {entry["path"] for entry in entries}

    assert "conf/phase2_multi_run/latest_model_torque_shape_stability.yaml" in paths
    assert "conf/phase2_multi_run/flagella_count_failure_boundary_seed00.yaml" in paths
    assert "conf/phase2_multi_run/flagella_count_stability_narrow_seed00.yaml" in paths


def test_generic_multi_run_profile_exposes_metadata() -> None:
    entry = load_profile_entry(
        Path("conf/phase2_multi_run/latest_model_torque_shape_stability.yaml")
    )

    assert entry["kind"] == "generic_multi_run"
    assert entry["metadata"]["role"] == "sweep"
    assert entry["metadata"]["recommended_heatmap_profile"] is None


def test_latest_model_profile_base_overrides_only_intentional_diffs() -> None:
    profile = load_yaml(
        Path("conf/phase2_multi_run/latest_model_torque_shape_stability.yaml")
    )

    assert profile["base_overrides"] == {
        "time.duration_s": 1.0,
        "flagella.initial_helix_axis_from_rear_deg": 0,
        "seed.attach_seed": 0,
        "seed.phase_seed": 0,
    }
    assert profile["output"]["timestamp_subdir"] is False


def test_flagella_count_stability_candidates_expand_stiffness_axes() -> None:
    campaign = apply_campaign_cli_overrides(
        load_yaml(
            Path(
                "conf/phase2_multi_run/flagella_count_stability_candidates_seed00.yaml"
            )
        ),
        [],
    )

    conditions = build_campaign_conditions(campaign)

    assert len(conditions) == 18
    first = conditions[0]
    assert first["config_overrides"]["flagella"]["n_flagella"] == 4
    assert first["config_overrides"]["stiffness_scales"][
        "flag_spring"
    ] == pytest.approx(1.0)
    assert first["config_overrides"]["stiffness_scales"]["body"] == pytest.approx(1.0)
    assert any(
        condition["config_overrides"]["stiffness_scales"]["flag_spring"]
        == pytest.approx(2.0)
        and condition["config_overrides"]["stiffness_scales"]["body"]
        == pytest.approx(2.0)
        for condition in conditions
    )


def test_generic_multi_run_include_condition_ids_filters_in_requested_order() -> None:
    campaign = apply_campaign_cli_overrides(
        {
            "base_config": "conf/sim_swim_2010.yaml",
            "sweep": {
                "include_condition_ids": ["b2__a1", "b1__a0"],
                "axes": {
                    "b": {
                        "key": "seed.attach_seed",
                        "short_name": "b",
                        "values": [1, 2],
                        "ids": ["b1", "b2"],
                    },
                    "a": {
                        "key": "seed.phase_seed",
                        "short_name": "a",
                        "values": [0, 1],
                        "ids": ["a0", "a1"],
                    },
                },
            },
        },
        [],
    )

    conditions = build_campaign_conditions(campaign)

    assert [condition["condition_id"] for condition in conditions] == [
        "b2__a1",
        "b1__a0",
    ]
    assert conditions[0]["condition_index"] == 3
    assert conditions[0]["config_overrides"]["seed"]["attach_seed"] == 2
    assert conditions[0]["config_overrides"]["seed"]["phase_seed"] == 1


def test_generic_multi_run_include_condition_ids_rejects_unknown_id() -> None:
    campaign = apply_campaign_cli_overrides(
        {
            "base_config": "conf/sim_swim_2010.yaml",
            "sweep": {
                "include_condition_ids": ["missing"],
                "axes": {
                    "seed": {
                        "key": "seed.attach_seed",
                        "short_name": "seed",
                        "values": [0],
                        "ids": ["seed0"],
                    }
                },
            },
        },
        [],
    )

    with pytest.raises(ValueError, match="unknown condition IDs"):
        build_campaign_conditions(campaign)


def test_flagella_count_stability_narrow_seed00_profile_contract() -> None:
    campaign = apply_campaign_cli_overrides(
        load_yaml(
            Path("conf/phase2_multi_run/flagella_count_stability_narrow_seed00.yaml")
        ),
        [],
    )

    conditions = build_campaign_conditions(campaign)

    assert len(conditions) == 36
    assert campaign["base_overrides"]["seed"]["attach_seed"] == 0
    assert campaign["base_overrides"]["seed"]["phase_seed"] == 0
    assert campaign["base_overrides"]["motor"]["local_first_second_spring_scale"] == (
        pytest.approx(1.0)
    )
    assert campaign["plot"]["metrics"] == [
        "first_fail_t_s",
        "max_flag_bond_rel_err",
        "flag_bond_rel_err_max",
        "local_first_second_rel_err",
        "body_spring_max_stretch_ratio",
        "body_centerline_max_deviation_um",
    ]
    assert any(
        condition["config_overrides"]["flagella"]["n_flagella"] == 6
        and condition["config_overrides"]["stiffness_scales"]["flag_spring"]
        == pytest.approx(2.25)
        and condition["config_overrides"]["stiffness_scales"]["body"]
        == pytest.approx(2.5)
        for condition in conditions
    )


def test_flagella_count_stability_smoke_seed00_accepts_candidate_overrides() -> None:
    campaign = apply_campaign_cli_overrides(
        load_yaml(
            Path("conf/phase2_multi_run/flagella_count_stability_smoke_seed00.yaml")
        ),
        ["stiffness_scales.flag_spring=2.25", "stiffness_scales.body=1.5"],
    )

    conditions = build_campaign_conditions(campaign)

    assert [condition["condition_id"] for condition in conditions] == [
        "nf01",
        "nf02",
        "nf03",
    ]
    assert all(
        condition["config_overrides"]["stiffness_scales"]["flag_spring"]
        == pytest.approx(2.25)
        for condition in conditions
    )
    assert all(
        condition["config_overrides"]["stiffness_scales"]["body"] == pytest.approx(1.5)
        for condition in conditions
    )


def test_generic_multi_run_builds_conditions_and_cli_override() -> None:
    campaign = apply_campaign_cli_overrides(
        load_yaml(
            Path("conf/phase2_multi_run/latest_model_torque_shape_stability.yaml")
        ),
        ["sweep.axes.torque.values=[1.0e-20,2.0e-20]"],
    )

    conditions = build_campaign_conditions(campaign)

    assert [condition["condition_id"] for condition in conditions] == [
        "torque_1e-20",
        "torque_2e-20",
    ]
    assert conditions[0]["axis_values"]["torque"] == 1.0e-20
    assert conditions[1]["config_overrides"]["motor"]["torque_Nm"] == 2.0e-20


def test_generic_multi_run_manifests_record_model_profile(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign_path = tmp_path / "campaign.yaml"
    output_dir = tmp_path / "output"
    campaign_path.write_text(
        "\n".join(
            [
                "kind: generic_multi_run",
                "base_config: conf/sim_swim_2010.yaml",
                "base_overrides:",
                "  time.duration_s: 0.0001",
                "  motor.force_distribution: hook_coupled_body_reaction",
                "sweep:",
                "  axes:",
                "    torque:",
                "      key: motor.torque_Nm",
                "      values: [2.0e-20]",
                "output:",
                f"  base_dir: {output_dir}",
                "  timestamp_subdir: false",
                "  save_state_archive: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr("sim_swim.core.run_context._require_clean_git", lambda: None)

    run_campaign(["--campaign-config", str(campaign_path), "--sample-limit", "1"])

    root_manifest = json.loads(
        (output_dir / "manifest.json").read_text(encoding="utf-8")
    )
    run_manifest = json.loads(
        (output_dir / "run_manifest.json").read_text(encoding="utf-8")
    )
    for manifest in (root_manifest, run_manifest):
        assert manifest["source_config_path"] == "conf/sim_swim_2010.yaml"
        assert manifest["model_profile"] == {
            "year": 2010,
            "variant": "project",
            "resolution": "legacy_project",
            "implementation_status": "supported",
            "body_beads": 15,
            "flagellum_beads_per_filament": 11,
            "nominal_flagella_count": 3,
            "nominal_total_beads": 48,
        }
        assert manifest["time"]["duration_s"] == pytest.approx(0.1)
        assert manifest["time"]["dt_star"] == pytest.approx(1.0e-4)
        assert manifest["time"]["time_schema_source"] == "canonical"
    assert run_manifest["conditions"][0]["time"]["duration_s"] == pytest.approx(0.0001)
    assert run_manifest["conditions"][0]["time"]["time_schema_source"] == (
        "mixed_equivalent"
    )
    condition = run_manifest["conditions"][0]
    assert condition["dynamics"]["provenance"] == "paper_inspired_approximation"
    assert condition["dynamics"]["reaction_support_bead_counts"]
    assert condition["dynamics"]["reaction_fallback_used"] is False
    assert condition["geometry"]["actual"]["total_beads"] == 48
    assert condition["geometry"]["actual"]["body_diagonal_edges"] == 24


def test_issue113_seed_fixed_profile_builds_three_boundary_conditions() -> None:
    campaign = apply_campaign_cli_overrides(
        load_yaml(
            Path("conf/phase2_multi_run/flagella_count_failure_boundary_seed00.yaml")
        ),
        None,
    )

    conditions = build_campaign_conditions(campaign)

    assert [condition["condition_id"] for condition in conditions] == [
        "nf04",
        "nf05",
        "nf06",
    ]
    assert [condition["axis_values"]["n_flagella"] for condition in conditions] == [
        4,
        5,
        6,
    ]
    assert conditions[0]["config_overrides"]["seed"]["attach_seed"] == 0
    assert conditions[0]["config_overrides"]["seed"]["phase_seed"] == 0
    assert conditions[0]["config_overrides"]["motor"]["torque_Nm"] == 2.0e-20


def test_run_multi_run_wrapper_lists_generic_multi_run_kind(capsys) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/run_multi_run.py"),
        "phase2_run_multi_run_wrapper_generic_multi",
    )

    module.main(
        [
            "config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml",
            "list_kind=true",
        ]
    )

    assert capsys.readouterr().out.strip() == "generic_multi_run"


def test_run_sweep_wrapper_rejects_generic_multi_run_profile() -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/run_sweep.py"),
        "phase2_run_sweep_wrapper_reject_generic_multi",
    )

    try:
        module.main(
            ["config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml"]
        )
    except SystemExit as exc:
        assert "use run_multi_run.py" in str(exc)
    else:
        raise AssertionError("expected SystemExit for generic multi-run profile")


def test_generic_multi_run_plot_outputs_line_plots(tmp_path: Path) -> None:
    summary_csv = tmp_path / "summary.csv"
    summary_csv.write_text(
        "\n".join(
            [
                "condition_id,condition_label,axis_torque_label,axis_torque_index,first_fail_t_s,max_flag_bond_rel_err,hook_len_rel_err_max,axis_center_to_body_roll_ratio_mean",
                "torque_1p5e20,torque=1.5e-20,1.5e-20,0,0.5,0.3,0.1,120",
                "torque_2p0e20,torque=2.0e-20,2.0e-20,1,0.3,0.8,0.2,90",
                "torque_2p5e20,torque=2.5e-20,2.5e-20,2,0.1,1.5,0.4,60",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    module = _load_script(
        Path("scripts/02_phase2_analysis/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_generic_multi",
    )

    module.main(
        [
            "config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml",
            f"summary_csv={summary_csv}",
            f"output_dir={tmp_path / 'plots'}",
        ]
    )

    assert (tmp_path / "plots" / "plot_data.csv").is_file()
    assert (tmp_path / "plots" / "first_fail_t_s_vs_torque.png").is_file()
    assert (tmp_path / "plots" / "max_flag_bond_rel_err_vs_torque.png").is_file()


def test_generic_multi_run_plot_filters_extra_axes(tmp_path: Path) -> None:
    summary_csv = tmp_path / "summary.csv"
    summary_csv.write_text(
        "\n".join(
            [
                "condition_id,condition_label,axis_attach_seed_label,axis_attach_seed_index,axis_phase_seed_label,axis_phase_seed_value,axis_phase_seed_index,axis_n_flagella_label,axis_n_flagella_index,first_fail_t_s,hook_len_rel_err_max,max_flag_bond_rel_err,body_roll_net_abs_revolutions,axis_center_to_body_roll_ratio_mean",
                "as000__ps000__nf01,as=0 ps=0 nf=1,0,0,0,0,0,1,0,0.5,0.1,0.2,0.0,120",
                "as000__ps001__nf01,as=0 ps=1 nf=1,0,0,1,1,1,1,0,0.4,0.2,0.3,0.0,90",
                "as001__ps000__nf02,as=1 ps=0 nf=2,1,1,0,0,0,2,1,0.3,0.3,0.4,0.0,80",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    module = _load_script(
        Path("scripts/02_phase2_analysis/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_generic_multi_filter_axes",
    )

    module.main(
        [
            "config=conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml",
            f"summary_csv={summary_csv}",
            f"output_dir={tmp_path / 'plots'}",
        ]
    )

    plot_data = (tmp_path / "plots" / "plot_data.csv").read_text(encoding="utf-8")
    assert "as000__ps000__nf01" in plot_data
    assert "as001__ps000__nf02" in plot_data
    assert "as000__ps001__nf01" not in plot_data
    assert (tmp_path / "plots" / "first_fail_t_s_heatmap.png").is_file()


def test_generic_multi_run_plot_labels_follow_summary_axis_overrides() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/heatmaps/generic_multi_run.py"),
        "phase2_generic_multi_run_heatmap_labels",
    )
    axis = {
        "name": "torque",
        "labels": ["1.5e-20", "2e-20", "2.5e-20"],
    }
    rows = [
        {
            "axis_torque_index": "0",
            "axis_torque_label": "1e-20",
        },
        {
            "axis_torque_index": "1",
            "axis_torque_label": "2e-20",
        },
    ]

    assert module._axis_labels(rows, axis) == ["1e-20", "2e-20"]


def test_generic_multi_run_plot_accepts_run_dir(tmp_path: Path) -> None:
    summary_csv = tmp_path / "summary.csv"
    summary_csv.write_text(
        "\n".join(
            [
                "condition_id,condition_label,axis_torque_label,axis_torque_index,first_fail_t_s,max_flag_bond_rel_err,hook_len_rel_err_max,axis_center_to_body_roll_ratio_mean",
                "torque_1p5e20,torque=1.5e-20,1.5e-20,0,0.5,0.3,0.1,120",
                "torque_2p0e20,torque=2.0e-20,2.0e-20,1,0.3,0.8,0.2,90",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    module = _load_script(
        Path("scripts/02_phase2_analysis/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_generic_multi_run_dir",
    )

    module.main(
        [
            "config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml",
            f"run_dir={tmp_path}",
        ]
    )

    assert (tmp_path / "plots" / "plot_data.csv").is_file()
    assert (tmp_path / "plots" / "first_fail_t_s_vs_torque.png").is_file()


def test_generic_multi_run_plot_uses_fixed_output_base_dir(tmp_path: Path) -> None:
    run_dir = tmp_path / "fixed_run"
    run_dir.mkdir()
    profile = tmp_path / "fixed_profile.yaml"
    _write_generic_profile(profile, run_dir, timestamp_subdir=False)
    (run_dir / "summary.csv").write_text(
        "\n".join(
            [
                "condition_id,condition_label,axis_torque_label,axis_torque_index,first_fail_t_s",
                "torque_1p5e20,torque=1.5e-20,1.5e-20,0,0.5",
                "torque_2p0e20,torque=2.0e-20,2.0e-20,1,0.3",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    module = _load_script(
        Path("scripts/02_phase2_analysis/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_generic_multi_fixed_output",
    )

    module.main([f"config={profile}"])

    assert (run_dir / "plots" / "plot_data.csv").is_file()
    assert (run_dir / "plots" / "first_fail_t_s_vs_torque.png").is_file()


def test_generic_multi_run_plot_requires_run_dir_for_timestamped_output(
    tmp_path: Path,
) -> None:
    profile = tmp_path / "timestamped_profile.yaml"
    _write_generic_profile(
        profile,
        tmp_path / "timestamped_parent",
        timestamp_subdir=True,
    )
    module = _load_script(
        Path("scripts/02_phase2_analysis/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_generic_multi_timestamped_output",
    )

    with pytest.raises(ValueError, match="output.timestamp_subdir is true"):
        module.main([f"config={profile}"])


def test_plot_heatmap_lists_generic_multi_run_profiles(capsys) -> None:
    module = _load_script(
        Path("scripts/02_phase2_analysis/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_generic_multi_list",
    )

    module.main(["list_canonical_profiles=true"])

    output = capsys.readouterr().out
    assert "conf/phase2_multi_run/latest_model_torque_shape_stability.yaml" in output


def test_replay_load_inputs_uses_manifest_condition_order_and_output_dir(
    tmp_path: Path,
) -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_generic_multi_inputs",
    )
    input_dir = tmp_path / "replay"
    input_dir.mkdir()
    (input_dir / "summary.csv").write_text(
        "\n".join(
            [
                "condition_id,condition_label",
                "torque_2p0e20,torque=2.0e-20",
                "torque_1p5e20,torque=1.5e-20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    external_root = tmp_path / "campaign_outputs"
    for condition_id in ("torque_1p5e20", "torque_2p0e20"):
        condition_dir = external_root / condition_id
        condition_dir.mkdir(parents=True)
        (condition_dir / "state_archive.npz").write_bytes(b"")
    manifest = {
        "base_config": "conf/sim_swim_2010.yaml",
        "condition_order": ["torque_1p5e20", "torque_2p0e20"],
        "conditions": [
            {
                "condition_id": "torque_1p5e20",
                "output_dir": str(external_root / "torque_1p5e20"),
                "config_overrides": {},
            },
            {
                "condition_id": "torque_2p0e20",
                "output_dir": str(external_root / "torque_2p0e20"),
                "config_overrides": {},
            },
        ],
    }
    (input_dir / "run_manifest.json").write_text(
        json.dumps(manifest),
        encoding="utf-8",
    )

    rows, records, base_cfg_path = module._load_inputs(input_dir)

    assert [row["condition_id"] for row in rows] == [
        "torque_1p5e20",
        "torque_2p0e20",
    ]
    assert records["torque_1p5e20"]["output_dir"] == str(
        external_root / "torque_1p5e20"
    )
    assert base_cfg_path == Path("conf/sim_swim_2010.yaml")


def test_replay_builds_refined_2015_geometry_from_campaign_record() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_refined_2015_geometry",
    )
    cfg = module._build_cfg(
        base_cfg_path=Path("conf/sim_swim_2015_paper.yaml"),
        condition_record={"config_overrides": {}},
        fps_out_3d=25.0,
    )

    from sim_swim.sim.core import Simulator

    simulator = Simulator(cfg)
    geometry = simulator.implementation_manifest()["geometry"]["actual"]
    assert geometry["total_beads"] == 120
    assert geometry["body_beads"] == 30
    assert geometry["flagellum_beads"] == [30, 30, 30]
    assert geometry["body_diagonal_edges"] == 0


def test_replay_applies_dotted_stage_a_condition_overrides() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_dotted_stage_a_overrides",
    )

    cfg = module._build_cfg(
        base_cfg_path=Path("conf/sim_swim_2015.yaml"),
        condition_record={
            "config_overrides": {
                "time.duration": {"value": 1.0, "unit": "tau"},
                "time.integration.dt_star": 1.0e-5,
                "motor.torque_Nm": 6.0e-19,
            }
        },
        fps_out_3d=25.0,
    )

    assert cfg.duration_star == pytest.approx(1.0)
    assert cfg.dt_star == pytest.approx(1.0e-5)
    assert cfg.motor.torque_Nm == pytest.approx(6.0e-19)


def test_replay_resamples_archived_states_for_display_only() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_state_resampling",
    )
    from sim_swim.sim.core import SimulationState

    states = [
        SimulationState(
            t=t,
            position_um=(t, 0.0, 0.0),
            quaternion=(0.0, 0.0, 0.0, 1.0),
            velocity_um_s=(1.0, 0.0, 0.0),
            omega_rad_s=(0.0, 0.0, 0.0),
            bead_positions_um=np.asarray([[t, 0.0, 0.0]]),
        )
        for t in (0.0, 1.0)
    ]

    replay_states = module._resample_states_for_replay(states, 5)

    assert [state.t for state in replay_states] == pytest.approx(
        [0.0, 0.25, 0.5, 0.75, 1.0]
    )
    assert replay_states[2].position_um == pytest.approx((0.5, 0.0, 0.0))
    assert np.allclose(replay_states[2].bead_positions_um, [[0.5, 0.0, 0.0]])


def test_replay_wrapper_accepts_config_run_dir_defaults(tmp_path: Path) -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_generic_multi_run_dir_args",
    )

    args = module._parse_args(
        [
            "config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml",
            f"run_dir={tmp_path}",
            "overwrite=true",
        ]
    )

    assert args.config == Path(
        "conf/phase2_multi_run/latest_model_torque_shape_stability.yaml"
    )
    assert args.input_dir == tmp_path
    assert args.output_dir == tmp_path / "replay"
    assert args.mode == "both"
    assert args.fps_out_3d == 10.0
    assert args.max_panels_per_grid == 9
    assert args.overwrite is True


def test_replay_wrapper_uses_fixed_output_base_dir(tmp_path: Path) -> None:
    run_dir = tmp_path / "fixed_run"
    profile = tmp_path / "fixed_profile.yaml"
    _write_generic_profile(profile, run_dir, timestamp_subdir=False)
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_generic_multi_fixed_output_args",
    )

    args = module._parse_args([f"config={profile}", "overwrite=true"])

    assert args.input_dir == run_dir
    assert args.output_dir == run_dir / "replay"
    assert args.mode == "both"
    assert args.fps_out_3d == 10.0
    assert args.max_panels_per_grid == 9


def test_replay_wrapper_reads_max_panels_per_grid_from_profile() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_max_panels_from_profile",
    )

    args = module._parse_args(
        ["config=conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml"]
    )

    assert args.max_panels_per_grid == 9


@pytest.mark.parametrize(
    ("n_conditions", "expected_shape"),
    [
        (1, (1, 1)),
        (2, (1, 2)),
        (3, (2, 2)),
        (4, (2, 2)),
        (5, (2, 3)),
        (9, (3, 3)),
        (10, (3, 4)),
        (36, (6, 6)),
    ],
)
def test_replay_auto_grid_shape_is_near_square(
    n_conditions: int,
    expected_shape: tuple[int, int],
) -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        f"phase2_replay_auto_grid_shape_{n_conditions}",
    )

    assert module._auto_grid_shape(n_conditions) == expected_shape


def test_replay_auto_grid_layout_preserves_condition_order() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_auto_grid_layout_order",
    )
    rows = [{"condition_id": f"cond_{index}"} for index in range(5)]

    n_rows, n_cols, positions = module._grid_layout_for_rows(rows)

    assert (n_rows, n_cols) == (2, 3)
    assert positions == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]


def test_replay_grid_layout_keeps_explicit_summary_positions() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_explicit_grid_layout",
    )
    rows = [
        {"condition_id": "a", "grid_row_index": "1", "grid_col_index": "0"},
        {"condition_id": "b", "grid_row_index": "0", "grid_col_index": "2"},
    ]

    n_rows, n_cols, positions = module._grid_layout_for_rows(rows)

    assert (n_rows, n_cols) == (2, 3)
    assert positions == [(1, 0), (0, 2)]


def test_replay_page_index_groups_preserve_condition_order() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_page_index_groups",
    )

    assert module._page_index_groups(36, 9) == [
        list(range(0, 9)),
        list(range(9, 18)),
        list(range(18, 27)),
        list(range(27, 36)),
    ]
    assert module._page_index_groups(3, 9) == [list(range(3))]


def test_replay_expands_flagella_count_condition_labels() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_expanded_condition_labels",
    )

    label = module._label_for_row(
        {"condition_label": "as=2, ps=1, nf=4", "condition_id": "nf04"}
    )

    assert label == "attach_seed=2, phase_seed=1, n_flagella=4"


def test_replay_marks_transient_first_fail_as_failure() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_transient_first_fail",
    )
    row = {
        "final_shape_pass_nonbody": "True",
        "first_fail_t_s": "0.3259",
        "first_fail_category_nonbody": "flag",
    }

    assert module._fail_label(row) == "FAIL flag@0.3259"
    assert module._row_passes_nonbody(row) is False


def test_replay_marks_no_first_fail_as_pass() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_no_first_fail",
    )
    row = {
        "final_shape_pass_nonbody": "True",
        "first_fail_t_s": "",
        "first_fail_category_nonbody": "",
    }

    assert module._fail_label(row) == "PASS"
    assert module._row_passes_nonbody(row) is True


def test_replay_uses_body_gate_for_generic_campaign_status() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_body_gate_status",
    )

    assert (
        module._fail_label(
            {
                "final_shape_pass_nonbody": "",
                "body_shape_pass": "True",
                "body_fail_category": "none",
            }
        )
        == "PASS"
    )
    assert (
        module._fail_label(
            {
                "final_shape_pass_nonbody": "True",
                "body_shape_pass": "False",
                "body_fail_category": "body_spring",
                "body_first_fail_t_s": "0.02725",
            }
        )
        == "FAIL body_spring@0.0272"
    )


def test_replay_marks_compact_campaign_status_pass() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_compact_campaign_pass",
    )

    assert module._fail_label({"status": "pass"}) == "PASS"


def test_replay_status_lines_include_motor_torque_after_time() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_status_motor_torque",
    )
    state = type("State", (), {"t": 0.125, "reverse_flagella": (), "flag_states": ()})()
    cfg = type(
        "Config",
        (),
        {
            "motor_torque_Nm": 2.0e-20,
            "motor": type("Motor", (), {"enable_switching": False})(),
        },
    )()

    assert module._replay_status_lines(state, cfg, "PASS") == [
        "RUN",
        "t = 0.125 s",
        "motor torque / flag = 2.00e-20 N m",
        "PASS",
    ]


def test_replay_status_lines_show_2015_tau_seconds_and_steps() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_2015_time_label",
    )
    state = type(
        "State", (), {"t": 1.0 / 1200.0, "reverse_flagella": (), "flag_states": ()}
    )()
    cfg = type(
        "Config",
        (),
        {
            "tau_s": 1.0 / 1200.0,
            "dt_star": 1.0e-5,
            "motor_torque_Nm": 1.2e-18,
            "model_profile": type("Profile", (), {"year": 2015})(),
            "motor": type("Motor", (), {"enable_switching": False})(),
        },
    )()

    assert module._replay_status_lines(state, cfg, "PASS") == [
        "RUN",
        "t = 1.000 τ (0.000833 s, 100,000 steps)",
        "motor torque / flag = 1.20e-18 N m",
        "PASS",
    ]


def test_replay_status_lines_show_reference_torque_2010_in_tau() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_2010_reference_torque_time_label",
    )
    state = type("State", (), {"t": 0.04, "reverse_flagella": (), "flag_states": ()})()
    cfg = type(
        "Config",
        (),
        {
            "tau_s": 0.04,
            "dt_star": 1.0e-3,
            "time_scale_policy": "reference_torque",
            "motor_torque_Nm": 2.5e-20,
            "model_profile": type("Profile", (), {"year": 2010})(),
            "motor": type("Motor", (), {"enable_switching": False})(),
        },
    )()

    assert module._replay_status_lines(state, cfg, "PASS") == [
        "RUN",
        "t = 1.000 τ (0.040000 s, 1,000 steps)",
        "motor torque / flag = 2.50e-20 N m",
        "PASS",
    ]


def test_replay_seed_offsets_are_deterministic() -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_seed_offsets",
    )
    rows = [
        {"axis_attach_seed_value": "1", "axis_phase_seed_value": "0"},
        {"axis_attach_seed_value": "0", "axis_phase_seed_value": "1"},
        {"axis_attach_seed_value": "0", "axis_phase_seed_value": "0"},
    ]

    offsets = module._seed_offsets(rows)

    assert list(offsets) == [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)]
    assert offsets[(0.0, 0.0)] == pytest.approx(-0.24)
    assert offsets[(1.0, 0.0)] == pytest.approx(0.24)
    assert module._seed_offsets([{}]) == {(0.0, 0.0): 0.0}


def test_replay_n_flagella_metrics_plot_has_explanatory_labels(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_explanatory_metrics_plot",
    )
    rows = []
    for n_flagella, attach_seed, phase_seed, fail_t in [
        (1, 0, 0, ""),
        (1, 0, 1, ""),
        (4, 0, 2, "0.8971"),
    ]:
        rows.append(
            {
                "duration_s": "1.0",
                "n_flagella": str(n_flagella),
                "axis_n_flagella_value": str(n_flagella),
                "axis_attach_seed_value": str(attach_seed),
                "axis_phase_seed_value": str(phase_seed),
                "final_shape_pass_nonbody": "True",
                "first_fail_t_s": fail_t,
                "first_fail_category_nonbody": "flag" if fail_t else "none",
                "max_flag_bond_rel_err": "1.08" if fail_t else "0.2",
                "axis_center_net_abs_revolutions_mean": "2.0",
                "axis_center_body_relative_net_abs_revolutions_mean": "1.8",
                "body_roll_net_abs_revolutions": "0.2",
                "axis_center_to_body_roll_ratio_mean": "10.0",
            }
        )
    original_close = module.plt.close
    monkeypatch.setattr(module.plt, "close", lambda _fig: None)

    out_path = module._plot_metrics(rows=rows, out_dir=tmp_path)
    fig = module.plt.gcf()

    assert out_path.is_file()
    assert len(fig.axes) == 6
    assert all(ax.get_xlabel() == "Number of flagella (n_flagella)" for ax in fig.axes)
    assert fig.axes[0].get_ylabel() == "Time to first QC failure or run end [s]"
    assert fig.axes[1].get_ylabel() == "Maximum flagellar bond relative error [-]"
    assert [text.get_text() for text in fig.legends[0].get_texts()] == [
        "PASS: no QC failure during run",
        "FAIL: QC threshold exceeded at least once",
    ]
    assert any("attach_seed=0" in text.get_text() for text in fig.axes[0].texts)
    assert any(line.get_linestyle() == "--" for line in fig.axes[1].lines)
    original_close(fig)


def test_replay_generic_campaign_with_n_flagella_column_uses_bar_plot(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/phase2_replay.py"),
        "phase2_replay_generic_campaign_metrics_plot",
    )
    rows = [
        {
            "condition_id": "torque_1p5e20",
            "condition_label": "torque=1.5e-20",
            "axis_torque_value": "1.5e-20",
            "n_flagella": "3",
        },
        {
            "condition_id": "torque_2p0e20",
            "condition_label": "torque=2.0e-20",
            "axis_torque_value": "2.0e-20",
            "n_flagella": "3",
        },
    ]
    calls: list[str] = []
    monkeypatch.setattr(
        module,
        "_plot_metrics_by_n_flagella",
        lambda **_kwargs: calls.append("n_flagella"),
    )
    monkeypatch.setattr(
        module,
        "_plot_metrics_as_bars",
        lambda **_kwargs: calls.append("bars"),
    )

    out_path = module._plot_metrics(rows=rows, out_dir=tmp_path)

    assert out_path == tmp_path / "shape_stability_metrics.png"
    assert calls == ["bars"]


def test_generic_multi_run_summary_fieldnames_include_body_shape_gate() -> None:
    fields = _summary_fieldnames(
        [
            {
                "condition_id": "nf04",
                "body_shape_pass": True,
                "body_fail_category": "none",
                "body_spring_max_stretch_ratio": 0.1,
                "body_bend_max_error_deg": 2.0,
                "body_centerline_max_deviation_um": 0.05,
                "body_triangle_area_ratio_min": 0.95,
            }
        ]
    )

    assert "body_shape_pass" in fields
    assert "body_fail_category" in fields
    assert "body_spring_max_stretch_ratio" in fields
    assert "body_bend_max_error_deg" in fields
    assert "body_centerline_max_deviation_um" in fields
    assert "body_triangle_area_ratio_min" in fields


def test_generic_heatmap_hydrates_compact_body_gate_and_axes(tmp_path: Path) -> None:
    module = _load_script(
        Path("src/sim_swim/analysis/heatmaps/generic_multi_run.py"),
        "generic_heatmap_compact_hydration",
    )
    condition_dir = tmp_path / "body2p0__phase0"
    condition_dir.mkdir()
    (condition_dir / "run_summary.json").write_text(
        json.dumps(
            {
                "gates": {
                    "shape_body": {"first_observed_fail_t_s": 0.0413},
                },
                "all_step_metrics": {
                    "body_spring_max_stretch_ratio": {"max": 1.25},
                    "body_triangle_area_ratio_min": {"min": 0.95, "max": 1.0},
                },
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "run_manifest.json").write_text(
        json.dumps(
            {
                "conditions": [
                    {
                        "condition_id": "body2p0__phase0",
                        "output_dir": str(condition_dir),
                        "axis_values": {"body_stiffness": 2.0},
                        "axis_labels": {"body_stiffness": "candidate_2p0"},
                        "axis_order": {"body_stiffness": 1},
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    rows = [{"condition_id": "body2p0__phase0", "body_shape_pass": "False"}]

    module._hydrate_compact_rows_from_manifest(rows, tmp_path)

    assert rows[0]["axis_body_stiffness_label"] == "candidate_2p0"
    assert rows[0]["axis_body_stiffness_index"] == "1"
    assert rows[0]["body_first_fail_t_s"] == "0.0413"
    assert rows[0]["body_spring_max_stretch_ratio"] == "1.25"
    assert rows[0]["body_triangle_area_ratio_min"] == "0.95"


def test_generic_compact_summary_retains_all_body_metrics(tmp_path: Path) -> None:
    (tmp_path / "run_summary.json").write_text(
        json.dumps(
            {
                "execution": {"status": "completed", "expected_total_steps": 100},
                "gates": {
                    "shape_nonbody": {},
                    "shape_body": {
                        "any_fail": True,
                        "first_failure_category": "body_spring",
                        "first_observed_fail_t_s": 0.0123,
                    },
                },
                "all_step_metrics": {
                    "body_spring_max_stretch_ratio": {"max": 1.2},
                    "body_bend_max_error_deg": {"max": 61.0},
                    "body_centerline_max_deviation_um": {"max": 2.1},
                    "body_triangle_area_ratio_min": {"min": 0.95, "max": 1.0},
                },
            }
        ),
        encoding="utf-8",
    )
    condition = {
        "condition_id": "body1p0",
        "condition_index": 0,
        "condition_label": "body=baseline",
        "axis_values": {"body_stiffness": 1.0},
        "axis_labels": {"body_stiffness": "baseline"},
        "axis_ids": {"body_stiffness": "body1p0"},
        "axis_order": {"body_stiffness": 0},
    }

    row = _condition_row(None, condition, tmp_path)  # type: ignore[arg-type]

    assert row["body_first_fail_t_s"] == 0.0123
    assert row["body_spring_max_stretch_ratio"] == 1.2
    assert row["body_bend_max_error_deg"] == 61.0
    assert row["body_centerline_max_deviation_um"] == 2.1
    assert row["body_triangle_area_ratio_min"] == 0.95
    assert row["axis_body_stiffness_label"] == "baseline"


def test_compact_summary_records_minimum_body_triangle_area_ratio() -> None:
    summary = OnlineRunSummary(expected_steps=2)
    initial = {
        "t_s": 0.0,
        "body_triangle_area_min": 4.0,
        "body_spring_max_stretch_ratio": 0.1,
        "body_bend_max_error_deg": 2.0,
        "body_centerline_max_deviation_um": 0.05,
    }
    summary.record_body(initial)
    summary.record_body({**initial, "t_s": 0.1, "body_triangle_area_min": 3.0})

    metric = summary.extrema["body_triangle_area_ratio_min"]
    assert metric["min"] == pytest.approx(0.75)
    assert metric["final"] == pytest.approx(0.75)
