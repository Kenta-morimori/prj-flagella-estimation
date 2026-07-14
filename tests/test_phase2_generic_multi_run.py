from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from sim_swim.analysis.cli_profiles import list_profile_entries, load_profile_entry
from sim_swim.analysis.multi_run_campaign import (
    apply_campaign_cli_overrides,
    build_campaign_conditions,
    load_yaml,
)
from sim_swim.analysis.sweeps.generic_multi_run import _summary_fieldnames


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
                "base_config: conf/sim_swim.yaml",
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
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
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
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
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
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
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
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
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
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_generic_multi_timestamped_output",
    )

    with pytest.raises(ValueError, match="output.timestamp_subdir is true"):
        module.main([f"config={profile}"])


def test_plot_heatmap_lists_generic_multi_run_profiles(capsys) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_generic_multi_list",
    )

    module.main(["list_canonical_profiles=true"])

    output = capsys.readouterr().out
    assert "conf/phase2_multi_run/latest_model_torque_shape_stability.yaml" in output


def test_replay_load_inputs_uses_manifest_condition_order_and_output_dir(
    tmp_path: Path,
) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
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
        "base_config": "conf/sim_swim.yaml",
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
    assert base_cfg_path == Path("conf/sim_swim.yaml")


def test_replay_wrapper_accepts_config_run_dir_defaults(tmp_path: Path) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
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
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
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
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
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
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
        f"phase2_replay_auto_grid_shape_{n_conditions}",
    )

    assert module._auto_grid_shape(n_conditions) == expected_shape


def test_replay_auto_grid_layout_preserves_condition_order() -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
        "phase2_replay_auto_grid_layout_order",
    )
    rows = [{"condition_id": f"cond_{index}"} for index in range(5)]

    n_rows, n_cols, positions = module._grid_layout_for_rows(rows)

    assert (n_rows, n_cols) == (2, 3)
    assert positions == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]


def test_replay_grid_layout_keeps_explicit_summary_positions() -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
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
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
        "phase2_replay_page_index_groups",
    )

    assert module._page_index_groups(36, 9) == [
        list(range(0, 9)),
        list(range(9, 18)),
        list(range(18, 27)),
        list(range(27, 36)),
    ]
    assert module._page_index_groups(3, 9) == [list(range(3))]


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
