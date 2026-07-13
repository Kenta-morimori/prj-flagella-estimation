from __future__ import annotations

import importlib.util
import json
from pathlib import Path

from sim_swim.analysis.cli_profiles import list_profile_entries, load_profile_entry
from sim_swim.analysis.multi_run_campaign import (
    apply_campaign_cli_overrides,
    build_campaign_conditions,
    load_yaml,
)


def _load_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_generic_multi_run_profile_is_listed() -> None:
    entries = list_profile_entries(role="sweep", canonical_only=True)
    paths = {entry["path"] for entry in entries}

    assert "conf/phase2_multi_run/latest_model_torque_shape_stability.yaml" in paths


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
    assert args.overwrite is True
