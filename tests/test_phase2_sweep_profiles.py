from __future__ import annotations

import builtins
import importlib.util
import json
from pathlib import Path

import pytest

from sim_swim.analysis.sweeps import shape_stability_grid
from sim_swim.analysis.cli_profiles import (
    args_from_profile,
    format_profile_description,
    format_profile_listing,
    key_value_args_to_cli_args,
    list_profile_entries,
    load_profile,
    load_profile_entry,
    split_config_key,
    sweep_aliases,
)


def _load_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_sweep_profile_converts_yaml_args_to_cli_args() -> None:
    profile = load_profile(Path("conf/phase2_sweeps/shape_stability_grid.yaml"))

    assert profile["kind"] == "shape_stability_grid"
    args = args_from_profile(profile)
    assert "--duration-s" in args
    assert args[args.index("--mode") + 1] == "preset"
    assert args[args.index("--config") + 1] == "conf/sim_swim.yaml"


def test_load_profile_entry_exposes_metadata() -> None:
    entry = load_profile_entry(Path("conf/phase2_sweeps/shape_stability_grid.yaml"))

    assert entry["path"] == "conf/phase2_sweeps/shape_stability_grid.yaml"
    assert entry["metadata"]["role"] == "sweep"
    assert entry["metadata"]["canonical"] is True
    assert (
        entry["metadata"]["recommended_heatmap_profile"]
        == "conf/phase2_sweeps/shape_stability_heatmap.yaml"
    )


def test_profile_catalog_lists_canonical_sweeps_only() -> None:
    entries = list_profile_entries(role="sweep", canonical_only=True)
    lines = format_profile_listing(entries)

    assert any("conf/phase2_sweeps/shape_stability_grid.yaml" in line for line in lines)
    assert not any("conf/phase2_sweeps/hook_overstretch.yaml" in line for line in lines)


def test_profile_description_lists_aliases_for_canonical_profile() -> None:
    entry = load_profile_entry(Path("conf/phase2_sweeps/shape_stability_grid.yaml"))
    lines = format_profile_description(entry, list_profile_entries())

    assert "status: canonical" in lines
    assert (
        "recommended_heatmap_profile: conf/phase2_sweeps/shape_stability_heatmap.yaml"
        in lines
    )
    assert "aliases: conf/phase2_sweeps/hook_overstretch.yaml" in lines


def test_split_config_key_extracts_key_value_profile_path() -> None:
    config, args = split_config_key(
        [
            "config=conf/phase2_sweeps/shape_stability_grid.yaml",
            "dry_run=true",
        ]
    )

    assert config == Path("conf/phase2_sweeps/shape_stability_grid.yaml")
    assert args == ["dry_run=true"]


def test_key_value_args_convert_to_argparse_options() -> None:
    args = key_value_args_to_cli_args(
        [
            "mode=first-second-grid",
            "flagella.initial_helix_axis_from_rear_deg=null",
            "time.duration_s=0.001",
            "motor.torque_Nm=0",
            "first_second_spring_scales=1",
            "dry_run=true",
        ],
        aliases=sweep_aliases("shape_stability_grid"),
    )

    assert args == [
        "--mode",
        "first-second-grid",
        "--initial-helix-axis-from-rear-deg",
        "null",
        "--duration-s",
        "0.001",
        "--torque-nm",
        "0",
        "--first-second-spring-scales",
        "1",
        "--dry-run",
    ]


def test_key_value_args_reject_false_boolean() -> None:
    with pytest.raises(ValueError, match="dry_run=false is not supported"):
        key_value_args_to_cli_args(["dry_run=false"])


def test_shape_stability_grid_initial_helix_axis_arg_defaults_to_posterior() -> None:
    args = shape_stability_grid._parse_args([])

    assert args.initial_helix_axis_from_rear_deg == pytest.approx(0.0)


def test_shape_stability_grid_initial_helix_axis_arg_accepts_null() -> None:
    args = shape_stability_grid._parse_args(
        ["--initial-helix-axis-from-rear-deg", "null"]
    )
    condition = shape_stability_grid.Condition("case", "preset", "case", {})
    overrides = shape_stability_grid._overrides_for_condition(args, condition)

    assert args.initial_helix_axis_from_rear_deg is None
    assert overrides["flagella"]["initial_helix_axis_from_rear_deg"] is None


def test_run_sweep_wrapper_lists_profile_kind(capsys) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/run_sweep.py"),
        "phase2_run_sweep_wrapper",
    )

    module.main(
        ["config=conf/phase2_sweeps/shape_stability_grid.yaml", "list_kind=true"]
    )

    assert capsys.readouterr().out.strip() == "shape_stability_grid"


def test_run_sweep_wrapper_keeps_hook_overstretch_alias(capsys) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/run_sweep.py"),
        "phase2_run_sweep_wrapper_alias",
    )

    module.main(["config=conf/phase2_sweeps/hook_overstretch.yaml", "list_kind=true"])

    assert capsys.readouterr().out.strip() == "hook_overstretch"


def test_run_sweep_wrapper_lists_canonical_profiles(capsys) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/run_sweep.py"),
        "phase2_run_sweep_wrapper_list_profiles",
    )

    module.main(["list_canonical_profiles=true"])

    output = capsys.readouterr().out
    assert "conf/phase2_sweeps/shape_stability_grid.yaml" in output
    assert "conf/phase2_sweeps/hook_overstretch.yaml" not in output
    assert "conf/phase2_sweeps/shape_stability_heatmap.yaml" not in output


def test_run_sweep_wrapper_describes_profile(capsys) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/run_sweep.py"),
        "phase2_run_sweep_wrapper_describe",
    )

    module.main(
        [
            "config=conf/phase2_sweeps/shape_stability_grid.yaml",
            "describe_profile=true",
        ]
    )

    output = capsys.readouterr().out
    assert "role: sweep" in output
    assert "status: canonical" in output
    assert "recommended_heatmap_profile:" in output


def test_run_sweep_wrapper_rejects_heatmap_profile() -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/run_sweep.py"),
        "phase2_run_sweep_wrapper_wrong_role",
    )

    with pytest.raises(SystemExit, match="has role 'heatmap'; use plot_heatmap.py"):
        module.main(["config=conf/phase2_sweeps/shape_stability_heatmap.yaml"])


def test_torque_distribution_profile_is_shape_stability_grid() -> None:
    profile = load_profile(Path("conf/phase2_sweeps/torque_distribution_grid.yaml"))

    assert profile["kind"] == "shape_stability_grid"
    args = args_from_profile(profile)
    assert args[args.index("--mode") + 1] == "torque-profile-grid"
    assert (
        args[args.index("--force-distributions") + 1]
        == "root_torque_segment_couples,root_torque_axis_projection"
    )
    assert args[args.index("--torque-distribution-profiles") + 1] == "diffusive,uniform"


def test_basal_freedom_profile_builds_issue103_conditions() -> None:
    profile = load_profile(Path("conf/phase2_sweeps/basal_freedom_diagnostic.yaml"))

    assert profile["kind"] == "shape_stability_grid"
    args = shape_stability_grid._parse_args(args_from_profile(profile))
    conditions = shape_stability_grid.build_conditions(args)

    assert [condition.condition_id for condition in conditions] == [
        "no_frame",
        "fp3",
        "ft1p5",
        "fp3_ft1p5_vector",
        "fp3_ft1p5_bearing",
    ]
    assert conditions[-1].scales["local_attach_frame_tangent_mode"] == ("basal_bearing")


def test_basal_freedom_position_only_profile_builds_issue103_followup_conditions() -> (
    None
):
    profile = load_profile(
        Path("conf/phase2_sweeps/basal_freedom_position_only_sweep.yaml")
    )

    assert profile["kind"] == "shape_stability_grid"
    args = shape_stability_grid._parse_args(args_from_profile(profile))
    conditions = shape_stability_grid.build_conditions(args)

    assert [condition.condition_id for condition in conditions] == [
        "no_frame",
        "fp1p25",
        "fp1p5",
        "fp2",
        "fp2p5",
        "fp3",
    ]
    assert all(
        condition.scales["local_attach_frame_tangent_scale"] == 1.0
        for condition in conditions
    )
    assert all(
        condition.scales["local_attach_frame_tangent_mode"] == "vector"
        for condition in conditions
    )


def test_shape_stability_grid_manifest_record_includes_campaign_fields() -> None:
    args = shape_stability_grid._parse_args(
        [
            "--mode",
            "attach-frame-grid",
            "--attach-frame-position-scales",
            "1,2",
            "--attach-frame-tangent-scales",
            "1,1.5",
        ]
    )
    condition = shape_stability_grid.Condition(
        "fp2_ft1p5",
        "attach-frame-grid",
        "attach frame condition",
        {
            "local_attach_frame_position_scale": 2.0,
            "local_attach_frame_tangent_scale": 1.5,
        },
    )

    record = shape_stability_grid._manifest_record(args, condition, condition_index=3)

    assert record["condition_index"] == 3
    assert record["condition_label"] == "attach frame condition"
    assert record["axis_values"]["attach_frame_position"] == 2.0
    assert record["axis_labels"]["attach_frame_tangent"] == "1.5"
    assert record["grid"]["row_axis"] == "attach_frame_tangent"
    assert record["grid"]["col_axis"] == "attach_frame_position"


def _write_replay_inputs(tmp_path: Path, condition_ids: list[str]) -> Path:
    input_dir = tmp_path / "replay"
    input_dir.mkdir()
    (input_dir / "summary.csv").write_text(
        "\n".join(
            [
                "condition_id,force_distribution,torque_distribution_profile,local_attach_frame_tangent_mode",
                *[
                    f"{condition_id},root_torque_segment_couples,diffusive,vector"
                    for condition_id in condition_ids
                ],
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    manifest = {
        "config": "conf/sim_swim.yaml",
        "conditions": [
            {"condition_id": condition_id, "config_overrides": {}}
            for condition_id in condition_ids
        ],
    }
    (input_dir / "run_manifest.json").write_text(
        json.dumps(manifest),
        encoding="utf-8",
    )
    for condition_id in condition_ids:
        condition_dir = input_dir / condition_id
        condition_dir.mkdir()
        (condition_dir / "state_archive.npz").write_bytes(b"")
    return input_dir


def test_issue97_replay_keeps_canonical_condition_order(tmp_path: Path) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
        "phase2_issue97_replay_order",
    )
    condition_ids = [
        "axis_projection_uniform_fp3_ft1p5",
        "segment_couples_uniform_fp3_ft1p5",
        "axis_projection_diffusive_fp3_ft1p5",
        "segment_couples_diffusive_fp3_ft1p5",
    ]
    input_dir = _write_replay_inputs(tmp_path, condition_ids)

    rows, _records, _base_cfg_path = module._load_inputs(input_dir)

    assert [row["condition_id"] for row in rows] == [
        "segment_couples_diffusive_fp3_ft1p5",
        "segment_couples_uniform_fp3_ft1p5",
        "axis_projection_diffusive_fp3_ft1p5",
        "axis_projection_uniform_fp3_ft1p5",
    ]


def test_issue103_replay_accepts_basal_freedom_conditions(tmp_path: Path) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
        "phase2_issue103_replay_order",
    )
    condition_ids = [
        "no_frame",
        "fp3",
        "ft1p5",
        "fp3_ft1p5_vector",
        "fp3_ft1p5_bearing",
    ]
    input_dir = _write_replay_inputs(tmp_path, condition_ids)

    rows, _records, _base_cfg_path = module._load_inputs(input_dir)

    assert [row["condition_id"] for row in rows] == condition_ids
    assert module._label_for_row(rows[-1]) == "fp3_ft1p5_bearing"


def test_issue103_replay_accepts_position_only_conditions(tmp_path: Path) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
        "phase2_issue103_position_only_replay_order",
    )
    condition_ids = [
        "no_frame",
        "fp1p25",
        "fp1p5",
        "fp2",
        "fp2p5",
        "fp3",
    ]
    input_dir = _write_replay_inputs(tmp_path, condition_ids)

    rows, _records, _base_cfg_path = module._load_inputs(input_dir)

    assert [row["condition_id"] for row in rows] == condition_ids
    assert module._label_for_row(rows[-1]) == "fp3"


def test_replay_module_import_does_not_require_cv2(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_import = builtins.__import__

    def guarded_import(
        name: str,
        globals: dict[str, object] | None = None,
        locals: dict[str, object] | None = None,
        fromlist: tuple[str, ...] = (),
        level: int = 0,
    ) -> object:
        if name == "cv2":
            raise AssertionError("cv2 should not be imported at module load time")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)

    _load_script(
        Path("scripts/01_simulate_swimming/render_shape_stability_grid_replay.py"),
        "phase2_replay_import_without_cv2",
    )


def test_shape_stability_grid_keeps_deprecated_torque_segment_profile_alias() -> None:
    args = shape_stability_grid._parse_args(
        ["--torque-segment-weight-profiles", "diffusive,uniform"]
    )

    assert args.torque_distribution_profiles == ["diffusive", "uniform"]


def test_plot_heatmap_wrapper_rejects_unknown_kind(tmp_path: Path) -> None:
    profile = tmp_path / "bad_heatmap.yaml"
    profile.write_text("kind: missing\n", encoding="utf-8")
    module = _load_script(
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper",
    )

    try:
        module.main(["--config", str(profile)])
    except SystemExit as exc:
        assert "Unknown heatmap kind" in str(exc)
    else:
        raise AssertionError("expected SystemExit for unknown heatmap kind")


def test_plot_heatmap_wrapper_lists_canonical_profiles(capsys) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_list_profiles",
    )

    module.main(["list_canonical_profiles=true"])

    output = capsys.readouterr().out
    assert "conf/phase2_sweeps/shape_stability_heatmap.yaml" in output
    assert "conf/phase2_sweeps/hook_overstretch_heatmap.yaml" not in output
    assert "conf/phase2_sweeps/shape_stability_grid.yaml" not in output


def test_plot_heatmap_wrapper_describes_alias_profile(capsys) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_describe",
    )

    module.main(
        [
            "config=conf/phase2_sweeps/hook_overstretch_heatmap.yaml",
            "describe_profile=true",
        ]
    )

    output = capsys.readouterr().out
    assert "role: heatmap" in output
    assert "status: alias" in output
    assert "alias_of: conf/phase2_sweeps/shape_stability_heatmap.yaml" in output


def test_plot_heatmap_wrapper_rejects_sweep_profile() -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_wrong_role",
    )

    with pytest.raises(SystemExit, match="has role 'sweep'; use run_sweep.py"):
        module.main(["config=conf/phase2_sweeps/shape_stability_grid.yaml"])


def test_heatmap_profiles_do_not_fix_output_dir() -> None:
    profile_paths = sorted(Path("conf/phase2_sweeps").glob("*heatmap.yaml"))

    assert profile_paths
    for profile_path in profile_paths:
        profile = load_profile(profile_path)
        raw_args = profile.get("args") or {}
        assert "output_dir" not in raw_args, profile_path


def test_plot_heatmap_wrapper_defaults_output_dir_next_to_summary(
    tmp_path: Path,
) -> None:
    profile = tmp_path / "heatmap.yaml"
    profile.write_text(
        "kind: shape_stability_grid\nargs:\n  mode: first-second-grid\n",
        encoding="utf-8",
    )
    summary_csv = tmp_path / "summary.csv"
    captured: dict[str, list[str]] = {}
    module = _load_script(
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_default_output",
    )
    module.HEATMAP_MAIN["shape_stability_grid"] = lambda args: captured.setdefault(
        "args", args
    )

    module.main(["config=" + str(profile), "summary_csv=" + str(summary_csv)])

    args = captured["args"]
    assert args[args.index("--output-dir") + 1] == str(tmp_path / "plots")


def test_plot_heatmap_wrapper_accepts_position_only_mode_override(
    tmp_path: Path,
) -> None:
    profile = tmp_path / "heatmap.yaml"
    profile.write_text(
        "kind: shape_stability_grid\nargs:\n  mode: attach-frame-grid\n",
        encoding="utf-8",
    )
    summary_csv = tmp_path / "summary.csv"
    captured: dict[str, list[str]] = {}
    module = _load_script(
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_position_only",
    )
    module.HEATMAP_MAIN["shape_stability_grid"] = lambda args: captured.setdefault(
        "args", args
    )

    module.main(
        [
            "config=" + str(profile),
            "mode=position-only-grid",
            "summary_csv=" + str(summary_csv),
        ]
    )

    args = captured["args"]
    mode_indices = [index for index, arg in enumerate(args) if arg == "--mode"]
    assert args[mode_indices[-1] + 1] == "position-only-grid"
    assert args[args.index("--output-dir") + 1] == str(tmp_path / "plots")


def test_plot_heatmap_wrapper_keeps_hook_overstretch_heatmap_alias(
    tmp_path: Path,
) -> None:
    summary_csv = tmp_path / "summary.csv"
    captured: dict[str, list[str]] = {}
    module = _load_script(
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_alias",
    )
    module.HEATMAP_MAIN["hook_overstretch"] = lambda args: captured.setdefault(
        "args", args
    )

    module.main(
        [
            "config=conf/phase2_sweeps/hook_overstretch_heatmap.yaml",
            "summary_csv=" + str(summary_csv),
        ]
    )

    args = captured["args"]
    assert args[args.index("--mode") + 1] == "attach-frame-grid"
    assert args[args.index("--output-dir") + 1] == str(tmp_path / "plots")


def test_plot_heatmap_wrapper_keeps_explicit_output_dir(tmp_path: Path) -> None:
    profile = tmp_path / "heatmap.yaml"
    profile_output_dir = tmp_path / "profile_plots"
    profile.write_text(
        "\n".join(
            [
                "kind: shape_stability_grid",
                "args:",
                "  mode: first-second-grid",
                f"  output_dir: {profile_output_dir}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    summary_csv = tmp_path / "summary.csv"
    captured: dict[str, list[str]] = {}
    module = _load_script(
        Path("scripts/01_simulate_swimming/plot_heatmap.py"),
        "phase2_plot_heatmap_wrapper_explicit_output",
    )
    module.HEATMAP_MAIN["shape_stability_grid"] = lambda args: captured.setdefault(
        "args", args
    )

    module.main(["config=" + str(profile), "summary_csv=" + str(summary_csv)])

    args = captured["args"]
    assert args.count("--output-dir") == 1
    assert args[args.index("--output-dir") + 1] == str(profile_output_dir)
