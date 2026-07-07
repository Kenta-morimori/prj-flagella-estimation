from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from sim_swim.analysis.cli_profiles import (
    args_from_profile,
    key_value_args_to_cli_args,
    load_profile,
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
