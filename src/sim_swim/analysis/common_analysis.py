"""Dispatch reusable Phase 2 analyses from the common dataset CLI."""

from __future__ import annotations

import argparse
from collections.abc import Callable
from pathlib import Path
from typing import Any
import sys

from sim_swim.analysis.cli_profiles import (
    args_from_profile,
    format_profile_description,
    format_profile_listing,
    key_value_args_to_cli_args,
    list_profile_entries,
    load_profile_entry,
    split_config_key,
    validate_profile_role,
)
from sim_swim.analysis.heatmaps import (
    dt_star_torque,
    generic_multi_run,
    hook_overstretch,
    local_scale_mode,
    motor_scale_collapse,
    shape_stability_grid,
)


HEATMAP_MAIN: dict[str, Callable[[list[str]], Any]] = {
    "motor_scale_collapse": motor_scale_collapse.main,
    "dt_star_torque": dt_star_torque.main,
    "local_scale_mode": local_scale_mode.main,
    "shape_stability_grid": shape_stability_grid.main,
    "hook_overstretch": hook_overstretch.main,
    "generic_multi_run": generic_multi_run.main,
}


def _heatmap_entries(*, canonical_only: bool = False) -> list[dict[str, object]]:
    entries = list_profile_entries(role="heatmap", canonical_only=canonical_only)
    entries.extend(
        entry
        for entry in list_profile_entries(role="sweep", canonical_only=canonical_only)
        if entry["kind"] == "generic_multi_run"
    )
    seen: set[str] = set()
    return [
        entry
        for entry in entries
        if not (str(entry["path"]) in seen or seen.add(str(entry["path"])))
    ]


def _has_option(args: list[str], option_name: str) -> bool:
    return any(arg == option_name or arg.startswith(f"{option_name}=") for arg in args)


def _first_option_value(args: list[str], option_name: str) -> str | None:
    for index, arg in enumerate(args):
        if arg == option_name and index + 1 < len(args):
            return args[index + 1]
        if arg.startswith(f"{option_name}="):
            return arg.split("=", 1)[1]
    return None


def _with_run_dir_defaults(args: list[str]) -> list[str]:
    run_dir = _first_option_value(args, "--run-dir")
    if run_dir is not None and not _has_option(args, "--summary-csv"):
        args = [*args, "--summary-csv", str(Path(run_dir) / "summary.csv")]
    if _has_option(args, "--output-dir"):
        return args
    if run_dir is not None:
        return [*args, "--output-dir", str(Path(run_dir) / "plots")]
    summary_csv = _first_option_value(args, "--summary-csv")
    if summary_csv is None:
        return args
    return [*args, "--output-dir", str(Path(summary_csv).parent / "plots")]


def run_heatmap(argv: list[str]) -> None:
    config_from_key, parser_argv = split_config_key(argv)
    parser = argparse.ArgumentParser(
        description="Plot a configured Phase 2 sweep heatmap."
    )
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--list-profiles", action="store_true")
    parser.add_argument("--list-canonical-profiles", action="store_true")
    parser.add_argument("--describe-profile", action="store_true")
    args, passthrough = parser.parse_known_args(key_value_args_to_cli_args(parser_argv))
    if config_from_key is not None and args.config is not None:
        parser.error("Use either config=PATH or --config PATH (not both)")
    if args.list_profiles or args.list_canonical_profiles:
        for line in format_profile_listing(
            _heatmap_entries(canonical_only=args.list_canonical_profiles)
        ):
            print(line)
        return
    config = config_from_key or args.config
    if config is None:
        parser.error("config=PATH or --config PATH is required")
    entry = load_profile_entry(config)
    kind = str(entry["kind"])
    if kind not in HEATMAP_MAIN:
        raise SystemExit(f"Unknown heatmap kind {kind!r}")
    if args.describe_profile:
        for line in format_profile_description(entry, _heatmap_entries()):
            print(line)
        return
    if kind != "generic_multi_run":
        validate_profile_role(entry, "heatmap")
    effective_args = (
        _with_run_dir_defaults(
            ["--campaign-config", str(config)] + key_value_args_to_cli_args(passthrough)
        )
        if kind == "generic_multi_run"
        else _with_run_dir_defaults(args_from_profile(entry) + passthrough)
    )
    HEATMAP_MAIN[kind](effective_args)


def dispatch(analysis_kind: str, argv: list[str]) -> None:
    if analysis_kind == "heatmap":
        run_heatmap(argv)
        return
    if analysis_kind == "run-summary":
        from sim_swim.analysis.run_summary import write_run_summary

        parser = argparse.ArgumentParser()
        parser.add_argument("--input-dir", type=Path, required=True)
        parser.add_argument("--overwrite", action="store_true")
        args = parser.parse_args(argv)
        print(write_run_summary(args.input_dir, overwrite=args.overwrite))
        return
    if analysis_kind == "2010-torque-dt":
        from sim_swim.analysis.torque_dt_stability_campaign import summarize_campaign

        parser = argparse.ArgumentParser()
        parser.add_argument("--run-dir", type=Path, required=True)
        parser.add_argument(
            "--config",
            type=Path,
            default=Path(
                "conf/phase2_multi_run/2010_project_torque_dt_initial_screen.yaml"
            ),
        )
        args = parser.parse_args(argv)
        summarize_campaign(args.run_dir, config_path=args.config)
        print(args.run_dir / "qc_summary.json")
        return
    if analysis_kind == "issue61-2015-1tau":
        from sim_swim.analysis.issue61_2015_1tau import main

        main(argv)
        return
    if analysis_kind == "2010-fixed-performance":
        from sim_swim.analysis.torque_dt_stability_campaign import (
            render_fixed_real_time_qualitative_replay,
            summarize_campaign,
        )

        parser = argparse.ArgumentParser()
        parser.add_argument("--run-dir", type=Path, required=True)
        parser.add_argument(
            "--config",
            type=Path,
            default=Path(
                "conf/phase2_multi_run/2010_project_torque_fixed_real_time_0p5s_performance.yaml"
            ),
        )
        parser.add_argument("--allow-incomplete", action="store_true")
        parser.add_argument("--render-qualitative-replay", action="store_true")
        parser.add_argument("--fps-out-3d", type=float, default=20.0)
        parser.add_argument("--target-frame-count", type=int, default=101)
        args = parser.parse_args(argv)
        summarize_campaign(
            args.run_dir,
            config_path=args.config,
            allow_incomplete=args.allow_incomplete,
        )
        if args.render_qualitative_replay:
            print(
                render_fixed_real_time_qualitative_replay(
                    args.run_dir,
                    config_path=args.config,
                    allow_incomplete=args.allow_incomplete,
                    fps_out_3d=args.fps_out_3d,
                    target_frame_count=args.target_frame_count,
                )
            )
        print(args.run_dir / "performance_summary.csv")
        return
    if analysis_kind == "diagnose-158":
        from sim_swim.analysis.phase2_158_diagnostics import (
            DEFAULT_FAIL_CONDITION_IDS,
            Phase2158DiagnosticConfig,
            analyze_phase2_158_diagnostics,
            default_output_dir,
        )

        parser = argparse.ArgumentParser()
        parser.add_argument(
            "--campaign-config",
            type=Path,
            default=Path("conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml"),
        )
        parser.add_argument(
            "--input-dir",
            type=Path,
            default=Path("outputs/phase2_multi_run/flagella_count_duration_3s_r1"),
        )
        parser.add_argument("--output-dir", type=Path, default=None)
        parser.add_argument("--fail-condition-id", action="append", default=None)
        parser.add_argument("--first-fail-window-s", type=float, default=0.25)
        parser.add_argument("--overwrite", action="store_true")
        parser.add_argument("overrides", nargs="*")
        args = parser.parse_args(argv)
        output_dir = args.output_dir or default_output_dir()
        print(
            analyze_phase2_158_diagnostics(
                Phase2158DiagnosticConfig(
                    campaign_config=args.campaign_config,
                    input_dir=args.input_dir,
                    output_dir=output_dir,
                    fail_condition_ids=tuple(
                        args.fail_condition_id or DEFAULT_FAIL_CONDITION_IDS
                    ),
                    first_fail_window_s=args.first_fail_window_s,
                    overwrite=args.overwrite,
                    cli_overrides=tuple(args.overrides),
                )
            )
        )
        return
    module_name = {
        "2015-stage-a": "sim_swim.analysis.stage_a_2015_workflow",
        "spring-formulations": "sim_swim.analysis.spring_formulations_workflow",
        "task-d-2015": "sim_swim.analysis.task_d_2015_tau_linked",
        "partial-generic": "sim_swim.analysis.partial_generic_multi_run",
        "plan-2010-torque-dt": "sim_swim.analysis.torque_dt_stability",
        "plan-reference-torque": "sim_swim.analysis.reference_torque_comparison",
        "visualize-2010-torque-dt": "sim_swim.analysis.torque_dt_stability_visuals",
        "phase-seed-2d": "sim_swim.analysis.phase_seed_2d_replay",
        "distributions": "sim_swim.analysis.behavior_dataset_distributions",
    }.get(analysis_kind)
    if module_name is None:
        raise ValueError(f"Unsupported analysis kind: {analysis_kind}")
    module = __import__(module_name, fromlist=["main"])
    old_argv = sys.argv
    try:
        sys.argv = [module_name, *argv]
        module.main()
    finally:
        sys.argv = old_argv
