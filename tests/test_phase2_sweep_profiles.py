from __future__ import annotations

import importlib.util
from pathlib import Path

from sim_swim.analysis.cli_profiles import args_from_profile, load_profile


def _load_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_sweep_profile_converts_yaml_args_to_cli_args() -> None:
    profile = load_profile(Path("conf/phase2_sweeps/hook_overstretch.yaml"))

    assert profile["kind"] == "hook_overstretch"
    args = args_from_profile(profile)
    assert "--duration-s" in args
    assert args[args.index("--mode") + 1] == "preset"
    assert args[args.index("--config") + 1] == "conf/sim_swim.yaml"


def test_run_sweep_wrapper_lists_profile_kind(capsys) -> None:
    module = _load_script(
        Path("scripts/01_simulate_swimming/run_sweep.py"),
        "phase2_run_sweep_wrapper",
    )

    module.main(["--config", "conf/phase2_sweeps/hook_overstretch.yaml", "--list-kind"])

    assert capsys.readouterr().out.strip() == "hook_overstretch"


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
