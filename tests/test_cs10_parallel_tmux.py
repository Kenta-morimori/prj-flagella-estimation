from __future__ import annotations

import importlib.util
import json
from datetime import datetime
from pathlib import Path
import subprocess
import sys
from zoneinfo import ZoneInfo

import pytest


ROOT = Path(__file__).resolve().parents[1]
ISSUE203 = ROOT / "conf/phase2_parallel/issue203_uniform_torque_profile/job.yaml"


def _load_script(name: str):
    path = ROOT / "scripts/cs10/parallel_tmux.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_start_records_nas_output_root_and_typed_command(
    tmp_path: Path, monkeypatch
) -> None:
    helper = _load_script("cs10_parallel_tmux_start")
    control, output_root = tmp_path / "control", tmp_path / "parallel_job"
    monkeypatch.setattr(helper, "_require_tmux", lambda: "/usr/bin/tmux")
    monkeypatch.setattr(helper, "_runtime_python", lambda: Path("/runtime/python"))
    monkeypatch.setattr(
        helper, "_git_info", lambda: {"commit": "abc", "branch": "test"}
    )
    monkeypatch.setattr(helper, "_require_output_base", lambda: tmp_path / "nas")
    monkeypatch.setattr(
        helper,
        "_output_paths",
        lambda label, output_base: (control, output_root),
    )

    commands: list[list[str]] = []

    def fake_run(command, **kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(
            command, 1 if "has-session" in command else 0
        )

    monkeypatch.setattr(helper.subprocess, "run", fake_run)
    record = helper.start(ISSUE203, session="issue203", label="issue203_uniform")

    assert record["output_root"] == str(output_root)
    assert record["execution"]["max_workers"] == 8
    assert "--output-root" in record["command"]
    assert (control / "launch.json").is_file()
    assert str(output_root) in (control / "launch.sh").read_text(encoding="utf-8")
    assert any(command[1:3] == ["new-session", "-d"] for command in commands)


def test_output_paths_keep_control_local_and_put_job_on_nas(
    monkeypatch, tmp_path: Path
) -> None:
    helper = _load_script("cs10_parallel_tmux_output_paths")
    monkeypatch.setattr(helper, "REPOSITORY_ROOT", tmp_path / "repository")
    monkeypatch.setattr(
        helper,
        "_now",
        lambda: datetime(2026, 8, 28, 12, 34, 56, tzinfo=ZoneInfo("Asia/Tokyo")),
    )
    monkeypatch.setattr(
        helper, "uuid4", lambda: type("Id", (), {"hex": "abc123def456"})()
    )

    control, output_root = helper._output_paths("example", tmp_path / "nas/outputs")

    assert control == (
        tmp_path / "repository/outputs/2026-08-28/123456/cs10_parallel/example"
    )
    assert output_root == (
        tmp_path / "nas/outputs/2026-08-28/123456/parallel/example__abc123def456"
    )


def test_require_output_base_creates_missing_project_output_directory(
    monkeypatch, tmp_path: Path
) -> None:
    helper = _load_script("cs10_parallel_tmux_create_nas")
    output_base = tmp_path / "nas/project/outputs"
    monkeypatch.setattr(helper, "CS10_OUTPUT_BASE", output_base)
    monkeypatch.setattr(helper, "_require_nas_mount", lambda: None)

    assert helper._require_output_base() == output_base
    assert output_base.is_dir()


def test_require_output_base_reports_directory_creation_failure(
    monkeypatch, tmp_path: Path
) -> None:
    helper = _load_script("cs10_parallel_tmux_create_failure")
    blocked_parent = tmp_path / "blocked"
    blocked_parent.write_text("not a directory", encoding="utf-8")
    monkeypatch.setattr(helper, "CS10_OUTPUT_BASE", blocked_parent / "outputs")
    monkeypatch.setattr(helper, "_require_nas_mount", lambda: None)

    with pytest.raises(
        RuntimeError, match="could not create cs10 NAS output directory"
    ):
        helper._require_output_base()


def test_require_output_base_checks_write_access(monkeypatch, tmp_path: Path) -> None:
    helper = _load_script("cs10_parallel_tmux_writable_nas")
    output_base = tmp_path / "nas"
    output_base.mkdir()
    monkeypatch.setattr(helper, "CS10_OUTPUT_BASE", output_base)
    monkeypatch.setattr(helper, "_require_nas_mount", lambda: None)

    assert helper._require_output_base() == output_base


def test_require_output_base_reports_write_failure(monkeypatch, tmp_path: Path) -> None:
    helper = _load_script("cs10_parallel_tmux_unwritable_nas")
    output_base = tmp_path / "nas"
    output_base.mkdir()
    monkeypatch.setattr(helper, "CS10_OUTPUT_BASE", output_base)
    monkeypatch.setattr(helper, "_require_nas_mount", lambda: None)

    def fail_temporary_file(**kwargs):
        raise OSError("read-only filesystem")

    monkeypatch.setattr(helper.tempfile, "TemporaryFile", fail_temporary_file)

    with pytest.raises(RuntimeError, match="NAS output directory is not writable"):
        helper._require_output_base()


def test_require_nas_mount_rejects_unmounted_path(monkeypatch, tmp_path: Path) -> None:
    helper = _load_script("cs10_parallel_tmux_unmounted_nas")
    mount_root = tmp_path / "work01"
    mount_root.mkdir()
    monkeypatch.setattr(helper, "CS10_NAS_MOUNT_ROOT", mount_root)

    with pytest.raises(RuntimeError, match="cs10 NAS mount is unavailable"):
        helper._require_nas_mount()


def test_start_does_not_create_tmux_when_nas_preflight_fails(monkeypatch) -> None:
    helper = _load_script("cs10_parallel_tmux_nas_preflight")
    monkeypatch.setattr(helper, "_require_tmux", lambda: "/usr/bin/tmux")
    monkeypatch.setattr(helper, "_runtime_python", lambda: Path("/runtime/python"))
    monkeypatch.setattr(
        helper, "_git_info", lambda: {"commit": "abc", "branch": "test"}
    )

    def fail_output_preflight() -> Path:
        raise RuntimeError("NAS is not writable")

    monkeypatch.setattr(helper, "_require_output_base", fail_output_preflight)
    commands: list[list[str]] = []

    def fake_run(command, **kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(helper.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError, match="NAS is not writable"):
        helper.start(ISSUE203, session="issue203", label="issue203_uniform")

    assert commands == []


def test_runtime_python_includes_repository_source_in_preflight(
    monkeypatch, tmp_path: Path
) -> None:
    helper = _load_script("cs10_parallel_tmux_runtime")
    source = tmp_path / "src"
    runtime = tmp_path / ".venv-cs10/bin/python"
    runtime.parent.mkdir(parents=True)
    runtime.touch()
    monkeypatch.setattr(helper, "SOURCE_ROOT", source)
    monkeypatch.setattr(helper, "REPOSITORY_ROOT", tmp_path)
    monkeypatch.setenv("PYTHONPATH", "/existing/path")
    captured: dict[str, object] = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["env"] = kwargs["env"]
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(helper.subprocess, "run", fake_run)

    assert helper._runtime_python() == runtime
    assert captured["command"] == [
        str(runtime),
        "-c",
        "import matplotlib, numpy, yaml, sim_swim",
    ]
    assert captured["env"]["PYTHONPATH"] == f"{source}:/existing/path"


def test_status_follows_canonical_condition_symlinks(tmp_path: Path) -> None:
    helper = _load_script("cs10_parallel_tmux_status")
    control, output_root = tmp_path / "control", tmp_path / "parallel_job"
    control.mkdir()
    (output_root / "campaign/conditions").mkdir(parents=True)
    source = tmp_path / "child/condition"
    source.mkdir(parents=True)
    (source / "run_summary.json").write_text("{}", encoding="utf-8")
    (output_root / "campaign/conditions/condition_a").symlink_to(source)
    (control / "launch.json").write_text(
        json.dumps({"output_root": str(output_root)}), encoding="utf-8"
    )
    (output_root / "job_manifest.json").write_text(
        json.dumps({"status": "succeeded", "configs": [{"status": "succeeded"}]}),
        encoding="utf-8",
    )
    (output_root / "campaign/campaign_completion.json").write_text(
        json.dumps({"status": "completed", "exit_code": 0}), encoding="utf-8"
    )

    result = helper.status(control)

    assert result["job"]["succeeded"] == 1
    assert result["campaign"]["run_summary_count"] == 1
