from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys


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


def test_start_records_deterministic_output_root_and_typed_command(
    tmp_path: Path, monkeypatch
) -> None:
    helper = _load_script("cs10_parallel_tmux_start")
    control, output_root = tmp_path / "control", tmp_path / "parallel_job"
    monkeypatch.setattr(helper, "_require_tmux", lambda: "/usr/bin/tmux")
    monkeypatch.setattr(helper, "_runtime_python", lambda: Path("/runtime/python"))
    monkeypatch.setattr(
        helper, "_git_info", lambda: {"commit": "abc", "branch": "test"}
    )
    monkeypatch.setattr(helper, "_output_paths", lambda label: (control, output_root))

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
