from __future__ import annotations

import json
from pathlib import Path

import pytest

from flagella_estimation.core.run_context import init_run
from flagella_estimation.core.utils.git import GitInfo
from flagella_estimation.core.utils.time import RunTime


def test_init_run_creates_log_and_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "flagella_estimation.core.run_context.require_committed_state",
        lambda: GitInfo(commit="a" * 40, commit_short="a" * 7, branch="main"),
    )
    monkeypatch.setattr(
        "flagella_estimation.core.run_context.now_jst",
        lambda: RunTime(date="2026-01-21", time="120000"),
    )

    ctx = init_run(base_dir=tmp_path, input_info={"data_path": "dummy.mp4"})

    assert ctx.out.root.is_dir()
    assert ctx.out.tracking_dir.is_dir()

    log_path = ctx.out.root / "run.log"
    manifest_path = ctx.out.root / "manifest.json"
    assert log_path.is_file()
    assert manifest_path.is_file()

    log_text = log_path.read_text(encoding="utf-8")
    assert "Run start" in log_text

    m = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert m["git"]["commit"] == "a" * 40
    assert m["run_time"]["date"] == "2026-01-21"
    assert Path(m["outputs"]["tracking_dir"]).name == "tracking"
