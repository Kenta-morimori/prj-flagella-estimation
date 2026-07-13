from __future__ import annotations

import json
from pathlib import Path

import pytest

from sim_swim.core.run_context import init_run

pytestmark = pytest.mark.light


def test_init_run_creates_log_and_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("sim_swim.core.run_context._require_clean_git", lambda: None)
    monkeypatch.setattr(
        "sim_swim.core.run_context._now_jst",
        lambda: ("2026-02-18", "210000"),
    )
    ctx = init_run(base_dir=tmp_path, input_info={"config": "conf/sim_swim.yaml"})

    assert ctx.out.root.is_dir()
    assert ctx.out.sim_dir.is_dir()
    assert ctx.out.render_dir.is_dir()
    assert ctx.out.render2d_dir.is_dir()

    log_path = ctx.out.root / "run.log"
    manifest_path = ctx.out.root / "manifest.json"
    assert log_path.is_file()
    assert manifest_path.is_file()

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["input"]["config"] == "conf/sim_swim.yaml"
    assert Path(manifest["outputs"]["sim_dir"]).name == "sim"
    assert manifest["run_time"]["timestamp_subdir"] is True


def test_init_run_can_use_base_dir_as_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("sim_swim.core.run_context._require_clean_git", lambda: None)
    monkeypatch.setattr(
        "sim_swim.core.run_context._now_jst",
        lambda: ("2026-02-18", "210000"),
    )

    run_dir = tmp_path / "run"
    ctx = init_run(
        base_dir=run_dir,
        input_info={"config": "conf/sim_swim.yaml"},
        timestamp_subdir=False,
    )

    assert ctx.out.root == run_dir
    assert (run_dir / "manifest.json").is_file()
    assert (run_dir / "sim").is_dir()
    assert not (run_dir / "2026-02-18").exists()


def test_init_run_fixed_root_requires_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("sim_swim.core.run_context._require_clean_git", lambda: None)

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "manifest.json").write_text("{}", encoding="utf-8")

    with pytest.raises(FileExistsError, match="overwrite=true"):
        init_run(
            base_dir=run_dir,
            input_info={"config": "conf/sim_swim.yaml"},
            timestamp_subdir=False,
        )


def test_init_run_fixed_root_overwrite_replaces_existing_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr("sim_swim.core.run_context._require_clean_git", lambda: None)

    run_dir = tmp_path / "run"
    run_dir.mkdir()
    stale = run_dir / "stale.txt"
    stale.write_text("old", encoding="utf-8")

    ctx = init_run(
        base_dir=run_dir,
        input_info={"config": "conf/sim_swim.yaml"},
        timestamp_subdir=False,
        overwrite=True,
    )

    assert ctx.out.root == run_dir
    assert not stale.exists()
    assert (run_dir / "manifest.json").is_file()


def test_init_run_aborts_when_dirty_and_leaves_no_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "sim_swim.core.run_context._run",
        lambda cmd: " M src/sim_swim/core/run_context.py",
    )

    with pytest.raises(RuntimeError) as exc_info:
        init_run(base_dir=tmp_path, input_info={"config": "conf/sim_swim.yaml"})

    msg = str(exc_info.value)
    assert "未コミット変更があるため実行を中止" in msg
    assert "git status --porcelain" in msg
    assert "commit/stash" in msg
    assert not any(tmp_path.iterdir())
