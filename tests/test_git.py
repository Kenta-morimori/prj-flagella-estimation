from __future__ import annotations

import pytest

from flagella_estimation.core.utils import git as gitmod
from flagella_estimation.core.utils.git import GitNotReadyError


def test_require_committed_state_ok(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(cmd: list[str]) -> str:
        if cmd == ["git", "rev-parse", "--is-inside-work-tree"]:
            return "true"
        if cmd == ["git", "status", "--porcelain"]:
            return ""
        if cmd == ["git", "rev-parse", "HEAD"]:
            return "a" * 40
        if cmd == ["git", "rev-parse", "--abbrev-ref", "HEAD"]:
            return "main"
        if cmd == ["git", "rev-parse", "--short", "HEAD"]:
            return "a" * 7
        raise AssertionError(f"unexpected cmd: {cmd}")

    monkeypatch.setattr(gitmod, "_run", fake_run)
    info = gitmod.require_committed_state()
    assert info.commit == "a" * 40
    assert info.commit_short == "a" * 7
    assert info.branch == "main"


def test_require_committed_state_dirty(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_run(cmd: list[str]) -> str:
        if cmd == ["git", "rev-parse", "--is-inside-work-tree"]:
            return "true"
        if cmd == ["git", "status", "--porcelain"]:
            return " M src/x.py"
        return ""

    monkeypatch.setattr(gitmod, "_run", fake_run)
    with pytest.raises(GitNotReadyError):
        gitmod.require_committed_state()
