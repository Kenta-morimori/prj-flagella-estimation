from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def _load_module() -> object:
    path = Path("tools/codex/check_issue_script_lifecycle.py")
    spec = importlib.util.spec_from_file_location("issue_script_lifecycle", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_checker_rejects_active_reference(monkeypatch: pytest.MonkeyPatch) -> None:
    module = _load_module()
    monkeypatch.setattr(module, "deleted_scripts", lambda base: ["scripts/old.py"])
    monkeypatch.setattr(
        module, "active_references", lambda path: ["tests/test_old.py:1:scripts/old.py"]
    )
    with pytest.raises(SystemExit, match="active consumers"):
        module.main(["--base", "base"])


def test_checker_allows_historical_only_reference(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _load_module()
    monkeypatch.setattr(module, "deleted_scripts", lambda base: ["scripts/old.py"])
    monkeypatch.setattr(module, "active_references", lambda path: [])
    module.main(["--base", "base"])
