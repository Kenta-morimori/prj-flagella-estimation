from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def _load_replay_cli() -> object:
    path = Path("scripts/03_dataset_building/replay_dataset.py")
    spec = importlib.util.spec_from_file_location("common_replay_cli", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    ("source_option", "view", "expected"),
    [
        ("--dataset-dir", "2d", "--no-render-3d"),
        ("--dataset-dir", "3d", "--no-render-2d"),
        ("--run-dir", "2d", "2d"),
        ("--run-dir", "3d+2d", "3d+2d"),
    ],
)
def test_replay_cli_propagates_view(
    monkeypatch: pytest.MonkeyPatch, source_option: str, view: str, expected: str
) -> None:
    module = _load_replay_cli()
    captured: list[str] = []
    if source_option == "--dataset-dir":
        monkeypatch.setattr(
            module, "replay_dataset", lambda argv: captured.extend(argv)
        )
    else:
        monkeypatch.setattr(module, "replay_run", lambda argv: captured.extend(argv))

    module.main([source_option, "/tmp/input", "--view", view])

    assert expected in captured
    if source_option == "--run-dir":
        assert captured[captured.index("--view") + 1] == view
