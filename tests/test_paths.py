from __future__ import annotations

from pathlib import Path

import pytest

from flagella_estimation.core.io.paths import make_output_dirs


def test_make_output_dirs_creates_expected_structure(tmp_path: Path) -> None:
    out = make_output_dirs(tmp_path, "2026-01-21", "120000")
    assert out.root == tmp_path / "2026-01-21" / "120000"
    assert out.tracking_dir == out.root / "tracking"
    assert out.tracking_dir.is_dir()


def test_make_output_dirs_is_not_overwritten(tmp_path: Path) -> None:
    make_output_dirs(tmp_path, "2026-01-21", "120000")
    with pytest.raises(FileExistsError):
        make_output_dirs(tmp_path, "2026-01-21", "120000")
