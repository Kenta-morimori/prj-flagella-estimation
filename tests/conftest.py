from __future__ import annotations

import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@pytest.fixture
def fake_h264_writer(monkeypatch: pytest.MonkeyPatch) -> None:
    """Use a file-producing writer for tests that do not test FFmpeg itself."""

    class Writer:
        def __init__(self, path: Path, **_: object) -> None:
            self.path = path
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.touch()

        def isOpened(self) -> bool:
            return True

        def write(self, frame: object) -> None:
            del frame

        def release(self) -> None:
            return None

    monkeypatch.setattr("sim_swim.render.video_writer._FFmpegVideoWriter", Writer)
