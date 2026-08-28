from pathlib import Path

import numpy as np

from sim_swim.render.video_writer import _FFmpegVideoWriter


def test_ffmpeg_writer_pads_odd_dimensions_to_yuv420p(
    monkeypatch, tmp_path: Path
) -> None:
    """Odd-size replay frames must not make the H.264 subprocess exit early."""

    captured: list[str] = []

    def fake_popen(command, **_kwargs):
        captured.extend(command)

        class RunningProcess:
            def __init__(self) -> None:
                import io

                self.stdin = io.BytesIO()

            def poll(self):
                return None

            def wait(self):
                return 0

        return RunningProcess()

    monkeypatch.setattr("sim_swim.render.video_writer.resolve_ffmpeg", lambda: "ffmpeg")
    monkeypatch.setattr("sim_swim.render.video_writer.subprocess.Popen", fake_popen)
    writer = _FFmpegVideoWriter(tmp_path / "odd.mp4", fps=25.0, frame_size=(1439, 960))
    writer.write(np.zeros((960, 1439, 3), dtype=np.uint8))
    writer.release()

    assert "pad=ceil(iw/2)*2:ceil(ih/2)*2" in captured
