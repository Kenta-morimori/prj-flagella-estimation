from pathlib import Path

import pytest

from sim_swim.analysis.video_transcode import build_ffmpeg_command, transcode_tree


def test_h264_transcode_command_is_quicktime_compatible(tmp_path: Path) -> None:
    command = build_ffmpeg_command(
        "ffmpeg", tmp_path / "input.mp4", tmp_path / "output.mp4"
    )

    assert ["-c:v", "libx264"] == command[
        command.index("-c:v") : command.index("-c:v") + 2
    ]
    assert ["-pix_fmt", "yuv420p"] == command[
        command.index("-pix_fmt") : command.index("-pix_fmt") + 2
    ]
    assert ["-movflags", "+faststart"] == command[
        command.index("-movflags") : command.index("-movflags") + 2
    ]


def test_transcode_rejects_output_below_input(tmp_path: Path) -> None:
    source = tmp_path / "source"
    source.mkdir()
    (source / "movie.mp4").touch()

    with pytest.raises(ValueError, match="outside input_dir"):
        transcode_tree(source, source / "h264", dry_run=True)
