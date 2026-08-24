from pathlib import Path

import pytest

from sim_swim.analysis.video_transcode import (
    _update_video_entries,
    build_ffmpeg_command,
    transcode_tree,
)


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


def test_promotion_updates_nested_replay_video_manifest(tmp_path: Path) -> None:
    video = tmp_path / "grid.mp4"
    manifest = {"render_video": {"grid": {"pages": [{"video": {"path": str(video)}}]}}}

    _update_video_entries(
        manifest, {str(video.resolve()): {"source": str(video), "kind": "promotion"}}
    )

    video_entry = manifest["render_video"]["grid"]["pages"][0]["video"]
    assert video_entry["selected_codec"] == "libx264"
    assert video_entry["attempted_codecs"] == ["ffmpeg:libx264"]


def test_promotion_updates_stale_path_by_unique_video_name(tmp_path: Path) -> None:
    video = tmp_path / "grid.mp4"
    manifest = {"video": {"path": "/stale/cs10/grid.mp4"}}

    _update_video_entries(manifest, {str(video.resolve()): {"source": str(video)}})

    assert manifest["video"]["path"] == str(video)
    assert manifest["video"]["selected_codec"] == "libx264"
