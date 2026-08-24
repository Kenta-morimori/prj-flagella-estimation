"""Reproducible H.264 transcode for rendered MP4 artifact trees."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
from typing import Any


def _probe(path: Path, ffprobe: str) -> dict[str, Any]:
    completed = subprocess.run(
        [
            ffprobe,
            "-v",
            "error",
            "-show_entries",
            "format=duration,format_name:stream=codec_name,codec_tag_string",
            "-of",
            "json",
            str(path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(completed.stdout)


def _codec(probe: dict[str, Any]) -> str | None:
    streams = probe.get("streams", [])
    return str(streams[0].get("codec_name")) if streams else None


def _duration_s(probe: dict[str, Any]) -> float:
    try:
        return float(probe["format"]["duration"])
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError("ffprobe did not report a usable video duration") from exc


def build_ffmpeg_command(ffmpeg: str, source: Path, destination: Path) -> list[str]:
    return [
        ffmpeg,
        "-y",
        "-i",
        str(source),
        "-map",
        "0:v:0",
        "-an",
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-movflags",
        "+faststart",
        str(destination),
    ]


def transcode_tree(
    input_dir: Path,
    output_dir: Path,
    *,
    ffmpeg: str = "ffmpeg",
    ffprobe: str = "ffprobe",
    overwrite: bool = False,
    dry_run: bool = False,
) -> Path:
    input_dir, output_dir = input_dir.resolve(), output_dir.resolve()
    if not input_dir.is_dir():
        raise FileNotFoundError(input_dir)
    if output_dir == input_dir or output_dir.is_relative_to(input_dir):
        raise ValueError("output_dir must be outside input_dir")
    if output_dir.exists() and not overwrite:
        raise FileExistsError(output_dir)
    if output_dir.exists() and overwrite:
        shutil.rmtree(output_dir)
    videos = sorted(input_dir.rglob("*.mp4"))
    if not videos:
        raise FileNotFoundError(f"No MP4 files below {input_dir}")
    records: list[dict[str, Any]] = []
    for source in videos:
        destination = output_dir / source.relative_to(input_dir)
        record: dict[str, Any] = {
            "source": str(source),
            "destination": str(destination),
            "command": build_ffmpeg_command(ffmpeg, source, destination),
        }
        if not dry_run:
            destination.parent.mkdir(parents=True, exist_ok=True)
            before = _probe(source, ffprobe)
            subprocess.run(
                record["command"], check=True, capture_output=True, text=True
            )
            after = _probe(destination, ffprobe)
            if _codec(after) != "h264":
                raise RuntimeError(f"Expected H.264 output: {destination}")
            if abs(_duration_s(before) - _duration_s(after)) > 0.05:
                raise RuntimeError(f"Duration changed during transcode: {destination}")
            record.update({"input_probe": before, "output_probe": after})
        records.append(record)
    if dry_run:
        return output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "transcode_manifest.json").write_text(
        json.dumps(
            {
                "kind": "h264_mp4_transcode",
                "input_dir": str(input_dir),
                "output_dir": str(output_dir),
                "video_count": len(records),
                "records": records,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return output_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--ffmpeg", default="ffmpeg")
    parser.add_argument("--ffprobe", default="ffprobe")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    print(
        transcode_tree(
            args.input_dir,
            args.output_dir,
            ffmpeg=args.ffmpeg,
            ffprobe=args.ffprobe,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
        )
    )
