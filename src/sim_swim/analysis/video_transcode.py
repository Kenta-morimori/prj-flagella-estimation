"""Reproducible H.264 transcode for rendered MP4 artifact trees."""

from __future__ import annotations

import argparse
import json
import os
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


def _update_video_entries(value: Any, promoted: dict[str, dict[str, Any]]) -> None:
    if isinstance(value, list):
        for item in value:
            _update_video_entries(item, promoted)
    elif isinstance(value, dict):
        path = value.get("path")
        key = str(Path(path).resolve()) if isinstance(path, str) else ""
        record = promoted.get(key)
        if record is None and isinstance(path, str):
            matches = [
                candidate
                for candidate_key, candidate in promoted.items()
                if Path(candidate_key).name == Path(path).name
            ]
            if len(matches) == 1:
                record = matches[0]
        if record is not None:
            value["path"] = str(record["source"])
            value["selected_codec"] = "libx264"
            value["attempted_codecs"] = ["ffmpeg:libx264"]
        for item in value.values():
            _update_video_entries(item, promoted)


def promote_tree_in_place(
    input_dir: Path,
    staging_dir: Path,
    *,
    manifest_json: Path | None = None,
    ffmpeg: str = "ffmpeg",
    ffprobe: str = "ffprobe",
    overwrite: bool = False,
) -> Path:
    """Replace MP4s with verified H.264 copies and update an optional manifest."""
    input_dir, staging_dir = input_dir.resolve(), staging_dir.resolve()
    transcode_tree(
        input_dir,
        staging_dir,
        ffmpeg=ffmpeg,
        ffprobe=ffprobe,
        overwrite=overwrite,
    )
    transcode_manifest = json.loads(
        (staging_dir / "transcode_manifest.json").read_text(encoding="utf-8")
    )
    promoted: dict[str, dict[str, Any]] = {}
    for record in transcode_manifest["records"]:
        source = Path(str(record["source"]))
        destination = Path(str(record["destination"]))
        os.replace(destination, source)
        promoted[str(source.resolve())] = record
    shutil.rmtree(staging_dir)
    if manifest_json is not None:
        data = json.loads(manifest_json.read_text(encoding="utf-8"))
        _update_video_entries(data, promoted)
        manifest_json.write_text(
            json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
    promotion_manifest = input_dir / "h264_promotion_manifest.json"
    promotion_manifest.write_text(
        json.dumps(
            {
                "kind": "h264_mp4_promotion",
                "input_dir": str(input_dir),
                "video_count": len(promoted),
                "records": list(promoted.values()),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return promotion_manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--manifest-json", type=Path)
    parser.add_argument("--ffmpeg", default="ffmpeg")
    parser.add_argument("--ffprobe", default="ffprobe")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.replace:
        if args.dry_run:
            parser.error("--replace cannot be combined with --dry-run")
        print(
            promote_tree_in_place(
                args.input_dir,
                args.output_dir,
                manifest_json=args.manifest_json,
                ffmpeg=args.ffmpeg,
                ffprobe=args.ffprobe,
                overwrite=args.overwrite,
            )
        )
    else:
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
