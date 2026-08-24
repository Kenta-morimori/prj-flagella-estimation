"""Canonical on-disk layout migration for the two Issue #203 references."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any


_CONDITION_ID = re.compile(r"^as(?P<attach>\d{3})__ps(?P<phase>\d{3})__nf(?P<n>\d{2})$")


@dataclass(frozen=True)
class MigrationResult:
    root: Path
    profile: str
    condition_count: int
    composite_count: int
    applied: bool


def canonical_artifact_name(condition_id: str) -> str:
    """Return the display/path name without changing the logical condition ID."""
    match = _CONDITION_ID.fullmatch(condition_id)
    if not match:
        raise ValueError(f"Unsupported Issue #203 condition_id: {condition_id}")
    return "nf{n}_as{attach}_ps{phase}".format(**match.groupdict())


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def _probe_video(path: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [
            "ffprobe",
            "-v",
            "error",
            "-count_frames",
            "-show_entries",
            "stream=codec_name,pix_fmt,nb_read_frames",
            "-of",
            "json",
            str(path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    streams = json.loads(completed.stdout).get("streams", [])
    if len(streams) != 1:
        raise ValueError(f"Expected exactly one video stream: {path}")
    stream = streams[0]
    if stream.get("codec_name") != "h264" or stream.get("pix_fmt") != "yuv420p":
        raise ValueError(f"Expected H.264/yuv420p composite: {path}")
    frames = int(stream.get("nb_read_frames") or 0)
    if frames <= 0:
        raise ValueError(f"Expected decodable composite frames: {path}")
    return {"codec": "h264", "pix_fmt": "yuv420p", "frame_count": frames}


def _condition_records(root: Path) -> list[dict[str, Any]]:
    manifest = _read_json(root / "run_manifest.json")
    records = list(manifest.get("conditions") or [])
    if not records:
        raise ValueError(f"Missing conditions in {root / 'run_manifest.json'}")
    names = [canonical_artifact_name(str(record["condition_id"])) for record in records]
    if len(names) != len(set(names)):
        raise ValueError("Canonical artifact names are not unique")
    return records


def _validate_condition_sources(root: Path, records: list[dict[str, Any]]) -> None:
    for record in records:
        condition_id = str(record["condition_id"])
        source = root / condition_id
        destination = root / canonical_artifact_name(condition_id)
        if not source.is_dir() and not destination.is_dir():
            raise FileNotFoundError(f"Missing condition artifact: {source}")
        if source.is_dir() and destination.is_dir():
            raise FileExistsError(f"Both source and destination exist: {condition_id}")


def _update_run_manifests(root: Path, records: list[dict[str, Any]]) -> None:
    expected = {
        str(record["condition_id"]): canonical_artifact_name(
            str(record["condition_id"])
        )
        for record in records
    }
    for name in ("run_manifest.json", "manifest.json"):
        path = root / name
        data = _read_json(path)
        for record in data.get("conditions") or []:
            condition_id = str(record["condition_id"])
            if condition_id in expected:
                record["output_dir"] = expected[condition_id]
        _write_json(path, data)


def _rename_conditions(root: Path, records: list[dict[str, Any]]) -> None:
    for record in records:
        condition_id = str(record["condition_id"])
        source = root / condition_id
        destination = root / canonical_artifact_name(condition_id)
        if source.is_dir():
            source.rename(destination)
    _update_run_manifests(root, records)


def _composite_source(root: Path, profile: str) -> Path:
    return (
        root
        / "analysis"
        / ("composite_h264" if profile == "uniform" else "issue203_composite")
    )


def _composite_condition_ids(records: list[dict[str, Any]]) -> list[str]:
    return [
        str(record["condition_id"])
        for record in records
        if int((record.get("axis_values") or {}).get("n_flagella", -1)) <= 3
    ]


def _validate_composites(
    root: Path, profile: str, records: list[dict[str, Any]]
) -> list[tuple[str, dict[str, Any]]]:
    source = _composite_source(root, profile)
    expected = _composite_condition_ids(records)
    if not source.is_dir():
        raise FileNotFoundError(source)
    results: list[tuple[str, dict[str, Any]]] = []
    for condition_id in expected:
        path = source / f"{condition_id}_composite.mp4"
        if not path.is_file():
            raise FileNotFoundError(path)
        results.append((condition_id, _probe_video(path)))
    if len(list(source.glob("*_composite.mp4"))) != len(expected):
        raise ValueError(f"Unexpected composite video count: {source}")
    return results


def _uniform_composite_manifests(
    root: Path, condition_ids: list[str]
) -> dict[str, dict[str, Any]]:
    legacy = root / "analysis" / "composite"
    manifests: dict[str, dict[str, Any]] = {}
    for condition_id in condition_ids:
        path = legacy / f"{condition_id}_composite_manifest.json"
        if not path.is_file():
            raise FileNotFoundError(path)
        manifests[condition_id] = _read_json(path)
    return manifests


def _move_composites(
    root: Path,
    profile: str,
    results: list[tuple[str, dict[str, Any]]],
) -> None:
    analysis = root / "analysis"
    source = _composite_source(root, profile)
    target = analysis / "composite"
    condition_ids = [condition_id for condition_id, _ in results]
    legacy_manifests = (
        _uniform_composite_manifests(root, condition_ids)
        if profile == "uniform"
        else {}
    )
    if target.is_dir() and target != source:
        if profile != "uniform":
            raise FileExistsError(target)
        legacy = analysis / ".issue203_legacy_composite_mp4v"
        if legacy.exists():
            raise FileExistsError(legacy)
        target.rename(legacy)
    source.rename(target)
    (target / "transcode_manifest.json").unlink(missing_ok=True)
    composite_entries: list[dict[str, Any]] = []
    for condition_id, probe in results:
        name = canonical_artifact_name(condition_id)
        old_video = target / f"{condition_id}_composite.mp4"
        video = target / f"{name}_composite.mp4"
        old_video.rename(video)
        old_manifest = target / f"{condition_id}_composite_manifest.json"
        manifest = target / f"{name}_composite_manifest.json"
        data = legacy_manifests.get(condition_id)
        if data is None and old_manifest.is_file():
            data = _read_json(old_manifest)
        if data is None:
            raise FileNotFoundError(old_manifest)
        old_manifest.unlink(missing_ok=True)
        data.update(
            {
                "movie": str(video),
                "selected_codec": "libx264",
                "attempted_codecs": ["ffmpeg:libx264"],
                "frame_count": probe["frame_count"],
            }
        )
        _write_json(manifest, data)
        composite_entries.append(
            {
                "condition_id": condition_id,
                "artifact_name": name,
                "movie": str(video),
                **probe,
            }
        )
    _write_json(
        target / "manifest.json",
        {
            "kind": "phase2_issue203_composite_collection",
            "profile": profile,
            "artifact_name_format": "nfNN_asNNN_psNNN",
            "condition_count": len(composite_entries),
            "conditions": composite_entries,
        },
    )
    if profile == "uniform":
        shutil.rmtree(analysis / ".issue203_legacy_composite_mp4v")


def migrate_reference(
    root: Path, profile: str, *, apply: bool = False
) -> MigrationResult:
    if profile not in {"uniform", "diffusive"}:
        raise ValueError("profile must be uniform or diffusive")
    root = root.resolve()
    records = _condition_records(root)
    _validate_condition_sources(root, records)
    results = _validate_composites(root, profile, records)
    if apply:
        _rename_conditions(root, records)
        _move_composites(root, profile, results)
    return MigrationResult(
        root=root,
        profile=profile,
        condition_count=len(records),
        composite_count=len(results),
        applied=apply,
    )


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uniform-root", type=Path, required=True)
    parser.add_argument("--diffusive-root", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args(argv)
    results = [
        migrate_reference(args.uniform_root, "uniform", apply=args.apply),
        migrate_reference(args.diffusive_root, "diffusive", apply=args.apply),
    ]
    print(
        json.dumps(
            [result.__dict__ | {"root": str(result.root)} for result in results],
            ensure_ascii=False,
            indent=2,
        )
    )
