#!/usr/bin/env python3
"""Stage any completed hydrodynamics multi-run as a portable reference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import tempfile
from typing import Any


ROOT_FILES = ("run_manifest.json", "summary.csv", "campaign_completion.json")
CONDITION_FILES = ("run_summary.json", "state_archive.npz", "hydro_archive.npz")


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _condition_dir(root: Path, record: dict[str, Any]) -> Path:
    source = Path(str(record.get("output_dir", "")))
    for candidate in (source, root / source.name, root / "conditions" / source.name):
        if candidate.is_dir():
            return candidate
    return root / "conditions" / str(record["condition_id"])


def stage(campaign_dir: Path, reference_dir: Path, *, overwrite: bool = False) -> Path:
    manifest_path = campaign_dir / "run_manifest.json"
    manifest = _json(manifest_path)
    records = list(manifest.get("conditions", []))
    if not records:
        raise ValueError("run_manifest.json has no conditions")
    required_root = [manifest_path]
    required_root.extend(
        campaign_dir / name
        for name in ROOT_FILES[1:]
        if (campaign_dir / name).is_file()
    )
    missing = [str(path) for path in required_root if not path.is_file()]
    for record in records:
        source = _condition_dir(campaign_dir, record)
        missing.extend(
            str(source / name)
            for name in CONDITION_FILES
            if not (source / name).is_file()
        )
    if missing:
        raise FileNotFoundError(
            "missing hydrodynamics reference inputs: " + ", ".join(missing)
        )
    if reference_dir.exists() and not overwrite:
        raise FileExistsError(reference_dir)
    reference_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{reference_dir.name}.staging-", dir=reference_dir.parent
        )
    )
    try:
        for source in required_root:
            shutil.copy2(source, staging / source.name)
        for record in records:
            source = _condition_dir(campaign_dir, record)
            target = staging / "conditions" / str(record["condition_id"])
            target.mkdir(parents=True)
            for name in CONDITION_FILES:
                shutil.copy2(source / name, target / name)
        (staging / "reference_manifest.json").write_text(
            json.dumps(
                {
                    "kind": "phase2_hydrodynamics_reference",
                    "source_campaign": str(campaign_dir),
                    "condition_count": len(records),
                    "copied_artifacts": [
                        *[path.name for path in required_root],
                        "conditions/*/{run_summary.json,state_archive.npz,hydro_archive.npz}",
                    ],
                    "raw_hydrodynamics_archives_included": True,
                },
                ensure_ascii=False,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        if reference_dir.exists():
            shutil.rmtree(reference_dir)
        staging.replace(reference_dir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return reference_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    print(stage(args.campaign_dir, args.reference_dir, overwrite=args.overwrite))


if __name__ == "__main__":
    main()
