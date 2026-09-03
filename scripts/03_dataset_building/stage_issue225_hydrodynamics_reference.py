#!/usr/bin/env python3
"""Stage the complete 2-second hydrodynamics campaign as a reusable reference."""

from __future__ import annotations
import argparse
import json
import shutil
import tempfile
from pathlib import Path
from typing import Any

ROOT_FILES = ("run_manifest.json", "summary.csv", "campaign_completion.json", "run.log")


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _condition(root: Path, record: dict[str, Any]) -> Path:
    raw = Path(str(record.get("output_dir", "")))
    for p in (raw, root / raw.name, root / "conditions" / raw.name):
        if p.is_dir():
            return p
    return root / "conditions" / str(record["condition_id"])


def stage(campaign_dir: Path, reference_dir: Path, *, overwrite: bool = False) -> Path:
    manifest = _json(campaign_dir / "run_manifest.json")
    records = manifest.get("conditions", [])
    if len(records) != 36:
        raise ValueError(f"expected 36 conditions, found {len(records)}")
    missing = [
        str(campaign_dir / name)
        for name in ROOT_FILES
        if not (campaign_dir / name).is_file()
    ]
    for r in records:
        d = _condition(campaign_dir, r)
        missing += [
            str(d / name)
            for name in ("run_summary.json", "state_archive.npz", "hydro_archive.npz")
            if not (d / name).is_file()
        ]
    if missing:
        raise FileNotFoundError(
            "missing required hydrodynamics artifacts: " + ", ".join(missing)
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
        for name in ROOT_FILES:
            shutil.copy2(campaign_dir / name, staging / name)
        for r in records:
            src = _condition(campaign_dir, r)
            dst = staging / "conditions" / str(r["condition_id"])
            dst.mkdir(parents=True)
            for name in ("run_summary.json", "state_archive.npz", "hydro_archive.npz"):
                shutil.copy2(src / name, dst / name)
        (staging / "reference_manifest.json").write_text(
            json.dumps(
                {
                    "kind": "phase2_issue225_2s_hydrodynamics_reference",
                    "source_campaign": str(campaign_dir),
                    "condition_count": 36,
                    "copied_artifacts": [
                        *ROOT_FILES,
                        "conditions/*/{run_summary.json,state_archive.npz,hydro_archive.npz}",
                    ],
                    "raw_hydrodynamics_archives_included": True,
                },
                ensure_ascii=False,
                indent=2,
            )
            + "\n"
        )
        if reference_dir.exists():
            shutil.rmtree(reference_dir)
        staging.replace(reference_dir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return reference_dir


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign-dir", type=Path, required=True)
    p.add_argument("--reference-dir", type=Path, required=True)
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args(argv)
    print(stage(a.campaign_dir, a.reference_dir, overwrite=a.overwrite))


if __name__ == "__main__":
    main()
