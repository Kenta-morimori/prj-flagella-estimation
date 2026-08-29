#!/usr/bin/env python3
"""Stage the compact Issue #215 campaign artifacts into a local reference."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
from typing import Any


ROOT_FILES = ("run_manifest.json", "summary.csv", "campaign_completion.json", "run.log")


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _condition_dir(root: Path, record: dict[str, Any]) -> Path:
    candidate = Path(str(record.get("output_dir", "")))
    for path in (
        candidate,
        root / candidate.name,
        root / "conditions" / candidate.name,
    ):
        if path.is_dir():
            return path
    return root / "conditions" / str(record["condition_id"])


def _strict_status(summary: dict[str, Any]) -> tuple[bool, str]:
    if summary.get("execution", {}).get("status") != "completed":
        return False, "execution_not_completed"
    for name in ("finite", "shape_nonbody", "shape_body"):
        gate = summary.get("gates", {}).get(name, {})
        if gate.get("status") == "available" and (
            gate.get("any_fail") or not gate.get("final_pass")
        ):
            return False, f"{name}_qc_failed"
    return True, "pass"


def stage(
    campaign_dir: Path,
    reference_dir: Path,
    *,
    raw_campaign_root: str | None = None,
    overwrite: bool = False,
) -> Path:
    manifest_path = campaign_dir / "run_manifest.json"
    analysis = campaign_dir / "analysis" / "motion_features"
    if (
        not manifest_path.is_file()
        or not (campaign_dir / "summary.csv").is_file()
        or not analysis.is_dir()
    ):
        raise FileNotFoundError(
            "campaign requires run_manifest.json, summary.csv, and analysis/motion_features"
        )
    manifest = _read_json(manifest_path)
    conditions = list(manifest.get("conditions", []))
    if len(conditions) != 36:
        raise ValueError(f"expected 36 conditions, found {len(conditions)}")
    if reference_dir.exists():
        if not overwrite:
            raise FileExistsError(reference_dir)
        shutil.rmtree(reference_dir)
    reference_dir.mkdir(parents=True)
    for name in ROOT_FILES:
        if (source := campaign_dir / name).is_file():
            shutil.copy2(source, reference_dir / name)
    shutil.copytree(analysis, reference_dir / "analysis" / "motion_features")
    qc_rows: list[dict[str, Any]] = []
    for record in conditions:
        condition_id = str(record["condition_id"])
        source = _condition_dir(campaign_dir, record) / "run_summary.json"
        if not source.is_file():
            raise FileNotFoundError(f"{condition_id}: missing run_summary.json")
        strict_pass, strict_reason = _strict_status(_read_json(source))
        target = reference_dir / "conditions" / condition_id / "run_summary.json"
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        axis = record.get("axis_values", {})
        qc_rows.append(
            {
                "condition_id": condition_id,
                "n_flagella": axis.get("n_flagella"),
                "attach_seed": axis.get("attach_seed"),
                "phase_seed": axis.get("phase_seed"),
                "strict_pass": strict_pass,
                "strict_reason": strict_reason,
            }
        )
    with (reference_dir / "qc_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qc_rows[0]))
        writer.writeheader()
        writer.writerows(qc_rows)
    provenance = {
        "kind": "phase2_issue215_compact_reference",
        "raw_campaign_root": raw_campaign_root or str(campaign_dir),
        "raw_archive_transferred": False,
        "condition_count": len(conditions),
        "run_summary_count": len(qc_rows),
        "strict_pass_count": sum(bool(row["strict_pass"]) for row in qc_rows),
        "strict_nonpass_count": sum(not bool(row["strict_pass"]) for row in qc_rows),
        "git": manifest.get("git"),
        "copied_artifacts": [
            *ROOT_FILES,
            "conditions/*/run_summary.json",
            "analysis/motion_features",
            "qc_summary.csv",
        ],
    }
    (reference_dir / "reference_manifest.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return reference_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument(
        "--raw-campaign-root",
        help="Original NAS campaign path; recorded as provenance without copying raw artifacts.",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    print(
        stage(
            args.campaign_dir,
            args.reference_dir,
            raw_campaign_root=args.raw_campaign_root,
            overwrite=args.overwrite,
        )
    )


if __name__ == "__main__":
    main()
