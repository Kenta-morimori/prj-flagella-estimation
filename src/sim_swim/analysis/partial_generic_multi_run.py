"""Export completed conditions from an interrupted generic multi-run campaign.

The normal generic runner writes ``summary.csv`` and ``run_manifest.json`` only
after every requested condition has finished.  This helper creates a separate,
explicitly partial analysis input without modifying that live campaign root.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

from sim_swim.analysis.multi_run_campaign import (
    apply_campaign_cli_overrides,
    build_campaign_conditions,
    campaign_axes_metadata,
    load_yaml,
)
from sim_swim.analysis.sweeps.generic_multi_run import (
    _condition_row,
    _manifest_condition_record,
    _summary_fieldnames,
)
from sim_swim.sim.params import SimulationConfig


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _completed_summary(path: Path, expected_duration_s: float) -> bool:
    if not path.is_file():
        return False
    execution = dict(_read_json(path).get("execution", {}) or {})
    if execution.get("status") != "completed":
        return False
    observed = execution.get("observed_final_t_s")
    expected = execution.get("expected_final_step_summary_t_s", expected_duration_s)
    try:
        return abs(float(observed) - float(expected)) <= 1.0e-9
    except (TypeError, ValueError):
        return False


def export_completed_campaign(
    *, campaign_config: Path, run_dir: Path, output_dir: Path, overwrite: bool
) -> Path:
    """Write an analysis-only manifest containing only completed source runs."""
    if output_dir.exists() and not overwrite:
        raise FileExistsError(f"Output already exists: {output_dir}; pass --overwrite")
    if output_dir.exists():
        import shutil

        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    campaign = apply_campaign_cli_overrides(load_yaml(campaign_config), [])
    conditions = build_campaign_conditions(campaign)
    base_config_path = Path(str(campaign["base_config"]))
    base_cfg = load_yaml(base_config_path)

    completed: list[tuple[dict[str, Any], dict[str, Any]]] = []
    excluded: list[dict[str, str]] = []
    for condition in conditions:
        condition_dir = run_dir / condition["condition_id"]
        cfg = SimulationConfig.from_dict(base_cfg).with_overrides(
            condition["config_overrides"]
        )
        summary_path = condition_dir / "run_summary.json"
        archive_path = condition_dir / "state_archive.npz"
        trajectory_path = condition_dir / "trajectory.csv"
        hydro_path = condition_dir / "hydro_archive.npz"
        if not _completed_summary(summary_path, cfg.time.duration_s):
            excluded.append(
                {"condition_id": condition["condition_id"], "reason": "not_completed"}
            )
            continue
        if (
            not archive_path.is_file()
            or not trajectory_path.is_file()
            or (cfg.hydrodynamics.enabled and not hydro_path.is_file())
        ):
            excluded.append(
                {"condition_id": condition["condition_id"], "reason": "missing_archive"}
            )
            continue
        completed.append(
            (
                _condition_row(cfg, condition, condition_dir),
                _manifest_condition_record(
                    run_dir,
                    condition,
                    time_manifest=cfg.time_manifest(),
                    hydrodynamics_enabled=cfg.hydrodynamics.enabled,
                ),
            )
        )

    if not completed:
        raise ValueError(f"No completed conditions found under {run_dir}")
    completed.sort(
        key=lambda item: (
            int(item[1]["axis_values"].get("n_flagella", 0)),
            int(item[1]["axis_values"].get("phase_seed", 0)),
            int(item[1]["axis_values"].get("attach_seed", 0)),
        )
    )
    rows = [row for row, _ in completed]
    records = [record for _, record in completed]
    summary_path = output_dir / "summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_summary_fieldnames(rows))
        writer.writeheader()
        writer.writerows(rows)

    manifest = {
        "kind": "generic_multi_run_partial_analysis",
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "partial": True,
        "campaign_config": str(campaign_config),
        "base_config": str(base_config_path),
        "output_root": str(run_dir),
        "summary_csv": str(summary_path),
        "axes": campaign_axes_metadata(campaign),
        "condition_order": [record["condition_id"] for record in records],
        "conditions": records,
        "excluded_conditions": excluded,
    }
    manifest_path = output_dir / "run_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "kind": "partial_generic_multi_run_export",
                "source_run_dir": str(run_dir),
                "run_manifest_json": str(manifest_path),
                "completed_condition_count": len(records),
                "excluded_condition_count": len(excluded),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (output_dir / "run.log").write_text(
        f"completed_condition_count={len(records)}\nexcluded_condition_count={len(excluded)}\n",
        encoding="utf-8",
    )
    return output_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-config", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    print(
        export_completed_campaign(
            campaign_config=args.campaign_config,
            run_dir=args.run_dir,
            output_dir=args.output_dir,
            overwrite=args.overwrite,
        )
    )


if __name__ == "__main__":
    main()
