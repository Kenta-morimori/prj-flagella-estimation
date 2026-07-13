#!/usr/bin/env python3
"""Run a generic Phase 2 multi-condition campaign."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

from sim_swim.analysis.flagella_count_behavior import (
    save_state_archive,
    write_trajectory_csv,
)
from sim_swim.analysis.multi_run_campaign import (
    apply_campaign_cli_overrides,
    build_campaign_conditions,
    campaign_axes_metadata,
    load_yaml,
    summary_axis_fields,
)
from sim_swim.analysis.sweeps.shape_stability_grid import (
    SUMMARY_FIELDS,
    Condition,
    _axis_center_phase_summary,
    _read_step_rows,
    _summary_row,
)
from sim_swim.core.run_context import init_run
from sim_swim.sim.core import Simulator
from sim_swim.sim.helix_retention_gate import summarize_single_flagellum_helix_retention
from sim_swim.sim.params import SimulationConfig


def _parse_args(argv: list[str] | None = None) -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-config", type=Path, required=True)
    parser.add_argument("--sample-limit", type=int, default=None)
    parser.add_argument("--progress-interval", type=int, default=None)
    parser.add_argument("--save-state-archive", type=str, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args, passthrough = parser.parse_known_args(argv)
    return args, passthrough


def _parse_bool(value: Any, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off"}:
        return False
    raise ValueError(f"Invalid boolean value: {value}")


def _condition_row(
    cfg: SimulationConfig,
    condition: dict[str, Any],
    condition_dir: Path,
) -> dict[str, Any]:
    rows = _read_step_rows(condition_dir / "step_summary.csv")
    axis_center_summary = _axis_center_phase_summary(
        _read_step_rows(condition_dir / "flag_helix_axis_diagnostics.csv")
    )
    helix_summary = summarize_single_flagellum_helix_retention(
        rows,
        min_steps=min(50, max(len(rows) - 1, 1)),
        min_net_abs_spin_revolutions=0.0,
    )
    if not rows:
        raise RuntimeError(f"No rows found in {condition_dir / 'step_summary.csv'}")
    shape_condition = Condition(
        condition_id=condition["condition_id"],
        mode="generic-multi-run",
        description=condition["condition_label"],
        scales={},
    )
    row = _summary_row(
        cfg,
        shape_condition,
        condition_dir,
        rows[-1],
        rows,
        helix_summary,
        axis_center_summary,
    )
    row.update(summary_axis_fields(condition))
    return row


def _manifest_condition_record(root: Path, condition: dict[str, Any]) -> dict[str, Any]:
    return {
        "condition_id": condition["condition_id"],
        "condition_index": condition["condition_index"],
        "condition_label": condition["condition_label"],
        "output_dir": str(root / condition["condition_id"]),
        "config_overrides": condition["config_overrides"],
        "axis_values": condition["axis_values"],
        "axis_labels": condition["axis_labels"],
        "axis_ids": condition["axis_ids"],
        "axis_order": condition["axis_order"],
        "grid": {
            "row_axis": condition.get("grid_row_axis"),
            "row_index": condition.get("grid_row_index"),
            "row_label": condition.get("grid_row_label"),
            "col_axis": condition.get("grid_col_axis"),
            "col_index": condition.get("grid_col_index"),
            "col_label": condition.get("grid_col_label"),
        },
    }


def _summary_fieldnames(rows: list[dict[str, Any]]) -> list[str]:
    base = list(SUMMARY_FIELDS)
    extras: list[str] = []
    seen = set(base)
    for row in rows:
        for key in row:
            if key not in seen:
                extras.append(key)
                seen.add(key)
    return [*base, *extras]


def run_campaign(argv: list[str] | None = None) -> Path:
    args, passthrough = _parse_args(argv)
    raw_campaign = load_yaml(args.campaign_config)
    campaign = apply_campaign_cli_overrides(raw_campaign, passthrough)
    conditions = build_campaign_conditions(campaign)
    if args.sample_limit is not None:
        conditions = conditions[: args.sample_limit]

    if args.dry_run:
        for condition in conditions:
            print(condition["condition_id"])
        return Path()

    output_base_dir = Path(
        str(campaign.get("output", {}).get("base_dir") or "outputs/phase2_multi_run")
    )
    output_cfg = dict(campaign.get("output", {}) or {})
    timestamp_subdir = _parse_bool(
        output_cfg.get("timestamp_subdir"),
        default=True,
    )
    overwrite = args.overwrite or bool(output_cfg.get("overwrite", False))
    ctx = init_run(
        base_dir=output_base_dir,
        timestamp_subdir=timestamp_subdir,
        overwrite=overwrite,
        input_info={
            "campaign_config": str(args.campaign_config),
            "overrides": passthrough,
            "sample_limit": args.sample_limit,
        },
    )
    logger = ctx.logger
    progress_interval = (
        args.progress_interval
        if args.progress_interval is not None
        else int(campaign.get("output", {}).get("progress_interval", 1000))
    )
    save_state_archive_enabled = _parse_bool(
        args.save_state_archive,
        default=bool(campaign.get("output", {}).get("save_state_archive", True)),
    )
    logger.info("Loaded generic multi-run config: %s", args.campaign_config)
    logger.info("Campaign conditions=%d", len(conditions))

    base_config_path = Path(str(campaign["base_config"]))
    base_cfg = load_yaml(base_config_path)
    rows: list[dict[str, Any]] = []
    total = len(conditions)
    for index, condition in enumerate(conditions, start=1):
        logger.info(
            "[%d/%d] generic_multi_run %s", index, total, condition["condition_id"]
        )
        condition_dir = ctx.out.root / condition["condition_id"]
        condition_dir.mkdir(parents=True, exist_ok=False)
        cfg = SimulationConfig.from_dict(base_cfg).with_overrides(
            condition["config_overrides"]
        )
        simulator = Simulator(cfg)
        states = simulator.run(
            cfg.time.duration_s,
            step_summary_dir=condition_dir,
            stop_on_shape_fail=False,
            progress_interval=progress_interval,
        )
        if save_state_archive_enabled:
            save_state_archive(condition_dir / "state_archive.npz", states)
            write_trajectory_csv(condition_dir / "trajectory.csv", states)
        rows.append(_condition_row(cfg, condition, condition_dir))

    summary_path = ctx.out.root / "summary.csv"
    fieldnames = _summary_fieldnames(rows)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    manifest = {
        "kind": "generic_multi_run",
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "campaign_config": str(args.campaign_config),
        "base_config": str(base_config_path),
        "summary_csv": str(summary_path),
        "git": {
            "commit": ctx.git.commit,
            "commit_short": ctx.git.commit_short,
            "branch": ctx.git.branch,
            "is_clean": ctx.git.is_clean,
        },
        "output_root": str(ctx.out.root),
        "save_state_archive": save_state_archive_enabled,
        "replay": dict(campaign.get("replay", {}) or {}),
        "plot": dict(campaign.get("plot", {}) or {}),
        "axes": campaign_axes_metadata(campaign),
        "condition_order": [condition["condition_id"] for condition in conditions],
        "conditions": [
            _manifest_condition_record(ctx.out.root, condition)
            for condition in conditions
        ],
    }
    run_manifest_path = ctx.out.root / "run_manifest.json"
    run_manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    manifest_path = ctx.out.root / "manifest.json"
    if manifest_path.is_file():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        existing["input"]["campaign"] = {
            "base_config": str(base_config_path),
            "replay": manifest["replay"],
            "plot": manifest["plot"],
        }
        existing["outputs"]["summary_csv"] = str(summary_path)
        existing["outputs"]["run_manifest_json"] = str(run_manifest_path)
        manifest_path.write_text(
            json.dumps(existing, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    logger.info("Wrote summary: %s", summary_path)
    logger.info("Wrote run manifest: %s", run_manifest_path)
    print(summary_path)
    return summary_path


def main(argv: list[str] | None = None) -> None:
    run_campaign(argv)


if __name__ == "__main__":
    main()
