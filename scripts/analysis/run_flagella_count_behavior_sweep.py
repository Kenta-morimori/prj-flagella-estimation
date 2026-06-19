#!/usr/bin/env python3
"""Run Phase 2.8 flagella-count behavior simulations as one batch."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo

import yaml
from tqdm.auto import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.flagella_count_behavior import (
    apply_analysis_cli_overrides,
    normalize_base_overrides,
    save_state_archive,
    write_trajectory_csv,
)
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _write_yaml(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


def _merge_nested(dst: dict[str, Any], src: dict[str, Any]) -> dict[str, Any]:
    merged = dict(dst)
    for key, value in src.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_nested(merged[key], value)
        else:
            merged[key] = value
    return merged


def _git_info() -> dict[str, Any]:
    def run(cmd: list[str]) -> str:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True).strip()

    try:
        commit = run(["git", "rev-parse", "HEAD"])
        return {
            "commit": commit,
            "commit_short": run(["git", "rev-parse", "--short", "HEAD"]),
            "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
            "is_clean": run(["git", "status", "--porcelain"]) == "",
        }
    except Exception:
        return {
            "commit": "unknown",
            "commit_short": "unknown",
            "branch": "unknown",
            "is_clean": False,
        }


def _now_jst() -> str:
    return datetime.now(ZoneInfo("Asia/Tokyo")).isoformat()


def _sample_id(n_flagella: int, seed: int) -> str:
    return f"nf{int(n_flagella):02d}_seed{int(seed):03d}"


def build_conditions(config: dict[str, Any]) -> list[dict[str, Any]]:
    sweep = config.get("sweep", {}) or {}
    n_flagella_values = [int(v) for v in sweep.get("n_flagella", [])]
    seeds = [int(v) for v in sweep.get("seeds", [])]
    conditions: list[dict[str, Any]] = []
    for n_flagella in n_flagella_values:
        for seed in seeds:
            sample_id = _sample_id(n_flagella, seed)
            conditions.append(
                {
                    "sample_id": sample_id,
                    "condition_tag": f"n_flagella={n_flagella},seed={seed}",
                    "n_flagella": n_flagella,
                    "seed": seed,
                }
            )
    return conditions


def _build_sample_config(
    *,
    analysis_config: dict[str, Any],
    base_config: dict[str, Any],
    condition: dict[str, Any],
) -> SimulationConfig:
    base_overrides = normalize_base_overrides(
        analysis_config.get("base_overrides", {}) or {}
    )
    sample_overrides = {
        "flagella": {"n_flagella": int(condition["n_flagella"])},
        "seed": {"global_seed": int(condition["seed"])},
    }
    overrides = _merge_nested(base_overrides, sample_overrides)
    return SimulationConfig.from_dict(base_config).with_overrides(overrides)


def _setup_sample_logger(log_path: Path) -> logging.Logger:
    logger = logging.getLogger(f"flagella_count_behavior.{log_path.parent.name}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(log_path, encoding="utf-8")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def _run_sample(
    *,
    cfg: SimulationConfig,
    condition: dict[str, Any],
    sample_dir: Path,
    raw_dir: Path,
    log_path: Path,
    stop_on_shape_fail: bool,
    progress_interval: int | None,
) -> dict[str, Any]:
    logger = _setup_sample_logger(log_path)
    logger.info("Sample start: %s", condition["sample_id"])
    logger.info("condition=%s", condition)
    logger.info("effective_config=%s", cfg)
    sim = Simulator(cfg)
    states = sim.run(
        cfg.time.duration_s,
        logger=logger,
        progress_interval=progress_interval,
        step_summary_dir=raw_dir,
        stop_on_shape_fail=stop_on_shape_fail,
    )
    archive_path = raw_dir / "state_archive.npz"
    trajectory_path = raw_dir / "trajectory.csv"
    save_state_archive(archive_path, states)
    write_trajectory_csv(trajectory_path, states)
    logger.info("Sample end: %s states=%d", condition["sample_id"], len(states))
    return {
        "status": "completed",
        "state_count": len(states),
        "outputs": {
            "raw_dir": str(raw_dir),
            "step_summary_csv": str(raw_dir / "step_summary.csv"),
            "flag_helix_axis_diagnostics_csv": str(
                raw_dir / "flag_helix_axis_diagnostics.csv"
            ),
            "initial_geometry_summary_json": str(
                raw_dir / "initial_geometry_summary.json"
            ),
            "trajectory_csv": str(trajectory_path),
            "state_archive_npz": str(archive_path),
            "run_log": str(log_path),
        },
    }


def run_batch(
    *,
    analysis_config_path: Path,
    dry_run: bool,
    overwrite: bool,
    stop_on_shape_fail: bool,
    sample_limit: int | None,
    progress_interval: int | None,
    cli_overrides: list[str] | None = None,
) -> Path:
    raw_analysis_config = _load_yaml(analysis_config_path)
    cli_overrides = cli_overrides or []
    analysis_config = apply_analysis_cli_overrides(
        raw_analysis_config,
        cli_overrides,
    )
    base_config_path = Path(analysis_config.get("base_config", "conf/sim_swim.yaml"))
    base_config = _load_yaml(base_config_path)

    dataset_id = str(analysis_config["dataset_id"])
    run_batch_id = str(analysis_config.get("run_batch_id", dataset_id))
    output_cfg = analysis_config.get("output", {}) or {}
    run_batch_dir = Path(output_cfg["run_batch_dir"])
    configs_dir = run_batch_dir / "configs"
    samples_root = run_batch_dir / "samples"

    conditions = build_conditions(analysis_config)
    if sample_limit is not None:
        conditions = conditions[: int(sample_limit)]

    run_batch_dir.mkdir(parents=True, exist_ok=True)
    configs_dir.mkdir(parents=True, exist_ok=True)
    samples_root.mkdir(parents=True, exist_ok=True)
    effective_analysis_config_path = run_batch_dir / "analysis_config_used.yaml"
    _write_yaml(effective_analysis_config_path, analysis_config)

    manifest_samples: list[dict[str, Any]] = []
    for condition in tqdm(conditions, desc="flagella sweep", unit="sample"):
        sample_id = str(condition["sample_id"])
        sample_dir = samples_root / sample_id
        raw_dir = sample_dir / "raw"
        log_path = sample_dir / "run.log"
        sample_config_path = configs_dir / f"{sample_id}.yaml"
        cfg = _build_sample_config(
            analysis_config=analysis_config,
            base_config=base_config,
            condition=condition,
        )
        _write_yaml(sample_config_path, asdict(cfg))

        sample_record: dict[str, Any] = {
            **condition,
            "dataset_id": dataset_id,
            "run_batch_id": run_batch_id,
            "config_path": str(sample_config_path),
            "sample_dir": str(sample_dir),
            "raw_dir": str(raw_dir),
            "run_log": str(log_path),
            "duration_s": float(cfg.time.duration_s),
            "dt_star": float(cfg.dt_star),
            "torque_Nm": float(cfg.motor.torque_Nm),
            "force_distribution": cfg.motor.force_distribution,
            "cli_overrides": list(cli_overrides),
        }

        if dry_run:
            sample_record["status"] = "planned"
            manifest_samples.append(sample_record)
            continue

        if sample_dir.exists() and overwrite:
            shutil.rmtree(sample_dir)
        sample_dir.mkdir(parents=True, exist_ok=True)
        raw_dir.mkdir(parents=True, exist_ok=True)

        step_summary = raw_dir / "step_summary.csv"
        if step_summary.exists() and not overwrite:
            sample_record["status"] = "skipped_existing"
            sample_record["outputs"] = {
                "raw_dir": str(raw_dir),
                "step_summary_csv": str(step_summary),
                "flag_helix_axis_diagnostics_csv": str(
                    raw_dir / "flag_helix_axis_diagnostics.csv"
                ),
                "initial_geometry_summary_json": str(
                    raw_dir / "initial_geometry_summary.json"
                ),
                "trajectory_csv": str(raw_dir / "trajectory.csv"),
                "state_archive_npz": str(raw_dir / "state_archive.npz"),
                "run_log": str(log_path),
            }
            manifest_samples.append(sample_record)
            continue

        try:
            sample_record.update(
                _run_sample(
                    cfg=cfg,
                    condition=condition,
                    sample_dir=sample_dir,
                    raw_dir=raw_dir,
                    log_path=log_path,
                    stop_on_shape_fail=stop_on_shape_fail,
                    progress_interval=progress_interval,
                )
            )
        except Exception as exc:
            sample_record["status"] = "failed"
            sample_record["error"] = repr(exc)
        manifest_samples.append(sample_record)

    manifest = {
        "run_batch_id": run_batch_id,
        "dataset_id": dataset_id,
        "created_at": _now_jst(),
        "analysis_config": str(analysis_config_path),
        "effective_analysis_config": analysis_config,
        "effective_analysis_config_yaml": str(effective_analysis_config_path),
        "base_config": str(base_config_path),
        "feature_schema": str(
            analysis_config.get(
                "feature_schema",
                "conf/analysis/flagella_count_behavior_features.yaml",
            )
        ),
        "base_overrides": analysis_config.get("base_overrides", {}) or {},
        "cli_overrides": list(cli_overrides),
        "sweep": analysis_config.get("sweep", {}) or {},
        "output": {
            "run_batch_dir": str(run_batch_dir),
            "dataset_dir": str(output_cfg.get("dataset_dir", "")),
            "configs_dir": str(configs_dir),
            "samples_dir": str(samples_root),
        },
        "dry_run": bool(dry_run),
        "git": _git_info(),
        "samples": manifest_samples,
    }
    manifest_path = run_batch_dir / "run_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("conf/analysis/flagella_count_behavior_dataset.yaml"),
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--stop-on-shape-fail", action="store_true")
    parser.add_argument("--sample-limit", type=int, default=None)
    parser.add_argument("--progress-interval", type=int, default=None)
    parser.add_argument(
        "overrides",
        nargs="*",
        help=(
            "Optional analysis/simulation overrides as key=value "
            "(e.g. dataset_id=my_dataset time.duration_s=0.25)"
        ),
    )
    args = parser.parse_args()

    manifest_path = run_batch(
        analysis_config_path=args.config,
        dry_run=bool(args.dry_run),
        overwrite=bool(args.overwrite),
        stop_on_shape_fail=bool(args.stop_on_shape_fail),
        sample_limit=args.sample_limit,
        progress_interval=args.progress_interval,
        cli_overrides=args.overrides,
    )
    print(manifest_path)


if __name__ == "__main__":
    main()
