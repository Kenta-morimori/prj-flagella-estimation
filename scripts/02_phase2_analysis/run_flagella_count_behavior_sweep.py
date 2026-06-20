#!/usr/bin/env python3
"""Run Phase 2.8 flagella-count behavior simulations as one batch."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime
import hashlib
import json
import logging
from pathlib import Path
import shutil
import subprocess
import sys
import time
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


def _canonical_data(data: dict[str, Any]) -> dict[str, Any]:
    dumped = yaml.safe_dump(data, sort_keys=True)
    return yaml.safe_load(dumped) or {}


def _config_fingerprint(data: dict[str, Any]) -> str:
    canonical = json.dumps(
        _canonical_data(data),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _existing_sample_config_path(
    *,
    sample_config_used_path: Path,
    batch_config_path: Path,
) -> Path | None:
    if sample_config_used_path.is_file():
        return sample_config_used_path
    if batch_config_path.is_file():
        return batch_config_path
    return None


def _ensure_existing_config_matches(
    *,
    existing_config_path: Path | None,
    expected_config: dict[str, Any],
    sample_id: str,
    step_summary: Path,
) -> None:
    if existing_config_path is None:
        raise RuntimeError(
            f"Existing raw output for {sample_id} cannot be reused because "
            f"{step_summary} exists but no previous sample config was found. "
            "Use --overwrite or choose a new run_batch_id/output.run_batch_dir."
        )

    existing_config = _load_yaml(existing_config_path)
    if _config_fingerprint(existing_config) != _config_fingerprint(expected_config):
        raise RuntimeError(
            f"Existing raw output for {sample_id} was created with a different "
            f"sample config ({existing_config_path}). Use --overwrite or choose "
            "a new run_batch_id/output.run_batch_dir."
        )


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


def _sample_id_split(n_flagella: int, attach_seed: int, phase_seed: int) -> str:
    return f"nf{int(n_flagella):02d}_as{int(attach_seed):03d}_ps{int(phase_seed):03d}"


def build_conditions(config: dict[str, Any]) -> list[dict[str, Any]]:
    sweep = config.get("sweep", {}) or {}
    n_flagella_values = [int(v) for v in sweep.get("n_flagella", [])]
    conditions: list[dict[str, Any]] = []
    if "attach_seeds" in sweep or "phase_seeds" in sweep:
        fallback_seeds = [int(v) for v in sweep.get("seeds", [0])]
        attach_seeds = [int(v) for v in sweep.get("attach_seeds", fallback_seeds)]
        phase_seeds = [int(v) for v in sweep.get("phase_seeds", fallback_seeds)]
        for n_flagella in n_flagella_values:
            for attach_seed in attach_seeds:
                for phase_seed in phase_seeds:
                    sample_id = _sample_id_split(n_flagella, attach_seed, phase_seed)
                    conditions.append(
                        {
                            "sample_id": sample_id,
                            "condition_tag": (
                                f"n_flagella={n_flagella},"
                                f"attach_seed={attach_seed},"
                                f"phase_seed={phase_seed}"
                            ),
                            "n_flagella": n_flagella,
                            "seed": attach_seed,
                            "attach_seed": attach_seed,
                            "phase_seed": phase_seed,
                        }
                    )
    else:
        seeds = [int(v) for v in sweep.get("seeds", [])]
        for n_flagella in n_flagella_values:
            for seed in seeds:
                sample_id = _sample_id(n_flagella, seed)
                conditions.append(
                    {
                        "sample_id": sample_id,
                        "condition_tag": f"n_flagella={n_flagella},seed={seed}",
                        "n_flagella": n_flagella,
                        "seed": seed,
                        "attach_seed": seed,
                        "phase_seed": seed,
                    }
                )
    return conditions


def _runner_config(config: dict[str, Any]) -> dict[str, Any]:
    raw = dict(config.get("runner", {}) or {})
    step_summary_stride = max(1, int(raw.get("step_summary_stride", 1)))
    state_stride = max(1, int(raw.get("state_stride", 1)))
    flush_interval_steps = max(1, int(raw.get("flush_interval_steps", 100)))
    sample_order = str(raw.get("sample_order", "grouped"))
    if sample_order not in {"grouped", "interleave_n_flagella"}:
        raise ValueError(
            "runner.sample_order must be 'grouped' or 'interleave_n_flagella'"
        )
    return {
        "step_summary_stride": step_summary_stride,
        "state_stride": state_stride,
        "flush_interval_steps": flush_interval_steps,
        "sample_order": sample_order,
    }


def order_conditions(
    conditions: list[dict[str, Any]],
    *,
    sample_order: str,
) -> list[dict[str, Any]]:
    if sample_order == "grouped":
        return list(conditions)
    if sample_order != "interleave_n_flagella":
        raise ValueError(
            "runner.sample_order must be 'grouped' or 'interleave_n_flagella'"
        )
    return sorted(
        conditions,
        key=lambda condition: (
            int(condition["attach_seed"]),
            int(condition["phase_seed"]),
            int(condition["n_flagella"]),
        ),
    )


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
        "seed": {
            "global_seed": int(condition["seed"]),
            "attach_seed": int(condition["attach_seed"]),
            "phase_seed": int(condition["phase_seed"]),
        },
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
    step_summary_stride: int,
    state_stride: int,
    flush_interval_steps: int,
) -> dict[str, Any]:
    logger = _setup_sample_logger(log_path)
    logger.info("Sample start: %s", condition["sample_id"])
    logger.info("condition=%s", condition)
    logger.info("effective_config=%s", cfg)
    logger.info(
        (
            "runner_options=step_summary_stride:%d,state_stride:%d,"
            "flush_interval_steps:%d"
        ),
        step_summary_stride,
        state_stride,
        flush_interval_steps,
    )
    sim = Simulator(cfg)
    simulation_start = time.perf_counter()
    states = sim.run(
        cfg.time.duration_s,
        logger=logger,
        progress_interval=progress_interval,
        step_summary_dir=raw_dir,
        stop_on_shape_fail=stop_on_shape_fail,
        step_summary_stride=step_summary_stride,
        state_stride=state_stride,
        flush_interval_steps=flush_interval_steps,
    )
    simulation_elapsed_s = time.perf_counter() - simulation_start
    archive_path = raw_dir / "state_archive.npz"
    trajectory_path = raw_dir / "trajectory.csv"
    archive_start = time.perf_counter()
    save_state_archive(archive_path, states)
    archive_elapsed_s = time.perf_counter() - archive_start
    trajectory_start = time.perf_counter()
    write_trajectory_csv(trajectory_path, states)
    trajectory_elapsed_s = time.perf_counter() - trajectory_start
    logger.info("Sample end: %s states=%d", condition["sample_id"], len(states))
    return {
        "status": "completed",
        "state_count": len(states),
        "timing": {
            "simulation_elapsed_s": simulation_elapsed_s,
            "archive_elapsed_s": archive_elapsed_s,
            "trajectory_elapsed_s": trajectory_elapsed_s,
        },
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
    runner_cfg = _runner_config(analysis_config)

    dataset_id = str(analysis_config["dataset_id"])
    run_batch_id = str(analysis_config.get("run_batch_id", dataset_id))
    output_cfg = analysis_config.get("output", {}) or {}
    run_batch_dir = Path(output_cfg["run_batch_dir"])
    configs_dir = run_batch_dir / "configs"
    samples_root = run_batch_dir / "samples"

    conditions = build_conditions(analysis_config)
    conditions = order_conditions(
        conditions,
        sample_order=str(runner_cfg["sample_order"]),
    )
    if sample_limit is not None:
        conditions = conditions[: int(sample_limit)]

    run_batch_dir.mkdir(parents=True, exist_ok=True)
    configs_dir.mkdir(parents=True, exist_ok=True)
    samples_root.mkdir(parents=True, exist_ok=True)
    effective_analysis_config_path = run_batch_dir / "analysis_config_used.yaml"

    manifest_samples: list[dict[str, Any]] = []
    for condition in tqdm(conditions, desc="flagella sweep", unit="sample"):
        sample_id = str(condition["sample_id"])
        sample_dir = samples_root / sample_id
        raw_dir = sample_dir / "raw"
        log_path = sample_dir / "run.log"
        sample_config_path = configs_dir / f"{sample_id}.yaml"
        sample_config_used_path = sample_dir / "sample_config_used.yaml"
        cfg = _build_sample_config(
            analysis_config=analysis_config,
            base_config=base_config,
            condition=condition,
        )
        sample_config = asdict(cfg)
        sample_config_fingerprint = _config_fingerprint(sample_config)

        sample_record: dict[str, Any] = {
            **condition,
            "dataset_id": dataset_id,
            "run_batch_id": run_batch_id,
            "config_path": str(sample_config_path),
            "sample_config_used_path": str(sample_config_used_path),
            "sample_config_fingerprint": sample_config_fingerprint,
            "sample_dir": str(sample_dir),
            "raw_dir": str(raw_dir),
            "run_log": str(log_path),
            "duration_s": float(cfg.time.duration_s),
            "dt_star": float(cfg.dt_star),
            "torque_Nm": float(cfg.motor.torque_Nm),
            "force_distribution": cfg.motor.force_distribution,
            "attach_seed": (
                int(cfg.seed.attach_seed)
                if cfg.seed.attach_seed is not None
                else int(cfg.seed.global_seed)
            ),
            "phase_seed": (
                int(cfg.seed.phase_seed)
                if cfg.seed.phase_seed is not None
                else int(cfg.seed.global_seed)
            ),
            "cli_overrides": list(cli_overrides),
        }

        if dry_run:
            step_summary = raw_dir / "step_summary.csv"
            if step_summary.exists() and not overwrite:
                existing_config_path = _existing_sample_config_path(
                    sample_config_used_path=sample_config_used_path,
                    batch_config_path=sample_config_path,
                )
                _ensure_existing_config_matches(
                    existing_config_path=existing_config_path,
                    expected_config=sample_config,
                    sample_id=sample_id,
                    step_summary=step_summary,
                )
                if not sample_config_path.is_file():
                    _write_yaml(sample_config_path, sample_config)
            else:
                _write_yaml(sample_config_path, sample_config)
            sample_record["status"] = "planned"
            manifest_samples.append(sample_record)
            continue

        if sample_dir.exists() and overwrite:
            shutil.rmtree(sample_dir)
        sample_dir.mkdir(parents=True, exist_ok=True)
        raw_dir.mkdir(parents=True, exist_ok=True)

        step_summary = raw_dir / "step_summary.csv"
        if step_summary.exists() and not overwrite:
            existing_config_path = _existing_sample_config_path(
                sample_config_used_path=sample_config_used_path,
                batch_config_path=sample_config_path,
            )
            _ensure_existing_config_matches(
                existing_config_path=existing_config_path,
                expected_config=sample_config,
                sample_id=sample_id,
                step_summary=step_summary,
            )
            if not sample_config_path.is_file():
                _write_yaml(sample_config_path, sample_config)
            if not sample_config_used_path.is_file():
                _write_yaml(sample_config_used_path, sample_config)
            sample_record["status"] = "skipped_existing"
            sample_record["outputs"] = {
                "raw_dir": str(raw_dir),
                "sample_config_used_yaml": str(sample_config_used_path),
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

        _write_yaml(sample_config_path, sample_config)
        _write_yaml(sample_config_used_path, sample_config)
        try:
            sample_record["started_at"] = _now_jst()
            sample_start = time.perf_counter()
            sample_record.update(
                _run_sample(
                    cfg=cfg,
                    condition=condition,
                    sample_dir=sample_dir,
                    raw_dir=raw_dir,
                    log_path=log_path,
                    stop_on_shape_fail=stop_on_shape_fail,
                    progress_interval=progress_interval,
                    step_summary_stride=int(runner_cfg["step_summary_stride"]),
                    state_stride=int(runner_cfg["state_stride"]),
                    flush_interval_steps=int(runner_cfg["flush_interval_steps"]),
                )
            )
        except Exception as exc:
            sample_record["status"] = "failed"
            sample_record["error"] = repr(exc)
        finally:
            if "started_at" in sample_record:
                sample_record["ended_at"] = _now_jst()
                sample_record["elapsed_s"] = time.perf_counter() - sample_start
        manifest_samples.append(sample_record)

    _write_yaml(effective_analysis_config_path, analysis_config)
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
                "conf/phase2_analysis/flagella_count_behavior_features.yaml",
            )
        ),
        "base_overrides": analysis_config.get("base_overrides", {}) or {},
        "cli_overrides": list(cli_overrides),
        "sweep": analysis_config.get("sweep", {}) or {},
        "runner": runner_cfg,
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
        default=Path("conf/phase2_analysis/flagella_count_behavior_dataset.yaml"),
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
