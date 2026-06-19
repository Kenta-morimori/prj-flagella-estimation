#!/usr/bin/env python3
"""Render 3D/2D outputs from a flagella-count behavior raw sample archive."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import subprocess
from typing import Any

import yaml

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.flagella_count_behavior import (
    archive_metadata,
    load_state_archive,
    write_trajectory_csv,
)
from sim_swim.render.project2d import project_states
from sim_swim.render.render3d import save_swim_movie
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _git_info() -> dict[str, Any]:
    def run(cmd: list[str]) -> str:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True).strip()

    try:
        return {
            "commit": run(["git", "rev-parse", "HEAD"]),
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


def _setup_logger(log_path: Path) -> logging.Logger:
    logger = logging.getLogger(f"flagella_count_behavior.replay.{log_path.parent.name}")
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
    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    logger.addHandler(stream)
    return logger


def _infer_config_path(sample_dir: Path) -> Path:
    return sample_dir.parent.parent / "configs" / f"{sample_dir.name}.yaml"


def _infer_archive_path(sample_dir: Path) -> Path:
    return sample_dir / "raw" / "state_archive.npz"


def _build_render_sampling_overrides(
    *,
    out_all_steps_3d: bool | None,
    fps_out_3d: float | None,
    fps_out_2d: float | None,
) -> dict[str, Any]:
    sampling: dict[str, Any] = {}
    if out_all_steps_3d is not None:
        sampling["out_all_steps_3d"] = bool(out_all_steps_3d)
    if fps_out_3d is not None:
        if fps_out_3d <= 0.0:
            raise ValueError("--fps-out-3d must be positive")
        sampling["fps_out_3d"] = float(fps_out_3d)
    if fps_out_2d is not None:
        if fps_out_2d <= 0.0:
            raise ValueError("--fps-out-2d must be positive")
        sampling["fps_out_2d"] = float(fps_out_2d)
    if not sampling:
        return {}
    return {"output_sampling": sampling}


def _render_sampling_summary(cfg: SimulationConfig) -> dict[str, Any]:
    return {
        "out_all_steps_3d": cfg.output_sampling.out_all_steps_3d,
        "fps_out_3d": cfg.output_sampling.fps_out_3d,
        "fps_out_2d": cfg.output_sampling.fps_out_2d,
    }


def render_sample(
    *,
    sample_dir: Path,
    output_dir: Path,
    config_path: Path | None,
    archive_path: Path | None,
    out_all_steps_3d: bool | None = False,
    fps_out_3d: float | None = None,
    fps_out_2d: float | None = None,
) -> Path:
    sample_dir = sample_dir.resolve()
    config_path = (config_path or _infer_config_path(sample_dir)).resolve()
    archive_path = (archive_path or _infer_archive_path(sample_dir)).resolve()

    raw_cfg = _load_yaml(config_path)
    cfg = SimulationConfig.from_dict(raw_cfg)
    render_sampling_overrides = _build_render_sampling_overrides(
        out_all_steps_3d=out_all_steps_3d,
        fps_out_3d=fps_out_3d,
        fps_out_2d=fps_out_2d,
    )
    if render_sampling_overrides:
        cfg = cfg.with_overrides(render_sampling_overrides)
    states = load_state_archive(archive_path)

    simulator = Simulator(cfg)
    render_dir = output_dir / "render"
    render2d_dir = output_dir / "render2d"
    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "run.log"
    logger = _setup_logger(log_path)
    logger.info("Replay start: sample_dir=%s", sample_dir)
    logger.info("config_path=%s", config_path)
    logger.info("archive_path=%s", archive_path)
    logger.info("state_count=%d", len(states))
    logger.info("render_sampling=%s", _render_sampling_summary(cfg))
    logger.info(
        "render_sampling_overrides=%s",
        render_sampling_overrides.get("output_sampling", {}),
    )

    trajectory_path = output_dir / "trajectory.csv"
    write_trajectory_csv(trajectory_path, states)
    save_swim_movie(states, cfg, simulator.rig, render_dir)
    project_states(states, cfg, simulator.rig, render2d_dir)

    manifest = {
        "git": _git_info(),
        "input": {
            "sample_dir": str(sample_dir),
            "config_path": str(config_path),
            "archive_path": str(archive_path),
        },
        "archive": archive_metadata(
            sample_id=sample_dir.name,
            config_path=str(config_path),
        ),
        "render_sampling": _render_sampling_summary(cfg),
        "render_sampling_overrides": render_sampling_overrides.get(
            "output_sampling", {}
        ),
        "outputs": {
            "root": str(output_dir),
            "trajectory_csv": str(trajectory_path),
            "render_dir": str(render_dir),
            "render2d_dir": str(render2d_dir),
            "log": str(log_path),
        },
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    logger.info("Replay end")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--archive", type=Path, default=None)
    sampling = parser.add_mutually_exclusive_group()
    sampling.add_argument(
        "--out-all-steps-3d",
        dest="out_all_steps_3d",
        action="store_true",
        help="Render every archived 3D state instead of sampling by fps.",
    )
    sampling.add_argument(
        "--no-out-all-steps-3d",
        dest="out_all_steps_3d",
        action="store_false",
        help="Sample 3D replay frames by --fps-out-3d.",
    )
    parser.set_defaults(out_all_steps_3d=False)
    parser.add_argument(
        "--fps-out-3d",
        type=float,
        default=None,
        help="3D replay frame rate when --no-out-all-steps-3d is active.",
    )
    parser.add_argument(
        "--fps-out-2d",
        type=float,
        default=None,
        help="2D replay frame rate.",
    )
    args = parser.parse_args()

    sample_dir = args.sample_dir
    output_dir = args.output_dir or (sample_dir / "replay")
    render_sample(
        sample_dir=sample_dir,
        output_dir=output_dir,
        config_path=args.config,
        archive_path=args.archive,
        out_all_steps_3d=args.out_all_steps_3d,
        fps_out_3d=args.fps_out_3d,
        fps_out_2d=args.fps_out_2d,
    )
    print(output_dir)


if __name__ == "__main__":
    main()
