#!/usr/bin/env python3
"""Replay Phase 3 common clip datasets for visual review."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from flagella_estimation.phase3.replay import (
    ReplayConfig,
    render_3d_2d_grid_mp4,
    render_contact_sheet,
)


def _optional_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no"}:
        return False
    raise argparse.ArgumentTypeError("expected true/false")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset_dir", type=Path)
    parser.add_argument("output_dir", type=Path, nargs="?")
    parser.add_argument("--n-flagella", type=int, default=None)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--clip-index", type=int, default=None)
    parser.add_argument("--time-band", default=None)
    parser.add_argument("--qc-label", default=None)
    parser.add_argument("--training-candidate", type=_optional_bool, default=None)
    parser.add_argument("--max-clips", type=int, default=12)
    parser.add_argument(
        "--frames-per-clip",
        type=int,
        default=4,
        help="contact-sheet mode only; MP4 grid always uses all clip frames",
    )
    parser.add_argument(
        "--mode",
        choices=("mp4-grid", "contact-sheet"),
        default="mp4-grid",
        help="visual replay artifact to create",
    )
    parser.add_argument(
        "--panel-layout",
        choices=("3d+2d", "3d", "2d"),
        default="3d+2d",
        help="MP4 grid panel layout",
    )
    args = parser.parse_args(argv)
    cfg = ReplayConfig(
        dataset_dir=args.dataset_dir,
        output_dir=args.output_dir,
        n_flagella=args.n_flagella,
        run_id=args.run_id,
        clip_index=args.clip_index,
        time_band=args.time_band,
        qc_label=args.qc_label,
        training_candidate=args.training_candidate,
        max_clips=args.max_clips,
        frames_per_clip=args.frames_per_clip,
        panel_layout=args.panel_layout,
    )
    output_dir = (
        render_contact_sheet(cfg)
        if args.mode == "contact-sheet"
        else render_3d_2d_grid_mp4(cfg)
    )
    print(output_dir)


if __name__ == "__main__":
    main()
