"""Render full-duration, body-only 2D replay grids for phase-seed studies."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
import platform
import re
import subprocess
from typing import Any
from zoneinfo import ZoneInfo

import cv2
import numpy as np

from flagella_estimation.phase3.render import select_frames
from sim_swim.analysis.flagella_count_behavior import load_state_archive
from sim_swim.render.body2d import BodyCapsuleRenderConfig, render_body_capsule_frame
from sim_swim.render.grid_movie import GridLayout, compose_grid_frame, write_mp4_grid


_CONDITION_ID = re.compile(r"^as(?P<attach>\d+)__ps(?P<phase>\d+)__nf(?P<n>\d+)$")


@dataclass(frozen=True)
class PhaseSeedReplayConfig:
    run_dir: Path
    output_dir: Path
    duration_s: float = 2.0
    fps: float = 25.0
    image_size_px: int = 96
    pixel_size_um: float = 0.1
    panel_width_px: int = 256
    panel_height_px: int = 182
    overwrite: bool = False


@dataclass(frozen=True)
class ReplayCondition:
    condition_id: str
    attach_seed: int
    phase_seed: int
    n_flagella: int
    states: list[Any]
    first_fail_t_s: float | None


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _first_fail_t_s(summary: dict[str, Any]) -> float | None:
    observed: list[float] = []
    for gate_name in ("shape_nonbody", "shape_body"):
        value = dict((summary.get("gates") or {}).get(gate_name) or {}).get(
            "first_observed_fail_t_s"
        )
        if value is not None:
            observed.append(float(value))
    return min(observed) if observed else None


def _sample_full_duration(
    states: list[Any], *, duration_s: float, fps: float
) -> list[Any]:
    sampled = select_frames(states, fps)
    result = [state for state in sampled if float(state.t) < duration_s - 1.0e-9]
    expected_count = int(round(duration_s * fps))
    if len(result) != expected_count:
        raise ValueError(
            "archive does not provide the requested full-duration replay: "
            f"expected_frames={expected_count}, observed_frames={len(result)}, "
            f"duration_s={duration_s:g}, fps={fps:g}"
        )
    return result


def _discover_conditions(cfg: PhaseSeedReplayConfig) -> list[ReplayCondition]:
    discovered: list[ReplayCondition] = []
    for condition_dir in sorted(cfg.run_dir.iterdir()):
        if not condition_dir.is_dir():
            continue
        match = _CONDITION_ID.match(condition_dir.name)
        if match is None or int(match["attach"]) != 0:
            continue
        archive_path = condition_dir / "state_archive.npz"
        summary_path = condition_dir / "run_summary.json"
        if not archive_path.is_file() or not summary_path.is_file():
            continue
        summary = _read_json(summary_path)
        execution = dict(summary.get("execution") or {})
        if execution.get("status") != "completed":
            continue
        discovered.append(
            ReplayCondition(
                condition_id=condition_dir.name,
                attach_seed=int(match["attach"]),
                phase_seed=int(match["phase"]),
                n_flagella=int(match["n"]),
                states=_sample_full_duration(
                    load_state_archive(archive_path),
                    duration_s=cfg.duration_s,
                    fps=cfg.fps,
                ),
                first_fail_t_s=_first_fail_t_s(summary),
            )
        )
    discovered.sort(key=lambda item: (item.n_flagella, item.phase_seed))
    expected = {
        (n_flagella, phase_seed)
        for n_flagella in range(1, 5)
        for phase_seed in range(3)
    }
    actual = {(item.n_flagella, item.phase_seed) for item in discovered}
    if actual != expected:
        missing = sorted(expected - actual)
        unexpected = sorted(actual - expected)
        raise ValueError(
            "phase-seed replay requires completed attach_seed=0 conditions for "
            f"n=1..4 and ps=0..2; missing={missing}, unexpected={unexpected}"
        )
    return discovered


def _render_panel(
    condition: ReplayCondition,
    state: Any,
    *,
    panel_width_px: int,
    panel_height_px: int,
    body_render_cfg: BodyCapsuleRenderConfig,
) -> np.ndarray:
    title_height = 34
    has_failed = (
        condition.first_fail_t_s is not None
        and float(state.t) >= condition.first_fail_t_s
    )
    panel = np.full((panel_height_px, panel_width_px, 3), 255, dtype=np.uint8)
    title_color = (225, 235, 255) if has_failed else (245, 245, 245)
    panel[:title_height] = title_color
    status = "QC FAIL" if has_failed else "strict pass"
    status_color = (0, 0, 190) if has_failed else (35, 90, 35)
    cv2.putText(
        panel,
        f"n={condition.n_flagella}  ps={condition.phase_seed}  t={float(state.t):.2f}s",
        (7, 14),
        cv2.FONT_HERSHEY_SIMPLEX,
        0.42,
        (35, 35, 35),
        1,
        cv2.LINE_AA,
    )
    cv2.putText(
        panel,
        status,
        (7, 29),
        cv2.FONT_HERSHEY_SIMPLEX,
        0.38,
        status_color,
        1,
        cv2.LINE_AA,
    )
    frame, _ = render_body_capsule_frame(state, body_render_cfg)
    panel[title_height:] = cv2.cvtColor(
        cv2.resize(
            frame,
            (panel_width_px, panel_height_px - title_height),
            interpolation=cv2.INTER_NEAREST,
        ),
        cv2.COLOR_GRAY2BGR,
    )
    return panel


def _git_info() -> dict[str, str | None]:
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    return {"commit": commit}


def render_phase_seed_2d_replay(cfg: PhaseSeedReplayConfig) -> Path:
    """Write the fixed 4-by-3, 2D-only full-duration comparison replay."""
    if cfg.fps <= 0.0 or cfg.duration_s <= 0.0:
        raise ValueError("duration_s and fps must be positive")
    if cfg.output_dir.exists():
        if not cfg.overwrite:
            raise FileExistsError(f"Output already exists: {cfg.output_dir}")
        import shutil

        shutil.rmtree(cfg.output_dir)
    conditions = _discover_conditions(cfg)
    cfg.output_dir.mkdir(parents=True)
    body_render_cfg = BodyCapsuleRenderConfig(
        image_size_px=cfg.image_size_px,
        pixel_size_um=cfg.pixel_size_um,
        tracking_center=True,
    )
    layout = GridLayout(
        rows=4,
        cols=3,
        cell_width_px=cfg.panel_width_px,
        cell_height_px=cfg.panel_height_px,
    )
    video_path = cfg.output_dir / "phase_seed_2d_grid.mp4"

    def frames() -> Any:
        for frame_index in range(int(round(cfg.duration_s * cfg.fps))):
            yield compose_grid_frame(
                [
                    _render_panel(
                        condition,
                        condition.states[frame_index],
                        panel_width_px=cfg.panel_width_px,
                        panel_height_px=cfg.panel_height_px,
                        body_render_cfg=body_render_cfg,
                    )
                    for condition in conditions
                ],
                layout,
            )

    result = write_mp4_grid(video_path, frames=frames(), fps=cfg.fps)
    final_frame = compose_grid_frame(
        [
            _render_panel(
                condition,
                condition.states[-1],
                panel_width_px=cfg.panel_width_px,
                panel_height_px=cfg.panel_height_px,
                body_render_cfg=body_render_cfg,
            )
            for condition in conditions
        ],
        layout,
    )
    final_png_path = cfg.output_dir / "phase_seed_2d_grid_final.png"
    if not cv2.imwrite(str(final_png_path), final_frame):
        raise RuntimeError(f"Failed to write {final_png_path}")
    created_at = datetime.now(ZoneInfo("Asia/Tokyo")).isoformat()
    manifest = {
        "pipeline_name": "phase_seed_full_duration_2d_replay",
        "created_at": created_at,
        "input_run_dir": str(cfg.run_dir),
        "duration_s": cfg.duration_s,
        "fps": cfg.fps,
        "frame_count": int(round(cfg.duration_s * cfg.fps)),
        "layout": {
            "rows": 4,
            "columns": 3,
            "row": "n_flagella",
            "column": "phase_seed",
        },
        "render": {
            "mode": "body_capsule_orthographic_v1",
            "image_size_px": cfg.image_size_px,
            "pixel_size_um": cfg.pixel_size_um,
            "tracking_center": True,
            "flagella_rendered": False,
        },
        "conditions": [
            {
                "condition_id": item.condition_id,
                "n_flagella": item.n_flagella,
                "attach_seed": item.attach_seed,
                "phase_seed": item.phase_seed,
                "first_fail_t_s": item.first_fail_t_s,
            }
            for item in conditions
        ],
        "outputs": {
            "mp4": str(video_path),
            "final_png": str(final_png_path),
            **result.to_manifest(),
        },
        "git": _git_info(),
        "environment": {"python": platform.python_version(), "opencv": cv2.__version__},
    }
    (cfg.output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (cfg.output_dir / "run.log").write_text(
        "\n".join(
            [
                f"created_at={created_at}",
                f"input_run_dir={cfg.run_dir}",
                f"output_dir={cfg.output_dir}",
                f"duration_s={cfg.duration_s:g}",
                f"fps={cfg.fps:g}",
                f"condition_count={len(conditions)}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return cfg.output_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--duration-s", type=float, default=2.0)
    parser.add_argument("--fps", type=float, default=25.0)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    print(
        render_phase_seed_2d_replay(
            PhaseSeedReplayConfig(
                run_dir=args.run_dir,
                output_dir=args.output_dir,
                duration_s=args.duration_s,
                fps=args.fps,
                overwrite=args.overwrite,
            )
        )
    )
