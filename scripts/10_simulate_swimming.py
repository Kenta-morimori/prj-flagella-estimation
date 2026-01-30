from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import typer
import yaml

from flagella_estimation.core.run_context import init_run
from flagella_sim.render.project2d import heading_from_quat, project_states
from flagella_sim.render.render3d import save_swim_movie
from flagella_sim.sim.core import Simulator
from flagella_sim.sim.params import SimulationConfig, merge_overrides

app = typer.Typer(
    add_completion=False,
    help="Simulate swimming and render projections (MVP skeleton).",
)


def _load_config(path: Path) -> Dict[str, Any]:
    raw_text = path.read_text(encoding="utf-8") if path.exists() else ""
    return yaml.safe_load(raw_text) or {}


def _to_nested_overrides(items: Optional[List[str]]) -> Dict[str, Any]:
    return merge_overrides({}, items)


@app.command()
def main(
    config: Path = typer.Option(Path("conf/sim_swim.yaml"), "--config", "-c"),
    overrides: List[str] = typer.Argument(
        None,
        help="Optional overrides as key=value (e.g., flagella.n_flagella=6)",
    ),
) -> None:
    """Phase2 用のシミュレーション＆投影エントリ。"""

    raw_cfg = _load_config(config)
    override_dict = _to_nested_overrides(overrides)
    cfg = SimulationConfig.from_dict(raw_cfg).with_overrides(override_dict)

    output_base = raw_cfg.get("output", {}).get("base_dir", "outputs")
    ctx = init_run(
        base_dir=output_base,
        input_info={"config": str(config), "overrides": overrides or []},
    )
    logger = ctx.logger
    logger.info("Loaded simulation config: %s", cfg)
    logger.info("Overrides: %s", overrides)

    duration_s = float(raw_cfg.get("time", {}).get("duration_s", 0.1))
    simulator = Simulator(cfg)
    states = simulator.run(duration_s)

    # 保存: 3D軌跡（全ステップ）
    traj_path = ctx.out.sim_dir / "trajectory.csv"
    traj_path.parent.mkdir(parents=True, exist_ok=True)
    traj_df = pd.DataFrame(
        [
            {
                "t": st.t,
                "x": st.position_um[0],
                "y": st.position_um[1],
                "z": st.position_um[2],
                "qx": st.quaternion[0],
                "qy": st.quaternion[1],
                "qz": st.quaternion[2],
                "qw": st.quaternion[3],
                "vx": st.velocity_um_s[0],
                "vy": st.velocity_um_s[1],
                "vz": st.velocity_um_s[2],
                "wx": st.omega_rad_s[0],
                "wy": st.omega_rad_s[1],
                "wz": st.omega_rad_s[2],
            }
            for st in states
        ]
    )
    traj_df.to_csv(traj_path, index=False)
    logger.info("Saved trajectory to %s", traj_path)

    # レンダリング
    save_swim_movie(states, cfg, simulator.rig, ctx.out.render_dir)
    project_states(states, cfg, simulator.rig, ctx.out.render2d_dir)

    # tracking出力（投影座標）
    track_rows = []
    px_per_um = 1.0 / cfg.render.pixel_size_um
    body_major_px = cfg.body.length_total_um * px_per_um
    body_minor_px = cfg.body.diameter_um * px_per_um
    for frame_idx, st in enumerate(states):
        cx = cfg.render.image_size_px / 2 + st.position_um[0] * px_per_um
        cy = cfg.render.image_size_px / 2 + st.position_um[1] * px_per_um
        theta = heading_from_quat(st.quaternion)
        track_rows.append(
            {
                "frame": frame_idx,
                "track_id": 0,
                "cx": cx,
                "cy": cy,
                "theta": theta,
                "major": body_major_px,
                "minor": body_minor_px,
                "vx": st.velocity_um_s[0] * px_per_um,
                "vy": st.velocity_um_s[1] * px_per_um,
                "is_valid": True,
            }
        )
    track_df = pd.DataFrame(track_rows)
    track_path = ctx.out.tracking_dir / "track.csv"
    track_df.to_csv(track_path, index=False)
    butt_path = ctx.out.tracking_dir / "butt.json"
    butt_path.write_text("{}", encoding="utf-8")
    # overlay用 mp4（2D投影を流用）
    overlay_src = ctx.out.render2d_dir / "projection.mp4"
    overlay_dst = ctx.out.tracking_dir / "overlay.mp4"
    overlay_dst.write_bytes(overlay_src.read_bytes())

    # manifest に出力一覧を追記
    manifest_path = ctx.out.root / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except Exception:
        manifest = {}
    outputs = manifest.get("outputs", {})
    outputs.update(
        {
            "trajectory_csv": str(traj_path),
            "render3d": str(ctx.out.render_dir),
            "render2d": str(ctx.out.render2d_dir),
        }
    )
    manifest["outputs"] = outputs
    manifest["files"] = [
        str(traj_path.relative_to(ctx.out.root)),
        "render/swim3d.mp4",
        "render/swim3d_final.png",
        "render2d/projection.mp4",
    ]
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    logger.info("Manifest updated with simulation outputs")


if __name__ == "__main__":
    app()
