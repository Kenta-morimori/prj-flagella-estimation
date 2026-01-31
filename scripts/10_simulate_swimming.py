from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import typer
import yaml

from flagella_estimation.core.run_context import init_run
from flagella_sim.render.project2d import project_states
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
    duration_s: Optional[float] = typer.Option(
        None, help="Override time.duration_s (seconds)"
    ),
    fps_out: Optional[float] = typer.Option(
        None, help="Override time.fps_out (frames per second)"
    ),
    render_flagella: bool = typer.Option(
        None,
        help="Enable flagella rendering (overrides render.render_flagella)",
    ),
    overrides: List[str] = typer.Argument(
        None,
        help="Optional overrides as key=value (e.g., flagella.n_flagella=6)",
    ),
) -> None:
    """Phase2 用のシミュレーション＆投影エントリ。"""

    raw_cfg = _load_config(config)
    override_dict = _to_nested_overrides(overrides)
    if duration_s is not None:
        override_dict.setdefault("time", {})["duration_s"] = duration_s
    if fps_out is not None:
        override_dict.setdefault("time", {})["fps_out"] = fps_out
    if render_flagella is not None:
        override_dict.setdefault("render", {})["render_flagella"] = render_flagella
    cfg = SimulationConfig.from_dict(raw_cfg).with_overrides(override_dict)

    output_base = raw_cfg.get("output", {}).get("base_dir", "outputs")
    ctx = init_run(
        base_dir=output_base,
        input_info={"config": str(config), "overrides": overrides or []},
    )
    logger = ctx.logger
    logger.info("Loaded simulation config (effective): %s", cfg)
    logger.info("Overrides: %s", override_dict if override_dict else "None")

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
