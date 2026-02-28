from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
import typer
import yaml

from sim_swim.core.run_context import init_run
from sim_swim.render.project2d import project_states
from sim_swim.render.render3d import save_swim_movie
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig, merge_overrides

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
        None, help="Override output_sampling.fps_out_2d (frames per second)"
    ),
    render_flagella: Optional[bool] = typer.Option(
        None,
        help="Enable flagella rendering in 3D (overrides render.render_flagella)",
    ),
    render_flagella_2d: Optional[bool] = typer.Option(
        None,
        help="Enable flagella rendering in 2D (overrides render.render_flagella_2d)",
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
        override_dict.setdefault("output_sampling", {})["fps_out_2d"] = fps_out
    if render_flagella is not None:
        override_dict.setdefault("render", {})["render_flagella"] = render_flagella
    if render_flagella_2d is not None:
        override_dict.setdefault("render", {})["render_flagella_2d"] = (
            render_flagella_2d
        )
    cfg = SimulationConfig.from_dict(raw_cfg).with_overrides(override_dict)

    output_base = cfg.output.base_dir
    ctx = init_run(
        base_dir=output_base,
        input_info={"config": str(config), "overrides": overrides or []},
    )
    logger = ctx.logger
    logger.info("Loaded simulation config (effective): %s", cfg)
    logger.info("Overrides: %s", override_dict if override_dict else "None")
    cfg.validate_time_scaling()
    logger.info("[time-scale] τ is fixed to 1.0 in stable integration mode")

    if cfg.use_eta_b3_torque:
        logger.info(
            "[time-scale] input motor.torque_Nm=-1; using internal mode (tau=1)"
        )
        logger.info(
            (
                "[time-scale] b=%.6e um, b_m=%.6e m, η=%.6e Pa·s, "
                "T_scale=η b^3=%.6e N·m, T_motor=%.6e N·m, τ=%.6e s"
            ),
            cfg.scale.b_um,
            cfg.b_m,
            cfg.viscosity_Pa_s,
            cfg.torque_scale_Nm,
            cfg.motor_torque_Nm,
            cfg.tau_s,
        )
    elif cfg.is_motor_off_torque:
        logger.info("[time-scale] input motor.torque_Nm=0; motor OFF mode (tau=1)")
        logger.info(
            (
                "[time-scale] b=%.6e um, b_m=%.6e m, η=%.6e Pa·s, "
                "T_scale=η b^3=%.6e N·m, T_motor=0, τ=%.6e s"
            ),
            cfg.scale.b_um,
            cfg.b_m,
            cfg.viscosity_Pa_s,
            cfg.torque_scale_Nm,
            cfg.tau_s,
        )
    else:
        logger.info(
            "[time-scale] using physical motor.torque_Nm: %.6e N·m",
            cfg.input_torque_Nm,
        )
        logger.info(
            (
                "[time-scale] b=%.6e um, b_m=%.6e m, η=%.6e Pa·s, "
                "T_scale=|T_motor|=%.6e N·m, T_motor=%.6e N·m, τ=%.6e s"
            ),
            cfg.scale.b_um,
            cfg.b_m,
            cfg.viscosity_Pa_s,
            cfg.torque_scale_Nm,
            cfg.motor_torque_Nm,
            cfg.tau_s,
        )
    logger.info(
        "[time-step ] Δt_internal=%.6e s, τ=%.6e s, Δt*=dt_star=%.6e (=Δt/τ)",
        cfg.dt_s,
        cfg.tau_s,
        cfg.dt_star,
    )
    logger.info("[time-step ] output_dt_s(config)=%.6e s", cfg.output_dt_s)
    logger.info(
        (
            "[duration  ] t_end=duration_s=%.6e s, "
            "t_end*=duration_star=%.6e (=t_end/τ), steps=%d"
        ),
        cfg.time.duration_s,
        cfg.duration_star,
        cfg.total_steps,
    )

    n_layers = cfg.compute_body_n_layers()
    l_over_b = cfg.body.length_total_um / max(cfg.scale.b_um, 1e-12)
    logger.info(
        "Body discretization: b_um=%.6f, L_um=%.6f, L_over_b=%.6f, dz_over_b=%.6f, n_layers=%d",
        cfg.scale.b_um,
        cfg.body.length_total_um,
        l_over_b,
        cfg.body.prism.dz_over_b,
        n_layers,
    )

    for name, ok, actual, expected in cfg.paper_reference_checks():
        if ok:
            logger.info(
                "Paper check OK: %s (actual=%s, expected=%s)", name, actual, expected
            )
        else:
            logger.warning(
                "Paper check mismatch: %s (actual=%s, expected=%s)",
                name,
                actual,
                expected,
            )

    sim_duration_s = float(cfg.time.duration_s)
    simulator = Simulator(cfg)
    states = simulator.run(
        sim_duration_s, logger=logger, step_summary_dir=ctx.out.sim_dir
    )

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
            "step_summary_csv": str(ctx.out.sim_dir / "step_summary.csv"),
            "render3d": str(ctx.out.render_dir),
            "render2d": str(ctx.out.render2d_dir),
        }
    )
    manifest["outputs"] = outputs
    manifest["files"] = [
        str(traj_path.relative_to(ctx.out.root)),
        "sim/step_summary.csv",
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
