from __future__ import annotations

from copy import deepcopy
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

_RUNTIME_DEFAULTS: Dict[str, Any] = {
    "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
    "body": {
        "prism": {"n_prism": 3, "dz_over_b": 0.5, "radius_over_b": 0.5, "axis": "x"},
        "length_total_um": 2.0,
    },
    "flagella": {
        "n_flagella": 3,
        "placement_mode": "uniform",
        "discretization": {"ds_over_b": 0.58},
        "bond_L_over_b": 0.58,
        "length_over_b": 5.8,
        "helix_init": {"radius_over_b": 0.2, "pitch_over_b": 1.0},
    },
    "potentials": {
        "spring": {"H_over_T_over_b": 10.0, "s": 0.1},
        "bend": {"kb_over_T": 20.0},
        "torsion": {"kt_over_T": 10.0},
        "spring_spring_repulsion": {
            "A_ss_over_T": 1.0,
            "a_ss_over_b": 0.2,
            "cutoff_over_b": 0.2,
        },
    },
    "hook": {"enabled": True, "threshold_deg": 90.0, "kb_over_T": 20.0},
    "run_tumble": {
        "run_tau": 1200.0,
        "tumble_tau": 800.0,
        "semicoiled_tau": 400.0,
        "curly1_tau": 400.0,
    },
    "time": {"duration_s": 0.1, "dt_s": 2.5e-7},
    "output_sampling": {"out_all_steps_3d": True, "fps_out_2d": 25.0},
    "render": {
        "image_size_px": 256,
        "pixel_size_um": 0.203,
        "flagella_linewidth_px": 2.0,
        "render_flagella": True,
        "save_frames_3d": True,
        "follow_camera_3d": True,
        "view_range_um": 10.0,
        "timestamp_3d": True,
        "timestamp_fmt": "t = {t:.3f} s",
        "label_flagella": True,
        "follow_camera_2d": False,
        "save_frames_2d": True,
    },
    "brownian": {"enabled": False, "temperature_K": 298.0},
    "seed": {"global_seed": 0},
    "output": {"base_dir": "outputs"},
}


def _load_config(path: Path) -> Dict[str, Any]:
    raw_text = path.read_text(encoding="utf-8") if path.exists() else ""
    return yaml.safe_load(raw_text) or {}


def _deep_merge(dst: Dict[str, Any], src: Dict[str, Any]) -> Dict[str, Any]:
    for key, value in src.items():
        if isinstance(value, dict) and isinstance(dst.get(key), dict):
            _deep_merge(dst[key], value)
        else:
            dst[key] = value
    return dst


def _expand_phase25_config(raw_cfg: Dict[str, Any]) -> Dict[str, Any]:
    cfg = _deep_merge(deepcopy(_RUNTIME_DEFAULTS), raw_cfg)

    scale_raw = cfg.setdefault("scale", {})
    b_m_raw = scale_raw.get("b_m")
    if b_m_raw is not None:
        scale_raw["b_um"] = float(b_m_raw) * 1e6

    bead_raw = cfg.get("bead", {}) or {}
    if "a_star" in bead_raw and "bead_radius_a_over_b" not in scale_raw:
        scale_raw["bead_radius_a_over_b"] = float(bead_raw["a_star"])

    spring_raw = cfg.get("spring", {}) or {}
    elastic_raw = cfg.get("elastic", {}) or {}
    rep_raw = cfg.get("repulsion", {}) or {}
    hook_raw = cfg.setdefault("hook", {})
    potentials_raw = cfg.setdefault("potentials", {})
    spring_p = potentials_raw.setdefault("spring", {})
    bend_p = potentials_raw.setdefault("bend", {})
    torsion_p = potentials_raw.setdefault("torsion", {})
    rep_p = potentials_raw.setdefault("spring_spring_repulsion", {})
    flag_raw = cfg.setdefault("flagella", {})
    flag_dis = flag_raw.setdefault("discretization", {})

    if "L_star" in spring_raw:
        l_star = float(spring_raw["L_star"])
        flag_dis["ds_over_b"] = l_star
        flag_raw["bond_L_over_b"] = l_star
    if "H_star" in spring_raw:
        spring_p["H_over_T_over_b"] = float(spring_raw["H_star"])
    if "s" in spring_raw:
        spring_p["s"] = float(spring_raw["s"])

    if "kb_star" in elastic_raw:
        bend_p["kb_over_T"] = float(elastic_raw["kb_star"])
        hook_raw.setdefault("kb_over_T", float(elastic_raw["kb_star"]))
    if "kt_star" in elastic_raw:
        torsion_p["kt_over_T"] = float(elastic_raw["kt_star"])

    if "A_ss_star" in rep_raw:
        rep_p["A_ss_over_T"] = float(rep_raw["A_ss_star"])
    if "a_ss_star" in rep_raw:
        rep_p["a_ss_over_b"] = float(rep_raw["a_ss_star"])
    if "cutoff_star" in rep_raw:
        rep_p["cutoff_over_b"] = float(rep_raw["cutoff_star"])

    if "theta_threshold_deg" in hook_raw and "threshold_deg" not in hook_raw:
        hook_raw["threshold_deg"] = float(hook_raw["theta_threshold_deg"])

    integrator_raw = cfg.get("integrator", {}) or {}
    time_raw = cfg.setdefault("time", {})
    if "dt_s" not in time_raw and "dt_star" in integrator_raw:
        b_m = float(scale_raw.get("b_m", float(scale_raw.get("b_um", 1.0)) * 1e-6))
        eta = float((cfg.get("fluid", {}) or {}).get("viscosity_Pa_s", 1e-3))
        torque = float((cfg.get("motor", {}) or {}).get("torque_Nm", 4e-18))
        tau_s = eta * (b_m**3) / max(abs(torque), 1e-30)
        time_raw["dt_s"] = float(integrator_raw["dt_star"]) * tau_s

    return cfg


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
        help="Enable flagella rendering (overrides render.render_flagella)",
    ),
    overrides: List[str] = typer.Argument(
        None,
        help="Optional overrides as key=value (e.g., flagella.n_flagella=6)",
    ),
) -> None:
    """Phase2 用のシミュレーション＆投影エントリ。"""

    raw_cfg = _expand_phase25_config(_load_config(config))
    override_dict = _to_nested_overrides(overrides)
    if duration_s is not None:
        override_dict.setdefault("time", {})["duration_s"] = duration_s
    if fps_out is not None:
        override_dict.setdefault("output_sampling", {})["fps_out_2d"] = fps_out
    if render_flagella is not None:
        override_dict.setdefault("render", {})["render_flagella"] = render_flagella
    cfg = SimulationConfig.from_dict(raw_cfg).with_overrides(override_dict)

    output_base = cfg.output.base_dir
    ctx = init_run(
        base_dir=output_base,
        input_info={"config": str(config), "overrides": overrides or []},
    )
    logger = ctx.logger
    logger.info("Loaded simulation config (effective): %s", cfg)
    logger.info("Overrides: %s", override_dict if override_dict else "None")

    logger.info(
        (
            "[time-scale] b=%.6e um, b_m=%.6e m, eta=%.6e Pa*s, "
            "torque=%.6e N*m, tau_s=%.6e s"
        ),
        cfg.scale.b_um,
        cfg.b_m,
        cfg.viscosity_Pa_s,
        cfg.motor.torque_Nm,
        cfg.tau_s,
    )
    logger.info(
        "[time-step ] dt_s=%.6e s, dt_star=%.6e (=dt_s/tau_s)",
        cfg.dt_s,
        cfg.dt_star,
    )
    logger.info(
        "[duration  ] duration_s=%.6e s, duration_star=%.6e, total_steps=%d",
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

    sim_duration_s = float(cfg.time.duration_s)
    simulator = Simulator(cfg)
    states = simulator.run(sim_duration_s, logger=logger)

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
        "render/movie_3d.mp4",
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
