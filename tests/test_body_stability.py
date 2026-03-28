from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _triplet_angle_rad(positions_m: np.ndarray, triplet: np.ndarray) -> float:
    i, j, k = (int(triplet[0]), int(triplet[1]), int(triplet[2]))
    u = positions_m[i] - positions_m[j]
    v = positions_m[k] - positions_m[j]
    nu = max(float(np.linalg.norm(u)), 1e-18)
    nv = max(float(np.linalg.norm(v)), 1e-18)
    c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return float(np.arccos(c))


def _body_spring_rows(sim: Simulator) -> np.ndarray:
    pairs = sim.model.spring_pairs
    bi = sim.model.bead_is_body[pairs[:, 0]]
    bj = sim.model.bead_is_body[pairs[:, 1]]
    return np.where(bi & bj)[0]


def _body_bending_rows(sim: Simulator) -> np.ndarray:
    return np.where(sim.model.bending_flag_ids < 0)[0]


def _body_axis(sim: Simulator, positions_m: np.ndarray) -> np.ndarray:
    first = np.mean(positions_m[sim.model.body_layer_indices[0]], axis=0)
    last = np.mean(positions_m[sim.model.body_layer_indices[-1]], axis=0)
    axis = last - first
    n = max(float(np.linalg.norm(axis)), 1e-18)
    return axis / n


def _body_com(sim: Simulator, positions_m: np.ndarray) -> np.ndarray:
    return np.mean(positions_m[sim.model.body_indices], axis=0)


def _body_shape_metrics(sim: Simulator, positions_m: np.ndarray) -> tuple[float, float]:
    body_spring_rows = _body_spring_rows(sim)
    spring_pairs = sim.model.spring_pairs[body_spring_rows]
    spring_l0 = sim.model.spring_rest_lengths_m[body_spring_rows]
    spring_l = np.linalg.norm(
        positions_m[spring_pairs[:, 1]] - positions_m[spring_pairs[:, 0]],
        axis=1,
    )
    stretch = np.abs(spring_l - spring_l0) / np.maximum(spring_l0, 1e-30)

    body_bending_rows = _body_bending_rows(sim)
    theta = np.asarray(
        [
            _triplet_angle_rad(positions_m, sim.model.bending_triplets[row])
            for row in body_bending_rows
        ],
        dtype=float,
    )
    theta0, _ = sim.engine._state_angles_rad()
    theta0_body = theta0[body_bending_rows]
    err_deg = np.abs(theta - theta0_body) * (180.0 / np.pi)
    return float(np.max(stretch)), float(np.max(err_deg))


def _triangle_area(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    return 0.5 * float(np.linalg.norm(np.cross(b - a, c - a)))


def _body_triangle_area_range(
    sim: Simulator, positions_m: np.ndarray
) -> tuple[float, float]:
    areas = []
    for layer in sim.model.body_layer_indices:
        pts = positions_m[layer.astype(int, copy=False)]
        areas.append(_triangle_area(pts[0], pts[1], pts[2]))
    return float(min(areas)), float(max(areas))


def _body_centerline_max_deviation_um(sim: Simulator, positions_m: np.ndarray) -> float:
    centers = np.asarray(
        [
            np.mean(positions_m[layer.astype(int, copy=False)], axis=0)
            for layer in sim.model.body_layer_indices
        ],
        dtype=float,
    )
    p0 = centers[0]
    p1 = centers[-1]
    d = p1 - p0
    d_norm = max(float(np.linalg.norm(d)), 1e-18)
    u = d / d_norm
    proj = p0[None, :] + ((centers - p0[None, :]) @ u)[:, None] * u[None, :]
    dev = np.linalg.norm(centers - proj, axis=1)
    return float(np.max(dev) * 1.0e6)


def _make_cfg() -> SimulationConfig:
    return SimulationConfig.from_dict(
        {
            "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
            "body": {
                "prism": {
                    "n_prism": 3,
                    "dz_over_b": 0.5,
                    "radius_over_b": 0.5,
                    "axis": "x",
                },
                "length_total_um": 2.0,
            },
            "flagella": {
                "n_flagella": 0,
                "placement_mode": "uniform",
                "discretization": {"ds_over_b": 0.58},
                "bond_L_over_b": 0.58,
                "length_over_b": 2.32,
                "helix_init": {"radius_over_b": 0.25, "pitch_over_b": 2.5},
            },
            "fluid": {"viscosity_Pa_s": 1.0e-3},
            "motor": {
                "torque_Nm": 0.0,
                "reverse_n_flagella": 0,
                "enable_switching": False,
            },
            "potentials": {
                "spring": {"H_over_T_over_b": 10.0, "s": 0.1},
                "bend": {
                    "kb_over_T": 20.0,
                    "theta0_deg": {
                        "normal": 142.0,
                        "semicoiled": 90.0,
                        "curly1": 105.0,
                    },
                },
                "torsion": {
                    "kt_over_T": 10.0,
                    "phi0_deg": {
                        "normal": -60.0,
                        "semicoiled": 65.0,
                        "curly1": 120.0,
                    },
                },
                "spring_spring_repulsion": {
                    "A_ss_over_T": 1.0,
                    "a_ss_over_b": 0.2,
                    "cutoff_over_b": 0.2,
                },
            },
            "hook": {"enabled": False, "threshold_deg": 90.0, "kb_over_T": 20.0},
            "projection": {"enable_body_rigid_projection": False},
            "run_tumble": {
                "run_tau": 20.0,
                "tumble_tau": 8.0,
                "semicoiled_tau": 4.0,
                "curly1_tau": 4.0,
            },
            "time": {"duration_s": 0.008, "dt_s": 1.0e-3},
            "output_sampling": {"out_all_steps_3d": True, "fps_out_2d": 25.0},
            "brownian": {
                "enabled": False,
                "temperature_K": 298.0,
                "method": "cholesky",
                "jitter": 1.0e-20,
            },
            "render": {
                "image_size_px": 128,
                "pixel_size_um": 0.203,
                "flagella_linewidth_px": 2.0,
                "render_flagella": False,
                "render_flagella_2d": False,
                "save_frames_3d": False,
                "follow_camera_3d": True,
                "view_range_um": 8.0,
                "timestamp_3d": False,
                "label_flagella": False,
                "follow_camera_2d": False,
                "center_body_in_2d": True,
                "save_frames_2d": False,
            },
            "seed": {"global_seed": 0},
            "output": {"base_dir": "outputs"},
        }
    )


def test_body_stability_static_equilibrium_projection_off(tmp_path: Path) -> None:
    cfg = _make_cfg()
    sim = Simulator(cfg)

    init_pos = sim.model.positions_m.copy()
    init_com = _body_com(sim, init_pos)
    sim.run(cfg.time.duration_s, step_summary_dir=tmp_path / "sim")
    final_pos = sim.model.positions_m.copy()

    assert np.isfinite(final_pos).all()

    max_stretch, max_bend_err_deg = _body_shape_metrics(sim, final_pos)
    assert max_stretch < 2.0e-3
    assert max_bend_err_deg < 1.0

    com_drift_um = float(np.linalg.norm(_body_com(sim, final_pos) - init_com) * 1.0e6)
    assert com_drift_um < 2.0e-3

    body_diag = tmp_path / "sim" / "body_constraint_diagnostics.csv"
    assert body_diag.is_file()
    with body_diag.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows
    required_cols = {
        "step",
        "t_s",
        "body_spring_max_stretch_ratio",
        "body_spring_mean_stretch_ratio",
        "body_bend_max_error_deg",
        "body_bend_mean_error_deg",
        "body_triangle_area_min",
        "body_triangle_area_max",
        "body_centerline_max_deviation_um",
        "com_x_um",
        "com_y_um",
        "com_z_um",
    }
    assert required_cols.issubset(set(rows[0].keys()))


def test_body_uniform_force_translation_projection_off() -> None:
    levels = {
        "low": 2.0e-19,
        "medium": 8.0e-19,
        "high": 2.0e-18,
    }

    disp_by_level: dict[str, float] = {}
    stable_found = False

    for name, force_x in levels.items():
        cfg = _make_cfg()
        sim = Simulator(cfg)

        body_idx = sim.model.body_indices.astype(int, copy=False)

        def _uniform_force(_pos: np.ndarray, _t_star: float) -> np.ndarray:
            out = np.zeros_like(sim.model.positions_m)
            out[body_idx, 0] = force_x
            return out

        sim.engine.set_external_force_callback(_uniform_force)
        init_pos = sim.model.positions_m.copy()
        init_com = _body_com(sim, init_pos)
        init_area_min, _ = _body_triangle_area_range(sim, init_pos)

        for _ in range(80):
            sim.engine.step(cfg.dt_star)

        final_pos = sim.model.positions_m.copy()
        final_com = _body_com(sim, final_pos)

        disp_um = float(np.linalg.norm(final_com - init_com) * 1.0e6)
        disp_by_level[name] = disp_um

        area_min, _ = _body_triangle_area_range(sim, final_pos)
        centerline_dev_um = _body_centerline_max_deviation_um(sim, final_pos)
        assert np.isfinite(final_pos).all()
        if area_min > 0.6 * init_area_min and centerline_dev_um < 0.8:
            stable_found = True

    assert disp_by_level["low"] > 1.0e-4
    assert disp_by_level["medium"] > 1.0e-4
    assert disp_by_level["high"] > 1.0e-4
    ratio = max(disp_by_level.values()) / max(min(disp_by_level.values()), 1e-12)
    assert ratio > 1.2
    assert stable_found


def test_body_rotation_couple_projection_off() -> None:
    cfg = _make_cfg()
    sim = Simulator(cfg)

    first_layer = sim.model.body_layer_indices[0].astype(int, copy=False)
    last_layer = sim.model.body_layer_indices[-1].astype(int, copy=False)
    f_each = 1.0e-19

    def _couple_force(_pos: np.ndarray, _t_star: float) -> np.ndarray:
        out = np.zeros_like(sim.model.positions_m)
        out[first_layer, 1] = f_each
        out[last_layer, 1] = -f_each
        return out

    sim.engine.set_external_force_callback(_couple_force)

    init_pos = sim.model.positions_m.copy()
    init_axis = _body_axis(sim, init_pos)
    init_com = _body_com(sim, init_pos)
    init_area_min, _ = _body_triangle_area_range(sim, init_pos)

    for _ in range(120):
        sim.engine.step(cfg.dt_star)

    final_pos = sim.model.positions_m.copy()
    final_axis = _body_axis(sim, final_pos)
    final_com = _body_com(sim, final_pos)

    angle_deg = math.degrees(
        math.acos(float(np.clip(np.dot(init_axis, final_axis), -1.0, 1.0)))
    )
    com_drift_um = float(np.linalg.norm(final_com - init_com) * 1.0e6)
    area_min, _ = _body_triangle_area_range(sim, final_pos)
    centerline_dev_um = _body_centerline_max_deviation_um(sim, final_pos)

    assert np.isfinite(final_pos).all()
    assert angle_deg > 0.1
    assert com_drift_um < 1.5
    assert area_min > 0.5 * init_area_min
    assert centerline_dev_um < 1.2
