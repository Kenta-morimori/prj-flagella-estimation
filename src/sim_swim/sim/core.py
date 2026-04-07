"""bead-spring + RPY のシミュレーションラッパ。"""

from __future__ import annotations

from dataclasses import dataclass
import json
import logging
import math
from pathlib import Path
import time
from typing import Any, List, Tuple

import numpy as np

from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.debug_summary import (
    BodyConstraintDiagnosticsRecorder,
    BodyConstraintLocalDiagnosticsRecorder,
    StepSummaryRecorder,
)
from sim_swim.sim.flagella_geometry import FlagellaRig
from sim_swim.sim.params import SimulationConfig

M_TO_UM = 1e6


def _quat_normalize(q: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(q)
    if norm == 0:
        return np.array([0.0, 0.0, 0.0, 1.0], dtype=float)
    return q / norm


def _quat_multiply(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    """Hamilton積でクォータニオンを合成する。"""
    x1, y1, z1, w1 = q1
    x2, y2, z2, w2 = q2
    return np.array(
        [
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
            w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        ],
        dtype=float,
    )


def _rotate_vec(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """クォータニオンでベクトルを回転する。"""
    v_q = np.array([v[0], v[1], v[2], 0.0], dtype=float)
    q_conj = np.array([-q[0], -q[1], -q[2], q[3]], dtype=float)
    return _quat_multiply(_quat_multiply(q, v_q), q_conj)[:3]


def _quat_to_rotmat(q: np.ndarray) -> np.ndarray:
    """クォータニオンから回転行列を計算する。"""
    x, y, z, w = q
    return np.array(
        [
            [1 - 2 * (y**2 + z**2), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x**2 + z**2), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x**2 + y**2)],
        ],
        dtype=float,
    )


def _quat_from_two_vectors(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a_n = a / max(np.linalg.norm(a), 1e-12)
    b_n = b / max(np.linalg.norm(b), 1e-12)
    dot = float(np.clip(np.dot(a_n, b_n), -1.0, 1.0))

    if dot < -0.999999:
        axis = np.array([0.0, 1.0, 0.0], dtype=float)
        if abs(a_n[1]) > 0.9:
            axis = np.array([0.0, 0.0, 1.0], dtype=float)
        axis = axis - np.dot(axis, a_n) * a_n
        axis = axis / max(np.linalg.norm(axis), 1e-12)
        return np.array([axis[0], axis[1], axis[2], 0.0], dtype=float)

    cross = np.cross(a_n, b_n)
    q = np.array([cross[0], cross[1], cross[2], 1.0 + dot], dtype=float)
    return _quat_normalize(q)


def _omega_from_quats(q_prev: np.ndarray, q_curr: np.ndarray, dt: float) -> np.ndarray:
    if dt <= 0:
        return np.zeros(3, dtype=float)
    q_prev_conj = np.array([-q_prev[0], -q_prev[1], -q_prev[2], q_prev[3]], dtype=float)
    dq = _quat_multiply(q_curr, q_prev_conj)
    dq = _quat_normalize(dq)

    if dq[3] < 0.0:
        dq = -dq

    v = dq[:3]
    w = float(np.clip(dq[3], -1.0, 1.0))
    nv = np.linalg.norm(v)
    if nv < 1e-12:
        return np.zeros(3, dtype=float)

    axis = v / nv
    angle = 2.0 * math.atan2(nv, w)
    return axis * (angle / dt)


def _triplet_angle_deg(positions_m: np.ndarray, triplet: np.ndarray) -> float:
    i, j, k = (int(triplet[0]), int(triplet[1]), int(triplet[2]))
    u = positions_m[i] - positions_m[j]
    v = positions_m[k] - positions_m[j]
    nu = max(float(np.linalg.norm(u)), 1e-18)
    nv = max(float(np.linalg.norm(v)), 1e-18)
    c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return float(np.rad2deg(np.arccos(c)))


def _dihedral_angle_deg(positions_m: np.ndarray, quad: np.ndarray) -> float:
    a_i, b_i, c_i, d_i = (int(quad[0]), int(quad[1]), int(quad[2]), int(quad[3]))
    a = positions_m[a_i]
    b = positions_m[b_i]
    c = positions_m[c_i]
    d = positions_m[d_i]

    b0 = a - b
    b1 = c - b
    b2 = d - c
    b1n = b1 / max(float(np.linalg.norm(b1)), 1e-18)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1n, v), w))
    return float(np.rad2deg(np.arctan2(y, x)))


def _wrap_deg(deg: float) -> float:
    return float((deg + 180.0) % 360.0 - 180.0)


def _estimate_helix_radius_pitch_over_b(
    points_m: np.ndarray,
    b_m: float,
) -> tuple[float, float]:
    centered = points_m - np.mean(points_m, axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    axis = vh[0]
    axis = axis / max(float(np.linalg.norm(axis)), 1e-18)
    if float(np.dot(axis, points_m[-1] - points_m[0])) < 0.0:
        axis = -axis

    ref = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(ref, axis))) > 0.9:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)
    e1 = np.cross(axis, ref)
    e1 = e1 / max(float(np.linalg.norm(e1)), 1e-18)
    e2 = np.cross(axis, e1)

    origin = points_m[0]
    rel = points_m - origin
    s = rel @ axis
    u = rel @ e1
    v = rel @ e2

    mat = np.column_stack([u, v, np.ones_like(u)])
    rhs = -(u * u + v * v)
    coef, *_ = np.linalg.lstsq(mat, rhs, rcond=None)
    a, b, c = coef
    cx = -0.5 * a
    cy = -0.5 * b
    radius = np.sqrt(max(cx * cx + cy * cy - c, 0.0))

    phase = np.unwrap(np.arctan2(v - cy, u - cx))
    if phase.size < 2 or float(np.max(phase) - np.min(phase)) < 1e-9:
        return float(radius / max(b_m, 1e-18)), float("nan")
    slope, _ = np.polyfit(phase, s, 1)
    pitch = abs(float(slope)) * 2.0 * np.pi
    return float(radius / max(b_m, 1e-18)), float(pitch / max(b_m, 1e-18))


@dataclass
class SimulationState:
    """1時刻分の観測互換データ。"""

    t: float
    position_um: Tuple[float, float, float]
    quaternion: Tuple[float, float, float, float]
    velocity_um_s: Tuple[float, float, float]
    omega_rad_s: Tuple[float, float, float]
    bead_positions_um: np.ndarray
    flag_states: Tuple[int, ...] = ()
    reverse_flagella: Tuple[int, ...] = ()


class Simulator:
    """時間ループと出力変換を担う薄いラッパ。"""

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.config.validate_time_scaling()
        self.model = ModelBuilder(config).build()
        self.engine = DynamicsEngine(self.model, config)
        self.rig = FlagellaRig(
            body_layer_indices=[arr.copy() for arr in self.model.body_layer_indices],
            body_ring_edges=self.model.body_ring_edges.copy(),
            body_vertical_edges=self.model.body_vertical_edges.copy(),
            flagella_indices=[arr.copy() for arr in self.model.flagella_indices],
        )
        self.initial_geometry_summary = self._build_initial_geometry_summary()

    def _build_initial_geometry_summary(self) -> dict[str, Any]:
        pos = self.model.positions_m
        b_m = max(self.config.b_m, 1e-30)

        first_layer = np.mean(pos[self.model.body_layer_indices[0]], axis=0)
        last_layer = np.mean(pos[self.model.body_layer_indices[-1]], axis=0)
        body_axis = last_layer - first_layer
        body_axis = body_axis / max(float(np.linalg.norm(body_axis)), 1e-18)
        rear_dir = -body_axis

        bend_target = float(self.config.potentials.bend.theta0_deg["normal"])
        torsion_target = float(self.config.potentials.torsion.phi0_deg["normal"])

        summary: dict[str, Any] = {
            "flagella": {
                "init_mode": str(self.config.flagella.init_mode),
                "stub_mode": str(self.config.flagella.stub_mode),
                "n_flagella": int(self.config.flagella.n_flagella),
                "n_beads_per_flagellum": int(
                    self.config.flagella.n_beads_per_flagellum
                ),
                "bond_L_over_b": float(self.config.flagella.bond_L_over_b),
                "target_theta0_deg_normal": bend_target,
                "target_phi0_deg_normal": torsion_target,
                "source_of_truth": (
                    [
                        "bond_L_over_b",
                        "potentials.bend.theta0_deg.normal",
                        "potentials.torsion.phi0_deg.normal",
                        "n_beads_per_flagellum",
                        "n_flagella",
                    ]
                    if self.config.flagella.init_mode == "paper_table1"
                    else [
                        "helix_init.radius_over_b",
                        "helix_init.pitch_over_b",
                        "bond_L_over_b",
                        "length_over_b",
                        "n_flagella",
                    ]
                ),
            },
            "per_flagellum": [],
        }

        for f_id, idxs in enumerate(self.model.flagella_indices):
            idx = idxs.astype(int, copy=False)
            pts = pos[idx]
            bonds = np.linalg.norm(pts[1:] - pts[:-1], axis=1)
            contour_len_um = float(np.sum(bonds) * M_TO_UM)
            end_to_end_len_um = float(np.linalg.norm(pts[-1] - pts[0]) * M_TO_UM)
            bond_mean_um = float(np.mean(bonds) * M_TO_UM)
            bond_min_um = float(np.min(bonds) * M_TO_UM)
            bond_max_um = float(np.max(bonds) * M_TO_UM)
            radius_over_b, pitch_over_b = _estimate_helix_radius_pitch_over_b(
                pts,
                b_m,
            )

            attach_first_len_um = float("nan")
            hook_angle_deg = float("nan")
            for triplet in self.model.hook_triplets:
                a, first, second = (int(triplet[0]), int(triplet[1]), int(triplet[2]))
                if first == int(idx[0]) and second == int(idx[1]):
                    attach_first_len_um = float(
                        np.linalg.norm(pos[first] - pos[a]) * M_TO_UM
                    )
                    hook_angle_deg = _triplet_angle_deg(pos, triplet)
                    break

            bend_rows = np.where(self.model.bending_flag_ids == f_id)[0]
            if bend_rows.size > 0:
                bend_deg = np.asarray(
                    [
                        _triplet_angle_deg(pos, self.model.bending_triplets[row])
                        for row in bend_rows
                    ],
                    dtype=float,
                )
                bend_err = np.abs(bend_deg - bend_target)
            else:
                bend_deg = np.asarray([float("nan")], dtype=float)
                bend_err = np.asarray([float("nan")], dtype=float)

            torsion_rows = np.where(self.model.torsion_flag_ids == f_id)[0]
            if torsion_rows.size > 0:
                torsion_deg = np.asarray(
                    [
                        _dihedral_angle_deg(pos, self.model.torsion_quads[row])
                        for row in torsion_rows
                    ],
                    dtype=float,
                )
                torsion_err = np.asarray(
                    [abs(_wrap_deg(float(v - torsion_target))) for v in torsion_deg],
                    dtype=float,
                )
            else:
                torsion_deg = np.asarray([float("nan")], dtype=float)
                torsion_err = np.asarray([float("nan")], dtype=float)

            tangent0 = pts[1] - pts[0]
            tangent0 = tangent0 / max(float(np.linalg.norm(tangent0)), 1e-18)
            tangent_vs_rear_deg = float(
                np.rad2deg(
                    np.arccos(float(np.clip(np.dot(tangent0, rear_dir), -1.0, 1.0)))
                )
            )

            summary["per_flagellum"].append(
                {
                    "flag_id": int(f_id),
                    "bead_count": int(idx.size),
                    "contour_length_um": contour_len_um,
                    "end_to_end_length_um": end_to_end_len_um,
                    "attach_first_length_um": attach_first_len_um,
                    "bond_length_mean_um": bond_mean_um,
                    "bond_length_min_um": bond_min_um,
                    "bond_length_max_um": bond_max_um,
                    "derived_helix_radius_over_b": radius_over_b,
                    "derived_helix_pitch_over_b": pitch_over_b,
                    "initial_hook_angle_deg": hook_angle_deg,
                    "initial_bend_angle_mean_deg": float(np.mean(bend_deg)),
                    "initial_bend_angle_max_deg": float(np.max(bend_deg)),
                    "initial_bend_err_mean_deg": float(np.mean(bend_err)),
                    "initial_bend_err_max_deg": float(np.max(bend_err)),
                    "initial_torsion_angle_mean_deg": (
                        float(
                            np.rad2deg(
                                np.arctan2(
                                    np.mean(np.sin(np.deg2rad(torsion_deg))),
                                    np.mean(np.cos(np.deg2rad(torsion_deg))),
                                )
                            )
                        )
                        if np.isfinite(torsion_deg).all()
                        else float("nan")
                    ),
                    "initial_torsion_angle_max_deg": float(np.max(torsion_deg)),
                    "initial_torsion_err_mean_deg": float(np.mean(torsion_err)),
                    "initial_torsion_err_max_deg": float(np.max(torsion_err)),
                    "initial_tangent_vs_rear_direction_angle_deg": tangent_vs_rear_deg,
                }
            )

        return summary

    def _observe(self, t: float, prev: SimulationState | None) -> SimulationState:
        body_pts = self.model.positions_m[self.model.body_indices]
        center_m = body_pts.mean(axis=0)

        first_layer = self.model.positions_m[self.model.body_layer_indices[0]].mean(
            axis=0
        )
        last_layer = self.model.positions_m[self.model.body_layer_indices[-1]].mean(
            axis=0
        )
        axis = last_layer - first_layer
        if np.linalg.norm(axis) < 1e-15:
            axis = np.array([1.0, 0.0, 0.0], dtype=float)
        q = _quat_from_two_vectors(np.array([1.0, 0.0, 0.0], dtype=float), axis)

        if prev is None:
            vel_um_s = np.zeros(3, dtype=float)
            omega = np.zeros(3, dtype=float)
        else:
            dt = max(t - prev.t, 1e-12)
            prev_center_m = np.array(prev.position_um, dtype=float) / M_TO_UM
            vel_um_s = ((center_m - prev_center_m) / dt) * M_TO_UM
            prev_q = np.array(prev.quaternion, dtype=float)
            omega = _omega_from_quats(prev_q, q, dt)

        return SimulationState(
            t=round(t, 9),
            position_um=tuple((center_m * M_TO_UM).tolist()),
            quaternion=tuple(q.tolist()),
            velocity_um_s=tuple(vel_um_s.tolist()),
            omega_rad_s=tuple(omega.tolist()),
            bead_positions_um=self.model.positions_m.copy() * M_TO_UM,
            flag_states=tuple(int(s) for s in self.model.flag_states.tolist()),
            reverse_flagella=tuple(
                int(i) for i in self.model.reverse_flagella.tolist()
            ),
        )

    def run(
        self,
        duration_s: float,
        logger: logging.Logger | None = None,
        progress_interval: int | None = None,
        step_summary_dir: Path | None = None,
    ) -> List[SimulationState]:
        """与えた時間だけシミュレーションして状態列を返す。

        Args:
            duration_s: シミュレーション時間 [s]。
            logger: 進捗を出力するロガー（任意）。
            progress_interval: 進捗ログのステップ間隔。None なら自動設定。
        """

        tau_s = self.config.tau_s
        dt_star = max(self.config.dt_star, 1e-12)
        duration_star = max(float(duration_s), 0.0) / max(tau_s, 1e-30)
        total_steps = max(1, int(math.ceil(duration_star / dt_star)))

        states: List[SimulationState] = []
        wall_start = time.perf_counter()

        if progress_interval is None:
            progress_interval = 1000 if total_steps >= 10000 else 100
        progress_interval = max(1, int(progress_interval))

        if logger is not None:
            logger.info(
                (
                    "Simulation loop start: total_steps=%d, dt_star=%.6e, "
                    "dt_s=%.6e s, duration_star=%.6e, progress_interval=%d"
                ),
                total_steps,
                dt_star,
                self.config.dt_s,
                duration_star,
                progress_interval,
            )

        prev = None
        states.append(self._observe(0.0, prev))
        debug_recorder = (
            StepSummaryRecorder(self.model, self.config, step_summary_dir)
            if step_summary_dir is not None
            else None
        )
        if step_summary_dir is not None:
            step_summary_dir.mkdir(parents=True, exist_ok=True)
            initial_summary_path = step_summary_dir / "initial_geometry_summary.json"
            initial_summary_path.write_text(
                json.dumps(self.initial_geometry_summary, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
            if logger is not None:
                logger.info(
                    "Saved initial geometry summary to %s", initial_summary_path
                )
        body_diag_recorder = (
            BodyConstraintDiagnosticsRecorder(self.model, self.config, step_summary_dir)
            if step_summary_dir is not None and duration_s <= 0.01
            else None
        )
        body_local_diag_recorder = (
            BodyConstraintLocalDiagnosticsRecorder(self.model, step_summary_dir)
            if (step_summary_dir is not None and duration_s <= 0.01)
            else None
        )

        for step in range(total_steps):
            t_star_before = self.engine.t_star
            step_diag = self.engine.step(dt_star)

            if debug_recorder is not None:
                debug_recorder.record(step=step, t_star=t_star_before, diag=step_diag)
            if body_diag_recorder is not None:
                body_diag_recorder.record(
                    step=step,
                    t_s=t_star_before * tau_s,
                    diag=step_diag,
                )
            if body_local_diag_recorder is not None:
                body_local_diag_recorder.record(
                    step=step,
                    t_s=t_star_before * tau_s,
                    pos_after=step_diag.positions_after_m,
                )

            t_now = self.engine.t_star * tau_s
            prev = states[-1]
            states.append(self._observe(t_now, prev))

            completed = step + 1
            if logger is not None and (
                completed % progress_interval == 0 or completed == total_steps
            ):
                logger.info(
                    (
                        "Simulation progress: step=%d/%d (%.1f%%), "
                        "t_star=%.6f, t_s=%.6f s"
                    ),
                    completed,
                    total_steps,
                    (completed / total_steps) * 100.0,
                    self.engine.t_star,
                    self.engine.t_star * tau_s,
                )

        if logger is not None:
            logger.info(
                "Simulation loop end: total_steps=%d, elapsed=%.2f s",
                total_steps,
                time.perf_counter() - wall_start,
            )

        if debug_recorder is not None:
            step_csv = debug_recorder.write_csv()
            if logger is not None:
                logger.info("Saved step summary to %s", step_csv)
        if body_diag_recorder is not None:
            body_csv = body_diag_recorder.write_csv()
            if logger is not None:
                logger.info("Saved body diagnostics to %s", body_csv)
        if body_local_diag_recorder is not None:
            body_local_csv = body_local_diag_recorder.write_csv()
            if logger is not None:
                logger.info("Saved body local diagnostics to %s", body_local_csv)

        return states
