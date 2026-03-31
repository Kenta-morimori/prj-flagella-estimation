"""シミュレーション診断CSVの生成。"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

import numpy as np

from sim_swim.dynamics.engine import StepDiagnostics
from sim_swim.model.types import PolymorphState, SimModel
from sim_swim.sim.params import SimulationConfig

M_TO_UM = 1.0e6

STEP_SUMMARY_COLUMNS = [
    "step",
    "t_star",
    "dt_star",
    "t_s",
    "dt_s",
    "tau_s",
    "dt_internal_s",
    "pos_all_finite",
    "any_nan",
    "any_inf",
    "mean_disp_um",
    "max_disp_um",
    "bond_count_body_body",
    "bond_len_mean_body_body_um",
    "flag_intra_count",
    "bond_len_mean_flag_intra_um",
    "flag_intra_len_mean_um",
    "flag_intra_len_min_um",
    "flag_intra_len_max_um",
    "hook_count",
    "bond_len_mean_body_flag_um",
    "hook_len_mean_um",
    "hook_len_min_um",
    "hook_len_max_um",
    "hook_len_mean_over_b",
    "hook_len_min_over_b",
    "hook_len_max_over_b",
    "hook_len_rel_err_mean",
    "hook_len_rel_err_max",
    "hook_angle_mean_deg",
    "hook_angle_min_deg",
    "hook_angle_max_deg",
    "hook_angle_err_mean_deg",
    "hook_angle_err_max_deg",
    "flag_bend_err_mean_deg",
    "flag_bend_err_max_deg",
    "flag_torsion_err_mean_deg",
    "flag_torsion_err_max_deg",
    "flag_bond_len_mean_over_b",
    "flag_bond_len_min_over_b",
    "flag_bond_len_max_over_b",
    "flag_bond_rel_err_mean",
    "flag_bond_rel_err_max",
    "F_total_mean_body",
    "F_total_mean_flag",
    "F_total_mean_all",
    "F_motor_mean_body",
    "F_motor_mean_flag",
    "F_spring_mean_body",
    "F_spring_mean_flag",
    "F_bend_mean_body",
    "F_bend_mean_flag",
    "F_torsion_mean_body",
    "F_torsion_mean_flag",
    "F_repulsion_mean_body",
    "F_repulsion_mean_flag",
    "F_hook_mean_body",
    "F_hook_mean_flag",
    "torque_for_forces_Nm",
    "motor_torque_Nm",
    "torsion_fd_eps_m",
    "torsion_fd_eps_over_b",
    "motor_degenerate_axis_count",
    "motor_split_rank_deficient_count",
    "motor_bond_length_clipped_count",
    "flag_state_changed",
    "projection_body_rigid_enabled",
    "projection_hook_length_enabled",
    "projection_basal_link_direction_enabled",
    "projection_flagella_chain_length_enabled",
    "projection_flagella_template_enabled",
    "brownian_enabled",
    "brownian_disp_mean_um",
]

BODY_CONSTRAINT_DIAGNOSTICS_COLUMNS = [
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
    "F_spring_mean_body",
    "F_bend_mean_body",
    "F_torsion_mean_body",
    "F_repulsion_mean_body",
    "F_total_mean_body",
]

BODY_CONSTRAINT_LOCAL_DIAGNOSTICS_COLUMNS = [
    "step",
    "t_s",
    "spring_pair_i",
    "spring_pair_j",
    "rest_um",
    "dist_um",
    "stretch_ratio",
    "triplet_i",
    "triplet_j",
    "triplet_k",
    "angle_deg",
    "theta0_deg",
    "angle_error_deg",
    "layer_idx",
    "face_idx",
    "triangle_area",
]

RAD_TO_DEG = 180.0 / np.pi


def _wrap_angle(rad: np.ndarray | float) -> np.ndarray | float:
    return (np.asarray(rad) + np.pi) % (2.0 * np.pi) - np.pi


def _triplet_angles_rad(positions_m: np.ndarray, triplets: np.ndarray) -> np.ndarray:
    if triplets.size == 0:
        return np.zeros((0,), dtype=float)
    i = triplets[:, 0]
    j = triplets[:, 1]
    k = triplets[:, 2]
    u = positions_m[i] - positions_m[j]
    v = positions_m[k] - positions_m[j]
    nu = np.linalg.norm(u, axis=1)
    nv = np.linalg.norm(v, axis=1)
    denom = np.maximum(nu * nv, 1e-18)
    cos_t = np.sum(u * v, axis=1) / denom
    cos_t = np.clip(cos_t, -1.0, 1.0)
    return np.arccos(cos_t)


def _torsion_angles_rad(positions_m: np.ndarray, quads: np.ndarray) -> np.ndarray:
    if quads.size == 0:
        return np.zeros((0,), dtype=float)
    out = np.zeros((quads.shape[0],), dtype=float)
    for idx, (a_i, b_i, c_i, d_i) in enumerate(quads):
        a = positions_m[int(a_i)]
        b = positions_m[int(b_i)]
        c = positions_m[int(c_i)]
        d = positions_m[int(d_i)]

        b0 = a - b
        b1 = c - b
        b2 = d - c
        b1n = b1 / max(float(np.linalg.norm(b1)), 1e-18)
        v = b0 - np.dot(b0, b1n) * b1n
        w = b2 - np.dot(b2, b1n) * b1n
        x = float(np.dot(v, w))
        y = float(np.dot(np.cross(b1n, v), w))
        out[idx] = float(np.arctan2(y, x))
    return out


def _mean_norm(forces: np.ndarray, mask: np.ndarray) -> float:
    if forces.size == 0:
        return 0.0
    selected = forces[mask]
    if selected.size == 0:
        return 0.0
    return float(np.mean(np.linalg.norm(selected, axis=1)))


def _pair_distance_stats_um(
    positions_m: np.ndarray, spring_pairs: np.ndarray, pair_rows: np.ndarray
) -> tuple[int, float, float, float]:
    if pair_rows.size == 0:
        return 0, float("nan"), float("nan"), float("nan")
    pairs = spring_pairs[pair_rows]
    diff = positions_m[pairs[:, 0]] - positions_m[pairs[:, 1]]
    dist_um = np.linalg.norm(diff, axis=1) * M_TO_UM
    return (
        int(pair_rows.size),
        float(np.mean(dist_um)),
        float(np.min(dist_um)),
        float(np.max(dist_um)),
    )


def _pair_rel_error_stats(
    positions_m: np.ndarray,
    spring_pairs: np.ndarray,
    spring_rests_m: np.ndarray,
    pair_rows: np.ndarray,
) -> tuple[float, float]:
    if pair_rows.size == 0:
        return float("nan"), float("nan")
    pairs = spring_pairs[pair_rows]
    rests = np.maximum(spring_rests_m[pair_rows], 1e-30)
    lens = np.linalg.norm(positions_m[pairs[:, 1]] - positions_m[pairs[:, 0]], axis=1)
    rel = np.abs(lens - rests) / rests
    return float(np.mean(rel)), float(np.max(rel))


def _triangle_area(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    return 0.5 * float(np.linalg.norm(np.cross(b - a, c - a)))


def _point_to_line_distance(
    points: np.ndarray, p0: np.ndarray, p1: np.ndarray
) -> np.ndarray:
    d = p1 - p0
    d_norm = float(np.linalg.norm(d))
    if d_norm <= 1e-18:
        return np.linalg.norm(points - p0[None, :], axis=1)
    u = d / d_norm
    proj = p0[None, :] + ((points - p0[None, :]) @ u)[:, None] * u[None, :]
    return np.linalg.norm(points - proj, axis=1)


@dataclass
class BodyConstraintDiagnosticsRecorder:
    """body-only安定性評価用の診断CSVを保存する。"""

    model: SimModel
    cfg: SimulationConfig
    out_dir: Path
    _csv_fp: TextIO | None = field(init=False, default=None, repr=False)
    _writer: csv.DictWriter | None = field(init=False, default=None, repr=False)

    def __post_init__(self) -> None:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.body_diag_path = self.out_dir / "body_constraint_diagnostics.csv"
        self._csv_fp = self.body_diag_path.open("w", encoding="utf-8", newline="")
        self._writer = csv.DictWriter(
            self._csv_fp, fieldnames=BODY_CONSTRAINT_DIAGNOSTICS_COLUMNS
        )
        self._writer.writeheader()

        self.body_indices = self.model.body_indices.astype(int, copy=False)
        self.body_mask = self.model.bead_is_body.astype(bool)
        self.spring_pairs = self.model.spring_pairs.astype(int, copy=False)
        self.spring_rests = self.model.spring_rest_lengths_m.astype(float, copy=False)
        if self.spring_pairs.size == 0:
            self.body_spring_rows = np.zeros((0,), dtype=int)
        else:
            bi = self.model.bead_is_body[self.spring_pairs[:, 0]]
            bj = self.model.bead_is_body[self.spring_pairs[:, 1]]
            self.body_spring_rows = np.where(bi & bj)[0]
        self.body_bending_rows = np.where(self.model.bending_flag_ids < 0)[0]
        self.body_theta0 = _triplet_angles_rad(
            self.model.positions_m,
            self.model.bending_triplets[self.body_bending_rows],
        )
        self.body_layers = [
            layer.astype(int, copy=False) for layer in self.model.body_layer_indices
        ]

    def record(self, step: int, t_s: float, diag: StepDiagnostics) -> None:
        if self._writer is None or self._csv_fp is None:
            raise RuntimeError("body diagnostics writer is not initialized")

        pos_after = diag.positions_after_m

        if self.body_spring_rows.size > 0:
            pairs = self.spring_pairs[self.body_spring_rows]
            rests = np.maximum(self.spring_rests[self.body_spring_rows], 1e-30)
            lens = np.linalg.norm(
                pos_after[pairs[:, 1]] - pos_after[pairs[:, 0]], axis=1
            )
            stretch = np.abs(lens - rests) / rests
            body_spring_max_stretch_ratio = float(np.max(stretch))
            body_spring_mean_stretch_ratio = float(np.mean(stretch))
        else:
            body_spring_max_stretch_ratio = 0.0
            body_spring_mean_stretch_ratio = 0.0

        if self.body_bending_rows.size > 0:
            theta = _triplet_angles_rad(
                pos_after,
                self.model.bending_triplets[self.body_bending_rows],
            )
            err_deg = np.abs(theta - self.body_theta0) * RAD_TO_DEG
            body_bend_max_error_deg = float(np.max(err_deg))
            body_bend_mean_error_deg = float(np.mean(err_deg))
        else:
            body_bend_max_error_deg = 0.0
            body_bend_mean_error_deg = 0.0

        tri_areas = []
        centroids = []
        for layer in self.body_layers:
            pts = pos_after[layer]
            if pts.shape[0] >= 3:
                tri_areas.append(_triangle_area(pts[0], pts[1], pts[2]))
            centroids.append(np.mean(pts, axis=0))
        if tri_areas:
            body_triangle_area_min = float(np.min(tri_areas))
            body_triangle_area_max = float(np.max(tri_areas))
        else:
            body_triangle_area_min = 0.0
            body_triangle_area_max = 0.0

        centers = np.asarray(centroids, dtype=float)
        if centers.shape[0] >= 2:
            d = _point_to_line_distance(centers, centers[0], centers[-1])
            body_centerline_max_deviation_um = float(np.max(d)) * M_TO_UM
        else:
            body_centerline_max_deviation_um = 0.0

        body_pos = pos_after[self.body_indices]
        com = np.mean(body_pos, axis=0) * M_TO_UM

        row = {
            "step": int(step),
            "t_s": float(t_s),
            "body_spring_max_stretch_ratio": body_spring_max_stretch_ratio,
            "body_spring_mean_stretch_ratio": body_spring_mean_stretch_ratio,
            "body_bend_max_error_deg": body_bend_max_error_deg,
            "body_bend_mean_error_deg": body_bend_mean_error_deg,
            "body_triangle_area_min": body_triangle_area_min,
            "body_triangle_area_max": body_triangle_area_max,
            "body_centerline_max_deviation_um": body_centerline_max_deviation_um,
            "com_x_um": float(com[0]),
            "com_y_um": float(com[1]),
            "com_z_um": float(com[2]),
            "F_spring_mean_body": _mean_norm(diag.spring_forces, self.body_mask),
            "F_bend_mean_body": _mean_norm(diag.bend_forces, self.body_mask),
            "F_torsion_mean_body": _mean_norm(diag.torsion_forces, self.body_mask),
            "F_repulsion_mean_body": _mean_norm(
                diag.repulsion_forces,
                self.body_mask,
            ),
            "F_total_mean_body": _mean_norm(diag.total_forces, self.body_mask),
        }
        self._writer.writerow(row)
        self._csv_fp.flush()

    def write_csv(self) -> Path:
        if self._csv_fp is not None:
            self._csv_fp.close()
            self._csv_fp = None
            self._writer = None
        return self.body_diag_path


@dataclass
class BodyConstraintLocalDiagnosticsRecorder:
    """body局所崩れ要因の診断CSVを保存する。"""

    model: SimModel
    out_dir: Path
    _csv_fp: TextIO | None = field(init=False, default=None, repr=False)
    _writer: csv.DictWriter | None = field(init=False, default=None, repr=False)

    def __post_init__(self) -> None:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.body_local_diag_path = (
            self.out_dir / "body_constraint_local_diagnostics.csv"
        )
        self._csv_fp = self.body_local_diag_path.open("w", encoding="utf-8", newline="")
        self._writer = csv.DictWriter(
            self._csv_fp,
            fieldnames=BODY_CONSTRAINT_LOCAL_DIAGNOSTICS_COLUMNS,
        )
        self._writer.writeheader()

        spring_pairs = self.model.spring_pairs.astype(int, copy=False)
        spring_rests = self.model.spring_rest_lengths_m.astype(float, copy=False)
        if spring_pairs.size == 0:
            self.body_spring_rows = np.zeros((0,), dtype=int)
        else:
            bi = self.model.bead_is_body[spring_pairs[:, 0]]
            bj = self.model.bead_is_body[spring_pairs[:, 1]]
            self.body_spring_rows = np.where(bi & bj)[0]
        self.body_spring_pairs = spring_pairs
        self.body_spring_rests = spring_rests

        self.body_bending_rows = np.where(self.model.bending_flag_ids < 0)[0]
        self.body_triplets = self.model.bending_triplets.astype(int, copy=False)
        self.body_theta0_deg = np.rad2deg(
            _triplet_angles_rad(
                self.model.positions_m,
                self.body_triplets[self.body_bending_rows],
            )
        )
        self.body_layers = [
            layer.astype(int, copy=False) for layer in self.model.body_layer_indices
        ]

    def record(self, step: int, t_s: float, pos_after: np.ndarray) -> None:
        if self._writer is None or self._csv_fp is None:
            raise RuntimeError("body local diagnostics writer is not initialized")

        if self.body_spring_rows.size > 0:
            pairs = self.body_spring_pairs[self.body_spring_rows]
            rests = np.maximum(self.body_spring_rests[self.body_spring_rows], 1e-30)
            dist_m = np.linalg.norm(
                pos_after[pairs[:, 1]] - pos_after[pairs[:, 0]], axis=1
            )
            stretch = np.abs(dist_m - rests) / rests
            for (i_raw, j_raw), rest_m, dist_val, stretch_val in zip(
                pairs,
                rests,
                dist_m,
                stretch,
            ):
                self._writer.writerow(
                    {
                        "step": int(step),
                        "t_s": float(t_s),
                        "spring_pair_i": int(i_raw),
                        "spring_pair_j": int(j_raw),
                        "rest_um": float(rest_m * M_TO_UM),
                        "dist_um": float(dist_val * M_TO_UM),
                        "stretch_ratio": float(stretch_val),
                        "triplet_i": "",
                        "triplet_j": "",
                        "triplet_k": "",
                        "angle_deg": "",
                        "theta0_deg": "",
                        "angle_error_deg": "",
                        "layer_idx": "",
                        "face_idx": "",
                        "triangle_area": "",
                    }
                )

        if self.body_bending_rows.size > 0:
            triplets = self.body_triplets[self.body_bending_rows]
            theta_deg = np.rad2deg(_triplet_angles_rad(pos_after, triplets))
            err_deg = np.abs(theta_deg - self.body_theta0_deg)
            for (i_raw, j_raw, k_raw), theta_val, theta0_val, err_val in zip(
                triplets,
                theta_deg,
                self.body_theta0_deg,
                err_deg,
            ):
                self._writer.writerow(
                    {
                        "step": int(step),
                        "t_s": float(t_s),
                        "spring_pair_i": "",
                        "spring_pair_j": "",
                        "rest_um": "",
                        "dist_um": "",
                        "stretch_ratio": "",
                        "triplet_i": int(i_raw),
                        "triplet_j": int(j_raw),
                        "triplet_k": int(k_raw),
                        "angle_deg": float(theta_val),
                        "theta0_deg": float(theta0_val),
                        "angle_error_deg": float(err_val),
                        "layer_idx": "",
                        "face_idx": "",
                        "triangle_area": "",
                    }
                )

        for layer_idx, layer in enumerate(self.body_layers):
            pts = pos_after[layer]
            if pts.shape[0] < 3:
                continue
            area = _triangle_area(pts[0], pts[1], pts[2])
            self._writer.writerow(
                {
                    "step": int(step),
                    "t_s": float(t_s),
                    "spring_pair_i": "",
                    "spring_pair_j": "",
                    "rest_um": "",
                    "dist_um": "",
                    "stretch_ratio": "",
                    "triplet_i": "",
                    "triplet_j": "",
                    "triplet_k": "",
                    "angle_deg": "",
                    "theta0_deg": "",
                    "angle_error_deg": "",
                    "layer_idx": int(layer_idx),
                    "face_idx": 0,
                    "triangle_area": float(area),
                }
            )

        self._csv_fp.flush()

    def write_csv(self) -> Path:
        if self._csv_fp is not None:
            self._csv_fp.close()
            self._csv_fp = None
            self._writer = None
        return self.body_local_diag_path


@dataclass
class StepSummaryRecorder:
    """ステップごとの診断指標をCSVへ保存する。"""

    model: SimModel
    cfg: SimulationConfig
    out_dir: Path
    _csv_fp: TextIO | None = field(init=False, default=None, repr=False)
    _writer: csv.DictWriter | None = field(init=False, default=None, repr=False)

    def __post_init__(self) -> None:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.step_summary_path = self.out_dir / "step_summary.csv"
        self._csv_fp = self.step_summary_path.open("w", encoding="utf-8", newline="")
        self._writer = csv.DictWriter(self._csv_fp, fieldnames=STEP_SUMMARY_COLUMNS)
        self._writer.writeheader()

        self.body_mask = self.model.bead_is_body.astype(bool)
        self.flag_mask = ~self.body_mask
        self.spring_pairs = self.model.spring_pairs.astype(int, copy=False)
        self.spring_rests_m = self.model.spring_rest_lengths_m.astype(float, copy=False)
        self.flag_bending_rows = np.where(self.model.bending_flag_ids >= 0)[0]
        self.flag_torsion_rows = np.where(self.model.torsion_flag_ids >= 0)[0]
        self.bend_map = self.cfg.potentials.bend.theta0_deg or {
            "normal": 142.0,
            "semicoiled": 90.0,
            "curly1": 105.0,
        }
        self.torsion_map = self.cfg.potentials.torsion.phi0_deg or {
            "normal": -60.0,
            "semicoiled": 65.0,
            "curly1": 120.0,
        }
        self.prev_flag_states = self.model.flag_states.copy()

        if self.spring_pairs.size == 0:
            self.body_body_rows = np.zeros((0,), dtype=int)
            self.flag_intra_rows = np.zeros((0,), dtype=int)
            self.body_flag_rows = np.zeros((0,), dtype=int)
            return

        i = self.spring_pairs[:, 0]
        j = self.spring_pairs[:, 1]
        fi = self.model.bead_flag_ids[i]
        fj = self.model.bead_flag_ids[j]
        bi = self.body_mask[i]
        bj = self.body_mask[j]

        self.body_body_rows = np.where(bi & bj)[0]
        self.body_flag_rows = np.where(np.logical_xor(bi, bj))[0]
        self.flag_intra_rows = np.where((~bi) & (~bj) & (fi == fj) & (fi >= 0))[0]

    def record(self, step: int, t_star: float, diag: StepDiagnostics) -> None:
        pos_after = diag.positions_after_m
        disp_um = np.linalg.norm(pos_after - diag.positions_before_m, axis=1) * M_TO_UM
        any_nan = bool(np.isnan(pos_after).any())
        any_inf = bool(np.isinf(pos_after).any())
        pos_all_finite = bool(np.isfinite(pos_after).all())

        brownian_disp_mean_um = float("nan")
        if diag.brownian_enabled:
            brownian_disp = np.linalg.norm(diag.brownian_disp_m, axis=1) * M_TO_UM
            brownian_disp_mean_um = float(np.mean(brownian_disp))

        (
            bond_count_body_body,
            bond_len_mean_body_body_um,
            _bond_len_min_body_body_um,
            _bond_len_max_body_body_um,
        ) = _pair_distance_stats_um(pos_after, self.spring_pairs, self.body_body_rows)
        (
            flag_intra_count,
            flag_intra_len_mean_um,
            flag_intra_len_min_um,
            flag_intra_len_max_um,
        ) = _pair_distance_stats_um(pos_after, self.spring_pairs, self.flag_intra_rows)
        (
            hook_count,
            hook_len_mean_um,
            hook_len_min_um,
            hook_len_max_um,
        ) = _pair_distance_stats_um(pos_after, self.spring_pairs, self.body_flag_rows)
        hook_len_mean_over_b = (
            hook_len_mean_um / max(self.cfg.scale.b_um, 1e-12)
            if hook_count > 0
            else float("nan")
        )
        hook_len_min_over_b = (
            hook_len_min_um / max(self.cfg.scale.b_um, 1e-12)
            if hook_count > 0
            else float("nan")
        )
        hook_len_max_over_b = (
            hook_len_max_um / max(self.cfg.scale.b_um, 1e-12)
            if hook_count > 0
            else float("nan")
        )
        hook_len_rel_err_mean, hook_len_rel_err_max = _pair_rel_error_stats(
            pos_after,
            self.spring_pairs,
            self.spring_rests_m,
            self.body_flag_rows,
        )

        hook_angle_mean_deg = float("nan")
        hook_angle_min_deg = float("nan")
        hook_angle_max_deg = float("nan")
        hook_angle_err_mean_deg = float("nan")
        hook_angle_err_max_deg = float("nan")
        if self.model.hook_triplets.size > 0:
            hook_angles_rad = _triplet_angles_rad(pos_after, self.model.hook_triplets)
            hook_angles_deg = hook_angles_rad * RAD_TO_DEG
            hook_angle_mean_deg = float(np.mean(hook_angles_deg))
            hook_angle_min_deg = float(np.min(hook_angles_deg))
            hook_angle_max_deg = float(np.max(hook_angles_deg))
            hook_angle_err_deg = np.maximum(
                float(self.cfg.hook.threshold_deg) - hook_angles_deg,
                0.0,
            )
            hook_angle_err_mean_deg = float(np.mean(hook_angle_err_deg))
            hook_angle_err_max_deg = float(np.max(hook_angle_err_deg))

        flag_bond_len_mean_over_b = (
            flag_intra_len_mean_um / max(self.cfg.scale.b_um, 1e-12)
            if flag_intra_count > 0
            else float("nan")
        )
        flag_bond_len_min_over_b = (
            flag_intra_len_min_um / max(self.cfg.scale.b_um, 1e-12)
            if flag_intra_count > 0
            else float("nan")
        )
        flag_bond_len_max_over_b = (
            flag_intra_len_max_um / max(self.cfg.scale.b_um, 1e-12)
            if flag_intra_count > 0
            else float("nan")
        )
        flag_bond_rel_err_mean, flag_bond_rel_err_max = _pair_rel_error_stats(
            pos_after,
            self.spring_pairs,
            self.spring_rests_m,
            self.flag_intra_rows,
        )
        flag_bend_err_mean_deg = float("nan")
        flag_bend_err_max_deg = float("nan")
        if self.flag_bending_rows.size > 0:
            rows = self.flag_bending_rows
            triplets = self.model.bending_triplets[rows]
            actual = _triplet_angles_rad(pos_after, triplets)
            target = np.zeros_like(actual)
            for idx, row in enumerate(rows):
                flag_id = int(self.model.bending_flag_ids[row])
                state = int(self.model.flag_states[flag_id])
                key = (
                    "normal"
                    if state == int(PolymorphState.NORMAL)
                    else "semicoiled"
                    if state == int(PolymorphState.SEMICOILED)
                    else "curly1"
                )
                target[idx] = np.deg2rad(float(self.bend_map[key]))
            err = np.abs(actual - target) * RAD_TO_DEG
            flag_bend_err_mean_deg = float(np.mean(err))
            flag_bend_err_max_deg = float(np.max(err))

        flag_torsion_err_mean_deg = float("nan")
        flag_torsion_err_max_deg = float("nan")
        if self.flag_torsion_rows.size > 0:
            rows = self.flag_torsion_rows
            quads = self.model.torsion_quads[rows]
            actual = _torsion_angles_rad(pos_after, quads)
            target = np.zeros_like(actual)
            for idx, row in enumerate(rows):
                flag_id = int(self.model.torsion_flag_ids[row])
                state = int(self.model.flag_states[flag_id])
                key = (
                    "normal"
                    if state == int(PolymorphState.NORMAL)
                    else "semicoiled"
                    if state == int(PolymorphState.SEMICOILED)
                    else "curly1"
                )
                target[idx] = np.deg2rad(float(self.torsion_map[key]))
            err = np.abs(_wrap_angle(actual - target)) * RAD_TO_DEG
            flag_torsion_err_mean_deg = float(np.mean(err))
            flag_torsion_err_max_deg = float(np.max(err))

        flag_state_changed = bool(
            self.model.flag_states.shape != self.prev_flag_states.shape
            or np.any(self.model.flag_states != self.prev_flag_states)
        )

        row: dict[str, float | int | bool] = {
            "step": int(step),
            "t_star": float(t_star),
            "dt_star": float(diag.dt_star),
            "t_s": float(t_star * self.cfg.tau_s),
            "dt_s": float(diag.dt_s),
            "tau_s": float(self.cfg.tau_s),
            "dt_internal_s": float(self.cfg.dt_s),
            "pos_all_finite": pos_all_finite,
            "any_nan": any_nan,
            "any_inf": any_inf,
            "mean_disp_um": float(np.mean(disp_um)),
            "max_disp_um": float(np.max(disp_um)),
            "bond_count_body_body": bond_count_body_body,
            "bond_len_mean_body_body_um": bond_len_mean_body_body_um,
            "flag_intra_count": flag_intra_count,
            "bond_len_mean_flag_intra_um": flag_intra_len_mean_um,
            "flag_intra_len_mean_um": flag_intra_len_mean_um,
            "flag_intra_len_min_um": flag_intra_len_min_um,
            "flag_intra_len_max_um": flag_intra_len_max_um,
            "hook_count": hook_count,
            "bond_len_mean_body_flag_um": hook_len_mean_um,
            "hook_len_mean_um": hook_len_mean_um,
            "hook_len_min_um": hook_len_min_um,
            "hook_len_max_um": hook_len_max_um,
            "hook_len_mean_over_b": hook_len_mean_over_b,
            "hook_len_min_over_b": hook_len_min_over_b,
            "hook_len_max_over_b": hook_len_max_over_b,
            "hook_len_rel_err_mean": hook_len_rel_err_mean,
            "hook_len_rel_err_max": hook_len_rel_err_max,
            "hook_angle_mean_deg": hook_angle_mean_deg,
            "hook_angle_min_deg": hook_angle_min_deg,
            "hook_angle_max_deg": hook_angle_max_deg,
            "hook_angle_err_mean_deg": hook_angle_err_mean_deg,
            "hook_angle_err_max_deg": hook_angle_err_max_deg,
            "flag_bend_err_mean_deg": flag_bend_err_mean_deg,
            "flag_bend_err_max_deg": flag_bend_err_max_deg,
            "flag_torsion_err_mean_deg": flag_torsion_err_mean_deg,
            "flag_torsion_err_max_deg": flag_torsion_err_max_deg,
            "flag_bond_len_mean_over_b": flag_bond_len_mean_over_b,
            "flag_bond_len_min_over_b": flag_bond_len_min_over_b,
            "flag_bond_len_max_over_b": flag_bond_len_max_over_b,
            "flag_bond_rel_err_mean": flag_bond_rel_err_mean,
            "flag_bond_rel_err_max": flag_bond_rel_err_max,
            "F_total_mean_body": _mean_norm(diag.total_forces, self.body_mask),
            "F_total_mean_flag": _mean_norm(diag.total_forces, self.flag_mask),
            "F_total_mean_all": float(
                np.mean(np.linalg.norm(diag.total_forces, axis=1))
            ),
            "F_motor_mean_body": _mean_norm(diag.motor_forces, self.body_mask),
            "F_motor_mean_flag": _mean_norm(diag.motor_forces, self.flag_mask),
            "F_spring_mean_body": _mean_norm(diag.spring_forces, self.body_mask),
            "F_spring_mean_flag": _mean_norm(diag.spring_forces, self.flag_mask),
            "F_bend_mean_body": _mean_norm(diag.bend_forces, self.body_mask),
            "F_bend_mean_flag": _mean_norm(diag.bend_forces, self.flag_mask),
            "F_torsion_mean_body": _mean_norm(diag.torsion_forces, self.body_mask),
            "F_torsion_mean_flag": _mean_norm(diag.torsion_forces, self.flag_mask),
            "F_repulsion_mean_body": _mean_norm(diag.repulsion_forces, self.body_mask),
            "F_repulsion_mean_flag": _mean_norm(diag.repulsion_forces, self.flag_mask),
            "F_hook_mean_body": _mean_norm(diag.hook_forces, self.body_mask),
            "F_hook_mean_flag": _mean_norm(diag.hook_forces, self.flag_mask),
            "torque_for_forces_Nm": float(self.cfg.torque_for_forces_Nm),
            "motor_torque_Nm": float(self.cfg.motor_torque_Nm),
            "torsion_fd_eps_m": float(diag.torsion_fd_eps_m),
            "torsion_fd_eps_over_b": float(
                diag.torsion_fd_eps_m / max(self.cfg.b_m, 1e-30)
            ),
            "motor_degenerate_axis_count": int(diag.motor_degenerate_axis_count),
            "motor_split_rank_deficient_count": int(
                diag.motor_split_rank_deficient_count
            ),
            "motor_bond_length_clipped_count": int(
                diag.motor_bond_length_clipped_count
            ),
            "flag_state_changed": flag_state_changed,
            "projection_body_rigid_enabled": bool(
                self.cfg.projection.enable_body_rigid_projection
            ),
            "projection_hook_length_enabled": bool(
                self.cfg.projection.enable_hook_length_projection
            ),
            "projection_basal_link_direction_enabled": bool(
                self.cfg.projection.enable_basal_link_direction_projection
            ),
            "projection_flagella_chain_length_enabled": bool(
                self.cfg.projection.enable_flagella_chain_length_projection
            ),
            "projection_flagella_template_enabled": bool(
                self.cfg.projection.enable_flagella_template_projection
            ),
            "brownian_enabled": bool(diag.brownian_enabled),
            "brownian_disp_mean_um": brownian_disp_mean_um,
        }
        if self._writer is None or self._csv_fp is None:
            raise RuntimeError("step summary writer is not initialized")
        self._writer.writerow(row)
        self._csv_fp.flush()
        self.prev_flag_states = self.model.flag_states.copy()

    def write_csv(self) -> Path:
        if self._csv_fp is not None:
            self._csv_fp.close()
            self._csv_fp = None
            self._writer = None
        return self.step_summary_path
