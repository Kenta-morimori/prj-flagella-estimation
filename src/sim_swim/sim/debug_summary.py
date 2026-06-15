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
    "finite_pass",
    "shape_pass_nonbody",
    "first_fail_category_nonbody",
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
    "local_attach_first_vs_body_axis_angle_deg",
    "local_attach_first_vs_body_axis_err_deg",
    "flag_bend_err_mean_deg",
    "flag_bend_err_max_deg",
    "flag_torsion_err_mean_deg",
    "flag_torsion_err_max_deg",
    "flag_bond_len_mean_over_b",
    "flag_bond_len_min_over_b",
    "flag_bond_len_max_over_b",
    "flag_bond_rel_err_mean",
    "flag_bond_rel_err_max",
    "body_displacement_um",
    "body_speed_um_s",
    "body_axis_step_angle_deg",
    "body_axis_cumulative_angle_deg",
    "body_axis_wobble_rms_deg",
    "body_angular_velocity_rms_rad_s",
    "bundle_axis_vs_body_axis_angle_deg",
    "bundle_axis_vs_rear_angle_deg",
    "bundle_rearward_projection",
    "bundle_tip_axis_dist_mean_um",
    "bundle_tip_axis_dist_max_um",
    "bundle_participation_ratio",
    "bundle_independent_flagella_count",
    "flag_tip_pair_dist_mean_um",
    "flag_tip_pair_dist_min_um",
    "flag_tip_pair_dist_max_um",
    "flag_flag_bead_pair_dist_min_um",
    "flag_flag_bead_pair_dist_mean_um",
    "flag_flag_close_pair_count",
    "flag_flag_repulsion_force_mean_N",
    "flag_flag_repulsion_force_max_N",
    "flag_flag_basal_repulsion_force_mean_N",
    "flag_flag_basal_repulsion_force_max_N",
    "F_total_mean_body",
    "F_total_mean_flag",
    "F_total_mean_all",
    "F_motor_mean_body",
    "F_motor_mean_flag",
    "motor_ra_len_um",
    "motor_rb_len_um",
    "motor_Ta_norm",
    "motor_Tb_norm",
    "motor_Fa_norm",
    "motor_Fb_norm",
    "motor_axis_vs_rear_direction_angle_deg",
    "local_twist_root_orientation_deg",
    "local_twist_tip_orientation_deg",
    "local_twist_abs_mean_deg",
    "local_twist_abs_max_deg",
    "local_twist_tip_activity_ratio",
    "flag_root_azimuth_deg",
    "flag_phase_deg",
    "flag_phase_rate_hz",
    "flag_body_phase_diff_deg",
    "flag_helix_spin_phase_deg",
    "flag_helix_spin_rate_hz",
    "flag_helix_spin_fit_r2",
    "motor_attach_force_norm",
    "motor_first_force_norm",
    "motor_second_force_norm",
    "motor_split_residual_norm",
    "motor_ta_dot_ra_abs",
    "motor_tb_dot_rb_abs",
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
    "F_body_equiv_load_mean",
    "F_body_equiv_load_max",
    "body_equiv_load_mode",
    "body_equiv_load_target_torque_Nm",
    "body_equiv_load_target_force_N",
    "body_equiv_attach_region_id",
    "local_attach_first_len_over_b",
    "local_attach_first_rel_err",
    "local_first_second_len_over_b",
    "local_first_second_rel_err",
    "local_second_third_len_over_b",
    "local_second_third_rel_err",
    "local_basal_bend_angle_deg",
    "local_basal_bend_err_deg",
    "local_first_torsion_angle_deg",
    "local_first_torsion_err_deg",
    "local_F_spring_attach_first",
    "local_F_spring_first_second",
    "local_F_spring_second_third",
    "local_F_bend_basal",
    "local_F_torsion_first",
    "local_F_motor_attach",
    "local_F_motor_first",
    "local_F_motor_second",
    "local_F_repulsion_basal_region",
    "flag_state_changed",
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

NONBODY_HOOK_REL_ERR_MAX_LIMIT = 1.0
NONBODY_HOOK_ANGLE_ERR_MAX_DEG_LIMIT = 30.0
NONBODY_FLAG_BOND_REL_ERR_MAX_LIMIT = 1.0
NONBODY_FLAG_BEND_ERR_MAX_DEG_LIMIT = 60.0
NONBODY_FLAG_TORSION_ERR_MAX_DEG_LIMIT = 120.0


def _wrap_angle(rad: np.ndarray | float) -> np.ndarray | float:
    return (np.asarray(rad) + np.pi) % (2.0 * np.pi) - np.pi


def _is_finite_number(value: float) -> bool:
    return bool(np.isfinite(float(value)))


def _body_frame_basis(
    positions_m: np.ndarray, model: SimModel
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if len(model.body_layer_indices) >= 2:
        first = model.body_layer_indices[0].astype(int, copy=False)
        last = model.body_layer_indices[-1].astype(int, copy=False)
        c_first = np.mean(positions_m[first], axis=0)
        c_last = np.mean(positions_m[last], axis=0)
        origin = c_first
        axis = c_last - c_first
        marker = positions_m[int(first[0])] - c_first if first.size > 0 else np.zeros(3)
    else:
        body_idx = model.body_indices.astype(int, copy=False)
        origin = np.mean(positions_m[body_idx], axis=0)
        axis = positions_m[int(body_idx[-1])] - positions_m[int(body_idx[0])]
        marker = positions_m[int(body_idx[0])] - origin

    axis_norm = float(np.linalg.norm(axis))
    if axis_norm <= 1e-18:
        axis = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        axis = axis / axis_norm

    marker = marker - float(np.dot(marker, axis)) * axis
    marker_norm = float(np.linalg.norm(marker))
    if marker_norm <= 1e-18:
        marker = np.array([0.0, 1.0, 0.0], dtype=float)
        marker = marker - float(np.dot(marker, axis)) * axis
        marker_norm = float(np.linalg.norm(marker))
    if marker_norm <= 1e-18:
        marker = np.cross(axis, np.array([0.0, 0.0, 1.0], dtype=float))
        marker_norm = float(np.linalg.norm(marker))
    e1 = marker / max(marker_norm, 1e-18)
    e2 = np.cross(axis, e1)
    e2 = e2 / max(float(np.linalg.norm(e2)), 1e-18)
    return origin, axis, e1, e2


def _to_body_frame(positions_m: np.ndarray, model: SimModel) -> np.ndarray:
    origin, axis, e1, e2 = _body_frame_basis(positions_m, model)
    rel = positions_m - origin
    return np.column_stack((rel @ axis, rel @ e1, rel @ e2))


def _estimate_helix_spin_offset_deg(points_body: np.ndarray) -> tuple[float, float]:
    if points_body.shape[0] < 5:
        return float("nan"), float("nan")

    centered = points_body - np.mean(points_body, axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    axis = vh[0]
    axis = axis / max(float(np.linalg.norm(axis)), 1e-18)
    if float(np.dot(axis, points_body[-1] - points_body[0])) < 0.0:
        axis = -axis

    ref = np.array([0.0, 1.0, 0.0], dtype=float)
    if abs(float(np.dot(ref, axis))) > 0.9:
        ref = np.array([0.0, 0.0, 1.0], dtype=float)
    e1 = ref - float(np.dot(ref, axis)) * axis
    e1 = e1 / max(float(np.linalg.norm(e1)), 1e-18)
    e2 = np.cross(axis, e1)
    e2 = e2 / max(float(np.linalg.norm(e2)), 1e-18)

    rel = points_body - points_body[0]
    s = rel @ axis
    u = rel @ e1
    v = rel @ e2
    mat = np.column_stack([u, v, np.ones_like(u)])
    rhs = -(u * u + v * v)
    coef, *_ = np.linalg.lstsq(mat, rhs, rcond=None)
    cx = -0.5 * float(coef[0])
    cy = -0.5 * float(coef[1])
    theta = np.unwrap(np.arctan2(v - cy, u - cx))
    if theta.size < 3 or not np.isfinite(theta).all():
        return float("nan"), float("nan")

    slope, intercept = np.polyfit(s, theta, 1)
    pred = slope * s + intercept
    ss_res = float(np.sum((theta - pred) ** 2))
    ss_tot = float(np.sum((theta - float(np.mean(theta))) ** 2))
    fit_r2 = 1.0 - ss_res / max(ss_tot, 1e-30)
    return float(np.rad2deg(intercept)), float(fit_r2)


def _check_nonbody_shape_pass(
    *,
    finite_pass: bool,
    has_hook_pair: bool,
    has_hook_angle: bool,
    has_flag_bond: bool,
    has_flag_bend: bool,
    has_flag_torsion: bool,
    local_attach_first_rel_err: float,
    hook_len_rel_err_max: float,
    hook_angle_err_max_deg: float,
    flag_bond_rel_err_max: float,
    flag_bend_err_max_deg: float,
    flag_torsion_err_max_deg: float,
) -> tuple[bool, str]:
    if not finite_pass:
        return False, "finite"

    hook_metrics: list[float] = []
    if has_hook_pair:
        hook_metrics.extend([local_attach_first_rel_err, hook_len_rel_err_max])
    if has_hook_angle:
        hook_metrics.append(hook_angle_err_max_deg)
    if any(not _is_finite_number(v) for v in hook_metrics):
        return False, "hook_nonfinite"
    if has_hook_pair and local_attach_first_rel_err > NONBODY_HOOK_REL_ERR_MAX_LIMIT:
        return False, "hook"
    if has_hook_pair and hook_len_rel_err_max > NONBODY_HOOK_REL_ERR_MAX_LIMIT:
        return False, "hook"
    if has_hook_angle and hook_angle_err_max_deg > NONBODY_HOOK_ANGLE_ERR_MAX_DEG_LIMIT:
        return False, "hook"

    flag_metrics: list[float] = []
    if has_flag_bond:
        flag_metrics.append(flag_bond_rel_err_max)
    if has_flag_bend:
        flag_metrics.append(flag_bend_err_max_deg)
    if has_flag_torsion:
        flag_metrics.append(flag_torsion_err_max_deg)
    if any(not _is_finite_number(v) for v in flag_metrics):
        return False, "flag_nonfinite"
    if has_flag_bond and flag_bond_rel_err_max > NONBODY_FLAG_BOND_REL_ERR_MAX_LIMIT:
        return False, "flag"
    if has_flag_bend and flag_bend_err_max_deg > NONBODY_FLAG_BEND_ERR_MAX_DEG_LIMIT:
        return False, "flag"
    if (
        has_flag_torsion
        and flag_torsion_err_max_deg > NONBODY_FLAG_TORSION_ERR_MAX_DEG_LIMIT
    ):
        return False, "flag"

    return True, "none"


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


def _mean_norm_indices(forces: np.ndarray, bead_indices: list[int]) -> float:
    if forces.size == 0:
        return float("nan")
    idx = [int(i) for i in bead_indices if int(i) >= 0]
    if not idx:
        return float("nan")
    selected = forces[np.asarray(idx, dtype=int)]
    return float(np.mean(np.linalg.norm(selected, axis=1)))


def _pair_row_lookup(spring_pairs: np.ndarray) -> dict[tuple[int, int], int]:
    out: dict[tuple[int, int], int] = {}
    for row, (i_raw, j_raw) in enumerate(spring_pairs):
        i = int(i_raw)
        j = int(j_raw)
        key = (i, j) if i < j else (j, i)
        out[key] = int(row)
    return out


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


def _unit_vector(vec: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vec))
    if norm <= 1e-18:
        return np.zeros(3, dtype=float)
    return vec / norm


def _angle_between_deg(a: np.ndarray, b: np.ndarray) -> float:
    a_u = _unit_vector(a)
    b_u = _unit_vector(b)
    if float(np.linalg.norm(a_u)) <= 1e-18 or float(np.linalg.norm(b_u)) <= 1e-18:
        return float("nan")
    return float(np.rad2deg(np.arccos(float(np.clip(np.dot(a_u, b_u), -1.0, 1.0)))))


def _body_center_axis_rear(
    positions_m: np.ndarray, model: SimModel
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    body_idx = model.body_indices.astype(int, copy=False)
    center = np.mean(positions_m[body_idx], axis=0)
    if len(model.body_layer_indices) >= 2:
        first = model.body_layer_indices[0].astype(int, copy=False)
        last = model.body_layer_indices[-1].astype(int, copy=False)
        c_first = np.mean(positions_m[first], axis=0)
        c_last = np.mean(positions_m[last], axis=0)
        axis = _unit_vector(c_last - c_first)
    else:
        axis = np.array([1.0, 0.0, 0.0], dtype=float)
    if float(np.linalg.norm(axis)) <= 1e-18:
        axis = np.array([1.0, 0.0, 0.0], dtype=float)
    rear = -axis
    return center, axis, rear


def _flag_tip_pair_stats_um(tips_m: np.ndarray) -> tuple[float, float, float]:
    if tips_m.shape[0] < 2:
        return float("nan"), float("nan"), float("nan")
    dists: list[float] = []
    for i in range(tips_m.shape[0]):
        for j in range(i + 1, tips_m.shape[0]):
            dists.append(float(np.linalg.norm(tips_m[i] - tips_m[j]) * M_TO_UM))
    arr = np.asarray(dists, dtype=float)
    return float(np.mean(arr)), float(np.min(arr)), float(np.max(arr))


def _attach_first_body_axis_metrics(
    positions_m: np.ndarray,
    hook_triplets: np.ndarray,
    body_axis: np.ndarray,
) -> dict[str, float]:
    if hook_triplets.size == 0:
        return {
            "local_attach_first_vs_body_axis_angle_deg": float("nan"),
            "local_attach_first_vs_body_axis_err_deg": float("nan"),
        }

    angles: list[float] = []
    for attach_raw, first_raw, _second_raw in hook_triplets.astype(int, copy=False):
        attach = int(attach_raw)
        first = int(first_raw)
        attach_first = positions_m[first] - positions_m[attach]
        angle = _angle_between_deg(attach_first, body_axis)
        if np.isfinite(angle):
            angles.append(float(angle))

    if not angles:
        return {
            "local_attach_first_vs_body_axis_angle_deg": float("nan"),
            "local_attach_first_vs_body_axis_err_deg": float("nan"),
        }
    arr = np.asarray(angles, dtype=float)
    return {
        "local_attach_first_vs_body_axis_angle_deg": float(np.mean(arr)),
        "local_attach_first_vs_body_axis_err_deg": float(np.max(np.abs(arr - 90.0))),
    }


def _bundle_metrics(
    positions_m: np.ndarray,
    model: SimModel,
    cfg: SimulationConfig,
    body_axis: np.ndarray,
    rear_dir: np.ndarray,
) -> dict[str, float | int]:
    n_flagella = len(model.flagella_indices)
    if n_flagella <= 0:
        return {
            "bundle_axis_vs_body_axis_angle_deg": float("nan"),
            "bundle_axis_vs_rear_angle_deg": float("nan"),
            "bundle_rearward_projection": float("nan"),
            "bundle_tip_axis_dist_mean_um": float("nan"),
            "bundle_tip_axis_dist_max_um": float("nan"),
            "bundle_participation_ratio": float("nan"),
            "bundle_independent_flagella_count": 0,
            "flag_tip_pair_dist_mean_um": float("nan"),
            "flag_tip_pair_dist_min_um": float("nan"),
            "flag_tip_pair_dist_max_um": float("nan"),
        }

    roots = np.asarray(
        [positions_m[int(idx[0])] for idx in model.flagella_indices],
        dtype=float,
    )
    tips = np.asarray(
        [positions_m[int(idx[-1])] for idx in model.flagella_indices],
        dtype=float,
    )
    root_center = np.mean(roots, axis=0)
    tip_center = np.mean(tips, axis=0)
    bundle_axis = _unit_vector(tip_center - root_center)
    if float(np.linalg.norm(bundle_axis)) <= 1e-18:
        bundle_axis_body_angle = float("nan")
        bundle_axis_rear_angle = float("nan")
        rearward_projection = float("nan")
        tip_axis_dist_um = np.full((n_flagella,), float("nan"), dtype=float)
    else:
        angle_forward = _angle_between_deg(bundle_axis, body_axis)
        angle_rear = _angle_between_deg(bundle_axis, -body_axis)
        bundle_axis_body_angle = float(np.nanmin([angle_forward, angle_rear]))
        bundle_axis_rear_angle = _angle_between_deg(bundle_axis, rear_dir)
        rearward_projection = float(np.dot(bundle_axis, rear_dir))
        p1 = root_center + bundle_axis
        tip_axis_dist_um = _point_to_line_distance(tips, root_center, p1) * M_TO_UM

    threshold_um = max(0.75 * float(cfg.scale.b_um), 1e-12)
    if np.isfinite(tip_axis_dist_um).all():
        participants = tip_axis_dist_um <= threshold_um
        participation_ratio = float(np.mean(participants))
        independent_count = int(np.sum(~participants))
        tip_axis_dist_mean_um = float(np.mean(tip_axis_dist_um))
        tip_axis_dist_max_um = float(np.max(tip_axis_dist_um))
    else:
        participation_ratio = float("nan")
        independent_count = n_flagella
        tip_axis_dist_mean_um = float("nan")
        tip_axis_dist_max_um = float("nan")

    pair_mean, pair_min, pair_max = _flag_tip_pair_stats_um(tips)
    return {
        "bundle_axis_vs_body_axis_angle_deg": bundle_axis_body_angle,
        "bundle_axis_vs_rear_angle_deg": bundle_axis_rear_angle,
        "bundle_rearward_projection": rearward_projection,
        "bundle_tip_axis_dist_mean_um": tip_axis_dist_mean_um,
        "bundle_tip_axis_dist_max_um": tip_axis_dist_max_um,
        "bundle_participation_ratio": participation_ratio,
        "bundle_independent_flagella_count": independent_count,
        "flag_tip_pair_dist_mean_um": pair_mean,
        "flag_tip_pair_dist_min_um": pair_min,
        "flag_tip_pair_dist_max_um": pair_max,
    }


def _flag_flag_repulsion_metrics(
    positions_m: np.ndarray,
    flag_bead_indices: np.ndarray,
    flag_bead_ids: np.ndarray,
    basal_flag_bead_indices: np.ndarray,
    cfg: SimulationConfig,
) -> dict[str, float | int]:
    if flag_bead_indices.size < 2:
        return {
            "flag_flag_bead_pair_dist_min_um": float("nan"),
            "flag_flag_bead_pair_dist_mean_um": float("nan"),
            "flag_flag_close_pair_count": 0,
            "flag_flag_repulsion_force_mean_N": 0.0,
            "flag_flag_repulsion_force_max_N": 0.0,
            "flag_flag_basal_repulsion_force_mean_N": 0.0,
            "flag_flag_basal_repulsion_force_max_N": 0.0,
        }

    cutoff_m = (
        float(cfg.potentials.spring_spring_repulsion.cutoff_over_b)
        * max(float(cfg.scale.b_um), 1e-12)
        * 1.0e-6
    )
    a_length_m = (
        float(cfg.potentials.spring_spring_repulsion.a_ss_over_b)
        * max(float(cfg.scale.b_um), 1e-12)
        * 1.0e-6
    )
    a_ss = float(cfg.potentials.spring_spring_repulsion.A_ss_over_T) * float(
        cfg.torque_for_forces_Nm
    )
    cutoff_eff = max(cutoff_m, 0.0)
    a_eff = max(a_length_m, 1e-12)

    beads = flag_bead_indices.astype(int, copy=False)
    flag_ids = flag_bead_ids.astype(int, copy=False)
    pts = positions_m[beads]
    diffs = pts[:, None, :] - pts[None, :, :]
    distances_full = np.linalg.norm(diffs, axis=2)
    different_flag = flag_ids[:, None] != flag_ids[None, :]
    upper = np.triu(np.ones(distances_full.shape, dtype=bool), k=1)
    mask = upper & different_flag
    distances = distances_full[mask]
    if distances.size == 0:
        return {
            "flag_flag_bead_pair_dist_min_um": float("nan"),
            "flag_flag_bead_pair_dist_mean_um": float("nan"),
            "flag_flag_close_pair_count": 0,
            "flag_flag_repulsion_force_mean_N": 0.0,
            "flag_flag_repulsion_force_max_N": 0.0,
            "flag_flag_basal_repulsion_force_mean_N": 0.0,
            "flag_flag_basal_repulsion_force_max_N": 0.0,
        }

    active_distances = distances[distances < cutoff_eff]
    active = (
        (a_ss / a_eff) * np.exp(-active_distances / a_eff)
        if active_distances.size > 0
        else np.zeros((0,), dtype=float)
    )
    basal = basal_flag_bead_indices.astype(int, copy=False)
    if basal.size > 0:
        basal_local = np.isin(beads, basal)
        basal_pair_mask = mask & (basal_local[:, None] | basal_local[None, :])
        basal_distances = distances_full[basal_pair_mask]
        active_basal_distances = basal_distances[basal_distances < cutoff_eff]
        if active_basal_distances.size > 0:
            basal_active = (a_ss / a_eff) * np.exp(-active_basal_distances / a_eff)
            basal_mean = float(np.mean(np.abs(basal_active)))
            basal_max = float(np.max(np.abs(basal_active)))
        else:
            basal_mean = 0.0
            basal_max = 0.0
    else:
        basal_mean = float("nan")
        basal_max = float("nan")

    return {
        "flag_flag_bead_pair_dist_min_um": float(np.min(distances) * M_TO_UM),
        "flag_flag_bead_pair_dist_mean_um": float(np.mean(distances) * M_TO_UM),
        "flag_flag_close_pair_count": int(active.size),
        "flag_flag_repulsion_force_mean_N": (
            float(np.mean(np.abs(active))) if active.size > 0 else 0.0
        ),
        "flag_flag_repulsion_force_max_N": (
            float(np.max(np.abs(active))) if active.size > 0 else 0.0
        ),
        "flag_flag_basal_repulsion_force_mean_N": basal_mean,
        "flag_flag_basal_repulsion_force_max_N": basal_max,
    }


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

        self.local_attach_idx = -1
        self.local_first_idx = -1
        self.local_second_idx = -1
        self.local_third_idx = -1
        self.local_flag_id = -1
        self.local_attach_first_row = -1
        self.local_first_second_row = -1
        self.local_second_third_row = -1
        self.local_basal_triplet = np.full((3,), -1, dtype=int)
        self.local_first_torsion_quad = np.full((4,), -1, dtype=int)
        self.prev_flag_root_azimuth_deg: float | None = None
        self.prev_flag_phase_deg: float | None = None
        self.prev_flag_phase_t_s: float | None = None
        self.prev_flag_helix_spin_offset_deg: float | None = None
        self.prev_flag_helix_spin_phase_deg: float | None = None
        self.prev_flag_helix_spin_t_s: float | None = None
        initial_center, initial_axis, _ = _body_center_axis_rear(
            self.model.positions_m,
            self.model,
        )
        self.initial_body_center_m = initial_center
        self.initial_body_axis = initial_axis
        self.prev_body_center_m: np.ndarray | None = None
        self.prev_body_axis: np.ndarray | None = None
        self.prev_body_t_s: float | None = None
        self.body_axis_cumulative_angle_deg = 0.0
        self.body_axis_wobble_square_sum = 0.0
        self.body_axis_wobble_count = 0
        self.body_omega_square_sum = 0.0
        self.body_omega_count = 0
        self.last_row: dict[str, float | int | bool] | None = None

        pair_rows = _pair_row_lookup(self.spring_pairs)

        def _find_pair_row(a: int, b: int) -> int:
            key = (a, b) if a < b else (b, a)
            return int(pair_rows.get(key, -1))

        if self.model.motor_triplets.size > 0:
            ib_raw, jf_raw, kf_raw = self.model.motor_triplets[0]
            self.local_attach_idx = int(ib_raw)
            self.local_first_idx = int(jf_raw)
            self.local_second_idx = int(kf_raw)
            self.local_flag_id = int(self.model.bead_flag_ids[self.local_first_idx])

            if self.local_flag_id >= 0:
                flag_chain = self.model.flagella_indices[self.local_flag_id].astype(
                    int, copy=False
                )
                pos = np.where(flag_chain == self.local_second_idx)[0]
                if pos.size > 0 and int(pos[0]) + 1 < int(flag_chain.size):
                    self.local_third_idx = int(flag_chain[int(pos[0]) + 1])

            self.local_attach_first_row = _find_pair_row(
                self.local_attach_idx, self.local_first_idx
            )
            self.local_first_second_row = _find_pair_row(
                self.local_first_idx, self.local_second_idx
            )
            if self.local_third_idx >= 0:
                self.local_second_third_row = _find_pair_row(
                    self.local_second_idx, self.local_third_idx
                )

            for triplet in self.model.hook_triplets.astype(int, copy=False):
                a, b, c = int(triplet[0]), int(triplet[1]), int(triplet[2])
                if (
                    a == self.local_attach_idx
                    and b == self.local_first_idx
                    and c == self.local_second_idx
                ):
                    self.local_basal_triplet = np.array([a, b, c], dtype=int)
                    break

            if self.local_flag_id >= 0 and self.flag_torsion_rows.size > 0:
                rows = self.flag_torsion_rows[
                    self.model.torsion_flag_ids[self.flag_torsion_rows]
                    == self.local_flag_id
                ]
                preferred_row = -1
                for row in rows:
                    quad = self.model.torsion_quads[int(row)].astype(int, copy=False)
                    if (
                        int(quad[0]) == self.local_first_idx
                        and int(quad[1]) == self.local_second_idx
                    ):
                        preferred_row = int(row)
                        break
                if preferred_row < 0 and rows.size > 0:
                    preferred_row = int(rows[0])
                if preferred_row >= 0:
                    self.local_first_torsion_quad = self.model.torsion_quads[
                        preferred_row
                    ].astype(int, copy=False)

        if self.spring_pairs.size == 0:
            self.body_body_rows = np.zeros((0,), dtype=int)
            self.flag_intra_rows = np.zeros((0,), dtype=int)
            self.body_flag_rows = np.zeros((0,), dtype=int)
            self.basal_flag_bead_indices = np.zeros((0,), dtype=int)
            self.flag_bead_indices = np.zeros((0,), dtype=int)
            self.flag_bead_ids = np.zeros((0,), dtype=int)
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

        basal_beads: list[int] = []
        flag_beads: list[int] = []
        flag_ids: list[int] = []
        for flag_id, flag_indices in enumerate(self.model.flagella_indices):
            chain = flag_indices.astype(int, copy=False)
            basal_beads.extend(chain[: min(2, chain.size)].tolist())
            flag_beads.extend(chain.tolist())
            flag_ids.extend([int(flag_id)] * int(chain.size))
        self.basal_flag_bead_indices = np.asarray(basal_beads, dtype=int)
        self.flag_bead_indices = np.asarray(flag_beads, dtype=int)
        self.flag_bead_ids = np.asarray(flag_ids, dtype=int)

    def _phase_reference_frame(
        self, positions_m: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        axis = np.zeros(3, dtype=float)
        if len(self.model.body_layer_indices) >= 2:
            first = self.model.body_layer_indices[0].astype(int, copy=False)
            last = self.model.body_layer_indices[-1].astype(int, copy=False)
            c_first = np.mean(positions_m[first], axis=0)
            c_last = np.mean(positions_m[last], axis=0)
            axis = c_last - c_first
        elif self.model.body_indices.size >= 2:
            i0 = int(self.model.body_indices[0])
            i1 = int(self.model.body_indices[-1])
            axis = positions_m[i1] - positions_m[i0]

        n_axis = float(np.linalg.norm(axis))
        if n_axis <= 1e-18:
            axis = np.array([1.0, 0.0, 0.0], dtype=float)
            n_axis = 1.0
        axis = axis / n_axis

        ref_candidates = (
            np.array([1.0, 0.0, 0.0], dtype=float),
            np.array([0.0, 1.0, 0.0], dtype=float),
            np.array([0.0, 0.0, 1.0], dtype=float),
        )
        e1 = np.zeros(3, dtype=float)
        for ref in ref_candidates:
            proj = ref - float(np.dot(ref, axis)) * axis
            n_proj = float(np.linalg.norm(proj))
            if n_proj > 1e-18:
                e1 = proj / n_proj
                break
        if float(np.linalg.norm(e1)) <= 1e-18:
            fallback = np.cross(axis, np.array([0.0, 1.0, 0.0], dtype=float))
            if float(np.linalg.norm(fallback)) <= 1e-18:
                fallback = np.cross(axis, np.array([0.0, 0.0, 1.0], dtype=float))
            e1 = fallback / max(float(np.linalg.norm(fallback)), 1e-18)
        e2 = np.cross(axis, e1)
        e2 = e2 / max(float(np.linalg.norm(e2)), 1e-18)
        return axis, e1, e2

    @staticmethod
    def _azimuth_deg(vec: np.ndarray, e1: np.ndarray, e2: np.ndarray) -> float:
        proj_x = float(np.dot(vec, e1))
        proj_y = float(np.dot(vec, e2))
        if abs(proj_x) <= 1e-18 and abs(proj_y) <= 1e-18:
            return float("nan")
        return float(np.rad2deg(np.arctan2(proj_y, proj_x)))

    @staticmethod
    def _wrap_deg(deg: float) -> float:
        return float((deg + 180.0) % 360.0 - 180.0)

    def record(self, step: int, t_star: float, diag: StepDiagnostics) -> None:
        pos_after = diag.positions_after_m
        disp_um = np.linalg.norm(pos_after - diag.positions_before_m, axis=1) * M_TO_UM
        any_nan = bool(np.isnan(pos_after).any())
        any_inf = bool(np.isinf(pos_after).any())
        pos_all_finite = bool(np.isfinite(pos_after).all())
        finite_pass = bool(pos_all_finite and (not any_nan) and (not any_inf))

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

        def _local_pair_metrics(row_idx: int) -> tuple[float, float]:
            if row_idx < 0:
                return float("nan"), float("nan")
            i = int(self.spring_pairs[row_idx, 0])
            j = int(self.spring_pairs[row_idx, 1])
            rest = max(float(self.spring_rests_m[row_idx]), 1e-30)
            dist = float(np.linalg.norm(pos_after[j] - pos_after[i]))
            return (
                float(dist / max(self.cfg.b_m, 1e-30)),
                float(abs(dist - rest) / rest),
            )

        local_attach_first_len_over_b, local_attach_first_rel_err = _local_pair_metrics(
            self.local_attach_first_row
        )
        local_first_second_len_over_b, local_first_second_rel_err = _local_pair_metrics(
            self.local_first_second_row
        )
        local_second_third_len_over_b, local_second_third_rel_err = _local_pair_metrics(
            self.local_second_third_row
        )

        local_basal_bend_angle_deg = float("nan")
        local_basal_bend_err_deg = float("nan")
        if np.all(self.local_basal_triplet >= 0):
            basal_angles = _triplet_angles_rad(
                pos_after, self.local_basal_triplet.reshape(1, 3)
            )
            local_basal_bend_angle_deg = float(basal_angles[0] * RAD_TO_DEG)
            local_basal_bend_err_deg = float(
                max(
                    float(self.cfg.hook.threshold_deg) - local_basal_bend_angle_deg, 0.0
                )
            )

        local_first_torsion_angle_deg = float("nan")
        local_first_torsion_err_deg = float("nan")
        if np.all(self.local_first_torsion_quad >= 0) and self.local_flag_id >= 0:
            tors = _torsion_angles_rad(
                pos_after, self.local_first_torsion_quad.reshape(1, 4)
            )
            local_first_torsion_angle_deg = float(tors[0] * RAD_TO_DEG)
            state = int(self.model.flag_states[self.local_flag_id])
            key = (
                "normal"
                if state == int(PolymorphState.NORMAL)
                else "semicoiled"
                if state == int(PolymorphState.SEMICOILED)
                else "curly1"
            )
            target = np.deg2rad(float(self.torsion_map[key]))
            local_first_torsion_err_deg = float(
                abs(float(_wrap_angle(float(tors[0]) - float(target)))) * RAD_TO_DEG
            )

        flag_root_azimuth_deg = float("nan")
        flag_phase_deg = float("nan")
        flag_phase_rate_hz = float("nan")
        flag_body_phase_diff_deg = float("nan")
        flag_helix_spin_phase_deg = float("nan")
        flag_helix_spin_rate_hz = float("nan")
        flag_helix_spin_fit_r2 = float("nan")
        if self.local_attach_idx >= 0 and self.local_first_idx >= 0:
            axis, e1, e2 = self._phase_reference_frame(pos_after)
            attach_vec = pos_after[self.local_attach_idx]
            first_vec = pos_after[self.local_first_idx]
            root_vec = first_vec - attach_vec
            root_proj = root_vec - float(np.dot(root_vec, axis)) * axis
            flag_root_azimuth_deg = self._azimuth_deg(root_proj, e1, e2)

            body_phase_deg = float("nan")
            if len(self.model.body_layer_indices) > 0:
                body_layer = self.model.body_layer_indices[0].astype(int, copy=False)
                body_center = np.mean(pos_after[body_layer], axis=0)
                if body_layer.size > 0:
                    body_marker = pos_after[int(body_layer[0])] - body_center
                    body_marker = body_marker - float(np.dot(body_marker, axis)) * axis
                    body_phase_deg = self._azimuth_deg(body_marker, e1, e2)

            if np.isfinite(flag_root_azimuth_deg):
                if (
                    self.prev_flag_root_azimuth_deg is None
                    or self.prev_flag_phase_deg is None
                ):
                    flag_phase_deg = flag_root_azimuth_deg
                    flag_phase_rate_hz = 0.0
                else:
                    delta = self._wrap_deg(
                        flag_root_azimuth_deg - self.prev_flag_root_azimuth_deg
                    )
                    flag_phase_deg = self.prev_flag_phase_deg + delta
                    prev_t_s = self.prev_flag_phase_t_s
                    if prev_t_s is not None:
                        dt_s = max(float(t_star * self.cfg.tau_s - prev_t_s), 1e-30)
                        flag_phase_rate_hz = (
                            (flag_phase_deg - self.prev_flag_phase_deg) / 360.0 / dt_s
                        )
                if np.isfinite(body_phase_deg):
                    flag_body_phase_diff_deg = self._wrap_deg(
                        flag_root_azimuth_deg - body_phase_deg
                    )
                self.prev_flag_root_azimuth_deg = flag_root_azimuth_deg
                self.prev_flag_phase_deg = flag_phase_deg
                self.prev_flag_phase_t_s = float(t_star * self.cfg.tau_s)

        if self.local_flag_id >= 0:
            flag_chain = self.model.flagella_indices[self.local_flag_id].astype(
                int, copy=False
            )
            body_frame_pos = _to_body_frame(pos_after, self.model)
            spin_offset_deg, flag_helix_spin_fit_r2 = _estimate_helix_spin_offset_deg(
                body_frame_pos[flag_chain]
            )
            if np.isfinite(spin_offset_deg):
                if (
                    self.prev_flag_helix_spin_offset_deg is None
                    or self.prev_flag_helix_spin_phase_deg is None
                ):
                    flag_helix_spin_phase_deg = spin_offset_deg
                    flag_helix_spin_rate_hz = 0.0
                else:
                    delta = self._wrap_deg(
                        spin_offset_deg - self.prev_flag_helix_spin_offset_deg
                    )
                    flag_helix_spin_phase_deg = (
                        self.prev_flag_helix_spin_phase_deg + delta
                    )
                    prev_t_s = self.prev_flag_helix_spin_t_s
                    if prev_t_s is not None:
                        dt_s = max(float(t_star * self.cfg.tau_s - prev_t_s), 1e-30)
                        flag_helix_spin_rate_hz = (
                            (
                                flag_helix_spin_phase_deg
                                - self.prev_flag_helix_spin_phase_deg
                            )
                            / 360.0
                            / dt_s
                        )
                self.prev_flag_helix_spin_offset_deg = spin_offset_deg
                self.prev_flag_helix_spin_phase_deg = flag_helix_spin_phase_deg
                self.prev_flag_helix_spin_t_s = float(t_star * self.cfg.tau_s)

        local_region_indices = [
            self.local_attach_idx,
            self.local_first_idx,
            self.local_second_idx,
            self.local_third_idx,
        ]
        local_torsion_indices = (
            self.local_first_torsion_quad.astype(int).tolist()
            if np.all(self.local_first_torsion_quad >= 0)
            else [self.local_first_idx, self.local_second_idx, self.local_third_idx]
        )
        shape_pass_nonbody, first_fail_category_nonbody = _check_nonbody_shape_pass(
            finite_pass=finite_pass,
            has_hook_pair=hook_count > 0,
            has_hook_angle=self.model.hook_triplets.size > 0,
            has_flag_bond=flag_intra_count > 0,
            has_flag_bend=self.flag_bending_rows.size > 0,
            has_flag_torsion=self.flag_torsion_rows.size > 0,
            local_attach_first_rel_err=local_attach_first_rel_err,
            hook_len_rel_err_max=hook_len_rel_err_max,
            hook_angle_err_max_deg=hook_angle_err_max_deg,
            flag_bond_rel_err_max=flag_bond_rel_err_max,
            flag_bend_err_max_deg=flag_bend_err_max_deg,
            flag_torsion_err_max_deg=flag_torsion_err_max_deg,
        )
        body_center_m, body_axis, rear_dir = _body_center_axis_rear(
            pos_after,
            self.model,
        )
        t_s = float(t_star * self.cfg.tau_s)
        body_displacement_um = float(
            np.linalg.norm(body_center_m - self.initial_body_center_m) * M_TO_UM
        )
        body_speed_um_s = float("nan")
        body_axis_step_angle_deg = 0.0
        if self.prev_body_center_m is not None and self.prev_body_t_s is not None:
            body_dt_s = max(t_s - self.prev_body_t_s, 1e-30)
            body_speed_um_s = float(
                np.linalg.norm(body_center_m - self.prev_body_center_m)
                * M_TO_UM
                / body_dt_s
            )
            body_axis_step_angle_deg = _angle_between_deg(
                self.prev_body_axis,
                body_axis,
            )
            if np.isfinite(body_axis_step_angle_deg):
                self.body_axis_cumulative_angle_deg += float(body_axis_step_angle_deg)
                omega_rad_s = np.deg2rad(float(body_axis_step_angle_deg)) / body_dt_s
                self.body_omega_square_sum += float(omega_rad_s * omega_rad_s)
                self.body_omega_count += 1
        body_wobble_deg = _angle_between_deg(self.initial_body_axis, body_axis)
        if np.isfinite(body_wobble_deg):
            self.body_axis_wobble_square_sum += float(body_wobble_deg**2)
            self.body_axis_wobble_count += 1
        body_axis_wobble_rms_deg = (
            float(
                np.sqrt(
                    self.body_axis_wobble_square_sum
                    / max(self.body_axis_wobble_count, 1)
                )
            )
            if self.body_axis_wobble_count > 0
            else float("nan")
        )
        body_angular_velocity_rms_rad_s = (
            float(np.sqrt(self.body_omega_square_sum / max(self.body_omega_count, 1)))
            if self.body_omega_count > 0
            else float("nan")
        )
        bundle_metrics = _bundle_metrics(
            pos_after,
            self.model,
            self.cfg,
            body_axis,
            rear_dir,
        )
        flag_flag_repulsion_metrics = _flag_flag_repulsion_metrics(
            pos_after,
            self.flag_bead_indices,
            self.flag_bead_ids,
            self.basal_flag_bead_indices,
            self.cfg,
        )
        attach_first_body_axis_metrics = _attach_first_body_axis_metrics(
            pos_after,
            self.model.hook_triplets,
            body_axis,
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
            "finite_pass": finite_pass,
            "shape_pass_nonbody": shape_pass_nonbody,
            "first_fail_category_nonbody": first_fail_category_nonbody,
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
            **attach_first_body_axis_metrics,
            "flag_bend_err_mean_deg": flag_bend_err_mean_deg,
            "flag_bend_err_max_deg": flag_bend_err_max_deg,
            "flag_torsion_err_mean_deg": flag_torsion_err_mean_deg,
            "flag_torsion_err_max_deg": flag_torsion_err_max_deg,
            "flag_bond_len_mean_over_b": flag_bond_len_mean_over_b,
            "flag_bond_len_min_over_b": flag_bond_len_min_over_b,
            "flag_bond_len_max_over_b": flag_bond_len_max_over_b,
            "flag_bond_rel_err_mean": flag_bond_rel_err_mean,
            "flag_bond_rel_err_max": flag_bond_rel_err_max,
            "body_displacement_um": body_displacement_um,
            "body_speed_um_s": body_speed_um_s,
            "body_axis_step_angle_deg": body_axis_step_angle_deg,
            "body_axis_cumulative_angle_deg": self.body_axis_cumulative_angle_deg,
            "body_axis_wobble_rms_deg": body_axis_wobble_rms_deg,
            "body_angular_velocity_rms_rad_s": body_angular_velocity_rms_rad_s,
            **bundle_metrics,
            **flag_flag_repulsion_metrics,
            "F_total_mean_body": _mean_norm(diag.total_forces, self.body_mask),
            "F_total_mean_flag": _mean_norm(diag.total_forces, self.flag_mask),
            "F_total_mean_all": float(
                np.mean(np.linalg.norm(diag.total_forces, axis=1))
            ),
            "F_motor_mean_body": _mean_norm(diag.motor_forces, self.body_mask),
            "F_motor_mean_flag": _mean_norm(diag.motor_forces, self.flag_mask),
            "motor_ra_len_um": float(diag.motor_ra_len_m * M_TO_UM),
            "motor_rb_len_um": float(diag.motor_rb_len_m * M_TO_UM),
            "motor_Ta_norm": float(diag.motor_Ta_norm),
            "motor_Tb_norm": float(diag.motor_Tb_norm),
            "motor_Fa_norm": float(diag.motor_Fa_norm),
            "motor_Fb_norm": float(diag.motor_Fb_norm),
            "motor_axis_vs_rear_direction_angle_deg": float(
                diag.motor_axis_vs_rear_direction_angle_deg
            ),
            "local_twist_root_orientation_deg": float(
                diag.local_twist_root_orientation_deg
            ),
            "local_twist_tip_orientation_deg": float(
                diag.local_twist_tip_orientation_deg
            ),
            "local_twist_abs_mean_deg": float(diag.local_twist_abs_mean_deg),
            "local_twist_abs_max_deg": float(diag.local_twist_abs_max_deg),
            "local_twist_tip_activity_ratio": float(
                diag.local_twist_tip_activity_ratio
            ),
            "flag_root_azimuth_deg": flag_root_azimuth_deg,
            "flag_phase_deg": flag_phase_deg,
            "flag_phase_rate_hz": flag_phase_rate_hz,
            "flag_body_phase_diff_deg": flag_body_phase_diff_deg,
            "flag_helix_spin_phase_deg": flag_helix_spin_phase_deg,
            "flag_helix_spin_rate_hz": flag_helix_spin_rate_hz,
            "flag_helix_spin_fit_r2": flag_helix_spin_fit_r2,
            "motor_attach_force_norm": float(diag.motor_attach_force_norm),
            "motor_first_force_norm": float(diag.motor_first_force_norm),
            "motor_second_force_norm": float(diag.motor_second_force_norm),
            "motor_split_residual_norm": float(diag.motor_split_residual_norm),
            "motor_ta_dot_ra_abs": float(diag.motor_ta_dot_ra_abs),
            "motor_tb_dot_rb_abs": float(diag.motor_tb_dot_rb_abs),
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
            "F_body_equiv_load_mean": float(diag.body_equiv_force_mean),
            "F_body_equiv_load_max": float(diag.body_equiv_force_max),
            "body_equiv_load_mode": str(diag.body_equiv_load_mode),
            "body_equiv_load_target_torque_Nm": float(
                diag.body_equiv_load_target_torque_Nm
            ),
            "body_equiv_load_target_force_N": float(
                diag.body_equiv_load_target_force_N
            ),
            "body_equiv_attach_region_id": int(diag.body_equiv_attach_region_id),
            "local_attach_first_len_over_b": local_attach_first_len_over_b,
            "local_attach_first_rel_err": local_attach_first_rel_err,
            "local_first_second_len_over_b": local_first_second_len_over_b,
            "local_first_second_rel_err": local_first_second_rel_err,
            "local_second_third_len_over_b": local_second_third_len_over_b,
            "local_second_third_rel_err": local_second_third_rel_err,
            "local_basal_bend_angle_deg": local_basal_bend_angle_deg,
            "local_basal_bend_err_deg": local_basal_bend_err_deg,
            "local_first_torsion_angle_deg": local_first_torsion_angle_deg,
            "local_first_torsion_err_deg": local_first_torsion_err_deg,
            "local_F_spring_attach_first": _mean_norm_indices(
                diag.spring_forces, [self.local_attach_idx, self.local_first_idx]
            ),
            "local_F_spring_first_second": _mean_norm_indices(
                diag.spring_forces, [self.local_first_idx, self.local_second_idx]
            ),
            "local_F_spring_second_third": _mean_norm_indices(
                diag.spring_forces, [self.local_second_idx, self.local_third_idx]
            ),
            "local_F_bend_basal": _mean_norm_indices(
                diag.hook_forces,
                [self.local_attach_idx, self.local_first_idx, self.local_second_idx],
            ),
            "local_F_torsion_first": _mean_norm_indices(
                diag.torsion_forces,
                local_torsion_indices,
            ),
            "local_F_motor_attach": float(diag.motor_attach_force_norm),
            "local_F_motor_first": float(diag.motor_first_force_norm),
            "local_F_motor_second": float(diag.motor_second_force_norm),
            "local_F_repulsion_basal_region": _mean_norm_indices(
                diag.repulsion_forces,
                local_region_indices,
            ),
            "flag_state_changed": flag_state_changed,
            "brownian_enabled": bool(diag.brownian_enabled),
            "brownian_disp_mean_um": brownian_disp_mean_um,
        }
        if self._writer is None or self._csv_fp is None:
            raise RuntimeError("step summary writer is not initialized")
        self._writer.writerow(row)
        self._csv_fp.flush()
        self.last_row = row
        self.prev_flag_states = self.model.flag_states.copy()
        self.prev_body_center_m = body_center_m.copy()
        self.prev_body_axis = body_axis.copy()
        self.prev_body_t_s = t_s

    def write_csv(self) -> Path:
        if self._csv_fp is not None:
            self._csv_fp.close()
            self._csv_fp = None
            self._writer = None
        return self.step_summary_path
