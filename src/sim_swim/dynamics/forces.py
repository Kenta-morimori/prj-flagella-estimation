"""ポテンシャル力とモータ力の計算。"""

from __future__ import annotations

import math

import numpy as np


def _safe_norm(v: np.ndarray, eps: float = 1e-18) -> float:
    return max(float(np.linalg.norm(v)), eps)


def _wrap_angle(rad: float) -> float:
    return (rad + math.pi) % (2.0 * math.pi) - math.pi


def compute_spring_forces(
    positions_m: np.ndarray,
    spring_pairs: np.ndarray,
    spring_rest_lengths_m: np.ndarray,
    h_const: float,
    s_limit_m: float,
    clamp_eps: float = 1e-6,
) -> np.ndarray:
    """Fraenkel型spring力を計算する。"""

    forces = np.zeros_like(positions_m)
    if spring_pairs.size == 0:
        return forces

    s = max(s_limit_m, 1e-15)
    h = float(h_const)

    for idx, (i, j) in enumerate(spring_pairs):
        ri = positions_m[int(i)]
        rj = positions_m[int(j)]
        dvec = ri - rj
        d = _safe_norm(dvec)
        x = d - float(spring_rest_lengths_m[idx])
        x = float(np.clip(x, -(1.0 - clamp_eps) * s, (1.0 - clamp_eps) * s))
        denom = 1.0 - (x * x) / (s * s)
        dU_dd = h * x / (denom * denom)
        fij = -dU_dd * (dvec / d)
        forces[int(i)] += fij
        forces[int(j)] -= fij

    return forces


def _bending_for_triplet(
    r_i: np.ndarray,
    r_j: np.ndarray,
    r_k: np.ndarray,
    kb: float,
    theta0_rad: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u = r_i - r_j
    v = r_k - r_j
    nu = _safe_norm(u)
    nv = _safe_norm(v)

    c = float(np.dot(u, v) / (nu * nv))
    c = float(np.clip(c, -1.0, 1.0))
    c0 = math.cos(theta0_rad)
    dU_dc = kb * (c - c0)

    dc_du = (v / (nu * nv)) - (c * u / (nu * nu))
    dc_dv = (u / (nu * nv)) - (c * v / (nv * nv))

    f_i = -dU_dc * dc_du
    f_k = -dU_dc * dc_dv
    f_j = -(f_i + f_k)
    return f_i, f_j, f_k


def compute_bending_forces(
    positions_m: np.ndarray,
    triplets: np.ndarray,
    theta0_rad: np.ndarray,
    kb: float,
) -> np.ndarray:
    """bending（三つ組）力を計算する。"""

    forces = np.zeros_like(positions_m)
    for t_idx, (i, j, k) in enumerate(triplets):
        f_i, f_j, f_k = _bending_for_triplet(
            positions_m[int(i)],
            positions_m[int(j)],
            positions_m[int(k)],
            kb=kb,
            theta0_rad=float(theta0_rad[t_idx]),
        )
        forces[int(i)] += f_i
        forces[int(j)] += f_j
        forces[int(k)] += f_k
    return forces


def _dihedral_angle(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> float:
    b0 = a - b
    b1 = c - b
    b2 = d - c

    b1n = b1 / _safe_norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n

    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1n, v), w))
    return math.atan2(y, x)


def _torsion_energy_quad(points: np.ndarray, phi0: float, kt: float) -> float:
    phi = _dihedral_angle(points[0], points[1], points[2], points[3])
    dphi = _wrap_angle(phi - phi0)
    return 0.5 * kt * dphi * dphi


def compute_torsion_forces(
    positions_m: np.ndarray,
    quads: np.ndarray,
    phi0_rad: np.ndarray,
    kt: float,
    fd_eps_m: float,
) -> np.ndarray:
    """torsion（四つ組）力を有限差分で計算する。"""

    forces = np.zeros_like(positions_m)
    eps = max(fd_eps_m, 1e-12)

    for q_idx, (i, j, k, ell) in enumerate(quads):
        ids = [int(i), int(j), int(k), int(ell)]
        base = positions_m[ids].copy()

        for local in range(4):
            for axis in range(3):
                plus = base.copy()
                minus = base.copy()
                plus[local, axis] += eps
                minus[local, axis] -= eps

                up = _torsion_energy_quad(plus, float(phi0_rad[q_idx]), kt)
                um = _torsion_energy_quad(minus, float(phi0_rad[q_idx]), kt)
                grad = (up - um) / (2.0 * eps)
                forces[ids[local], axis] += -grad

    return forces


def compute_hook_forces(
    positions_m: np.ndarray,
    hook_triplets: np.ndarray,
    kb_hook: float,
    threshold_deg: float,
) -> np.ndarray:
    """hook曲げ条件（閾値以下のみ有効）を計算する。"""

    forces = np.zeros_like(positions_m)
    threshold_rad = math.radians(threshold_deg)

    for i, j, k in hook_triplets:
        ri = positions_m[int(i)]
        rj = positions_m[int(j)]
        rk = positions_m[int(k)]

        u = ri - rj
        v = rk - rj
        nu = _safe_norm(u)
        nv = _safe_norm(v)
        c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
        theta = math.acos(c)

        if theta > threshold_rad:
            continue

        f_i, f_j, f_k = _bending_for_triplet(
            ri, rj, rk, kb=kb_hook, theta0_rad=math.pi / 2
        )
        forces[int(i)] += f_i
        forces[int(j)] += f_j
        forces[int(k)] += f_k

    return forces


def _closest_points_on_segments(
    p1: np.ndarray,
    q1: np.ndarray,
    p2: np.ndarray,
    q2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    u = q1 - p1
    v = q2 - p2
    w = p1 - p2
    a = float(np.dot(u, u))
    b = float(np.dot(u, v))
    c = float(np.dot(v, v))
    d = float(np.dot(u, w))
    e = float(np.dot(v, w))
    denom = a * c - b * b

    if denom <= 1e-20:
        # Segments are parallel or nearly parallel.
        # Fall back to checking the minimal distance among endpoint-to-segment pairs.
        #
        # Candidate 1: p1 to segment [p2, q2]
        c_len2 = max(c, 1e-20)
        t1 = float(np.dot(p1 - p2, v) / c_len2)
        t1 = float(np.clip(t1, 0.0, 1.0))
        pa1 = p1
        pb1 = p2 + t1 * v
        d1 = float(np.dot(pa1 - pb1, pa1 - pb1))

        # Candidate 2: q1 to segment [p2, q2]
        t2 = float(np.dot(q1 - p2, v) / c_len2)
        t2 = float(np.clip(t2, 0.0, 1.0))
        pa2 = q1
        pb2 = p2 + t2 * v
        d2 = float(np.dot(pa2 - pb2, pa2 - pb2))

        # Candidate 3: p2 to segment [p1, q1]
        a_len2 = max(a, 1e-20)
        s3 = float(np.dot(p2 - p1, u) / a_len2)
        s3 = float(np.clip(s3, 0.0, 1.0))
        pa3 = p1 + s3 * u
        pb3 = p2
        d3 = float(np.dot(pa3 - pb3, pa3 - pb3))

        # Candidate 4: q2 to segment [p1, q1]
        s4 = float(np.dot(q2 - p1, u) / a_len2)
        s4 = float(np.clip(s4, 0.0, 1.0))
        pa4 = p1 + s4 * u
        pb4 = q2
        d4 = float(np.dot(pa4 - pb4, pa4 - pb4))

        # Select the best candidate
        dists = [d1, d2, d3, d4]
        idx_min = int(np.argmin(dists))
        if idx_min == 0:
            s = 0.0
            t = t1
            pa, pb = pa1, pb1
        elif idx_min == 1:
            s = 1.0
            t = t2
            pa, pb = pa2, pb2
        elif idx_min == 2:
            s = s3
            t = 0.0
            pa, pb = pa3, pb3
        else:
            s = s4
            t = 1.0
            pa, pb = pa4, pb4
    else:
        s = np.clip((b * e - c * d) / denom, 0.0, 1.0)
        t = np.clip((a * e - b * d) / denom, 0.0, 1.0)
        pa = p1 + s * u
        pb = p2 + t * v
    dvec = pa - pb
    dist = _safe_norm(dvec)
    return pa, pb, s, t, dist


def compute_segment_repulsion_forces(
    positions_m: np.ndarray,
    spring_pairs: np.ndarray,
    segment_pair_indices: np.ndarray,
    a_ss: float,
    cutoff: float,
    a_length: float,
) -> np.ndarray:
    """spring-spring 反発力を近似計算する。"""

    forces = np.zeros_like(positions_m)
    if segment_pair_indices.size == 0:
        return forces

    cutoff_eff = max(cutoff, 0.0)
    a_eff = max(a_length, 1e-12)

    for seg_a, seg_b in segment_pair_indices:
        i, j = spring_pairs[int(seg_a)]
        k, ell = spring_pairs[int(seg_b)]

        p1 = positions_m[int(i)]
        q1 = positions_m[int(j)]
        p2 = positions_m[int(k)]
        q2 = positions_m[int(ell)]

        pa, pb, s, t, dist = _closest_points_on_segments(p1, q1, p2, q2)
        if dist >= cutoff_eff:
            continue

        dir_vec = (pa - pb) / _safe_norm(pa - pb)
        mag = (a_ss / a_eff) * math.exp(-dist / a_eff)
        f = mag * dir_vec

        forces[int(i)] += (1.0 - s) * f
        forces[int(j)] += s * f
        forces[int(k)] -= (1.0 - t) * f
        forces[int(ell)] -= t * f

    return forces


def compute_motor_forces(
    positions_m: np.ndarray,
    motor_triplets: np.ndarray,
    torque_per_flag: np.ndarray,
) -> np.ndarray:
    """総力ゼロ・指定軸トルクとなる最小ノルム力を算出する。"""

    forces = np.zeros_like(positions_m)
    if motor_triplets.size == 0:
        return forces

    for idx, (ib, jf, kf) in enumerate(motor_triplets):
        ids = [int(ib), int(jf), int(kf)]
        pts = positions_m[ids]

        axis_vec = pts[2] - pts[1]
        axis_norm = _safe_norm(axis_vec)
        if axis_norm <= 1e-18:
            continue
        e = axis_vec / axis_norm

        tau = float(torque_per_flag[idx])
        b = np.zeros(6, dtype=float)
        b[3:] = tau * e

        a = np.zeros((6, 9), dtype=float)
        a[0:3, 0:3] = np.eye(3)
        a[0:3, 3:6] = np.eye(3)
        a[0:3, 6:9] = np.eye(3)

        for row, p in enumerate(pts):
            cx = np.array(
                [
                    [0.0, -p[2], p[1]],
                    [p[2], 0.0, -p[0]],
                    [-p[1], p[0], 0.0],
                ],
                dtype=float,
            )
            col = 3 * row
            a[3:6, col : col + 3] = cx

        ata = a @ a.T
        try:
            x = a.T @ np.linalg.solve(ata, b)
        except np.linalg.LinAlgError:
            x = a.T @ np.linalg.pinv(ata) @ b

        for local, bead_idx in enumerate(ids):
            forces[bead_idx] += x[3 * local : 3 * local + 3]

    return forces
