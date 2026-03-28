"""overdamped dynamics エンジン。"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from sim_swim.dynamics.brownian import sample_brownian_displacement
from sim_swim.dynamics.forces import (
    MotorForceDiagnostics,
    compute_bead_steric_exclusion_forces,
    compute_bending_forces,
    compute_hook_forces,
    compute_motor_forces,
    compute_segment_repulsion_forces,
    compute_spring_forces,
    compute_torsion_forces,
)
from sim_swim.dynamics.hydro_rpy import compute_rpy_mobility
from sim_swim.model.types import PolymorphState, SimModel
from sim_swim.sim.params import SimulationConfig


def _triplet_angle_rad(r_i: np.ndarray, r_j: np.ndarray, r_k: np.ndarray) -> float:
    u = r_i - r_j
    v = r_k - r_j
    nu = max(float(np.linalg.norm(u)), 1e-18)
    nv = max(float(np.linalg.norm(v)), 1e-18)
    c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return math.acos(c)


def _dihedral_angle(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> float:
    b0 = a - b
    b1 = c - b
    b2 = d - c

    b1n = b1 / max(float(np.linalg.norm(b1)), 1e-18)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n

    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1n, v), w))
    return math.atan2(y, x)


def _wrap_angle(rad: float) -> float:
    return (rad + math.pi) % (2.0 * math.pi) - math.pi


@dataclass(frozen=True)
class StepDiagnostics:
    """1ステップ分の診断情報。"""

    dt_star: float
    dt_s: float
    positions_before_m: np.ndarray
    positions_after_m: np.ndarray
    spring_forces: np.ndarray
    bend_forces: np.ndarray
    torsion_forces: np.ndarray
    hook_forces: np.ndarray
    repulsion_forces: np.ndarray
    motor_forces: np.ndarray
    total_forces: np.ndarray
    brownian_enabled: bool
    brownian_disp_m: np.ndarray
    motor_degenerate_axis_count: int
    motor_split_rank_deficient_count: int
    motor_bond_length_clipped_count: int


class DynamicsEngine:
    """力計算と時間積分を担うクラス。"""

    def __init__(self, model: SimModel, cfg: SimulationConfig):
        self.model = model
        self.cfg = cfg
        self.t_star = 0.0
        self.rng = np.random.default_rng(cfg.seed.global_seed)

        torque = max(cfg.torque_for_forces_Nm, 1e-30)
        b_m = cfg.b_m
        self.spring_h = cfg.potentials.spring.H_over_T_over_b * torque / max(b_m, 1e-30)
        self.spring_s_m = cfg.potentials.spring.s * b_m
        self.k_bend = cfg.potentials.bend.kb_over_T * torque
        self.k_torsion = cfg.potentials.torsion.kt_over_T * torque
        self.k_hook = cfg.hook.kb_over_T * torque
        self.repulsion_A = cfg.potentials.spring_spring_repulsion.A_ss_over_T * torque
        self.repulsion_a_m = cfg.potentials.spring_spring_repulsion.a_ss_over_b * b_m
        self.repulsion_cutoff_m = (
            cfg.potentials.spring_spring_repulsion.cutoff_over_b * b_m
        )
        self.body_stiffness_scale = cfg.stiffness_scales.body
        self.flag_bend_stiffness_scale = cfg.stiffness_scales.flag_bend
        self.flag_torsion_stiffness_scale = cfg.stiffness_scales.flag_torsion
        self.constraint_projection_iters = 8
        self.enable_flagella_template_projection = (
            cfg.projection.enable_flagella_template_projection
        )
        self.enable_flagella_chain_length_projection_when_template_off = (
            cfg.projection.enable_flagella_chain_length_projection_when_template_off
        )
        self.enable_local_helix_when_template_off = bool(cfg.local_helix.enabled)
        self.local_helix_n_local = max(int(cfg.local_helix.n_local), 1)
        self.local_helix_k_radius = (
            float(cfg.local_helix.k_radius_over_torque) * torque / max(b_m, 1e-30)
        )
        self.local_helix_k_phase = float(cfg.local_helix.k_phase_over_torque) * torque
        k_radius_ref = torque / max(b_m, 1e-30)
        k_phase_ref = torque
        self.local_helix_alpha_radius = 0.2 * min(
            max(self.local_helix_k_radius / max(k_radius_ref, 1e-30), 0.0), 1.0
        )
        self.local_helix_alpha_phase = 0.2 * min(
            max(self.local_helix_k_phase / max(k_phase_ref, 1e-30), 0.0), 1.0
        )
        self.theta0_ref_rad, self.phi0_ref_rad = self._initial_reference_angles_rad()
        self.enable_basal_link_projection = bool(cfg.basal_link.enabled)
        self.basal_link_enforce_perpendicular_to_body_axis = bool(
            cfg.basal_link.enforce_perpendicular_to_body_axis
        )
        self.basal_link_projection_alpha = float(
            np.clip(cfg.basal_link.projection_alpha, 0.0, 1.0)
        )
        self.enable_steric_exclusion = bool(cfg.steric_exclusion.enabled)
        self.steric_epsilon = float(cfg.steric_exclusion.epsilon_over_torque) * torque
        self.steric_sigma_m = (
            float(cfg.steric_exclusion.sigma_over_2a) * 2.0 * self.model.bead_radius_m
        )
        self.steric_cutoff_m = (
            float(cfg.steric_exclusion.cutoff_over_sigma) * self.steric_sigma_m
        )
        spring_pairs = self.model.spring_pairs
        if spring_pairs.size == 0:
            self.body_spring_rows = np.zeros((0,), dtype=int)
            self.other_spring_rows = np.zeros((0,), dtype=int)
            self.hook_spring_rows = np.zeros((0,), dtype=int)
            self.flag_intra_spring_rows = np.zeros((0,), dtype=int)
        else:
            bi = self.model.bead_is_body[spring_pairs[:, 0]]
            bj = self.model.bead_is_body[spring_pairs[:, 1]]
            fi = self.model.bead_flag_ids[spring_pairs[:, 0]]
            fj = self.model.bead_flag_ids[spring_pairs[:, 1]]
            self.body_spring_rows = np.where(bi & bj)[0]
            self.other_spring_rows = np.where(~(bi & bj))[0]
            self.hook_spring_rows = np.where(np.logical_xor(bi, bj))[0]
            self.flag_intra_spring_rows = np.where(
                (~bi) & (~bj) & (fi == fj) & (fi >= 0)
            )[0]
        self.body_bending_rows = np.where(self.model.bending_flag_ids < 0)[0]
        self.flag_bending_rows = np.where(self.model.bending_flag_ids >= 0)[0]
        self.body_indices = self.model.body_indices.astype(int, copy=False)
        self.bead_is_body = self.model.bead_is_body.astype(bool, copy=False)
        self.flag_chain_edges: list[list[tuple[int, int, float]]] = []
        for idxs in self.model.flagella_indices:
            edges: list[tuple[int, int, float]] = []
            for j in range(1, int(idxs.shape[0])):
                i = int(idxs[j - 1])
                k = int(idxs[j])
                rest = float(
                    np.linalg.norm(
                        self.model.positions_m[k] - self.model.positions_m[i]
                    )
                )
                edges.append((i, k, rest))
            self.flag_chain_edges.append(edges)
        self.flag_first_bead_idx = np.full(
            (len(self.model.flagella_indices),), -1, dtype=int
        )
        self.flag_basal_rest_lengths_m = np.zeros(
            (len(self.model.flagella_indices),), dtype=float
        )
        self.flag_attach_body_idx = np.full(
            (len(self.model.flagella_indices),), -1, dtype=int
        )
        self.flag_hook_rest_lengths_m = np.zeros(
            (len(self.model.flagella_indices),), dtype=float
        )
        self.hook_anchor_mask = np.zeros((self.model.positions_m.shape[0],), dtype=bool)
        if self.hook_spring_rows.size > 0:
            for row in self.hook_spring_rows:
                i_raw, j_raw = self.model.spring_pairs[row]
                i = int(i_raw)
                j = int(j_raw)
                rest = float(self.model.spring_rest_lengths_m[row])
                if self.bead_is_body[i] and (not self.bead_is_body[j]):
                    self.hook_anchor_mask[j] = True
                    f_id = int(self.model.bead_flag_ids[j])
                    if f_id >= 0:
                        self.flag_attach_body_idx[f_id] = i
                        self.flag_first_bead_idx[f_id] = j
                        self.flag_basal_rest_lengths_m[f_id] = rest
                        self.flag_hook_rest_lengths_m[f_id] = rest
                elif self.bead_is_body[j] and (not self.bead_is_body[i]):
                    self.hook_anchor_mask[i] = True
                    f_id = int(self.model.bead_flag_ids[i])
                    if f_id >= 0:
                        self.flag_attach_body_idx[f_id] = j
                        self.flag_first_bead_idx[f_id] = i
                        self.flag_basal_rest_lengths_m[f_id] = rest
                        self.flag_hook_rest_lengths_m[f_id] = rest
        self.flag_template_local: list[np.ndarray] = []
        for f_id, idxs in enumerate(self.model.flagella_indices):
            idx = idxs.astype(int, copy=False)
            p0 = self.model.positions_m[idx[0]]
            attach = int(self.flag_attach_body_idx[f_id])
            if attach >= 0:
                axis = p0 - self.model.positions_m[attach]
            else:
                axis = self.model.positions_m[idx[-1]] - p0
            axis = axis / max(float(np.linalg.norm(axis)), 1e-18)

            t = self.model.positions_m[idx[min(1, idx.shape[0] - 1)]] - p0
            v = t - np.dot(t, axis) * axis
            if float(np.linalg.norm(v)) <= 1e-12:
                ref = np.array([1.0, 0.0, 0.0], dtype=float)
                if abs(float(np.dot(ref, axis))) > 0.9:
                    ref = np.array([0.0, 1.0, 0.0], dtype=float)
                v = np.cross(axis, ref)
            v = v / max(float(np.linalg.norm(v)), 1e-18)
            w = np.cross(axis, v)

            rel = self.model.positions_m[idx] - p0
            x = rel @ axis
            y = rel @ v
            z = rel @ w
            self.flag_template_local.append(np.column_stack([x, y, z]))
        self.local_helix_target_global_indices: list[np.ndarray] = []
        self.local_helix_radius_ref: list[float] = []
        self.local_helix_delta_phi_ref: list[np.ndarray] = []
        for f_id, idxs in enumerate(self.model.flagella_indices):
            idx = idxs.astype(int, copy=False)
            n_local = min(self.local_helix_n_local, int(idx.shape[0]))
            target = idx[:n_local].copy()
            self.local_helix_target_global_indices.append(target)
            if n_local == 0:
                self.local_helix_radius_ref.append(0.0)
                self.local_helix_delta_phi_ref.append(np.zeros((0,), dtype=float))
                continue

            attach = int(self.flag_attach_body_idx[f_id])
            if attach < 0:
                self.local_helix_radius_ref.append(0.0)
                self.local_helix_delta_phi_ref.append(
                    np.zeros((max(n_local - 1, 0),), dtype=float)
                )
                continue

            p0, _, v, w = self._build_flag_frame(
                self.model.positions_m,
                idx,
                attach,
            )
            rel = self.model.positions_m[target] - p0
            y = rel @ v
            z = rel @ w
            rho = np.sqrt(y * y + z * z)
            phi = np.arctan2(z, y)
            self.local_helix_radius_ref.append(float(np.mean(rho)))
            if n_local >= 2:
                delta = np.asarray(
                    [
                        _wrap_angle(float(phi[i + 1] - phi[i]))
                        for i in range(n_local - 1)
                    ],
                    dtype=float,
                )
            else:
                delta = np.zeros((0,), dtype=float)
            self.local_helix_delta_phi_ref.append(delta)
        self.body_bead_to_layer_idx = np.full(
            (self.model.positions_m.shape[0],), -1, dtype=int
        )
        for layer_idx, layer in enumerate(self.model.body_layer_indices):
            self.body_bead_to_layer_idx[layer.astype(int, copy=False)] = int(layer_idx)
        self.body_ref_center = np.mean(
            self.model.positions_m[self.body_indices], axis=0
        )
        self.body_ref_centered = (
            self.model.positions_m[self.body_indices] - self.body_ref_center
        )
        self.steric_pair_indices = self._build_steric_exclusion_pairs(
            exclude_same_flagellum=bool(cfg.steric_exclusion.exclude_same_flagellum),
            exclude_body_body=bool(cfg.steric_exclusion.exclude_body_body),
            exclude_hook_neighbors=bool(cfg.steric_exclusion.exclude_hook_neighbors),
        )

    def _body_axis_unit(self, positions_m: np.ndarray) -> np.ndarray:
        if not self.model.body_layer_indices:
            return np.array([1.0, 0.0, 0.0], dtype=float)
        first_layer = positions_m[self.model.body_layer_indices[0]].mean(axis=0)
        last_layer = positions_m[self.model.body_layer_indices[-1]].mean(axis=0)
        axis = last_layer - first_layer
        norm = float(np.linalg.norm(axis))
        if norm <= 1e-18:
            return np.array([1.0, 0.0, 0.0], dtype=float)
        return axis / norm

    def _build_steric_exclusion_pairs(
        self,
        exclude_same_flagellum: bool,
        exclude_body_body: bool,
        exclude_hook_neighbors: bool,
    ) -> np.ndarray:
        n_beads = int(self.model.positions_m.shape[0])
        if n_beads <= 1:
            return np.zeros((0, 2), dtype=int)

        hook_neighbor_pairs: set[tuple[int, int]] = set()
        if exclude_hook_neighbors:
            for attach, first, second in self.model.hook_triplets:
                i_attach = int(attach)
                i_first = int(first)
                i_second = int(second)
                hook_neighbor_pairs.add(
                    (min(i_attach, i_first), max(i_attach, i_first))
                )
                hook_neighbor_pairs.add(
                    (min(i_first, i_second), max(i_first, i_second))
                )
                hook_neighbor_pairs.add(
                    (min(i_attach, i_second), max(i_attach, i_second))
                )

        pairs: list[tuple[int, int]] = []
        for i in range(n_beads - 1):
            for j in range(i + 1, n_beads):
                bi = bool(self.model.bead_is_body[i])
                bj = bool(self.model.bead_is_body[j])
                fi = int(self.model.bead_flag_ids[i])
                fj = int(self.model.bead_flag_ids[j])

                if bi and bj:
                    if exclude_body_body:
                        continue
                    pairs.append((i, j))
                    continue

                if (not bi) and (not bj):
                    if fi == fj and exclude_same_flagellum:
                        continue
                    if fi == fj:
                        continue
                    pair = (i, j)
                    if pair in hook_neighbor_pairs:
                        continue
                    pairs.append(pair)
                    continue

                pair = (i, j)
                if pair in hook_neighbor_pairs:
                    continue
                pairs.append(pair)

        if not pairs:
            return np.zeros((0, 2), dtype=int)
        return np.asarray(pairs, dtype=int)

    def _build_flag_frame(
        self,
        positions_m: np.ndarray,
        idx: np.ndarray,
        attach: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        p0 = positions_m[int(idx[0])]
        if attach >= 0:
            axis = p0 - positions_m[attach]
        else:
            axis = positions_m[int(idx[-1])] - p0
        axis = axis / max(float(np.linalg.norm(axis)), 1e-18)

        t = positions_m[int(idx[min(1, idx.shape[0] - 1)])] - p0
        v = t - np.dot(t, axis) * axis
        if float(np.linalg.norm(v)) <= 1e-12:
            ref = np.array([1.0, 0.0, 0.0], dtype=float)
            if abs(float(np.dot(ref, axis))) > 0.9:
                ref = np.array([0.0, 1.0, 0.0], dtype=float)
            v = np.cross(axis, ref)
        v = v / max(float(np.linalg.norm(v)), 1e-18)
        w = np.cross(axis, v)
        return p0, axis, v, w

    def _initial_reference_angles_rad(self) -> tuple[np.ndarray, np.ndarray]:
        theta0 = np.zeros((self.model.bending_triplets.shape[0],), dtype=float)
        for idx, (i, j, k) in enumerate(self.model.bending_triplets):
            theta0[idx] = _triplet_angle_rad(
                self.model.positions_m[int(i)],
                self.model.positions_m[int(j)],
                self.model.positions_m[int(k)],
            )

        phi0 = np.zeros((self.model.torsion_quads.shape[0],), dtype=float)
        for idx, (i, j, k, ell) in enumerate(self.model.torsion_quads):
            phi0[idx] = _dihedral_angle(
                self.model.positions_m[int(i)],
                self.model.positions_m[int(j)],
                self.model.positions_m[int(k)],
                self.model.positions_m[int(ell)],
            )

        return theta0, phi0

    def _state_angles_rad(self) -> tuple[np.ndarray, np.ndarray]:
        bend_map = self.cfg.potentials.bend.theta0_deg or {
            "normal": 142.0,
            "semicoiled": 90.0,
            "curly1": 105.0,
        }
        torsion_map = self.cfg.potentials.torsion.phi0_deg or {
            "normal": -60.0,
            "semicoiled": 65.0,
            "curly1": 120.0,
        }

        theta0 = self.theta0_ref_rad.copy()
        phi0 = self.phi0_ref_rad.copy()

        for i, flag_id in enumerate(self.model.bending_flag_ids):
            if flag_id < 0:
                continue
            state = int(self.model.flag_states[int(flag_id)])
            key = (
                "normal"
                if state == int(PolymorphState.NORMAL)
                else "semicoiled"
                if state == int(PolymorphState.SEMICOILED)
                else "curly1"
            )
            theta0[i] = math.radians(float(bend_map[key]))

        for i, flag_id in enumerate(self.model.torsion_flag_ids):
            if flag_id < 0:
                continue
            state = int(self.model.flag_states[int(flag_id)])
            key = (
                "normal"
                if state == int(PolymorphState.NORMAL)
                else "semicoiled"
                if state == int(PolymorphState.SEMICOILED)
                else "curly1"
            )
            phi0[i] = _wrap_angle(math.radians(float(torsion_map[key])))

        return theta0, phi0

    def _update_run_tumble_state(self) -> None:
        if not self.cfg.motor.enable_switching:
            return

        rt = self.cfg.run_tumble

        run_tau = max(rt.run_tau, 0.0)
        tumble_tau = max(rt.tumble_tau, 0.0)
        semicoiled_tau = max(rt.semicoiled_tau, 0.0)
        curly1_tau = max(rt.curly1_tau, 0.0)

        self.model.flag_states[:] = int(PolymorphState.NORMAL)
        self.model.torque_signs[:] = 1.0

        if self.model.flag_states.size == 0:
            return

        cycle = max(run_tau + tumble_tau, 1e-12)
        phase_tau = self.t_star % cycle
        if phase_tau < run_tau:
            return

        tumble_phase_tau = phase_tau - run_tau
        reversed_flags = self.model.reverse_flagella
        if reversed_flags.size == 0:
            return

        self.model.torque_signs[reversed_flags] = -1.0
        if tumble_phase_tau < semicoiled_tau:
            self.model.flag_states[reversed_flags] = int(PolymorphState.SEMICOILED)
        elif tumble_phase_tau < (semicoiled_tau + curly1_tau):
            self.model.flag_states[reversed_flags] = int(PolymorphState.CURLY1)
        else:
            self.model.flag_states[reversed_flags] = int(PolymorphState.NORMAL)

    def _project_body_rigid(self, positions_m: np.ndarray) -> np.ndarray:
        if self.body_indices.size == 0:
            return positions_m

        body_curr = positions_m[self.body_indices]
        curr_center = np.mean(body_curr, axis=0)
        curr_centered = body_curr - curr_center

        h = self.body_ref_centered.T @ curr_centered
        u, _, vt = np.linalg.svd(h, full_matrices=False)
        r = vt.T @ u.T
        if np.linalg.det(r) < 0.0:
            vt[-1, :] *= -1.0
            r = vt.T @ u.T

        body_projected = self.body_ref_centered @ r.T + curr_center
        out = positions_m.copy()
        out[self.body_indices] = body_projected
        return out

    def _project_distance_pairs(
        self,
        positions_m: np.ndarray,
        pair_rows: np.ndarray,
        iterations: int,
        fixed_mask: np.ndarray | None = None,
    ) -> np.ndarray:
        if pair_rows.size == 0 or iterations <= 0:
            return positions_m

        out = positions_m.copy()
        pairs = self.model.spring_pairs[pair_rows]
        rests = self.model.spring_rest_lengths_m[pair_rows]
        for _ in range(iterations):
            for row, (i_raw, j_raw) in enumerate(pairs):
                i = int(i_raw)
                j = int(j_raw)
                d = out[j] - out[i]
                dist = float(np.linalg.norm(d))
                if dist <= 1e-18:
                    continue

                rest = float(rests[row])
                corr = ((dist - rest) / dist) * d
                fixed_i = self.bead_is_body[i] or (
                    fixed_mask is not None and fixed_mask[i]
                )
                fixed_j = self.bead_is_body[j] or (
                    fixed_mask is not None and fixed_mask[j]
                )
                wi = 0.0 if fixed_i else 1.0
                wj = 0.0 if fixed_j else 1.0
                wsum = wi + wj
                if wsum <= 0.0:
                    continue
                out[i] += (wi / wsum) * corr
                out[j] -= (wj / wsum) * corr
        return out

    def _project_flagella_chain_lengths(
        self, positions_m: np.ndarray, iterations: int
    ) -> np.ndarray:
        if iterations <= 0 or not self.flag_chain_edges:
            return positions_m

        out = positions_m.copy()
        for _ in range(iterations):
            for edges in self.flag_chain_edges:
                for i, j, rest in edges:
                    d = out[j] - out[i]
                    dist = float(np.linalg.norm(d))
                    if dist <= 1e-18:
                        continue
                    out[j] = out[i] + (rest / dist) * d
        return out

    def _project_basal_link_direction(self, positions_m: np.ndarray) -> np.ndarray:
        if self.flag_first_bead_idx.size == 0:
            return positions_m
        out = positions_m.copy()

        if self.basal_link_enforce_perpendicular_to_body_axis:
            body_axis = self._body_axis_unit(out)
            alpha = self.basal_link_projection_alpha
            for f_id in range(self.flag_first_bead_idx.shape[0]):
                first = int(self.flag_first_bead_idx[f_id])
                attach = int(self.flag_attach_body_idx[f_id])
                if first < 0 or attach < 0:
                    continue

                link = out[first] - out[attach]
                link_norm = max(float(np.linalg.norm(link)), 1e-18)
                current_dir = link / link_norm
                target_vec = current_dir - np.dot(current_dir, body_axis) * body_axis
                target_norm = float(np.linalg.norm(target_vec))
                if target_norm <= 1e-12:
                    fallback = np.array([1.0, 0.0, 0.0], dtype=float)
                    if abs(float(np.dot(fallback, body_axis))) > 0.9:
                        fallback = np.array([0.0, 1.0, 0.0], dtype=float)
                    target_vec = np.cross(body_axis, fallback)
                    target_norm = float(np.linalg.norm(target_vec))
                target_dir = target_vec / max(target_norm, 1e-18)

                blended = (1.0 - alpha) * current_dir + alpha * target_dir
                blended_norm = float(np.linalg.norm(blended))
                if blended_norm <= 1e-18:
                    continue
                rest = max(float(self.flag_hook_rest_lengths_m[f_id]), 1e-18)
                out[first] = out[attach] + rest * (blended / blended_norm)
            return out

        layer_centroids = [
            np.mean(out[layer.astype(int, copy=False)], axis=0)
            for layer in self.model.body_layer_indices
        ]

        for f_id in range(self.flag_first_bead_idx.shape[0]):
            first = int(self.flag_first_bead_idx[f_id])
            attach = int(self.flag_attach_body_idx[f_id])
            if first < 0 or attach < 0:
                continue
            layer_idx = int(self.body_bead_to_layer_idx[attach])
            if layer_idx < 0 or layer_idx >= len(layer_centroids):
                continue
            radial = out[attach] - layer_centroids[layer_idx]
            radial_norm = float(np.linalg.norm(radial))
            if radial_norm <= 1e-18:
                radial = out[first] - out[attach]
                radial_norm = float(np.linalg.norm(radial))
                if radial_norm <= 1e-18:
                    continue

            rest = max(float(self.flag_basal_rest_lengths_m[f_id]), 0.0)
            if rest <= 0.0:
                # Fallback: use current attach→first distance as rest length
                attach_to_first = out[first] - out[attach]
                attach_to_first_norm = float(np.linalg.norm(attach_to_first))
                if attach_to_first_norm <= 1e-18:
                    # Degenerate configuration: cannot determine a meaningful rest length
                    continue
                rest = attach_to_first_norm
            out[first] = out[attach] + (rest / radial_norm) * radial
        return out

    def _project_flagella_template(self, positions_m: np.ndarray) -> np.ndarray:
        if not self.model.flagella_indices:
            return positions_m
        out = positions_m.copy()
        for f_id, idxs in enumerate(self.model.flagella_indices):
            idx = idxs.astype(int, copy=False)
            if idx.size < 2:
                continue
            attach = int(self.flag_attach_body_idx[f_id])
            if attach < 0:
                continue

            p0, axis, v, w = self._build_flag_frame(out, idx, attach)

            local = self.flag_template_local[f_id]
            out[idx] = (
                p0
                + local[:, [0]] * axis[None, :]
                + local[:, [1]] * v[None, :]
                + local[:, [2]] * w[None, :]
            )
        return out

    def _project_local_helix_constraint(self, positions_m: np.ndarray) -> np.ndarray:
        if not self.model.flagella_indices:
            return positions_m

        out = positions_m.copy()
        for f_id, idxs in enumerate(self.model.flagella_indices):
            idx = idxs.astype(int, copy=False)
            target = self.local_helix_target_global_indices[f_id]
            if target.size == 0:
                continue

            attach = int(self.flag_attach_body_idx[f_id])
            if attach < 0:
                continue

            p0, axis, v, w = self._build_flag_frame(out, idx, attach)
            rel = out[target] - p0
            y = rel @ v
            z = rel @ w
            rho = np.sqrt(y * y + z * z)
            rho_safe = np.maximum(rho, 1e-18)
            radius_ref = max(float(self.local_helix_radius_ref[f_id]), 1e-18)
            radial_scale = np.clip((radius_ref - rho) / rho_safe, -1.0, 1.0)
            radial_dir = y[:, None] * v[None, :] + z[:, None] * w[None, :]
            out[target] += (
                self.local_helix_alpha_radius * radial_scale[:, None] * radial_dir
            )

            rel2 = out[target] - p0
            y2 = rel2 @ v
            z2 = rel2 @ w
            phi = np.arctan2(z2, y2)
            delta_ref = self.local_helix_delta_phi_ref[f_id]
            if phi.size < 2 or delta_ref.size != phi.size - 1:
                continue

            for j in range(phi.size - 1):
                err = _wrap_angle(float((phi[j + 1] - phi[j]) - delta_ref[j]))
                theta_corr = -self.local_helix_alpha_phase * err

                vec_j = y2[j] * v + z2[j] * w
                vec_k = y2[j + 1] * v + z2[j + 1] * w
                rho_j = max(float(np.linalg.norm(vec_j)), 1e-18)
                rho_k = max(float(np.linalg.norm(vec_k)), 1e-18)
                tan_j = np.cross(axis, vec_j / rho_j)
                tan_k = np.cross(axis, vec_k / rho_k)
                out[int(target[j])] += 0.5 * theta_corr * rho_j * tan_j
                out[int(target[j + 1])] -= 0.5 * theta_corr * rho_k * tan_k

        return out

    def _project_flagella_constraints(self, positions_m: np.ndarray) -> np.ndarray:
        out = positions_m
        if self.enable_flagella_template_projection:
            return self._project_flagella_template(out)

        if self.enable_flagella_chain_length_projection_when_template_off:
            out = self._project_flagella_chain_lengths(out, 1)
        if self.enable_local_helix_when_template_off:
            out = self._project_local_helix_constraint(out)
        return out

    def _project_hook_and_flag_bonds(self, positions_m: np.ndarray) -> np.ndarray:
        if self.constraint_projection_iters <= 0:
            return positions_m
        out = positions_m
        for _ in range(self.constraint_projection_iters):
            out = self._project_distance_pairs(out, self.hook_spring_rows, 1)
            if self.enable_basal_link_projection:
                out = self._project_basal_link_direction(out)
            out = self._project_flagella_constraints(out)
        return out

    def step(self, dt_star: float) -> StepDiagnostics:
        """1ステップ更新する。

        Args:
            dt_star: 無次元時間刻み（dt/tau）。
        """

        dt_star_eff = max(float(dt_star), 0.0)
        dt_s = dt_star_eff * self.cfg.tau_s
        self._update_run_tumble_state()

        theta0, phi0 = self._state_angles_rad()

        pos_before = self.model.positions_m.copy()
        pos = pos_before
        spring_forces = np.zeros_like(pos)
        if self.body_spring_rows.size > 0:
            spring_forces += compute_spring_forces(
                positions_m=pos,
                spring_pairs=self.model.spring_pairs[self.body_spring_rows],
                spring_rest_lengths_m=self.model.spring_rest_lengths_m[
                    self.body_spring_rows
                ],
                h_const=self.spring_h * self.body_stiffness_scale,
                s_limit_m=self.spring_s_m,
                clamp_eps=1e-3,
            )
        if self.other_spring_rows.size > 0:
            spring_forces += compute_spring_forces(
                positions_m=pos,
                spring_pairs=self.model.spring_pairs[self.other_spring_rows],
                spring_rest_lengths_m=self.model.spring_rest_lengths_m[
                    self.other_spring_rows
                ],
                h_const=self.spring_h,
                s_limit_m=self.spring_s_m,
                clamp_eps=1e-3,
            )
        bend_forces = np.zeros_like(pos)
        if self.body_bending_rows.size > 0:
            bend_forces += compute_bending_forces(
                positions_m=pos,
                triplets=self.model.bending_triplets[self.body_bending_rows],
                theta0_rad=theta0[self.body_bending_rows],
                kb=self.k_bend * self.body_stiffness_scale,
            )
        if self.flag_bending_rows.size > 0:
            bend_forces += compute_bending_forces(
                positions_m=pos,
                triplets=self.model.bending_triplets[self.flag_bending_rows],
                theta0_rad=theta0[self.flag_bending_rows],
                kb=self.k_bend * self.flag_bend_stiffness_scale,
            )
        torsion_forces = compute_torsion_forces(
            positions_m=pos,
            quads=self.model.torsion_quads,
            phi0_rad=phi0,
            kt=self.k_torsion * self.flag_torsion_stiffness_scale,
            fd_eps_m=max(self.model.b_m * 1e-1, 1e-7),
        )

        hook_forces = np.zeros_like(pos)
        if self.cfg.hook.enabled:
            hook_forces = compute_hook_forces(
                positions_m=pos,
                hook_triplets=self.model.hook_triplets,
                kb_hook=self.k_hook,
                threshold_deg=self.cfg.hook.threshold_deg,
            )

        repulsion_forces = compute_segment_repulsion_forces(
            positions_m=pos,
            spring_pairs=self.model.spring_pairs,
            segment_pair_indices=self.model.segment_pair_indices,
            a_ss=self.repulsion_A,
            cutoff=self.repulsion_cutoff_m,
            a_length=self.repulsion_a_m,
        )

        steric_forces = np.zeros_like(pos)
        if self.enable_steric_exclusion:
            steric_forces = compute_bead_steric_exclusion_forces(
                positions_m=pos,
                bead_pair_indices=self.steric_pair_indices,
                epsilon=self.steric_epsilon,
                sigma=self.steric_sigma_m,
                cutoff=self.steric_cutoff_m,
            )

        motor_forces = np.zeros_like(pos)
        motor_diag = MotorForceDiagnostics()
        if self.model.motor_triplets.shape[0] > 0:
            torque_per_flag = (
                self.cfg.motor_torque_Nm
                * self.model.torque_signs[: self.model.motor_triplets.shape[0]]
            )
            motor_forces, motor_diag = compute_motor_forces(
                positions_m=pos,
                motor_triplets=self.model.motor_triplets,
                torque_per_flag=torque_per_flag,
            )

        forces = (
            spring_forces
            + bend_forces
            + torsion_forces
            + hook_forces
            + repulsion_forces
            + motor_forces
            + steric_forces
        )

        mobility = compute_rpy_mobility(
            positions_m=pos_before,
            bead_radius_m=self.model.bead_radius_m,
            viscosity_Pa_s=self.cfg.fluid.viscosity_Pa_s,
        )

        drift = mobility @ forces.reshape(-1)
        xi = np.zeros_like(drift)
        if self.cfg.brownian.enabled:
            xi = sample_brownian_displacement(
                mobility=mobility,
                dt=dt_s,
                temperature_K=self.cfg.brownian.temperature_K,
                rng=self.rng,
                method=self.cfg.brownian.method,
                jitter=self.cfg.brownian.jitter,
            )

        brownian_disp = xi.reshape((-1, 3))
        pos_after = pos_before + (drift * dt_s + xi).reshape((-1, 3))
        pos_after = self._project_body_rigid(pos_after)
        pos_after = self._project_hook_and_flag_bonds(pos_after)
        self.model.positions_m = pos_after
        self.t_star += dt_star_eff
        return StepDiagnostics(
            dt_star=dt_star_eff,
            dt_s=dt_s,
            positions_before_m=pos_before,
            positions_after_m=pos_after,
            spring_forces=spring_forces,
            bend_forces=bend_forces,
            torsion_forces=torsion_forces,
            hook_forces=hook_forces,
            repulsion_forces=repulsion_forces,
            motor_forces=motor_forces,
            total_forces=forces,
            brownian_enabled=bool(self.cfg.brownian.enabled),
            brownian_disp_m=brownian_disp,
            motor_degenerate_axis_count=motor_diag.degenerate_axis_count,
            motor_split_rank_deficient_count=motor_diag.split_rank_deficient_count,
            motor_bond_length_clipped_count=motor_diag.bond_length_clipped_count,
        )
