"""overdamped dynamics エンジン。"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import numpy as np

from sim_swim.dynamics.brownian import sample_brownian_displacement
from sim_swim.dynamics.forces import (
    MotorForceDiagnostics,
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
    motor_ra_len_m: float
    motor_rb_len_m: float
    motor_Ta_norm: float
    motor_Tb_norm: float
    motor_Fa_norm: float
    motor_Fb_norm: float
    motor_axis_vs_rear_direction_angle_deg: float
    motor_attach_force_norm: float
    motor_first_force_norm: float
    motor_second_force_norm: float
    motor_split_residual_norm: float
    motor_ta_dot_ra_abs: float
    motor_tb_dot_rb_abs: float
    body_equiv_load_mode: str
    body_equiv_load_target_torque_Nm: float
    body_equiv_load_target_force_N: float
    body_equiv_attach_region_id: int
    body_equiv_force_mean: float
    body_equiv_force_max: float
    torsion_fd_eps_m: float


class DynamicsEngine:
    """力計算と時間積分を担うクラス。"""

    def __init__(self, model: SimModel, cfg: SimulationConfig):
        self.model = model
        self.cfg = cfg
        self.t_star = 0.0
        self.rng = np.random.default_rng(cfg.seed.global_seed)
        self.external_force_callback: (
            Callable[[np.ndarray, float], np.ndarray] | None
        ) = None

        torque = max(cfg.torque_for_forces_Nm, 1e-30)
        b_m = cfg.b_m
        self.spring_h = cfg.potentials.spring.H_over_T_over_b * torque / max(b_m, 1e-30)
        self.spring_s_m = cfg.potentials.spring.s * b_m
        self.k_bend = cfg.potentials.bend.kb_over_T * torque
        self.k_torsion = cfg.potentials.torsion.kt_over_T * torque
        self.torsion_fd_eps_m = max(cfg.potentials.torsion.fd_eps_over_b * b_m, 1e-12)
        self.k_hook = cfg.hook.kb_over_T * torque
        self.repulsion_A = cfg.potentials.spring_spring_repulsion.A_ss_over_T * torque
        self.repulsion_a_m = cfg.potentials.spring_spring_repulsion.a_ss_over_b * b_m
        self.repulsion_cutoff_m = (
            cfg.potentials.spring_spring_repulsion.cutoff_over_b * b_m
        )
        self.body_stiffness_scale = 50.0
        self.flag_bend_stiffness_scale = 300.0
        self.flag_torsion_stiffness_scale = 300.0
        self.theta0_ref_rad, self.phi0_ref_rad = self._initial_reference_angles_rad()
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
        self.body_spring_mask = np.zeros(
            (self.model.spring_pairs.shape[0],),
            dtype=bool,
        )
        self.body_spring_mask[self.body_spring_rows] = True
        if self.model.segment_pair_indices.size > 0:
            seg_pairs = self.model.segment_pair_indices
            body_body_seg = (
                self.body_spring_mask[seg_pairs[:, 0]]
                & self.body_spring_mask[seg_pairs[:, 1]]
            )
            self.segment_pair_indices_for_repulsion = seg_pairs[~body_body_seg]
        else:
            self.segment_pair_indices_for_repulsion = self.model.segment_pair_indices
        self.body_indices = self.model.body_indices.astype(int, copy=False)
        self.bead_is_body = self.model.bead_is_body.astype(bool, copy=False)

    def set_external_force_callback(
        self,
        callback: Callable[[np.ndarray, float], np.ndarray] | None,
    ) -> None:
        """各ステップの外力（N）を返すコールバックを設定する。"""

        self.external_force_callback = callback

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

    def _body_axis_unit(self, positions_m: np.ndarray) -> np.ndarray:
        if len(self.model.body_layer_indices) >= 2:
            first = self.model.body_layer_indices[0].astype(int, copy=False)
            last = self.model.body_layer_indices[-1].astype(int, copy=False)
            c_first = np.mean(positions_m[first], axis=0)
            c_last = np.mean(positions_m[last], axis=0)
            axis = c_last - c_first
            n_axis = float(np.linalg.norm(axis))
            if n_axis > 1e-18:
                return axis / n_axis

        if self.body_indices.shape[0] >= 2:
            i0 = int(self.body_indices[0])
            i1 = int(self.body_indices[-1])
            axis = positions_m[i1] - positions_m[i0]
            n_axis = float(np.linalg.norm(axis))
            if n_axis > 1e-18:
                return axis / n_axis

        axis_key = str(self.cfg.body.prism.axis).lower()
        if axis_key == "y":
            return np.array([0.0, 1.0, 0.0], dtype=float)
        if axis_key == "z":
            return np.array([0.0, 0.0, 1.0], dtype=float)
        return np.array([1.0, 0.0, 0.0], dtype=float)

    def _motor_axis_vs_rear_direction_angle_deg(self, positions_m: np.ndarray) -> float:
        if self.model.motor_triplets.size == 0:
            return float("nan")
        rear_dir = -self._body_axis_unit(positions_m)
        angles_deg: list[float] = []
        for ib, jf, kf in self.model.motor_triplets:
            r_b = positions_m[int(kf)] - positions_m[int(jf)]
            n_rb = float(np.linalg.norm(r_b))
            if n_rb <= 1e-18:
                continue
            axis = r_b / n_rb
            dot = float(np.clip(np.dot(axis, rear_dir), -1.0, 1.0))
            angles_deg.append(math.degrees(math.acos(dot)))
        if not angles_deg:
            return float("nan")
        return float(np.mean(np.asarray(angles_deg, dtype=float)))

    def _body_equiv_load_forces(self, positions_m: np.ndarray) -> np.ndarray:
        out = np.zeros_like(positions_m)
        cfg = self.cfg.body_equiv_load
        if (not cfg.enabled) or self.body_indices.size == 0:
            return out

        mode = str(cfg.mode).strip().lower()
        if mode in {"", "none", "off", "disabled"}:
            return out

        if len(self.model.body_layer_indices) == 0:
            return out
        rear_layer = self.model.body_layer_indices[0].astype(int, copy=False)
        if rear_layer.size < 2:
            return out

        axis = self._body_axis_unit(positions_m)
        rear_dir = -axis
        rear_center = np.mean(positions_m[rear_layer], axis=0)

        if mode == "pure_couple":
            torque = abs(float(cfg.target_torque_Nm))
            if torque <= 0.0:
                torque = abs(float(self.cfg.motor_torque_Nm))
            if torque <= 0.0:
                return out

            # Use the farthest rear-layer bead pair to avoid over-concentrated local twist.
            rear_pos = positions_m[rear_layer]
            pair_i = 0
            pair_j = 1
            pair_dist = -1.0
            for i in range(rear_layer.size):
                for j in range(i + 1, rear_layer.size):
                    d = float(np.linalg.norm(rear_pos[j] - rear_pos[i]))
                    if d > pair_dist:
                        pair_dist = d
                        pair_i = i
                        pair_j = j

            i0 = int(rear_layer[pair_i])
            i1 = int(rear_layer[pair_j])
            arm_vec = positions_m[i1] - positions_m[i0]
            arm_perp = arm_vec - float(np.dot(arm_vec, rear_dir)) * rear_dir
            arm = max(float(np.linalg.norm(arm_perp)), 1e-18)
            arm_hat = arm_perp / arm
            force_dir = np.cross(arm_hat, rear_dir)
            n_force = float(np.linalg.norm(force_dir))
            if n_force <= 1e-18:
                return out
            force_dir = force_dir / n_force
            f_mag = torque / arm
            out[i0] += force_dir * f_mag
            out[i1] -= force_dir * f_mag
            return out

        if mode == "attach_proxy_local":
            f_mag = abs(float(cfg.target_force_N))
            if f_mag <= 0.0:
                return out

            idx = int(cfg.attach_region_id) % int(rear_layer.size)
            i0 = int(rear_layer[idx])
            i1 = int(rear_layer[(idx + 1) % int(rear_layer.size)])
            i2 = int(rear_layer[(idx + 2) % int(rear_layer.size)])
            r0 = positions_m[i0] - rear_center
            tangential = np.cross(rear_dir, r0)
            n_t = float(np.linalg.norm(tangential))
            if n_t <= 1e-18:
                return out
            tangential = tangential / n_t

            out[i0] += tangential * f_mag
            out[i1] -= tangential * (0.5 * f_mag)
            out[i2] -= tangential * (0.5 * f_mag)
            return out

        if mode == "distributed_rear_load":
            f_mag = abs(float(cfg.target_force_N))
            if f_mag <= 0.0:
                return out
            for bead in rear_layer:
                i = int(bead)
                radial = positions_m[i] - rear_center
                tangential = np.cross(rear_dir, radial)
                n_t = float(np.linalg.norm(tangential))
                if n_t <= 1e-18:
                    continue
                out[i] += (tangential / n_t) * f_mag
            mean_force = np.mean(out[rear_layer], axis=0)
            out[rear_layer] -= mean_force[None, :]
            return out

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
            fd_eps_m=self.torsion_fd_eps_m,
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
            segment_pair_indices=self.segment_pair_indices_for_repulsion,
            a_ss=self.repulsion_A,
            cutoff=self.repulsion_cutoff_m,
            a_length=self.repulsion_a_m,
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
        motor_axis_vs_rear_direction_angle_deg = (
            self._motor_axis_vs_rear_direction_angle_deg(pos)
        )
        body_equiv_forces = self._body_equiv_load_forces(pos)

        forces = (
            spring_forces
            + bend_forces
            + torsion_forces
            + hook_forces
            + repulsion_forces
            + motor_forces
            + body_equiv_forces
        )
        if self.external_force_callback is not None:
            external_forces = self.external_force_callback(pos, self.t_star)
            if external_forces.shape != forces.shape:
                raise ValueError(
                    "External force callback must return shape "
                    f"{forces.shape}, got {external_forces.shape}."
                )
            forces = forces + external_forces

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
        self.model.positions_m = pos_after
        self.t_star += dt_star_eff
        body_equiv_norm = np.linalg.norm(body_equiv_forces, axis=1)
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
            motor_ra_len_m=float(motor_diag.ra_len_mean_m),
            motor_rb_len_m=float(motor_diag.rb_len_mean_m),
            motor_Ta_norm=float(motor_diag.Ta_norm_mean),
            motor_Tb_norm=float(motor_diag.Tb_norm_mean),
            motor_Fa_norm=float(motor_diag.Fa_norm_mean),
            motor_Fb_norm=float(motor_diag.Fb_norm_mean),
            motor_axis_vs_rear_direction_angle_deg=float(
                motor_axis_vs_rear_direction_angle_deg
            ),
            motor_attach_force_norm=float(motor_diag.attach_force_norm_mean),
            motor_first_force_norm=float(motor_diag.first_force_norm_mean),
            motor_second_force_norm=float(motor_diag.second_force_norm_mean),
            motor_split_residual_norm=float(motor_diag.split_residual_norm_mean),
            motor_ta_dot_ra_abs=float(motor_diag.ta_dot_ra_abs_mean),
            motor_tb_dot_rb_abs=float(motor_diag.tb_dot_rb_abs_mean),
            body_equiv_load_mode=str(self.cfg.body_equiv_load.mode),
            body_equiv_load_target_torque_Nm=float(
                self.cfg.body_equiv_load.target_torque_Nm
            ),
            body_equiv_load_target_force_N=float(
                self.cfg.body_equiv_load.target_force_N
            ),
            body_equiv_attach_region_id=int(self.cfg.body_equiv_load.attach_region_id),
            body_equiv_force_mean=float(np.mean(body_equiv_norm)),
            body_equiv_force_max=float(np.max(body_equiv_norm)),
            torsion_fd_eps_m=self.torsion_fd_eps_m,
        )
