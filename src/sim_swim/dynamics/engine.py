"""overdamped dynamics エンジン。"""

from __future__ import annotations

import math

import numpy as np

from sim_swim.dynamics.brownian import sample_brownian_displacement
from sim_swim.dynamics.forces import (
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


class DynamicsEngine:
    """力計算と時間積分を担うクラス。"""

    def __init__(self, model: SimModel, cfg: SimulationConfig):
        self.model = model
        self.cfg = cfg
        self.t = 0.0
        self.rng = np.random.default_rng(cfg.seed.global_seed)

        thermal = cfg.thermal_energy_J
        b_m = cfg.b_m
        self.spring_h = (
            cfg.potentials.spring.H_over_T_over_b * thermal / max(b_m, 1e-30)
        )
        self.spring_s_m = cfg.potentials.spring.s * b_m
        self.k_bend = cfg.potentials.bend.kb_over_T * thermal
        self.k_torsion = cfg.potentials.torsion.kt_over_T * thermal
        self.k_hook = cfg.hook.kb_over_T * thermal
        self.repulsion_A = cfg.potentials.spring_spring_repulsion.A_ss_over_T * thermal
        self.repulsion_a_m = cfg.potentials.spring_spring_repulsion.a_ss_over_b * b_m
        self.repulsion_cutoff_m = (
            cfg.potentials.spring_spring_repulsion.cutoff_over_b * b_m
        )

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

        theta0 = np.zeros((self.model.bending_triplets.shape[0],), dtype=float)
        phi0 = np.zeros((self.model.torsion_quads.shape[0],), dtype=float)

        for i, flag_id in enumerate(self.model.bending_flag_ids):
            if flag_id < 0:
                theta0[i] = math.pi
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
                phi0[i] = 0.0
                continue
            state = int(self.model.flag_states[int(flag_id)])
            key = (
                "normal"
                if state == int(PolymorphState.NORMAL)
                else "semicoiled"
                if state == int(PolymorphState.SEMICOILED)
                else "curly1"
            )
            phi0[i] = math.radians(float(torsion_map[key]))

        return theta0, phi0

    def _update_run_tumble_state(self) -> None:
        rt = self.cfg.run_tumble
        tau = self.cfg.tau_s

        run_s = max(rt.run_tau * tau, 0.0)
        tumble_s = max(rt.tumble_tau * tau, 0.0)
        semicoiled_s = max(rt.semicoiled_tau * tau, 0.0)
        curly1_s = max(rt.curly1_tau * tau, 0.0)

        self.model.flag_states[:] = int(PolymorphState.NORMAL)
        self.model.torque_signs[:] = 1.0

        if self.model.flag_states.size == 0:
            return

        cycle = max(run_s + tumble_s, 1e-12)
        phase = self.t % cycle
        if phase < run_s:
            return

        tumble_phase = phase - run_s
        reversed_flags = self.model.reverse_flagella
        if reversed_flags.size == 0:
            return

        self.model.torque_signs[reversed_flags] = -1.0
        if tumble_phase < semicoiled_s:
            self.model.flag_states[reversed_flags] = int(PolymorphState.SEMICOILED)
        elif tumble_phase < (semicoiled_s + curly1_s):
            self.model.flag_states[reversed_flags] = int(PolymorphState.CURLY1)
        else:
            self.model.flag_states[reversed_flags] = int(PolymorphState.NORMAL)

    def step(self, dt: float) -> None:
        """1ステップ更新する。

        Args:
            dt: 時間刻み [s]
        """

        dt_eff = max(float(dt), 0.0)
        self._update_run_tumble_state()

        theta0, phi0 = self._state_angles_rad()

        pos = self.model.positions_m
        forces = np.zeros_like(pos)

        forces += compute_spring_forces(
            positions_m=pos,
            spring_pairs=self.model.spring_pairs,
            spring_rest_lengths_m=self.model.spring_rest_lengths_m,
            h_const=self.spring_h,
            s_limit_m=self.spring_s_m,
        )
        forces += compute_bending_forces(
            positions_m=pos,
            triplets=self.model.bending_triplets,
            theta0_rad=theta0,
            kb=self.k_bend,
        )
        forces += compute_torsion_forces(
            positions_m=pos,
            quads=self.model.torsion_quads,
            phi0_rad=phi0,
            kt=self.k_torsion,
            fd_eps_m=max(self.model.b_m * 1e-4, 1e-12),
        )

        if self.cfg.hook.enabled:
            forces += compute_hook_forces(
                positions_m=pos,
                hook_triplets=self.model.hook_triplets,
                kb_hook=self.k_hook,
                threshold_deg=self.cfg.hook.threshold_deg,
            )

        forces += compute_segment_repulsion_forces(
            positions_m=pos,
            spring_pairs=self.model.spring_pairs,
            segment_pair_indices=self.model.segment_pair_indices,
            a_ss=self.repulsion_A,
            cutoff=self.repulsion_cutoff_m,
            a_length=self.repulsion_a_m,
        )

        if self.model.motor_triplets.shape[0] > 0:
            torque_per_flag = (
                self.cfg.motor.torque_Nm
                * self.model.torque_signs[: self.model.motor_triplets.shape[0]]
            )
            forces += compute_motor_forces(
                positions_m=pos,
                motor_triplets=self.model.motor_triplets,
                torque_per_flag=torque_per_flag,
            )

        mobility = compute_rpy_mobility(
            positions_m=pos,
            bead_radius_m=self.model.bead_radius_m,
            viscosity_Pa_s=self.cfg.fluid.viscosity_Pa_s,
        )

        drift = mobility @ forces.reshape(-1)
        xi = np.zeros_like(drift)
        if self.cfg.brownian.enabled:
            xi = sample_brownian_displacement(
                mobility=mobility,
                dt=dt_eff,
                temperature_K=self.cfg.brownian.temperature_K,
                rng=self.rng,
                method=self.cfg.brownian.method,
                jitter=self.cfg.brownian.jitter,
            )

        self.model.positions_m = pos + (drift * dt_eff + xi).reshape((-1, 3))
        self.t += dt_eff
