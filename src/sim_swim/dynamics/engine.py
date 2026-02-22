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

    def _state_angles_rad(self) -> tuple[np.ndarray, np.ndarray]:
        bend_map = self.cfg.potentials.bend.theta0_deg or {
            "normal": 25.0,
            "semicoiled": 55.0,
            "curly1": 75.0,
        }
        torsion_map = self.cfg.potentials.torsion.phi0_deg or {
            "normal": 15.0,
            "semicoiled": 95.0,
            "curly1": 145.0,
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
        self.model.flag_states[:] = int(PolymorphState.NORMAL)
        self.model.torque_signs[:] = 1.0

        if self.model.flag_states.size == 0:
            return

        cycle = max(rt.run_tau + rt.tumble_tau, 1e-9)
        phase = self.t % cycle
        if phase < rt.run_tau:
            return

        tumble_phase = phase - rt.run_tau
        reversed_flags = self.model.reverse_flagella
        if reversed_flags.size == 0:
            return

        self.model.torque_signs[reversed_flags] = -1.0
        if tumble_phase < rt.semicoiled_tau:
            self.model.flag_states[reversed_flags] = int(PolymorphState.SEMICOILED)
        elif tumble_phase < (rt.semicoiled_tau + rt.curly1_tau):
            self.model.flag_states[reversed_flags] = int(PolymorphState.CURLY1)
        else:
            self.model.flag_states[reversed_flags] = int(PolymorphState.NORMAL)

    def _spring_rest_lengths(self) -> np.ndarray:
        return np.where(
            self.model.spring_kinds == 0,
            self.model.body_bond_L_m,
            self.model.flag_bond_L_m,
        ).astype(float)

    def step(self, dt: float) -> None:
        """1ステップ更新する。

        Args:
            dt: 時間刻み [s]
        """

        dt_eff = max(float(dt), 0.0)
        self._update_run_tumble_state()

        spring_rest = self._spring_rest_lengths()
        theta0, phi0 = self._state_angles_rad()

        pos = self.model.positions_m
        forces = np.zeros_like(pos)

        forces += compute_spring_forces(
            positions_m=pos,
            spring_pairs=self.model.spring_pairs,
            spring_rest_lengths_m=spring_rest,
            h_const=self.cfg.potentials.spring.H,
            s_limit_m=self.cfg.potentials.spring.s_um * 1e-6,
        )
        forces += compute_bending_forces(
            positions_m=pos,
            triplets=self.model.bending_triplets,
            theta0_rad=theta0,
            kb=self.cfg.potentials.bend.kb,
        )
        forces += compute_torsion_forces(
            positions_m=pos,
            quads=self.model.torsion_quads,
            phi0_rad=phi0,
            kt=self.cfg.potentials.torsion.kt,
            fd_eps_m=max(self.model.ds_m * 1e-3, 1e-12),
        )

        if self.cfg.hook.enabled:
            forces += compute_hook_forces(
                positions_m=pos,
                hook_triplets=self.model.hook_triplets,
                kb_hook=self.cfg.hook.kb,
                threshold_deg=self.cfg.hook.threshold_deg,
            )

        forces += compute_segment_repulsion_forces(
            positions_m=pos,
            spring_pairs=self.model.spring_pairs,
            segment_pair_indices=self.model.segment_pair_indices,
            a_ss=self.cfg.potentials.spring_spring_repulsion.A_ss,
            cutoff=self.cfg.potentials.spring_spring_repulsion.cutoff_um * 1e-6,
            a_length=self.cfg.potentials.spring_spring_repulsion.a_ss_um * 1e-6,
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
