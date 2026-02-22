"""モデル幾何とトポロジの構築。"""

from __future__ import annotations

import math

import numpy as np

from sim_swim.model.types import PolymorphState, SimModel
from sim_swim.sim.params import SimulationConfig

UM_TO_M = 1e-6


def _build_helix_points_um(
    length_um: float,
    n_points: int,
    pitch_um: float,
    radius_um: float,
    phase: float,
) -> np.ndarray:
    s = np.linspace(0.0, max(length_um, 0.0), n_points)
    theta = 2.0 * math.pi * s / max(pitch_um, 1e-6) + phase
    y = radius_um * np.cos(theta)
    z = radius_um * np.sin(theta)
    y -= y[0]
    z -= z[0]
    return np.stack([s, y, z], axis=1)


def _segment_pairs_without_neighbors(spring_pairs: np.ndarray) -> np.ndarray:
    pairs: list[tuple[int, int]] = []
    m = spring_pairs.shape[0]
    for i in range(m):
        ai, aj = int(spring_pairs[i, 0]), int(spring_pairs[i, 1])
        for j in range(i + 1, m):
            bi, bj = int(spring_pairs[j, 0]), int(spring_pairs[j, 1])
            if len({ai, aj, bi, bj}) < 4:
                continue
            pairs.append((i, j))
    if not pairs:
        return np.zeros((0, 2), dtype=int)
    return np.asarray(pairs, dtype=int)


class ModelBuilder:
    """設定から bead-spring モデルを構築する。"""

    def __init__(self, cfg: SimulationConfig):
        self.cfg = cfg
        self.rng = np.random.default_rng(cfg.seed.global_seed)

    def build(self) -> SimModel:
        """シミュレーションモデルを構築して返す。"""

        cfg = self.cfg
        ds_um = max(cfg.discretization.ds_um, 1e-3)
        ds_m = ds_um * UM_TO_M

        n_body = max(2, int(math.floor(cfg.body.length_total_um / ds_um)) + 1)
        n_flag = max(2, int(math.floor(cfg.flagella.length_um / ds_um)) + 1)

        body_x_um = np.linspace(
            -cfg.body.length_total_um / 2.0,
            cfg.body.length_total_um / 2.0,
            n_body,
        )
        body_um = np.stack(
            [body_x_um, np.zeros_like(body_x_um), np.zeros_like(body_x_um)],
            axis=1,
        )

        body_attach_idx = n_body - 1
        body_radius_um = cfg.body.diameter_um / 2.0

        flagella_um: list[np.ndarray] = []
        base_offsets_um: list[np.ndarray] = []
        flagella_indices: list[np.ndarray] = []

        for i in range(cfg.flagella.n_flagella):
            if cfg.flagella.n_flagella <= 1:
                angle = 0.0
            else:
                angle = 2.0 * math.pi * i / cfg.flagella.n_flagella
            if cfg.flagella.placement_mode == "random":
                angle = float(self.rng.uniform(0.0, 2.0 * math.pi))

            radial = np.array([0.0, math.cos(angle), math.sin(angle)], dtype=float)
            attach_um = body_um[body_attach_idx] + radial * body_radius_um

            helix_local_um = _build_helix_points_um(
                length_um=cfg.flagella.length_um,
                n_points=n_flag,
                pitch_um=cfg.flagella.pitch_um,
                radius_um=cfg.flagella.radius_um,
                phase=angle,
            )
            flag_points_um = helix_local_um + attach_um

            start = n_body + sum(arr.shape[0] for arr in flagella_um)
            idx = np.arange(start, start + n_flag, dtype=int)
            flagella_indices.append(idx)
            flagella_um.append(flag_points_um)
            base_offsets_um.append(attach_um - body_um.mean(axis=0))

        all_um = (
            np.concatenate([body_um, *flagella_um], axis=0) if flagella_um else body_um
        )
        positions_m = all_um * UM_TO_M

        spring_pairs: list[tuple[int, int]] = []
        spring_kinds: list[int] = []
        bending_triplets: list[tuple[int, int, int]] = []
        bending_flag_ids: list[int] = []
        torsion_quads: list[tuple[int, int, int, int]] = []
        torsion_flag_ids: list[int] = []
        hook_triplets: list[tuple[int, int, int]] = []

        for i in range(n_body - 1):
            spring_pairs.append((i, i + 1))
            spring_kinds.append(0)
        for i in range(n_body - 2):
            bending_triplets.append((i, i + 1, i + 2))
            bending_flag_ids.append(-1)
        for i in range(n_body - 3):
            torsion_quads.append((i, i + 1, i + 2, i + 3))
            torsion_flag_ids.append(-1)

        for f_id, idx in enumerate(flagella_indices):
            for j in range(idx.shape[0] - 1):
                spring_pairs.append((int(idx[j]), int(idx[j + 1])))
                spring_kinds.append(1)
            for j in range(idx.shape[0] - 2):
                bending_triplets.append((int(idx[j]), int(idx[j + 1]), int(idx[j + 2])))
                bending_flag_ids.append(f_id)
            for j in range(idx.shape[0] - 3):
                torsion_quads.append(
                    (int(idx[j]), int(idx[j + 1]), int(idx[j + 2]), int(idx[j + 3]))
                )
                torsion_flag_ids.append(f_id)

            if idx.shape[0] >= 2:
                hook_triplets.append((body_attach_idx, int(idx[0]), int(idx[1])))

        spring_pairs_arr = np.asarray(spring_pairs, dtype=int)
        spring_kinds_arr = np.asarray(spring_kinds, dtype=int)
        bending_triplets_arr = np.asarray(bending_triplets, dtype=int)
        bending_flag_ids_arr = np.asarray(bending_flag_ids, dtype=int)
        torsion_quads_arr = np.asarray(torsion_quads, dtype=int)
        torsion_flag_ids_arr = np.asarray(torsion_flag_ids, dtype=int)
        hook_triplets_arr = np.asarray(hook_triplets, dtype=int)

        segment_pair_indices = _segment_pairs_without_neighbors(spring_pairs_arr)

        n_flagella = len(flagella_indices)
        reverse_n = int(np.clip(cfg.motor.reverse_n_flagella, 0, n_flagella))
        if reverse_n > 0:
            reverse_flagella = np.sort(
                self.rng.choice(np.arange(n_flagella), size=reverse_n, replace=False)
            )
        else:
            reverse_flagella = np.zeros((0,), dtype=int)

        helix_local_um = (
            _build_helix_points_um(
                length_um=cfg.flagella.length_um,
                n_points=n_flag,
                pitch_um=cfg.flagella.pitch_um,
                radius_um=cfg.flagella.radius_um,
                phase=0.0,
            )
            if n_flag > 0
            else np.zeros((0, 3), dtype=float)
        )

        return SimModel(
            positions_m=positions_m,
            body_indices=np.arange(n_body, dtype=int),
            flagella_indices=flagella_indices,
            spring_pairs=spring_pairs_arr,
            spring_kinds=spring_kinds_arr,
            bending_triplets=bending_triplets_arr,
            bending_flag_ids=bending_flag_ids_arr,
            torsion_quads=torsion_quads_arr,
            torsion_flag_ids=torsion_flag_ids_arr,
            hook_triplets=hook_triplets_arr,
            motor_triplets=hook_triplets_arr.copy(),
            segment_pair_indices=segment_pair_indices,
            body_bond_L_m=(cfg.body.bond_L_um or ds_um) * UM_TO_M,
            flag_bond_L_m=(cfg.flagella.bond_L_um or ds_um) * UM_TO_M,
            ds_m=ds_m,
            bead_radius_m=cfg.scale.bead_radius_a_um * UM_TO_M,
            reverse_flagella=reverse_flagella,
            flag_states=np.full((n_flagella,), int(PolymorphState.NORMAL), dtype=int),
            torque_signs=np.ones((n_flagella,), dtype=float),
            base_offsets_body_um=np.asarray(base_offsets_um, dtype=float)
            if base_offsets_um
            else np.zeros((0, 3), dtype=float),
            helix_local_um=helix_local_um,
        )
