"""モデル幾何とトポロジの構築。"""

from __future__ import annotations

import math

import numpy as np

from sim_swim.model.types import PolymorphState, SimModel
from sim_swim.sim.params import SimulationConfig

UM_TO_M = 1e-6


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


def _ensure_arr2(values: list[tuple[int, int]]) -> np.ndarray:
    if not values:
        return np.zeros((0, 2), dtype=int)
    return np.asarray(values, dtype=int)


def _ensure_arr3(values: list[tuple[int, int, int]]) -> np.ndarray:
    if not values:
        return np.zeros((0, 3), dtype=int)
    return np.asarray(values, dtype=int)


def _ensure_arr4(values: list[tuple[int, int, int, int]]) -> np.ndarray:
    if not values:
        return np.zeros((0, 4), dtype=int)
    return np.asarray(values, dtype=int)


class ModelBuilder:
    """設定から bead-spring モデルを構築する。"""

    def __init__(self, cfg: SimulationConfig):
        self.cfg = cfg
        self.rng = np.random.default_rng(cfg.seed.global_seed)

    def _body_prism_um(
        self,
    ) -> tuple[np.ndarray, list[np.ndarray], np.ndarray, np.ndarray]:
        cfg = self.cfg
        prism = cfg.body.prism

        n_prism = max(3, prism.n_prism)
        n_layers = max(2, prism.n_layers)
        dz_um = prism.dz_over_b * cfg.scale.b_um
        radius_um = prism.radius_over_b * cfg.scale.b_um

        x0 = -0.5 * (n_layers - 1) * dz_um

        body_points: list[np.ndarray] = []
        layer_indices: list[np.ndarray] = []
        ring_edges: list[tuple[int, int]] = []
        vertical_edges: list[tuple[int, int]] = []

        for layer_idx in range(n_layers):
            layer_start = len(body_points)
            if prism.axis != "x":
                raise ValueError(
                    f"Unsupported body.prism.axis='{prism.axis}'. v1は'x'のみ対応。"
                )

            x = x0 + layer_idx * dz_um
            for k in range(n_prism):
                angle = 2.0 * math.pi * k / n_prism
                y = radius_um * math.cos(angle)
                z = radius_um * math.sin(angle)
                body_points.append(np.array([x, y, z], dtype=float))

            idx = np.arange(layer_start, layer_start + n_prism, dtype=int)
            layer_indices.append(idx)

            for k in range(n_prism):
                i = int(idx[k])
                j = int(idx[(k + 1) % n_prism])
                ring_edges.append((i, j))

            if layer_idx > 0:
                prev = layer_indices[layer_idx - 1]
                for k in range(n_prism):
                    vertical_edges.append((int(prev[k]), int(idx[k])))

        return (
            np.asarray(body_points, dtype=float),
            layer_indices,
            _ensure_arr2(ring_edges),
            _ensure_arr2(vertical_edges),
        )

    def _flag_attach_indices(
        self, terminal_layer: np.ndarray, n_flagella: int
    ) -> np.ndarray:
        n_terminal = terminal_layer.shape[0]
        if n_flagella <= 0:
            return np.zeros((0,), dtype=int)
        slots = np.linspace(0.0, float(n_terminal), num=n_flagella, endpoint=False)
        picks = np.floor(slots).astype(int) % max(n_terminal, 1)
        return terminal_layer[picks]

    def build(self) -> SimModel:
        """シミュレーションモデルを構築して返す。"""

        cfg = self.cfg
        b_um = cfg.scale.b_um

        body_um, body_layers, ring_edges, vertical_edges = self._body_prism_um()
        n_body = body_um.shape[0]

        n_flagella = max(0, cfg.flagella.n_flagella)
        ds_flag_um = cfg.flagella.discretization.ds_over_b * b_um
        ds_flag_um = max(ds_flag_um, 1e-6)
        L_flag_um = cfg.flagella.length_over_b * b_um
        n_flag = max(2, int(math.floor(L_flag_um / ds_flag_um)) + 1)

        bond_L_flag_um = max(cfg.flagella.bond_L_over_b * b_um, 1e-6)

        attach_ids = self._flag_attach_indices(body_layers[-1], n_flagella)

        points_all = [body_um]
        flagella_indices: list[np.ndarray] = []
        spring_pairs: list[tuple[int, int]] = []
        spring_rest_lengths_m: list[float] = []
        bending_triplets: list[tuple[int, int, int]] = []
        bending_flag_ids: list[int] = []
        torsion_quads: list[tuple[int, int, int, int]] = []
        torsion_flag_ids: list[int] = []
        hook_triplets: list[tuple[int, int, int]] = []

        # Body edges as springs (ring + vertical)
        for i, j in ring_edges:
            spring_pairs.append((int(i), int(j)))
            dist_m = float(
                np.linalg.norm((body_um[int(i)] - body_um[int(j)]) * UM_TO_M)
            )
            spring_rest_lengths_m.append(dist_m)
        for i, j in vertical_edges:
            spring_pairs.append((int(i), int(j)))
            dist_m = float(
                np.linalg.norm((body_um[int(i)] - body_um[int(j)]) * UM_TO_M)
            )
            spring_rest_lengths_m.append(dist_m)

        # Body bending/torsion along each vertical chain
        n_prism = body_layers[0].shape[0]
        n_layers = len(body_layers)
        for k in range(n_prism):
            chain = [int(body_layers[layer_idx][k]) for layer_idx in range(n_layers)]
            for t in range(len(chain) - 2):
                bending_triplets.append((chain[t], chain[t + 1], chain[t + 2]))
                bending_flag_ids.append(-1)
            for t in range(len(chain) - 3):
                torsion_quads.append(
                    (chain[t], chain[t + 1], chain[t + 2], chain[t + 3])
                )
                torsion_flag_ids.append(-1)

        # Flagella
        start_index = n_body
        for f_id, attach_idx in enumerate(attach_ids):
            phase = 2.0 * math.pi * (f_id / max(n_flagella, 1))
            s = np.linspace(0.0, L_flag_um, n_flag)
            theta = (
                2.0
                * math.pi
                * s
                / max(cfg.flagella.helix_init.pitch_over_b * b_um, 1e-6)
                + phase
            )
            r_um = cfg.flagella.helix_init.radius_over_b * b_um

            x = s
            y = r_um * np.cos(theta) - r_um * math.cos(phase)
            z = r_um * np.sin(theta) - r_um * math.sin(phase)
            local = np.stack([x, y, z], axis=1)
            flag_points = local + body_um[int(attach_idx)]

            points_all.append(flag_points)
            idx = np.arange(start_index, start_index + n_flag, dtype=int)
            flagella_indices.append(idx)
            start_index += n_flag

            for j in range(idx.shape[0] - 1):
                spring_pairs.append((int(idx[j]), int(idx[j + 1])))
                spring_rest_lengths_m.append(bond_L_flag_um * UM_TO_M)
            for j in range(idx.shape[0] - 2):
                bending_triplets.append((int(idx[j]), int(idx[j + 1]), int(idx[j + 2])))
                bending_flag_ids.append(f_id)
            for j in range(idx.shape[0] - 3):
                torsion_quads.append(
                    (int(idx[j]), int(idx[j + 1]), int(idx[j + 2]), int(idx[j + 3]))
                )
                torsion_flag_ids.append(f_id)

            hook_triplets.append((int(attach_idx), int(idx[0]), int(idx[1])))

        positions_um = np.concatenate(points_all, axis=0)
        positions_m = positions_um * UM_TO_M

        spring_pairs_arr = _ensure_arr2(spring_pairs)
        spring_rest_arr = np.asarray(spring_rest_lengths_m, dtype=float)
        bending_triplets_arr = _ensure_arr3(bending_triplets)
        bending_flag_ids_arr = np.asarray(bending_flag_ids, dtype=int)
        torsion_quads_arr = _ensure_arr4(torsion_quads)
        torsion_flag_ids_arr = np.asarray(torsion_flag_ids, dtype=int)
        hook_triplets_arr = _ensure_arr3(hook_triplets)

        segment_pair_indices = _segment_pairs_without_neighbors(spring_pairs_arr)

        reverse_n = int(np.clip(cfg.motor.reverse_n_flagella, 0, n_flagella))
        if reverse_n > 0:
            reverse_flagella = np.sort(
                self.rng.choice(np.arange(n_flagella), size=reverse_n, replace=False)
            )
        else:
            reverse_flagella = np.zeros((0,), dtype=int)

        return SimModel(
            positions_m=positions_m,
            body_indices=np.arange(n_body, dtype=int),
            body_layer_indices=[arr.copy() for arr in body_layers],
            body_ring_edges=ring_edges,
            body_vertical_edges=vertical_edges,
            flagella_indices=flagella_indices,
            spring_pairs=spring_pairs_arr,
            spring_rest_lengths_m=spring_rest_arr,
            bending_triplets=bending_triplets_arr,
            bending_flag_ids=bending_flag_ids_arr,
            torsion_quads=torsion_quads_arr,
            torsion_flag_ids=torsion_flag_ids_arr,
            hook_triplets=hook_triplets_arr,
            motor_triplets=hook_triplets_arr.copy(),
            segment_pair_indices=segment_pair_indices,
            bead_radius_m=cfg.bead_radius_m,
            b_m=cfg.b_m,
            reverse_flagella=reverse_flagella,
            flag_states=np.full((n_flagella,), int(PolymorphState.NORMAL), dtype=int),
            torque_signs=np.ones((n_flagella,), dtype=float),
        )
