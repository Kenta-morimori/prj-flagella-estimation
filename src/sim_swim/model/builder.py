"""モデル幾何とトポロジの構築。"""

from __future__ import annotations

import math

import numpy as np

from sim_swim.model.types import PolymorphState, SimModel
from sim_swim.sim.params import SimulationConfig

UM_TO_M = 1e-6


def _solve_helix_dx_um(radius_um: float, pitch_um: float, bond_len_um: float) -> float:
    """隣接3D距離が bond_len_um になる軸方向刻み Δx を二分法で解く。"""

    r = max(float(radius_um), 0.0)
    p = max(float(pitch_um), 1e-12)
    target_len = max(float(bond_len_um), 1e-12)
    if r <= 0.0:
        return target_len

    def dist(dx: float) -> float:
        return math.sqrt(dx * dx + (2.0 * r * math.sin(math.pi * dx / p)) ** 2)

    lo = 0.0
    hi = target_len
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if dist(mid) < target_len:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _principal_axis(points: np.ndarray) -> np.ndarray:
    centered = points - np.mean(points, axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    axis = vh[0]
    norm = float(np.linalg.norm(axis))
    if norm <= 1e-12:
        return np.array([1.0, 0.0, 0.0], dtype=float)
    return axis / norm


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

        n_prism = prism.n_prism
        n_layers = max(2, cfg.compute_body_n_layers())
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
        if cfg.body.prism.n_prism != 3:
            raise ValueError("MVP: body.prism.n_prism must be 3")
        if not (0 <= cfg.flagella.n_flagella <= 3):
            raise ValueError("MVP: flagella.n_flagella must be in [0,3]")

        body_um, body_layers, ring_edges, vertical_edges = self._body_prism_um()
        n_body = body_um.shape[0]

        n_flagella = cfg.flagella.n_flagella
        bond_len_um = max(cfg.flagella.bond_L_over_b * b_um, 1e-6)
        L_flag_um = cfg.flagella.length_over_b * b_um
        n_flag = max(2, int(math.floor(L_flag_um / bond_len_um)) + 1)

        center_layer = body_layers[len(body_layers) // 2]
        attach_ids = self._flag_attach_indices(center_layer, n_flagella)
        hook_length_um = 0.25 * b_um
        body_axis = _principal_axis(body_um)
        rear_dir = -body_axis

        points_all = [body_um]
        flagella_indices: list[np.ndarray] = []
        spring_pairs: list[tuple[int, int]] = []
        bending_triplets: list[tuple[int, int, int]] = []
        bending_flag_ids: list[int] = []
        torsion_quads: list[tuple[int, int, int, int]] = []
        torsion_flag_ids: list[int] = []
        hook_triplets: list[tuple[int, int, int]] = []
        n_prism = body_layers[0].shape[0]
        n_layers = len(body_layers)

        # Body edges as springs (ring + vertical)
        for i, j in ring_edges:
            spring_pairs.append((int(i), int(j)))
        for i, j in vertical_edges:
            spring_pairs.append((int(i), int(j)))
        # Body side braces (diagonal) to suppress shear deformation between layers.
        for layer_idx in range(n_layers - 1):
            prev_layer = body_layers[layer_idx]
            curr_layer = body_layers[layer_idx + 1]
            for k in range(n_prism):
                i = int(prev_layer[k])
                j = int(curr_layer[(k + 1) % n_prism])
                spring_pairs.append((i, j))

        # Body bending/torsion along each vertical chain
        for layer in body_layers:
            a, b, c = int(layer[0]), int(layer[1]), int(layer[2])
            bending_triplets.append((b, a, c))
            bending_flag_ids.append(-1)
            bending_triplets.append((c, b, a))
            bending_flag_ids.append(-1)
            bending_triplets.append((a, c, b))
            bending_flag_ids.append(-1)

        for k in range(n_prism):
            chain = [int(body_layers[layer_idx][k]) for layer_idx in range(n_layers)]
            for t in range(len(chain) - 2):
                bending_triplets.append((chain[t], chain[t + 1], chain[t + 2]))
                bending_flag_ids.append(-1)

        # Flagella
        start_index = n_body
        for f_id, attach_idx in enumerate(attach_ids):
            phase = 2.0 * math.pi * (f_id / max(n_flagella, 1))
            pitch_um = max(cfg.flagella.helix_init.pitch_over_b * b_um, 1e-6)
            r_um = max(cfg.flagella.helix_init.radius_over_b * b_um, 0.0)
            dx_um = _solve_helix_dx_um(
                radius_um=r_um,
                pitch_um=pitch_um,
                bond_len_um=bond_len_um,
            )
            s = dx_um * np.arange(n_flag, dtype=float)
            theta = 2.0 * math.pi * s / pitch_um + phase

            x = s
            y = r_um * np.cos(theta) - r_um * math.cos(phase)
            z = -r_um * np.sin(theta) + r_um * math.sin(phase)
            attach_point = body_um[int(attach_idx)]
            axis_u = rear_dir
            hook_offset_um = hook_length_um * axis_u

            ref = np.array([1.0, 0.0, 0.0], dtype=float)
            if abs(float(np.dot(axis_u, ref))) > 0.9:
                ref = np.array([0.0, 0.0, 1.0], dtype=float)
            axis_v = np.cross(axis_u, ref)
            axis_v_norm = float(np.linalg.norm(axis_v))
            if axis_v_norm <= 1e-12:
                ref = np.array([0.0, 1.0, 0.0], dtype=float)
                axis_v = np.cross(axis_u, ref)
                axis_v_norm = float(np.linalg.norm(axis_v))
            axis_v /= max(axis_v_norm, 1e-12)
            axis_w = np.cross(axis_u, axis_v)

            local_points = np.column_stack([x, y, z])
            local_tangent0 = local_points[1] - local_points[0]
            local_tangent_norm = max(float(np.linalg.norm(local_tangent0)), 1e-12)
            local_u = local_tangent0 / local_tangent_norm
            local_ref = np.array([0.0, 0.0, 1.0], dtype=float)
            if abs(float(np.dot(local_u, local_ref))) > 0.9:
                local_ref = np.array([0.0, 1.0, 0.0], dtype=float)
            local_v = np.cross(local_ref, local_u)
            local_v /= max(float(np.linalg.norm(local_v)), 1e-12)
            local_w = np.cross(local_u, local_v)

            local_rot = np.column_stack([local_u, local_v, local_w])
            world_rot = np.column_stack([axis_u, axis_v, axis_w])
            rot = world_rot @ local_rot.T

            flag_points = attach_point + hook_offset_um + local_points @ rot.T
            tangent0 = flag_points[1] - flag_points[0]
            tangent0 /= max(float(np.linalg.norm(tangent0)), 1e-12)
            angle_deg = math.degrees(
                math.acos(float(np.clip(np.dot(tangent0, rear_dir), -1.0, 1.0)))
            )
            if angle_deg > 10.0 + 1e-8:
                raise ValueError(
                    "Flagellum base tangent is not aligned to rear direction:"
                    f" angle_deg={angle_deg:.6f}"
                )

            points_all.append(flag_points)
            idx = np.arange(start_index, start_index + n_flag, dtype=int)
            flagella_indices.append(idx)
            start_index += n_flag

            spring_pairs.append((int(attach_idx), int(idx[0])))

            for j in range(idx.shape[0] - 1):
                spring_pairs.append((int(idx[j]), int(idx[j + 1])))
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
        if spring_pairs_arr.size == 0:
            spring_rest_arr = np.zeros((0,), dtype=float)
        else:
            spring_rest_arr = np.linalg.norm(
                positions_m[spring_pairs_arr[:, 1]]
                - positions_m[spring_pairs_arr[:, 0]],
                axis=1,
            )
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

        n_beads = positions_m.shape[0]
        bead_is_body = np.zeros((n_beads,), dtype=bool)
        bead_is_body[:n_body] = True
        bead_flag_ids = np.full((n_beads,), -1, dtype=int)
        for f_id, idxs in enumerate(flagella_indices):
            bead_flag_ids[idxs] = int(f_id)

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
            bead_is_body=bead_is_body,
            bead_flag_ids=bead_flag_ids,
        )
