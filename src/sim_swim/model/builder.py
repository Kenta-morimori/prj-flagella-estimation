"""モデル幾何とトポロジの構築。"""

from __future__ import annotations

from itertools import combinations
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


def _generate_chain_from_bend_torsion(
    n_beads: int,
    bond_len_um: float,
    theta_deg: float,
    phi_deg: float,
) -> np.ndarray:
    """一定 bond / bend / torsion を満たす初期鎖を生成する。"""

    n = max(2, int(n_beads))
    bond = max(float(bond_len_um), 1e-12)
    # 入力theta_degはtriplet角（0..180）。接線の回転角はその補角を使う。
    theta = math.pi - math.radians(float(theta_deg))
    phi = math.radians(float(phi_deg))

    p = [np.zeros(3, dtype=float)]
    t_prev = np.array([1.0, 0.0, 0.0], dtype=float)
    p.append(p[-1] + bond * t_prev)
    if n == 2:
        return np.asarray(p, dtype=float)

    n0 = np.array([0.0, 1.0, 0.0], dtype=float)
    b0 = np.array([0.0, 0.0, 1.0], dtype=float)
    t_curr = math.cos(theta) * t_prev + math.sin(theta) * n0
    t_curr = t_curr / max(float(np.linalg.norm(t_curr)), 1e-12)
    p.append(p[-1] + bond * t_curr)

    binormal = b0.copy()
    for _ in range(2, n - 1):
        b_vec = np.cross(t_prev, t_curr)
        b_norm = float(np.linalg.norm(b_vec))
        if b_norm > 1e-12:
            binormal = b_vec / b_norm
        normal = np.cross(binormal, t_curr)
        normal = normal / max(float(np.linalg.norm(normal)), 1e-12)
        t_next = math.cos(theta) * t_curr + math.sin(theta) * (
            math.cos(phi) * normal + math.sin(phi) * binormal
        )
        t_next = t_next / max(float(np.linalg.norm(t_next)), 1e-12)
        p.append(p[-1] + bond * t_next)
        t_prev = t_curr
        t_curr = t_next

    return np.asarray(p, dtype=float)


def _principal_axis(points: np.ndarray) -> np.ndarray:
    centered = points - np.mean(points, axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    axis = vh[0]
    norm = float(np.linalg.norm(axis))
    if norm <= 1e-12:
        return np.array([1.0, 0.0, 0.0], dtype=float)
    return axis / norm


def _axis_basis(axis: np.ndarray, phase_hint: np.ndarray) -> np.ndarray:
    axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
    v = phase_hint - float(np.dot(phase_hint, axis)) * axis
    v_norm = float(np.linalg.norm(v))
    if v_norm <= 1e-12:
        candidates = (
            np.array([1.0, 0.0, 0.0], dtype=float),
            np.array([0.0, 1.0, 0.0], dtype=float),
            np.array([0.0, 0.0, 1.0], dtype=float),
        )
        for candidate in candidates:
            v = candidate - float(np.dot(candidate, axis)) * axis
            v_norm = float(np.linalg.norm(v))
            if v_norm > 1e-12:
                break
    v = v / max(v_norm, 1e-12)
    w = np.cross(axis, v)
    w = w / max(float(np.linalg.norm(w)), 1e-12)
    return np.column_stack([axis, v, w])


def _direction_from_rear_angle(
    rear_dir: np.ndarray,
    radial_unit: np.ndarray,
    angle_deg: float,
) -> np.ndarray:
    rear = rear_dir / max(float(np.linalg.norm(rear_dir)), 1e-12)
    radial = radial_unit / max(float(np.linalg.norm(radial_unit)), 1e-12)
    lateral = radial - float(np.dot(radial, rear)) * rear
    lateral_norm = float(np.linalg.norm(lateral))
    if lateral_norm <= 1e-12:
        lateral = np.cross(rear, np.array([0.0, 0.0, 1.0], dtype=float))
        lateral_norm = float(np.linalg.norm(lateral))
        if lateral_norm <= 1e-12:
            lateral = np.cross(rear, np.array([0.0, 1.0, 0.0], dtype=float))
            lateral_norm = float(np.linalg.norm(lateral))
    lateral = lateral / max(lateral_norm, 1e-12)
    angle = math.radians(float(angle_deg))
    direction = math.cos(angle) * rear + math.sin(angle) * lateral
    return direction / max(float(np.linalg.norm(direction)), 1e-12)


def _rotation_about_axis(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    axis = axis / max(float(np.linalg.norm(axis)), 1e-12)
    x, y, z = axis
    c = math.cos(float(angle_rad))
    s = math.sin(float(angle_rad))
    one_c = 1.0 - c
    return np.array(
        [
            [c + x * x * one_c, x * y * one_c - z * s, x * z * one_c + y * s],
            [y * x * one_c + z * s, c + y * y * one_c, y * z * one_c - x * s],
            [z * x * one_c - y * s, z * y * one_c + x * s, c + z * z * one_c],
        ],
        dtype=float,
    )


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
        self.attach_seed = int(
            cfg.seed.attach_seed
            if cfg.seed.attach_seed is not None
            else cfg.seed.global_seed
        )
        phase_seed = (
            cfg.seed.phase_seed
            if cfg.seed.phase_seed is not None
            else cfg.seed.global_seed
        )
        self.phase_rng = np.random.default_rng(int(phase_seed))

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

    def _collect_attach_candidates(self, body_layers: list[np.ndarray]) -> np.ndarray:
        """中心層とその前後層から候補点を集める。"""

        if not body_layers:
            return np.zeros((0,), dtype=int)

        center_idx = len(body_layers) // 2
        candidates: list[int] = []
        for layer_idx in (center_idx - 1, center_idx, center_idx + 1):
            if 0 <= layer_idx < len(body_layers):
                candidates.extend(int(i) for i in body_layers[layer_idx])
        return np.asarray(candidates, dtype=int)

    def _seeded_attach_indices(
        self,
        body_layers: list[np.ndarray],
        n_flagella: int,
    ) -> np.ndarray:
        if n_flagella <= 0:
            return np.zeros((0,), dtype=int)
        sequences = self._seeded_attach_sequences(body_layers, n_flagella)
        if not sequences:
            raise ValueError(
                "Not enough unique attach candidates for requested "
                "flagella count:"
                f" requested={n_flagella}, available=0"
            )
        return np.asarray(
            sequences[self.attach_seed % len(sequences)],
            dtype=int,
        )

    def _seeded_attach_sequences(
        self,
        body_layers: list[np.ndarray],
        n_flagella: int,
    ) -> list[tuple[int, ...]]:
        candidates = self._layered_attach_candidates(body_layers)
        if len(candidates) < n_flagella:
            raise ValueError(
                "Not enough unique attach candidates for requested "
                "flagella count:"
                f" requested={n_flagella}, available={len(candidates)}"
            )

        prefix = self._center_priority_attach_sequences(body_layers, n_flagella)
        prefix_keys = {tuple(sorted(seq)) for seq in prefix}
        unrestricted = [
            tuple(int(i) for i in seq)
            for seq in combinations(candidates, n_flagella)
            if tuple(sorted(int(i) for i in seq)) not in prefix_keys
        ]
        return prefix + unrestricted

    def _layered_attach_candidates(
        self,
        body_layers: list[np.ndarray],
    ) -> list[int]:
        if not body_layers:
            return []

        center_idx = len(body_layers) // 2
        layer_order = [center_idx, center_idx - 1, center_idx + 1]
        candidates: list[int] = []
        seen: set[int] = set()
        for layer_idx in layer_order:
            if not (0 <= layer_idx < len(body_layers)):
                continue
            for bead_idx in body_layers[layer_idx].astype(int, copy=False).tolist():
                if int(bead_idx) in seen:
                    continue
                candidates.append(int(bead_idx))
                seen.add(int(bead_idx))
        return candidates

    def _center_priority_attach_sequences(
        self,
        body_layers: list[np.ndarray],
        n_flagella: int,
    ) -> list[tuple[int, ...]]:
        if n_flagella <= 0:
            return [()]
        if not body_layers:
            return []

        center_idx = len(body_layers) // 2
        center = [
            int(i) for i in body_layers[center_idx].astype(int, copy=False).tolist()
        ]
        if n_flagella <= len(center):
            return self._rotating_prefix_sequences(center, n_flagella)

        outer = self._outer_attach_candidates_for_center_priority(body_layers)
        n_outer = n_flagella - len(center)
        if n_outer > len(outer):
            return []
        ordered_outer = self._rank_outer_attach_combinations(outer, n_outer)
        center_tuple = tuple(center)
        return [
            center_tuple + tuple(candidate["bead_idx"] for candidate in combo)
            for combo in ordered_outer
        ]

    def _rotating_prefix_sequences(
        self,
        beads: list[int],
        n_pick: int,
    ) -> list[tuple[int, ...]]:
        if n_pick <= 0:
            return [()]
        if n_pick == len(beads):
            return [tuple(beads)]
        return [
            tuple(beads[(start + offset) % len(beads)] for offset in range(n_pick))
            for start in range(len(beads))
        ]

    def _outer_attach_candidates_for_center_priority(
        self,
        body_layers: list[np.ndarray],
    ) -> list[dict[str, int]]:
        center_idx = len(body_layers) // 2
        out: list[dict[str, int]] = []
        # Body layers are generated from rear to front along the x-axis.
        for layer_priority, layer_idx in enumerate((center_idx - 1, center_idx + 1)):
            if not (0 <= layer_idx < len(body_layers)):
                continue
            for slot_idx, bead_idx in enumerate(
                body_layers[layer_idx].astype(int, copy=False).tolist()
            ):
                out.append(
                    {
                        "bead_idx": int(bead_idx),
                        "layer_priority": int(layer_priority),
                        "slot_idx": int(slot_idx),
                    }
                )
        return out

    def _rank_outer_attach_combinations(
        self,
        outer: list[dict[str, int]],
        n_pick: int,
    ) -> list[tuple[dict[str, int], ...]]:
        combos = list(combinations(outer, n_pick))

        def sort_key(combo: tuple[dict[str, int], ...]) -> tuple[object, ...]:
            slot_counts = {int(candidate["slot_idx"]) for candidate in combo}
            max_per_slot = max(
                sum(1 for candidate in combo if int(candidate["slot_idx"]) == slot_idx)
                for slot_idx in range(3)
            )
            rear_count = sum(
                1 for candidate in combo if int(candidate["layer_priority"]) == 0
            )
            layer_slot_order = tuple(
                (
                    int(candidate["layer_priority"]),
                    int(candidate["slot_idx"]),
                    int(candidate["bead_idx"]),
                )
                for candidate in combo
            )
            return (
                -len(slot_counts),
                max_per_slot,
                -rear_count,
                layer_slot_order,
            )

        return sorted(combos, key=sort_key)

    def build(self) -> SimModel:
        """シミュレーションモデルを構築して返す。"""

        cfg = self.cfg
        b_um = cfg.scale.b_um
        if cfg.body.prism.n_prism != 3:
            raise ValueError("MVP: body.prism.n_prism must be 3")
        if not (0 <= cfg.flagella.n_flagella <= 9):
            raise ValueError("MVP: flagella.n_flagella must be in [0,9]")
        if cfg.flagella.stub_mode not in (
            "minimal_basal_stub",
            "extended_basal_stub_5",
            "full_flagella",
        ):
            raise ValueError(
                "flagella.stub_mode must be "
                "'minimal_basal_stub', "
                "'extended_basal_stub_5' or 'full_flagella', "
                f"got: {cfg.flagella.stub_mode}"
            )

        body_um, body_layers, ring_edges, vertical_edges = self._body_prism_um()
        n_body = body_um.shape[0]

        n_flagella = cfg.flagella.n_flagella
        bond_len_um = max(cfg.flagella.bond_L_over_b * b_um, 1e-6)
        n_flag = max(2, int(cfg.flagella.n_beads_per_flagellum))

        # minimal_basal_stub モード: attach, first, second のみ (3 beads)
        if cfg.flagella.stub_mode == "minimal_basal_stub":
            n_flag = 3
        # 比較実験用: attach + 4-chain の 5 beads stub
        elif cfg.flagella.stub_mode == "extended_basal_stub_5":
            n_flag = 5

        placement_mode = str(cfg.flagella.placement_mode)
        center_layer = body_layers[len(body_layers) // 2]
        if placement_mode == "uniform" and n_flagella <= 3:
            attach_ids = self._flag_attach_indices(center_layer, n_flagella)
        elif placement_mode in ("uniform", "seeded_surface"):
            candidates = self._collect_attach_candidates(body_layers)
            unique_candidates = np.unique(candidates)
            if unique_candidates.shape[0] < n_flagella:
                raise ValueError(
                    "Not enough unique attach candidates for requested "
                    "flagella count:"
                    f" requested={n_flagella}, "
                    f"available={unique_candidates.shape[0]}"
                )
            if placement_mode == "uniform":
                attach_ids = self.rng.choice(
                    unique_candidates,
                    size=n_flagella,
                    replace=False,
                )
            else:
                attach_ids = self._seeded_attach_indices(
                    body_layers,
                    n_flagella,
                )
        else:
            raise ValueError(
                "Unsupported flagella.placement_mode:"
                f" {placement_mode}. Use 'uniform' or 'seeded_surface'."
            )
        initial_phase_mode = str(cfg.flagella.initial_phase_mode)
        if initial_phase_mode == "uniform":
            initial_phases = np.asarray(
                [
                    2.0 * math.pi * (f_id / max(n_flagella, 1))
                    for f_id in range(n_flagella)
                ],
                dtype=float,
            )
        elif initial_phase_mode == "seeded":
            initial_phases = self.phase_rng.uniform(
                0.0,
                2.0 * math.pi,
                size=n_flagella,
            )
        else:
            raise ValueError(
                "Unsupported flagella.initial_phase_mode:"
                f" {initial_phase_mode}. Use 'uniform' or 'seeded'."
            )
        hook_length_um = 0.25 * b_um
        body_axis = _principal_axis(body_um)
        rear_dir = -body_axis
        body_bead_to_layer = np.full((n_body,), -1, dtype=int)
        for layer_idx, layer in enumerate(body_layers):
            body_bead_to_layer[layer.astype(int, copy=False)] = int(layer_idx)

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
        # Body side braces (both diagonals) to suppress shear deformation between layers.
        for layer_idx in range(n_layers - 1):
            prev_layer = body_layers[layer_idx]
            curr_layer = body_layers[layer_idx + 1]
            for k in range(n_prism):
                i = int(prev_layer[k])
                j = int(curr_layer[(k + 1) % n_prism])
                spring_pairs.append((i, j))
                i2 = int(prev_layer[(k + 1) % n_prism])
                j2 = int(curr_layer[k])
                spring_pairs.append((i2, j2))

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
            phase = float(initial_phases[f_id])
            initial_helix_axis_from_rear_deg = (
                cfg.flagella.initial_helix_axis_from_rear_deg
            )
            # Posterior axis alignment derives a local-to-world basis from the
            # unphased shape, then applies phase in world space to avoid
            # canceling the seeded phase during alignment.
            shape_phase = 0.0 if initial_helix_axis_from_rear_deg is not None else phase
            init_mode = str(cfg.flagella.init_mode)
            if init_mode == "legacy_radius_pitch":
                pitch_um = max(cfg.flagella.helix_init.pitch_over_b * b_um, 1e-6)
                r_um = max(cfg.flagella.helix_init.radius_over_b * b_um, 0.0)
                dx_um = _solve_helix_dx_um(
                    radius_um=r_um,
                    pitch_um=pitch_um,
                    bond_len_um=bond_len_um,
                )
                s = dx_um * np.arange(n_flag, dtype=float)
                theta = 2.0 * math.pi * s / pitch_um + shape_phase

                x = s
                y = r_um * np.cos(theta) - r_um * math.cos(shape_phase)
                z = -r_um * np.sin(theta) + r_um * math.sin(shape_phase)
                local_points = np.column_stack([x, y, z])
            elif init_mode == "paper_table1":
                theta0_normal = float(cfg.potentials.bend.theta0_deg["normal"])
                phi0_normal = float(cfg.potentials.torsion.phi0_deg["normal"])
                local_points = _generate_chain_from_bend_torsion(
                    n_beads=n_flag,
                    bond_len_um=bond_len_um,
                    theta_deg=theta0_normal,
                    phi_deg=phi0_normal,
                )
                # 複数本で初期形状が完全一致しないよう、軸回り位相だけ回す。
                c = math.cos(shape_phase)
                s = math.sin(shape_phase)
                rot = np.array(
                    [
                        [1.0, 0.0, 0.0],
                        [0.0, c, -s],
                        [0.0, s, c],
                    ],
                    dtype=float,
                )
                local_points = local_points @ rot.T
            else:
                raise ValueError(
                    "Unsupported flagella.init_mode:"
                    f" {init_mode}. Use 'legacy_radius_pitch' or 'paper_table1'."
                )
            attach_point = body_um[int(attach_idx)]
            attach_layer_idx = int(body_bead_to_layer[int(attach_idx)])
            if attach_layer_idx < 0:
                raise ValueError(
                    "Failed to locate body layer for attach bead:"
                    f" attach_idx={int(attach_idx)}"
                )
            layer_centroid = np.mean(body_um[body_layers[attach_layer_idx]], axis=0)
            radial = attach_point - layer_centroid
            radial_norm = float(np.linalg.norm(radial))
            if radial_norm <= 1e-12:
                radial_unit = rear_dir
            else:
                radial_unit = radial / radial_norm
            hook_offset_um = hook_length_um * radial_unit
            if initial_helix_axis_from_rear_deg is None:
                axis_u = radial_unit

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
            else:
                local_helix_points = local_points[1:]
                local_axis = _principal_axis(local_helix_points)
                if (
                    float(
                        np.dot(
                            local_axis,
                            local_helix_points[-1] - local_helix_points[0],
                        )
                    )
                    < 0.0
                ):
                    local_axis = -local_axis
                local_phase_hint = local_helix_points[0] - np.mean(
                    local_helix_points,
                    axis=0,
                )
                target_axis = _direction_from_rear_angle(
                    rear_dir,
                    radial_unit,
                    initial_helix_axis_from_rear_deg,
                )
                local_rot = _axis_basis(local_axis, local_phase_hint)
                world_rot = _axis_basis(target_axis, radial_unit)
                rot = world_rot @ local_rot.T

            aligned_local_points = local_points @ rot.T
            if initial_helix_axis_from_rear_deg is not None:
                phase_rot = _rotation_about_axis(target_axis, phase)
                aligned_local_points = aligned_local_points @ phase_rot.T
            flag_points = attach_point + hook_offset_um + aligned_local_points
            tangent0 = flag_points[1] - flag_points[0]
            tangent0 /= max(float(np.linalg.norm(tangent0)), 1e-12)
            angle_deg = math.degrees(
                math.acos(float(np.clip(np.dot(tangent0, rear_dir), -1.0, 1.0)))
            )
            target_angle_deg = 90.0
            if (
                initial_helix_axis_from_rear_deg is None
                and abs(angle_deg - target_angle_deg) > 10.0 + 1e-8
            ):
                raise ValueError(
                    "Flagellum base tangent is not aligned to expected direction:"
                    f" angle_deg={angle_deg:.6f}, target_angle_deg={target_angle_deg:.6f}"
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
            flagella_attach_body_indices=np.asarray(attach_ids, dtype=int),
            flagella_initial_phases_rad=np.asarray(initial_phases, dtype=float),
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
