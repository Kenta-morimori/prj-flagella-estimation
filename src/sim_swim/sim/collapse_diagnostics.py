"""collapse onset diagnostics CSV writer."""

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
RAD_TO_DEG = 180.0 / np.pi

COLLAPSE_DIAGNOSTICS_COLUMNS = [
    "step",
    "t_s",
    "flagellum_id",
    "mean_radius_um",
    "std_radius_um",
    "mean_phase_diff_deg",
    "std_phase_diff_deg",
    "mean_bending_dev_deg",
    "mean_torsion_dev_deg",
    "hook_angle_deg",
    "min_interflagella_distance_um",
    "min_body_distance_um",
]

COLLAPSE_SUMMARY_COLUMNS = [
    "collapse_detected",
    "first_collapse_step",
    "first_collapse_t_s",
    "flagellum_id_at_min_distance",
    "global_min_interflagella_distance_um",
    "global_min_body_distance_um",
]


def _wrap_angle(rad: float) -> float:
    return (rad + np.pi) % (2.0 * np.pi) - np.pi


def _triplet_angle_rad(r_i: np.ndarray, r_j: np.ndarray, r_k: np.ndarray) -> float:
    u = r_i - r_j
    v = r_k - r_j
    nu = max(float(np.linalg.norm(u)), 1e-18)
    nv = max(float(np.linalg.norm(v)), 1e-18)
    c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return float(np.arccos(c))


def _torsion_angle_rad(
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
    return float(np.arctan2(y, x))


@dataclass
class CollapseDiagnosticsRecorder:
    """各ステップのcollapse診断をCSVへ保存する。"""

    model: SimModel
    cfg: SimulationConfig
    out_dir: Path
    _csv_fp: TextIO | None = field(init=False, default=None, repr=False)
    _writer: csv.DictWriter | None = field(init=False, default=None, repr=False)

    def __post_init__(self) -> None:
        collapse_cfg = self.cfg.diagnostics.collapse
        self.enabled = bool(collapse_cfg.enabled)
        if not self.enabled:
            self.diagnostics_path = self.out_dir.parent / "collapse_diagnostics.csv"
            self.summary_path = self.out_dir.parent / "collapse_summary.csv"
            return

        self.write_every_step = bool(collapse_cfg.write_every_step)
        self.max_points = max(int(collapse_cfg.max_flagella_points), 1)
        self.collapse_distance_um = float(collapse_cfg.collapse_distance_um)
        self.collapse_consecutive_steps = max(
            int(collapse_cfg.collapse_consecutive_steps),
            1,
        )

        self.diagnostics_path = self.out_dir.parent / "collapse_diagnostics.csv"
        self.summary_path = self.out_dir.parent / "collapse_summary.csv"
        self._csv_fp = self.diagnostics_path.open("w", encoding="utf-8", newline="")
        self._writer = csv.DictWriter(
            self._csv_fp,
            fieldnames=COLLAPSE_DIAGNOSTICS_COLUMNS,
        )
        self._writer.writeheader()

        self.flagella_indices = [
            idx.astype(int, copy=False) for idx in self.model.flagella_indices
        ]
        self.flag_point_indices = [
            idx[: min(idx.size, self.max_points)] for idx in self.flagella_indices
        ]
        self.body_indices = self.model.body_indices.astype(int, copy=False)

        self.attach_body_idx = np.full((len(self.flagella_indices),), -1, dtype=int)
        if self.model.spring_pairs.size > 0:
            body_mask = self.model.bead_is_body.astype(bool, copy=False)
            flag_ids = self.model.bead_flag_ids.astype(int, copy=False)
            for i_raw, j_raw in self.model.spring_pairs:
                i = int(i_raw)
                j = int(j_raw)
                if body_mask[i] and (not body_mask[j]):
                    f_id = int(flag_ids[j])
                    if f_id >= 0:
                        self.attach_body_idx[f_id] = i
                elif body_mask[j] and (not body_mask[i]):
                    f_id = int(flag_ids[i])
                    if f_id >= 0:
                        self.attach_body_idx[f_id] = j

        self.bend_rows_by_flag = [
            np.where(self.model.bending_flag_ids == f_id)[0]
            for f_id in range(len(self.flagella_indices))
        ]
        self.torsion_rows_by_flag = [
            np.where(self.model.torsion_flag_ids == f_id)[0]
            for f_id in range(len(self.flagella_indices))
        ]

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

        self._candidate_start_step = -1
        self._candidate_start_t = -1.0
        self._candidate_consecutive = 0
        self._collapse_detected = False
        self._first_collapse_step = -1
        self._first_collapse_t_s = -1.0

        self._global_min_inter_um = float("inf")
        self._global_min_body_um = float("inf")
        self._flag_id_at_min_distance = -1

    def _state_key(self, flag_id: int) -> str:
        state = int(self.model.flag_states[flag_id])
        if state == int(PolymorphState.NORMAL):
            return "normal"
        if state == int(PolymorphState.SEMICOILED):
            return "semicoiled"
        return "curly1"

    def _flag_frame(
        self,
        positions_m: np.ndarray,
        flag_id: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
        idx = self.flagella_indices[flag_id]
        if idx.size == 0:
            return None
        p0 = positions_m[int(idx[0])]
        attach = int(self.attach_body_idx[flag_id])
        if attach >= 0:
            axis = p0 - positions_m[attach]
        else:
            axis = positions_m[int(idx[-1])] - p0
        axis = axis / max(float(np.linalg.norm(axis)), 1e-18)

        t = positions_m[int(idx[min(1, idx.size - 1)])] - p0
        v = t - np.dot(t, axis) * axis
        if float(np.linalg.norm(v)) <= 1e-12:
            ref = np.array([1.0, 0.0, 0.0], dtype=float)
            if abs(float(np.dot(ref, axis))) > 0.9:
                ref = np.array([0.0, 1.0, 0.0], dtype=float)
            v = np.cross(axis, ref)
        v = v / max(float(np.linalg.norm(v)), 1e-18)
        w = np.cross(axis, v)
        return p0, axis, v, w

    def _body_axis(self, positions_m: np.ndarray) -> np.ndarray:
        first_layer = positions_m[self.model.body_layer_indices[0]].mean(axis=0)
        last_layer = positions_m[self.model.body_layer_indices[-1]].mean(axis=0)
        axis = last_layer - first_layer
        norm = float(np.linalg.norm(axis))
        if norm <= 1e-18:
            return np.array([1.0, 0.0, 0.0], dtype=float)
        return axis / norm

    def _should_write_step(self, step: int) -> bool:
        if self.cfg.time.duration_s <= 0.01:
            return True
        if self.write_every_step:
            return True
        return (step % 10) == 0

    def record(self, step: int, t_s: float, diag: StepDiagnostics) -> None:
        if not self.enabled:
            return
        if not self._should_write_step(step):
            return

        if self._writer is None or self._csv_fp is None:
            raise RuntimeError("collapse diagnostics writer is not initialized")

        pos = diag.positions_after_m
        body_axis = self._body_axis(pos)

        step_candidate = False
        for flag_id in range(len(self.flagella_indices)):
            frame = self._flag_frame(pos, flag_id)
            if frame is None:
                continue
            p0, axis, v, w = frame
            local_idx = self.flag_point_indices[flag_id]
            points = pos[local_idx]
            rel = points - p0
            y = rel @ v
            z = rel @ w
            rho_um = np.sqrt(y * y + z * z) * M_TO_UM
            phi = np.arctan2(z, y)

            mean_radius_um = float(np.mean(rho_um)) if rho_um.size > 0 else float("nan")
            std_radius_um = float(np.std(rho_um)) if rho_um.size > 0 else float("nan")
            if phi.size >= 2:
                dphi = np.asarray(
                    [
                        _wrap_angle(float(phi[i + 1] - phi[i]))
                        for i in range(phi.size - 1)
                    ],
                    dtype=float,
                )
                mean_phase_diff_deg = float(np.mean(dphi) * RAD_TO_DEG)
                std_phase_diff_deg = float(np.std(dphi) * RAD_TO_DEG)
            else:
                mean_phase_diff_deg = float("nan")
                std_phase_diff_deg = float("nan")

            bend_rows = self.bend_rows_by_flag[flag_id]
            if bend_rows.size > 0:
                state_key = self._state_key(flag_id)
                theta0 = np.deg2rad(float(self.bend_map[state_key]))
                vals = []
                for row in bend_rows:
                    i, j, k = self.model.bending_triplets[int(row)]
                    vals.append(
                        abs(
                            _triplet_angle_rad(
                                pos[int(i)],
                                pos[int(j)],
                                pos[int(k)],
                            )
                            - theta0
                        )
                        * RAD_TO_DEG
                    )
                mean_bending_dev_deg = float(np.mean(vals))
            else:
                mean_bending_dev_deg = float("nan")

            torsion_rows = self.torsion_rows_by_flag[flag_id]
            if torsion_rows.size > 0:
                state_key = self._state_key(flag_id)
                phi0 = np.deg2rad(float(self.torsion_map[state_key]))
                vals_t = []
                for row in torsion_rows:
                    a, b, c, d = self.model.torsion_quads[int(row)]
                    val = _torsion_angle_rad(
                        pos[int(a)],
                        pos[int(b)],
                        pos[int(c)],
                        pos[int(d)],
                    )
                    vals_t.append(abs(_wrap_angle(float(val - phi0))) * RAD_TO_DEG)
                mean_torsion_dev_deg = float(np.mean(vals_t))
            else:
                mean_torsion_dev_deg = float("nan")

            attach = int(self.attach_body_idx[flag_id])
            if attach >= 0:
                hook_vec = pos[int(self.flagella_indices[flag_id][0])] - pos[attach]
            else:
                hook_vec = axis
            hook_norm = max(float(np.linalg.norm(hook_vec)), 1e-18)
            hook_dir = hook_vec / hook_norm
            hook_angle_deg = float(
                np.degrees(
                    np.arccos(float(np.clip(np.dot(body_axis, hook_dir), -1.0, 1.0)))
                )
            )

            flag_all = pos[self.flagella_indices[flag_id]]
            others = [
                self.flagella_indices[idx]
                for idx in range(len(self.flagella_indices))
                if idx != flag_id
            ]
            if others:
                other_points = pos[np.concatenate(others)]
                diff = flag_all[:, None, :] - other_points[None, :, :]
                min_inter_um = float(np.min(np.linalg.norm(diff, axis=2)) * M_TO_UM)
            else:
                min_inter_um = float("inf")

            body_points = pos[self.body_indices]
            diff_body = flag_all[:, None, :] - body_points[None, :, :]
            min_body_um = float(np.min(np.linalg.norm(diff_body, axis=2)) * M_TO_UM)

            if min_inter_um < self.collapse_distance_um:
                step_candidate = True
            if min_inter_um < self._global_min_inter_um:
                self._global_min_inter_um = min_inter_um
                self._flag_id_at_min_distance = flag_id
            if min_body_um < self._global_min_body_um:
                self._global_min_body_um = min_body_um

            row = {
                "step": int(step),
                "t_s": float(t_s),
                "flagellum_id": int(flag_id),
                "mean_radius_um": mean_radius_um,
                "std_radius_um": std_radius_um,
                "mean_phase_diff_deg": mean_phase_diff_deg,
                "std_phase_diff_deg": std_phase_diff_deg,
                "mean_bending_dev_deg": mean_bending_dev_deg,
                "mean_torsion_dev_deg": mean_torsion_dev_deg,
                "hook_angle_deg": hook_angle_deg,
                "min_interflagella_distance_um": min_inter_um,
                "min_body_distance_um": min_body_um,
            }
            self._writer.writerow(row)

        self._csv_fp.flush()

        if step_candidate:
            if self._candidate_consecutive == 0:
                self._candidate_start_step = int(step)
                self._candidate_start_t = float(t_s)
            self._candidate_consecutive += 1
        else:
            self._candidate_start_step = -1
            self._candidate_start_t = -1.0
            self._candidate_consecutive = 0

        if (
            not self._collapse_detected
            and self._candidate_consecutive >= self.collapse_consecutive_steps
        ):
            self._collapse_detected = True
            self._first_collapse_step = self._candidate_start_step
            self._first_collapse_t_s = self._candidate_start_t

    def write_files(self) -> tuple[Path, Path]:
        if not self.enabled:
            return self.diagnostics_path, self.summary_path

        if self._csv_fp is not None:
            self._csv_fp.close()
            self._csv_fp = None
            self._writer = None

        summary = {
            "collapse_detected": bool(self._collapse_detected),
            "first_collapse_step": int(self._first_collapse_step),
            "first_collapse_t_s": float(self._first_collapse_t_s),
            "flagellum_id_at_min_distance": int(self._flag_id_at_min_distance),
            "global_min_interflagella_distance_um": float(self._global_min_inter_um),
            "global_min_body_distance_um": float(self._global_min_body_um),
        }
        if not np.isfinite(summary["global_min_interflagella_distance_um"]):
            summary["global_min_interflagella_distance_um"] = -1.0
        if not np.isfinite(summary["global_min_body_distance_um"]):
            summary["global_min_body_distance_um"] = -1.0
        if not self._collapse_detected:
            summary["first_collapse_step"] = -1
            summary["first_collapse_t_s"] = -1.0

        with self.summary_path.open("w", encoding="utf-8", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=COLLAPSE_SUMMARY_COLUMNS)
            writer.writeheader()
            writer.writerow(summary)

        return self.diagnostics_path, self.summary_path
