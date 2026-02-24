"""シミュレーション診断CSVの生成。"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from sim_swim.dynamics.engine import StepDiagnostics
from sim_swim.model.types import SimModel
from sim_swim.sim.params import SimulationConfig

M_TO_UM = 1.0e6

STEP_SUMMARY_COLUMNS = [
    "step",
    "t_star",
    "dt_star",
    "pos_all_finite",
    "mean_disp_um",
    "bond_len_mean_body_body_um",
    "bond_len_mean_flag_intra_um",
    "bond_len_mean_body_flag_um",
    "F_total_mean_body",
    "F_total_mean_flag",
    "F_total_mean_all",
    "brownian_enabled",
    "brownian_disp_mean_um",
]

STEP_SUMMARY_FULL_COLUMNS = [
    "step",
    "t_star",
    "dt_star",
    "t_s",
    "dt_s",
    "pos_all_finite",
    "any_nan",
    "any_inf",
    "mean_disp_um",
    "max_disp_um",
    "bond_len_mean_body_body_um",
    "bond_len_mean_flag_intra_um",
    "bond_len_mean_body_flag_um",
    "F_total_mean_body",
    "F_total_mean_flag",
    "F_total_mean_all",
    "F_motor_mean_body",
    "F_motor_mean_flag",
    "F_spring_mean_body",
    "F_spring_mean_flag",
    "F_bend_mean_body",
    "F_bend_mean_flag",
    "F_torsion_mean_body",
    "F_torsion_mean_flag",
    "F_repulsion_mean_body",
    "F_repulsion_mean_flag",
    "F_hook_mean_body",
    "F_hook_mean_flag",
    "motor_degenerate_axis_count",
    "motor_split_rank_deficient_count",
    "motor_bond_length_clipped_count",
    "brownian_enabled",
    "brownian_disp_mean_um",
]


def _mean_norm(forces: np.ndarray, mask: np.ndarray) -> float:
    if forces.size == 0:
        return 0.0
    selected = forces[mask]
    if selected.size == 0:
        return 0.0
    return float(np.mean(np.linalg.norm(selected, axis=1)))


def _mean_pair_distance_um(
    positions_m: np.ndarray, spring_pairs: np.ndarray, pair_rows: np.ndarray
) -> float:
    if pair_rows.size == 0:
        return float("nan")
    pairs = spring_pairs[pair_rows]
    diff = positions_m[pairs[:, 0]] - positions_m[pairs[:, 1]]
    dist_um = np.linalg.norm(diff, axis=1) * M_TO_UM
    return float(np.mean(dist_um))


@dataclass
class StepSummaryRecorder:
    """ステップごとの診断指標をCSVへ保存する。"""

    model: SimModel
    cfg: SimulationConfig
    out_dir: Path
    _rows_min: list[dict[str, float | int | bool]] = field(default_factory=list)
    _rows_full: list[dict[str, float | int | bool]] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.step_summary_path = self.out_dir / "step_summary.csv"
        self.step_summary_full_path = self.out_dir / "step_summary_full.csv"

        self.body_mask = self.model.bead_is_body.astype(bool)
        self.flag_mask = ~self.body_mask
        self.spring_pairs = self.model.spring_pairs.astype(int, copy=False)

        if self.spring_pairs.size == 0:
            self.body_body_rows = np.zeros((0,), dtype=int)
            self.flag_intra_rows = np.zeros((0,), dtype=int)
            self.body_flag_rows = np.zeros((0,), dtype=int)
            return

        i = self.spring_pairs[:, 0]
        j = self.spring_pairs[:, 1]
        fi = self.model.bead_flag_ids[i]
        fj = self.model.bead_flag_ids[j]
        bi = self.body_mask[i]
        bj = self.body_mask[j]

        self.body_body_rows = np.where(bi & bj)[0]
        self.body_flag_rows = np.where(np.logical_xor(bi, bj))[0]
        self.flag_intra_rows = np.where((~bi) & (~bj) & (fi == fj) & (fi >= 0))[0]

    def record(self, step: int, t_star: float, diag: StepDiagnostics) -> None:
        pos_after = diag.positions_after_m
        disp_um = np.linalg.norm(pos_after - diag.positions_before_m, axis=1) * M_TO_UM
        any_nan = bool(np.isnan(pos_after).any())
        any_inf = bool(np.isinf(pos_after).any())
        pos_all_finite = bool(np.isfinite(pos_after).all())

        brownian_disp_mean_um = float("nan")
        if diag.brownian_enabled:
            brownian_disp = np.linalg.norm(diag.brownian_disp_m, axis=1) * M_TO_UM
            brownian_disp_mean_um = float(np.mean(brownian_disp))

        row_full: dict[str, float | int | bool] = {
            "step": int(step),
            "t_star": float(t_star),
            "dt_star": float(diag.dt_star),
            "t_s": float(t_star * self.cfg.tau_s),
            "dt_s": float(diag.dt_s),
            "pos_all_finite": pos_all_finite,
            "any_nan": any_nan,
            "any_inf": any_inf,
            "mean_disp_um": float(np.mean(disp_um)),
            "max_disp_um": float(np.max(disp_um)),
            "bond_len_mean_body_body_um": _mean_pair_distance_um(
                pos_after, self.spring_pairs, self.body_body_rows
            ),
            "bond_len_mean_flag_intra_um": _mean_pair_distance_um(
                pos_after, self.spring_pairs, self.flag_intra_rows
            ),
            "bond_len_mean_body_flag_um": _mean_pair_distance_um(
                pos_after, self.spring_pairs, self.body_flag_rows
            ),
            "F_total_mean_body": _mean_norm(diag.total_forces, self.body_mask),
            "F_total_mean_flag": _mean_norm(diag.total_forces, self.flag_mask),
            "F_total_mean_all": float(
                np.mean(np.linalg.norm(diag.total_forces, axis=1))
            ),
            "F_motor_mean_body": _mean_norm(diag.motor_forces, self.body_mask),
            "F_motor_mean_flag": _mean_norm(diag.motor_forces, self.flag_mask),
            "F_spring_mean_body": _mean_norm(diag.spring_forces, self.body_mask),
            "F_spring_mean_flag": _mean_norm(diag.spring_forces, self.flag_mask),
            "F_bend_mean_body": _mean_norm(diag.bend_forces, self.body_mask),
            "F_bend_mean_flag": _mean_norm(diag.bend_forces, self.flag_mask),
            "F_torsion_mean_body": _mean_norm(diag.torsion_forces, self.body_mask),
            "F_torsion_mean_flag": _mean_norm(diag.torsion_forces, self.flag_mask),
            "F_repulsion_mean_body": _mean_norm(diag.repulsion_forces, self.body_mask),
            "F_repulsion_mean_flag": _mean_norm(diag.repulsion_forces, self.flag_mask),
            "F_hook_mean_body": _mean_norm(diag.hook_forces, self.body_mask),
            "F_hook_mean_flag": _mean_norm(diag.hook_forces, self.flag_mask),
            "motor_degenerate_axis_count": int(diag.motor_degenerate_axis_count),
            "motor_split_rank_deficient_count": int(
                diag.motor_split_rank_deficient_count
            ),
            "motor_bond_length_clipped_count": int(
                diag.motor_bond_length_clipped_count
            ),
            "brownian_enabled": bool(diag.brownian_enabled),
            "brownian_disp_mean_um": brownian_disp_mean_um,
        }

        row_min = {k: row_full[k] for k in STEP_SUMMARY_COLUMNS}
        self._rows_min.append(row_min)
        self._rows_full.append(row_full)

    def write_csv(self) -> tuple[Path, Path]:
        with self.step_summary_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=STEP_SUMMARY_COLUMNS)
            writer.writeheader()
            writer.writerows(self._rows_min)

        with self.step_summary_full_path.open("w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=STEP_SUMMARY_FULL_COLUMNS)
            writer.writeheader()
            writer.writerows(self._rows_full)

        return self.step_summary_path, self.step_summary_full_path
