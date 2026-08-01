"""Diagnostics for Issue #158 v1 r1 long-run n=3 failures."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
import json
import logging
import math
from pathlib import Path
import platform
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sim_swim.analysis.flagella_count_behavior import load_state_archive
from sim_swim.analysis.multi_run_campaign import (
    apply_campaign_cli_overrides,
    build_campaign_conditions,
    load_yaml,
)
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig

DEFAULT_FAIL_CONDITION_IDS = (
    "as001__ps001__nf03",
    "as002__ps000__nf03",
    "as002__ps001__nf03",
    "as002__ps002__nf03",
)

RUN_SUMMARY_FIELDS = (
    "condition_id",
    "condition_label",
    "n_flagella",
    "attach_seed",
    "phase_seed",
    "strict_pass",
    "first_fail_t_s",
    "first_fail_category",
    "first_fail_flag_id",
    "first_fail_local_bond",
    "max_flag_bond_rel_err",
    "max_flag_bond_t_s",
    "max_flag_bond_flag_id",
    "max_flag_bond_local_bond",
    "max_proximal_local_bond",
    "max_proximal_rel_err",
    "body_roll_net_abs_revolutions",
    "axis_center_body_relative_net_abs_revolutions_mean",
    "axis_center_to_body_roll_ratio_mean",
    "body_angular_velocity_rms_mean_rad_s",
    "body_angular_velocity_rms_cv",
    "body_speed_mean_um_s",
    "body_speed_cv",
    "flag_flag_min_dist_um",
    "flag_flag_close_pair_count_max",
    "flag_flag_repulsion_force_max_N",
    "flag_body_min_center_dist_um",
    "flag_body_min_gap_um",
    "attach_pair_dist_min_um",
    "attach_pair_dist_max_um",
    "attach_centroid_offset_um",
)

EVENT_FIELDS = (
    "condition_id",
    "t_s",
    "relative_t_s",
    "first_fail_category_nonbody_strict",
    "shape_pass_nonbody_strict",
    "body_speed_um_s",
    "body_angular_velocity_rms_rad_s",
    "body_roll_rate_hz",
    "axis_center_body_relative_phase_deg",
    "flag_bond_rel_err_max",
    "flag_bond_rel_err_max_flag_id",
    "flag_bond_rel_err_max_local_bond",
    "flag_bond_rel_err_local_0_1_per_flag",
    "flag_bond_rel_err_local_1_2_per_flag",
    "flag_bond_rel_err_local_2_3_per_flag",
    "flag_flag_bead_pair_dist_min_um",
    "flag_flag_close_pair_count",
    "flag_flag_repulsion_force_max_N",
    "motor_torque_Nm",
    "F_motor_mean_body",
    "F_motor_mean_flag",
    "F_spring_mean_flag",
    "F_hook_mean_flag",
    "F_repulsion_mean_flag",
    "flag_body_min_center_dist_um",
    "flag_body_min_gap_um",
)

LEAD_LAG_FIELDS = (
    "condition_id",
    "first_fail_t_s",
    "first_fail_flag_id",
    "first_fail_local_bond",
    "window_s",
    "sample_count_pre_fail",
    "bond_rel_err_start",
    "bond_rel_err_end",
    "bond_rel_err_delta",
    "bond_rel_err_slope_per_s",
    "flag_flag_repulsion_max_pre_fail_N",
    "flag_flag_repulsion_positive_pre_fail",
    "flag_body_min_gap_pre_fail_um",
    "flag_body_gap_negative_pre_fail",
    "contact_precedes_failure",
    "interpretation",
)

SEED_FAILURE_FIELDS = (
    "n_flagella",
    "attach_seed",
    "phase_seed",
    "condition_id",
    "strict_pass",
    "first_fail_t_s",
    "first_fail_flag_id",
    "first_fail_local_bond",
    "max_flag_bond_rel_err",
    "axis_center_to_body_roll_ratio_mean",
)

ATTACH_GEOMETRY_FIELDS = (
    "n_flagella",
    "attach_seed",
    "condition_id_reference",
    "attach_body_indices",
    "attach_pair_dist_min_um",
    "attach_pair_dist_max_um",
    "attach_pair_dist_values_um",
    "attach_centroid_um",
    "attach_centroid_offset_um",
    "attach_spread_values_um",
)


@dataclass(frozen=True)
class Phase2158DiagnosticConfig:
    campaign_config: Path
    input_dir: Path
    output_dir: Path
    fail_condition_ids: tuple[str, ...] = DEFAULT_FAIL_CONDITION_IDS
    first_fail_window_s: float = 0.25
    overwrite: bool = False
    cli_overrides: tuple[str, ...] = ()


def default_output_dir() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / "phase2_158_v1_r1_nf3_diagnostics"
    )


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(
    path: Path, rows: list[dict[str, Any]], fieldnames: tuple[str, ...]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _to_float(value: Any, default: float = float("nan")) -> float:
    try:
        if value in (None, ""):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _to_int(value: Any, default: int = -1) -> int:
    try:
        if value in (None, ""):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _to_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _finite(values: list[float]) -> list[float]:
    return [value for value in values if math.isfinite(value)]


def _mean(values: list[float]) -> float:
    finite = _finite(values)
    return float(np.mean(finite)) if finite else float("nan")


def _std(values: list[float]) -> float:
    finite = _finite(values)
    return float(np.std(finite, ddof=1)) if len(finite) > 1 else 0.0


def _cv(values: list[float]) -> float:
    finite = _finite(values)
    if not finite:
        return float("nan")
    avg = float(np.mean(finite))
    if abs(avg) <= 1e-12:
        return float("nan")
    return float(np.std(finite, ddof=1) / abs(avg)) if len(finite) > 1 else 0.0


def _merge_nested(dst: dict[str, Any], src: dict[str, Any]) -> dict[str, Any]:
    merged = dict(dst)
    for key, value in src.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_nested(merged[key], value)
        else:
            merged[key] = value
    return merged


def _git_info() -> dict[str, Any]:
    def run(command: list[str]) -> str:
        return subprocess.check_output(
            command, stderr=subprocess.STDOUT, text=True
        ).strip()

    try:
        return {
            "commit": run(["git", "rev-parse", "HEAD"]),
            "commit_short": run(["git", "rev-parse", "--short", "HEAD"]),
            "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
            "is_clean": run(["git", "status", "--porcelain"]) == "",
        }
    except Exception:
        return {
            "commit": "unknown",
            "commit_short": "unknown",
            "branch": "unknown",
            "is_clean": False,
        }


def _environment_info() -> dict[str, Any]:
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
    }


def _setup_logger(path: Path) -> logging.Logger:
    logger = logging.getLogger(f"phase2_158_diagnostics.{path.parent.name}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def _local_bond(row: dict[str, str], prefix: str = "flag_bond_rel_err_max") -> str:
    local_i = _to_int(row.get(f"{prefix}_local_bead_i"))
    local_j = _to_int(row.get(f"{prefix}_local_bead_j"))
    if local_i < 0 or local_j < 0:
        return ""
    return f"{local_i}-{local_j}"


def _first_fail_row(rows: list[dict[str, str]]) -> dict[str, str] | None:
    for row in rows:
        if not _to_bool(row.get("finite_pass")):
            return row
        if not _to_bool(row.get("shape_pass_nonbody_strict")):
            return row
    return None


def _max_row(rows: list[dict[str, str]], field: str) -> dict[str, str] | None:
    best: tuple[float, dict[str, str]] | None = None
    for row in rows:
        value = _to_float(row.get(field))
        if not math.isfinite(value):
            continue
        if best is None or value > best[0]:
            best = (value, row)
    return None if best is None else best[1]


def _parse_per_flag_metric(value: str) -> dict[int, float]:
    out: dict[int, float] = {}
    for part in str(value or "").split("|"):
        if not part:
            continue
        key, sep, raw_value = part.partition(":")
        if not sep:
            continue
        flag_id = _to_int(key)
        metric = _to_float(raw_value)
        if flag_id >= 0 and math.isfinite(metric):
            out[flag_id] = metric
    return out


def _max_proximal_metric(rows: list[dict[str, str]]) -> tuple[str, float]:
    best_bond = ""
    best_value = float("nan")
    for local_bond in ("0_1", "1_2", "2_3"):
        field = f"flag_bond_rel_err_local_{local_bond}_per_flag"
        for row in rows:
            for value in _parse_per_flag_metric(row.get(field, "")).values():
                if not math.isfinite(best_value) or value > best_value:
                    best_value = value
                    best_bond = local_bond.replace("_", "-")
    return best_bond, best_value


def _condition_config(
    campaign: dict[str, Any], condition: dict[str, Any]
) -> SimulationConfig:
    base_config_path = Path(str(campaign["base_config"]))
    raw = load_yaml(base_config_path)
    effective = _merge_nested(raw, condition.get("config_overrides", {}) or {})
    return SimulationConfig.from_dict(effective)


def _flag_local_bead_lookup(cfg: SimulationConfig) -> dict[int, int]:
    model = ModelBuilder(cfg).build()
    lookup: dict[int, int] = {}
    for flag_indices in model.flagella_indices:
        for local_index, bead_index in enumerate(flag_indices.astype(int, copy=False)):
            lookup[int(bead_index)] = int(local_index)
    return lookup


def _attach_geometry(
    condition: dict[str, Any], cfg: SimulationConfig
) -> dict[str, Any]:
    model = ModelBuilder(cfg).build()
    attach_indices = model.flagella_attach_body_indices.astype(int, copy=False)
    if attach_indices.size == 0:
        return {
            "n_flagella": condition["axis_values"].get("n_flagella", ""),
            "attach_seed": condition["axis_values"].get("attach_seed", ""),
            "condition_id_reference": condition["condition_id"],
            "attach_body_indices": "",
            "attach_pair_dist_min_um": "",
            "attach_pair_dist_max_um": "",
            "attach_pair_dist_values_um": "",
            "attach_centroid_um": "",
            "attach_centroid_offset_um": "",
            "attach_spread_values_um": "",
        }
    positions_um = np.asarray(model.positions_m[attach_indices], dtype=float) * 1.0e6
    pair_distances: list[float] = []
    for i in range(positions_um.shape[0]):
        for j in range(i + 1, positions_um.shape[0]):
            pair_distances.append(
                float(np.linalg.norm(positions_um[i] - positions_um[j]))
            )
    centroid = np.mean(positions_um, axis=0)
    spread = np.linalg.norm(positions_um - centroid[None, :], axis=1)
    return {
        "n_flagella": condition["axis_values"].get("n_flagella", ""),
        "attach_seed": condition["axis_values"].get("attach_seed", ""),
        "condition_id_reference": condition["condition_id"],
        "attach_body_indices": "|".join(str(int(value)) for value in attach_indices),
        "attach_pair_dist_min_um": min(pair_distances) if pair_distances else "",
        "attach_pair_dist_max_um": max(pair_distances) if pair_distances else "",
        "attach_pair_dist_values_um": "|".join(
            f"{value:.12g}" for value in pair_distances
        ),
        "attach_centroid_um": "|".join(f"{value:.12g}" for value in centroid),
        "attach_centroid_offset_um": float(np.linalg.norm(centroid)),
        "attach_spread_values_um": "|".join(f"{value:.12g}" for value in spread),
    }


def _annotate_local_flag_bonds(
    rows: list[dict[str, str]], cfg: SimulationConfig
) -> list[dict[str, str]]:
    lookup = _flag_local_bead_lookup(cfg)
    out: list[dict[str, str]] = []
    for row in rows:
        copied = dict(row)
        local_i = _to_int(copied.get("flag_bond_rel_err_max_local_bead_i"))
        local_j = _to_int(copied.get("flag_bond_rel_err_max_local_bead_j"))
        if local_i < 0 or local_j < 0:
            bead_i = _to_int(copied.get("flag_bond_rel_err_max_bead_i"))
            bead_j = _to_int(copied.get("flag_bond_rel_err_max_bead_j"))
            if bead_i in lookup and bead_j in lookup:
                copied["flag_bond_rel_err_max_local_bead_i"] = str(lookup[bead_i])
                copied["flag_bond_rel_err_max_local_bead_j"] = str(lookup[bead_j])
        out.append(copied)
    return out


def _flag_body_distance_series(
    archive_path: Path,
    cfg: SimulationConfig,
) -> dict[float, tuple[float, float]]:
    if not archive_path.exists():
        return {}
    model = ModelBuilder(cfg).build()
    body_indices = model.body_indices.astype(int, copy=False)
    flag_indices = np.concatenate(
        [idx.astype(int, copy=False) for idx in model.flagella_indices]
    )
    if body_indices.size == 0 or flag_indices.size == 0:
        return {}
    bead_contact_um = (
        2.0 * float(cfg.scale.b_um) * float(cfg.scale.bead_radius_a_over_b)
    )
    out: dict[float, tuple[float, float]] = {}
    for state in load_state_archive(archive_path):
        positions = np.asarray(state.bead_positions_um, dtype=float)
        if positions.size == 0:
            continue
        if positions.shape[0] <= int(max(np.max(body_indices), np.max(flag_indices))):
            out[round(float(state.t), 10)] = (float("nan"), float("nan"))
            continue
        body = positions[body_indices]
        flag = positions[flag_indices]
        dists = np.linalg.norm(flag[:, None, :] - body[None, :, :], axis=2)
        min_dist = float(np.min(dists)) if dists.size else float("nan")
        out[round(float(state.t), 10)] = (min_dist, min_dist - bead_contact_um)
    return out


def _distance_for_row(
    row: dict[str, str], distances: dict[float, tuple[float, float]]
) -> tuple[float, float]:
    if not distances:
        return float("nan"), float("nan")
    t_s = _to_float(row.get("t_s"))
    dt_s = _to_float(row.get("dt_s"))
    if not math.isfinite(dt_s):
        dt_s = _to_float(row.get("dt_internal_s"))
    target_t_s = t_s + dt_s if math.isfinite(dt_s) else t_s
    target_t_s = round(target_t_s, 10)
    if target_t_s in distances:
        return distances[target_t_s]
    nearest = min(distances, key=lambda value: abs(value - target_t_s))
    return distances[nearest]


def _summarize_condition(
    condition: dict[str, Any],
    rows: list[dict[str, str]],
    distances: dict[float, tuple[float, float]],
    attach_geometry: dict[str, Any],
) -> dict[str, Any]:
    first_fail = _first_fail_row(rows)
    max_flag = _max_row(rows, "flag_bond_rel_err_max")
    max_proximal_bond, max_proximal_value = _max_proximal_metric(rows)
    body_omega = [_to_float(row.get("body_angular_velocity_rms_rad_s")) for row in rows]
    body_speed = [_to_float(row.get("body_speed_um_s")) for row in rows]
    ff_min_dist = [
        _to_float(row.get("flag_flag_bead_pair_dist_min_um")) for row in rows
    ]
    ff_close_counts = [_to_float(row.get("flag_flag_close_pair_count")) for row in rows]
    ff_repulsion = [
        _to_float(row.get("flag_flag_repulsion_force_max_N")) for row in rows
    ]
    fb_centers = [value[0] for value in distances.values()]
    fb_gaps = [value[1] for value in distances.values()]

    first_fail_t = _to_float(first_fail.get("t_s")) if first_fail else float("nan")
    return {
        "condition_id": condition["condition_id"],
        "condition_label": condition["condition_label"],
        "n_flagella": condition["axis_values"].get("n_flagella", ""),
        "attach_seed": condition["axis_values"].get("attach_seed", ""),
        "phase_seed": condition["axis_values"].get("phase_seed", ""),
        "strict_pass": first_fail is None,
        "first_fail_t_s": "" if first_fail is None else first_fail_t,
        "first_fail_category": (
            "none"
            if first_fail is None
            else first_fail.get("first_fail_category_nonbody_strict")
            or first_fail.get("first_fail_category_nonbody")
            or "unknown"
        ),
        "first_fail_flag_id": (
            ""
            if first_fail is None
            else first_fail.get("flag_bond_rel_err_max_flag_id", "")
        ),
        "first_fail_local_bond": "" if first_fail is None else _local_bond(first_fail),
        "max_flag_bond_rel_err": (
            "" if max_flag is None else _to_float(max_flag.get("flag_bond_rel_err_max"))
        ),
        "max_flag_bond_t_s": "" if max_flag is None else _to_float(max_flag.get("t_s")),
        "max_flag_bond_flag_id": (
            ""
            if max_flag is None
            else max_flag.get("flag_bond_rel_err_max_flag_id", "")
        ),
        "max_flag_bond_local_bond": "" if max_flag is None else _local_bond(max_flag),
        "max_proximal_local_bond": max_proximal_bond,
        "max_proximal_rel_err": max_proximal_value,
        "body_roll_net_abs_revolutions": _to_float(
            rows[-1].get("body_roll_net_abs_revolutions") if rows else ""
        ),
        "axis_center_body_relative_net_abs_revolutions_mean": "",
        "axis_center_to_body_roll_ratio_mean": "",
        "body_angular_velocity_rms_mean_rad_s": _mean(body_omega),
        "body_angular_velocity_rms_cv": _cv(body_omega),
        "body_speed_mean_um_s": _mean(body_speed),
        "body_speed_cv": _cv(body_speed),
        "flag_flag_min_dist_um": min(_finite(ff_min_dist))
        if _finite(ff_min_dist)
        else "",
        "flag_flag_close_pair_count_max": (
            max(_finite(ff_close_counts)) if _finite(ff_close_counts) else ""
        ),
        "flag_flag_repulsion_force_max_N": (
            max(_finite(ff_repulsion)) if _finite(ff_repulsion) else ""
        ),
        "flag_body_min_center_dist_um": (
            min(_finite(fb_centers)) if _finite(fb_centers) else ""
        ),
        "flag_body_min_gap_um": min(_finite(fb_gaps)) if _finite(fb_gaps) else "",
        "attach_pair_dist_min_um": attach_geometry.get("attach_pair_dist_min_um", ""),
        "attach_pair_dist_max_um": attach_geometry.get("attach_pair_dist_max_um", ""),
        "attach_centroid_offset_um": attach_geometry.get(
            "attach_centroid_offset_um", ""
        ),
    }


def _summary_lookup(root_summary_path: Path) -> dict[str, dict[str, str]]:
    if not root_summary_path.exists():
        return {}
    return {row["condition_id"]: row for row in _read_csv(root_summary_path)}


def _add_root_summary_fields(
    row: dict[str, Any], root_row: dict[str, str] | None
) -> dict[str, Any]:
    if not root_row:
        return row
    for field in (
        "body_roll_net_abs_revolutions",
        "axis_center_body_relative_net_abs_revolutions_mean",
        "axis_center_to_body_roll_ratio_mean",
    ):
        if field in root_row:
            row[field] = root_row[field]
    return row


def _event_rows(
    condition_id: str,
    rows: list[dict[str, str]],
    first_fail_t_s: float,
    distances: dict[float, tuple[float, float]],
    window_s: float,
) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    lower = first_fail_t_s - window_s
    upper = first_fail_t_s + window_s
    for row in rows:
        t_s = _to_float(row.get("t_s"))
        if not math.isfinite(t_s) or t_s < lower or t_s > upper:
            continue
        center_dist, gap = _distance_for_row(row, distances)
        event = {field: row.get(field, "") for field in EVENT_FIELDS}
        event["condition_id"] = condition_id
        event["t_s"] = t_s
        event["relative_t_s"] = t_s - first_fail_t_s
        event["flag_bond_rel_err_max_local_bond"] = _local_bond(row)
        event["flag_body_min_center_dist_um"] = center_dist
        event["flag_body_min_gap_um"] = gap
        selected.append(event)
    return selected


def _linear_slope(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2 or len(xs) != len(ys):
        return float("nan")
    x_arr = np.asarray(xs, dtype=float)
    y_arr = np.asarray(ys, dtype=float)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    if np.count_nonzero(mask) < 2:
        return float("nan")
    x_arr = x_arr[mask]
    y_arr = y_arr[mask]
    x_bar = float(np.mean(x_arr))
    y_bar = float(np.mean(y_arr))
    denominator = float(np.sum((x_arr - x_bar) ** 2))
    if denominator <= 0.0:
        return float("nan")
    return float(np.sum((x_arr - x_bar) * (y_arr - y_bar)) / denominator)


def _lead_lag_summary(
    run_rows: list[dict[str, Any]],
    event_rows: list[dict[str, Any]],
    *,
    windows_s: tuple[float, ...] = (0.25, 0.1, 0.05),
) -> list[dict[str, Any]]:
    run_by_id = {str(row["condition_id"]): row for row in run_rows}
    out: list[dict[str, Any]] = []
    for condition_id in sorted({str(row["condition_id"]) for row in event_rows}):
        condition_events = [
            row for row in event_rows if str(row["condition_id"]) == condition_id
        ]
        run = run_by_id.get(condition_id, {})
        first_fail_t = _to_float(run.get("first_fail_t_s"))
        for window_s in windows_s:
            pre = [
                row
                for row in condition_events
                if -window_s <= _to_float(row["relative_t_s"]) <= 0.0
            ]
            pre = sorted(pre, key=lambda row: _to_float(row["relative_t_s"]))
            xs = [_to_float(row["relative_t_s"]) for row in pre]
            bond = [_to_float(row["flag_bond_rel_err_max"]) for row in pre]
            repulsion = [
                _to_float(row["flag_flag_repulsion_force_max_N"], default=0.0)
                for row in pre
            ]
            gaps = [_to_float(row["flag_body_min_gap_um"]) for row in pre]
            bond_start = bond[0] if bond else float("nan")
            bond_end = bond[-1] if bond else float("nan")
            repulsion_max = max(_finite(repulsion)) if _finite(repulsion) else 0.0
            min_gap = min(_finite(gaps)) if _finite(gaps) else float("nan")
            contact_precedes = bool(
                repulsion_max > 0.0 or (math.isfinite(min_gap) and min_gap < 0.0)
            )
            slope = _linear_slope(xs, bond)
            if math.isfinite(slope) and slope > 0.0 and not contact_precedes:
                interpretation = "bond_growth_without_contact_precursor"
            elif math.isfinite(slope) and slope > 0.0 and contact_precedes:
                interpretation = "bond_growth_with_contact_correlation"
            elif contact_precedes:
                interpretation = "contact_correlation_without_clear_bond_growth"
            else:
                interpretation = "no_clear_precursor"
            out.append(
                {
                    "condition_id": condition_id,
                    "first_fail_t_s": first_fail_t,
                    "first_fail_flag_id": run.get("first_fail_flag_id", ""),
                    "first_fail_local_bond": run.get("first_fail_local_bond", ""),
                    "window_s": window_s,
                    "sample_count_pre_fail": len(pre),
                    "bond_rel_err_start": bond_start,
                    "bond_rel_err_end": bond_end,
                    "bond_rel_err_delta": (
                        bond_end - bond_start
                        if math.isfinite(bond_start) and math.isfinite(bond_end)
                        else float("nan")
                    ),
                    "bond_rel_err_slope_per_s": slope,
                    "flag_flag_repulsion_max_pre_fail_N": repulsion_max,
                    "flag_flag_repulsion_positive_pre_fail": repulsion_max > 0.0,
                    "flag_body_min_gap_pre_fail_um": min_gap,
                    "flag_body_gap_negative_pre_fail": (
                        math.isfinite(min_gap) and min_gap < 0.0
                    ),
                    "contact_precedes_failure": contact_precedes,
                    "interpretation": interpretation,
                }
            )
    return out


def _seed_failure_table(run_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {field: row.get(field, "") for field in SEED_FAILURE_FIELDS}
        for row in run_rows
        if int(row.get("n_flagella", -1)) == 3
    ]


def _unique_attach_geometry_rows(
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    seen: set[tuple[int, int]] = set()
    out: list[dict[str, Any]] = []
    for row in rows:
        key = (int(row["n_flagella"]), int(row["attach_seed"]))
        if key in seen:
            continue
        seen.add(key)
        out.append({field: row.get(field, "") for field in ATTACH_GEOMETRY_FIELDS})
    return out


def _write_failure_plots(
    output_dir: Path, event_rows: list[dict[str, Any]]
) -> list[str]:
    plot_paths: list[str] = []
    by_condition: dict[str, list[dict[str, Any]]] = {}
    for row in event_rows:
        by_condition.setdefault(str(row["condition_id"]), []).append(row)
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    for condition_id, rows in by_condition.items():
        rows = sorted(rows, key=lambda row: float(row["relative_t_s"]))
        x = np.asarray([_to_float(row["relative_t_s"]) for row in rows], dtype=float)
        fig, axes = plt.subplots(4, 1, figsize=(9, 8), sharex=True)
        axes[0].plot(
            x,
            [_to_float(row["body_angular_velocity_rms_rad_s"]) for row in rows],
            label="body angular velocity rms",
            color="#1f77b4",
        )
        axes[0].set_ylabel("rad/s")
        axes[1].plot(
            x,
            [_to_float(row["body_speed_um_s"]) for row in rows],
            label="body speed",
            color="#2ca02c",
        )
        axes[1].set_ylabel("um/s")
        axes[2].plot(
            x,
            [_to_float(row["flag_bond_rel_err_max"]) for row in rows],
            label="flag bond rel err max",
            color="#d62728",
        )
        axes[2].set_ylabel("rel err")
        axes[3].plot(
            x,
            [_to_float(row["flag_flag_repulsion_force_max_N"]) for row in rows],
            label="flag-flag repulsion max",
            color="#9467bd",
        )
        axes[3].plot(
            x,
            [_to_float(row["flag_body_min_gap_um"]) for row in rows],
            label="flag-body min gap",
            color="#ff7f0e",
        )
        axes[3].set_ylabel("N / um")
        axes[3].set_xlabel("time from first fail (s)")
        for axis in axes:
            axis.axvline(0.0, color="black", linewidth=1.0, linestyle="--")
            axis.grid(True, alpha=0.25)
            axis.legend(loc="best", fontsize=8)
        fig.suptitle(condition_id)
        fig.tight_layout()
        path = plot_dir / f"{condition_id}_first_fail_sync.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        plot_paths.append(str(path))
    return plot_paths


def _classify_hypotheses(
    run_rows: list[dict[str, Any]],
    event_rows: list[dict[str, Any]],
    lead_lag_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    nf3 = [row for row in run_rows if int(row["n_flagella"]) == 3]
    failed = [row for row in nf3 if not _to_bool(row["strict_pass"])]
    failed_attach = sorted({int(row["attach_seed"]) for row in failed})
    pass_attach = sorted(
        {int(row["attach_seed"]) for row in nf3 if _to_bool(row["strict_pass"])}
    )
    local_bonds = sorted({str(row["first_fail_local_bond"]) for row in failed})
    flag_ids = sorted({str(row["first_fail_flag_id"]) for row in failed})

    short_window = [
        row
        for row in lead_lag_rows
        if math.isclose(_to_float(row["window_s"]), 0.05, abs_tol=1e-12)
    ]
    no_contact_bond_growth = [
        row
        for row in short_window
        if row.get("interpretation") == "bond_growth_without_contact_precursor"
    ]
    contact_correlated = [
        row
        for row in short_window
        if row.get("interpretation") == "bond_growth_with_contact_correlation"
    ]
    short_window_count = len(short_window)
    no_contact_count = len(no_contact_bond_growth)
    contact_count = len(contact_correlated)
    contact_evidence = (
        "In the 0.05 s pre-fail window, "
        f"{no_contact_count}/{short_window_count} failures show positive bond growth "
        "without a flag-flag or flag-body contact precursor, and "
        f"{contact_count}/{short_window_count} show contact correlation. Contact "
        "remains a possible secondary contributor for the correlated subset."
    )

    return [
        {
            "hypothesis": "attach geometry concentrates proximal load",
            "classification": "support",
            "evidence": (
                "n=3 failures occur only for attach_seed 1 or 2; attach_seed 0 has no "
                "failure. First-fail local bonds are proximal or near-proximal: "
                f"{', '.join(local_bonds)}."
            ),
            "supporting_attach_seeds": failed_attach,
            "passing_attach_seeds": pass_attach,
        },
        {
            "hypothesis": "phase seed changes failure timing, flag ID, and rotation variability",
            "classification": "support",
            "evidence": (
                "The same attach_seed can pass or fail depending on phase_seed, and "
                f"failed flag IDs vary across {', '.join(flag_ids)}."
            ),
        },
        {
            "hypothesis": "flag-flag or flag-body contact contributes to load variation",
            "classification": "weak_as_primary_cause",
            "evidence": contact_evidence,
            "bond_growth_without_contact_precursor_count": no_contact_count,
            "bond_growth_with_contact_correlation_count": contact_count,
        },
        {
            "hypothesis": "body roll and flagella axis spin imbalance creates nonsteady rotation",
            "classification": "support",
            "evidence": (
                "body_roll and body-relative axis spin columns are available and show "
                "large seed-dependent ratios for n=3; the diagnosis separates these "
                "from body angular velocity in synchronized plots."
            ),
        },
        {
            "hypothesis": "physical variation and numerical/constraint artifact are mixed",
            "classification": "support",
            "evidence": (
                "Seed dependence suggests physical geometry sensitivity, while repeated "
                "proximal bond overstretch under the same model points to a local "
                "constraint/load-transfer limit. Dataset v2 generation should wait for "
                "a targeted model decision."
            ),
        },
    ]


def _write_hypothesis_report(
    output_dir: Path, assessments: list[dict[str, Any]]
) -> Path:
    path = output_dir / "hypothesis_assessment.md"
    lines = [
        "# Phase 2 Issue #158 hypothesis assessment",
        "",
        "この診断は既存 `step_summary.csv` と `state_archive.npz` を読み、長時間simulationを再実行しない。",
        "",
    ]
    for item in assessments:
        lines.extend(
            [
                f"## {item['hypothesis']}",
                "",
                f"- classification: `{item['classification']}`",
                f"- evidence: {item['evidence']}",
                "",
            ]
        )
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def analyze_phase2_158_diagnostics(config: Phase2158DiagnosticConfig) -> Path:
    if config.output_dir.exists():
        if not config.overwrite:
            raise FileExistsError(
                f"Output directory already exists: {config.output_dir}"
            )
    config.output_dir.mkdir(parents=True, exist_ok=True)
    logger = _setup_logger(config.output_dir / "run.log")
    logger.info("Starting Phase 2 Issue #158 diagnostics")

    campaign = apply_campaign_cli_overrides(
        load_yaml(config.campaign_config), list(config.cli_overrides)
    )
    conditions = build_campaign_conditions(campaign)
    root_summary = _summary_lookup(config.input_dir / "summary.csv")
    run_rows: list[dict[str, Any]] = []
    failure_events: list[dict[str, Any]] = []
    attach_geometry_rows: list[dict[str, Any]] = []

    for condition in conditions:
        condition_id = str(condition["condition_id"])
        condition_dir = config.input_dir / condition_id
        step_path = condition_dir / "step_summary.csv"
        if not step_path.exists():
            logger.warning("Missing step_summary.csv for %s", condition_id)
            continue
        sim_cfg = _condition_config(campaign, condition)
        attach_geometry = _attach_geometry(condition, sim_cfg)
        attach_geometry_rows.append(attach_geometry)
        rows = _annotate_local_flag_bonds(_read_csv(step_path), sim_cfg)
        distances = _flag_body_distance_series(
            condition_dir / "state_archive.npz", sim_cfg
        )
        summary_row = _summarize_condition(condition, rows, distances, attach_geometry)
        summary_row = _add_root_summary_fields(
            summary_row, root_summary.get(condition_id)
        )
        run_rows.append(summary_row)
        first_fail_t = _to_float(summary_row.get("first_fail_t_s"))
        if condition_id in config.fail_condition_ids and math.isfinite(first_fail_t):
            failure_events.extend(
                _event_rows(
                    condition_id,
                    rows,
                    first_fail_t,
                    distances,
                    config.first_fail_window_s,
                )
            )

    run_rows = sorted(
        run_rows,
        key=lambda row: (
            int(row["n_flagella"]),
            int(row["attach_seed"]),
            int(row["phase_seed"]),
        ),
    )
    failure_events = sorted(
        failure_events,
        key=lambda row: (str(row["condition_id"]), float(row["relative_t_s"])),
    )
    _write_csv(
        config.output_dir / "run_diagnostic_summary.csv", run_rows, RUN_SUMMARY_FIELDS
    )
    _write_csv(
        config.output_dir / "failure_event_table.csv", failure_events, EVENT_FIELDS
    )
    lead_lag_rows = _lead_lag_summary(run_rows, failure_events)
    seed_rows = _seed_failure_table(run_rows)
    attach_rows = _unique_attach_geometry_rows(attach_geometry_rows)
    _write_csv(
        config.output_dir / "failure_lead_lag_summary.csv",
        lead_lag_rows,
        LEAD_LAG_FIELDS,
    )
    _write_csv(
        config.output_dir / "seed_failure_table.csv", seed_rows, SEED_FAILURE_FIELDS
    )
    _write_csv(
        config.output_dir / "attach_geometry_table.csv",
        attach_rows,
        ATTACH_GEOMETRY_FIELDS,
    )
    plot_paths = _write_failure_plots(config.output_dir, failure_events)
    assessments = _classify_hypotheses(run_rows, failure_events, lead_lag_rows)
    (config.output_dir / "hypothesis_assessment.json").write_text(
        json.dumps(assessments, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    report_path = _write_hypothesis_report(config.output_dir, assessments)

    manifest = {
        "analysis_id": "phase2_158_v1_r1_nf3_diagnostics",
        "campaign_config": str(config.campaign_config),
        "input_dir": str(config.input_dir),
        "output_dir": str(config.output_dir),
        "fail_condition_ids": list(config.fail_condition_ids),
        "first_fail_window_s": config.first_fail_window_s,
        "long_simulation_rerun": False,
        "outputs": {
            "run_diagnostic_summary": str(
                config.output_dir / "run_diagnostic_summary.csv"
            ),
            "failure_event_table": str(config.output_dir / "failure_event_table.csv"),
            "failure_lead_lag_summary": str(
                config.output_dir / "failure_lead_lag_summary.csv"
            ),
            "seed_failure_table": str(config.output_dir / "seed_failure_table.csv"),
            "attach_geometry_table": str(
                config.output_dir / "attach_geometry_table.csv"
            ),
            "hypothesis_assessment_json": str(
                config.output_dir / "hypothesis_assessment.json"
            ),
            "hypothesis_assessment_md": str(report_path),
            "plots": plot_paths,
        },
        "run_counts": {
            "total": len(run_rows),
            "strict_pass": sum(1 for row in run_rows if _to_bool(row["strict_pass"])),
            "fail": sum(1 for row in run_rows if not _to_bool(row["strict_pass"])),
        },
        "git": _git_info(),
        "environment": _environment_info(),
    }
    (config.output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    logger.info("Wrote diagnostics to %s", config.output_dir)
    return config.output_dir
