"""Paired ``diffusive`` versus ``uniform`` torque-profile analysis for #203."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

import matplotlib.pyplot as plt
import numpy as np
import yaml

from sim_swim.analysis.flagella_count_behavior import load_state_archive
from sim_swim.analysis.motion_feature_study import _body_axes
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


@dataclass(frozen=True)
class ComparisonConfig:
    uniform_run_dir: Path
    diffusive_run_dir: Path
    output_dir: Path
    allowed_n_flagella: tuple[int, ...] = (1, 2, 3)
    overwrite: bool = False


def load_config(path: Path) -> ComparisonConfig:
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if type(raw.get("overwrite", False)) is not bool:
        raise ValueError("overwrite must be a boolean")
    return ComparisonConfig(
        uniform_run_dir=Path(str(raw["uniform_run_dir"])),
        diffusive_run_dir=Path(str(raw["diffusive_run_dir"])),
        output_dir=Path(str(raw["output_dir"])),
        allowed_n_flagella=tuple(
            int(v) for v in raw.get("allowed_n_flagella", [1, 2, 3])
        ),
        overwrite=bool(raw.get("overwrite", False)),
    )


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["condition_id"]: row for row in csv.DictReader(handle)}


def _condition_dir(root: Path, record: dict[str, Any]) -> Path:
    value = Path(str(record["output_dir"]))
    return value if value.is_dir() else root / value.name


def _finite(value: str | float | None) -> float:
    try:
        parsed = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")
    return parsed if math.isfinite(parsed) else float("nan")


def _angle(a: np.ndarray, b: np.ndarray) -> float:
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        return float("nan")
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na <= 1e-12 or nb <= 1e-12:
        return float("nan")
    return math.degrees(math.acos(float(np.clip(abs(np.dot(a, b)) / (na * nb), -1, 1))))


def _axis_motion(axis: np.ndarray, times: np.ndarray) -> float:
    if len(axis) < 2:
        return float("nan")
    values = [_angle(axis[i - 1], axis[i]) for i in range(1, len(axis))]
    dt = np.diff(times)
    valid = np.asarray(values, dtype=float) / dt
    return float(np.nanmean(valid)) if np.isfinite(valid).any() else float("nan")


def _flag_axes(
    states: list[Any], flagella_indices: list[np.ndarray]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    axes: list[list[np.ndarray]] = []
    for state in states:
        one: list[np.ndarray] = []
        for indexes in flagella_indices:
            points = np.asarray(state.bead_positions_um, dtype=float)[
                np.asarray(indexes)[1:]
            ]
            centered = points - np.mean(points, axis=0)
            _, _, vh = np.linalg.svd(centered, full_matrices=False)
            vector = vh[0]
            if np.dot(vector, points[-1] - points[0]) < 0:
                vector *= -1
            one.append(vector / np.linalg.norm(vector))
        axes.append(one)
    array = np.asarray(axes, dtype=float)
    mean = np.full((len(states), 3), np.nan)
    alignment = np.full(len(states), np.nan)
    spread = np.full(len(states), np.nan)
    for index, vectors in enumerate(array):
        ref = vectors[0]
        aligned = np.asarray([v if np.dot(v, ref) >= 0 else -v for v in vectors])
        candidate = aligned.mean(axis=0)
        norm = np.linalg.norm(candidate)
        if norm > 1e-12:
            mean[index] = candidate / norm
        if len(vectors) > 1 and norm > 1e-12:
            alignment[index] = norm
            spread[index] = float(np.nanmean([_angle(v, mean[index]) for v in aligned]))
    return mean, alignment, spread


def _strict(summary: dict[str, Any]) -> tuple[bool, str]:
    execution = summary.get("execution", {})
    gates = summary.get("gates", {})
    if execution.get("status") != "completed":
        return False, "execution_not_completed"
    for name in ("finite", "shape_nonbody"):
        gate = gates.get(name, {})
        if (
            gate.get("status") != "available"
            or gate.get("any_fail")
            or not gate.get("final_pass")
        ):
            return False, f"{name}_qc_failed"
    body = gates.get("shape_body", {})
    if body.get("status") == "available" and body.get("any_fail"):
        return False, "shape_body_qc_failed"
    return True, "pass"


def _metrics(root: Path, record: dict[str, Any], row: dict[str, str]) -> dict[str, Any]:
    condition_dir = _condition_dir(root, record)
    summary_path = condition_dir / "run_summary.json"
    archive_path = condition_dir / "state_archive.npz"
    if not summary_path.is_file() or not archive_path.is_file():
        return {"availability": "missing_artifact"}
    summary = _read_json(summary_path)
    strict_pass, reason = _strict(summary)
    result: dict[str, Any] = {
        "availability": "available",
        "strict_pass": strict_pass,
        "strict_reason": reason,
        "first_fail_t_s": _finite(row.get("first_fail_t_s")),
        "first_fail_category": row.get("first_fail_category_nonbody", ""),
        "hook_len_rel_err_max": _finite(row.get("hook_len_rel_err_max")),
        "flag_bond_rel_err_max": _finite(row.get("flag_bond_rel_err_max")),
    }
    states = load_state_archive(archive_path)
    if len(states) < 2:
        result["availability"] = "insufficient_archive"
        return result
    t = np.asarray([state.t for state in states], dtype=float)
    position = np.asarray([state.position_um for state in states], dtype=float)
    displacement = float(np.linalg.norm(position[-1] - position[0]))
    duration = float(t[-1] - t[0])
    body = _body_axes(np.asarray([state.quaternion for state in states], dtype=float))
    cfg = SimulationConfig.from_dict(
        yaml.safe_load(
            Path(str(record.get("base_config", "conf/sim_swim_2010.yaml"))).read_text(
                encoding="utf-8"
            )
        )
    ).with_overrides(record.get("config_overrides", {}))
    mean_flag, alignment, spread = _flag_axes(
        states, Simulator(cfg).rig.flagella_indices
    )
    relation = np.asarray([_angle(a, b) for a, b in zip(body, mean_flag)])
    result.update(
        {
            "body_displacement_um": displacement,
            "body_mean_speed_um_s": displacement / duration
            if duration > 0
            else float("nan"),
            "body_axis_motion_rad_s": _axis_motion(body, t) * math.pi / 180.0,
            "mean_flagella_axis_motion_rad_s": _axis_motion(mean_flag, t)
            * math.pi
            / 180.0,
            "body_flagella_axis_angle_deg": float(np.nanmean(relation)),
            "flagella_axis_alignment": float(np.nanmean(alignment)),
            "flagella_axis_spread_deg": float(np.nanmean(spread)),
            "bundle_like_fraction": float(
                np.nanmean((alignment >= math.cos(math.radians(15))).astype(float))
            )
            if np.isfinite(alignment).any()
            else float("nan"),
        }
    )
    return result


METRICS = (
    "body_displacement_um",
    "body_mean_speed_um_s",
    "body_axis_motion_rad_s",
    "mean_flagella_axis_motion_rad_s",
    "body_flagella_axis_angle_deg",
    "flagella_axis_alignment",
    "flagella_axis_spread_deg",
    "bundle_like_fraction",
)


def _write(path: Path, rows: list[dict[str, Any]]) -> None:
    keys = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def analyze(cfg: ComparisonConfig) -> Path:
    if cfg.output_dir.exists():
        if not cfg.overwrite:
            raise FileExistsError(cfg.output_dir)
        import shutil

        shutil.rmtree(cfg.output_dir)
    cfg.output_dir.mkdir(parents=True)
    profiles = {"uniform": cfg.uniform_run_dir, "diffusive": cfg.diffusive_run_dir}
    manifests = {
        name: _read_json(path / "run_manifest.json") for name, path in profiles.items()
    }
    summaries = {name: _rows(path / "summary.csv") for name, path in profiles.items()}
    records = {
        name: {r["condition_id"]: r for r in manifest["conditions"]}
        for name, manifest in manifests.items()
    }
    paired: list[dict[str, Any]] = []
    for condition_id in sorted(set(records["uniform"]) | set(records["diffusive"])):
        source = records["uniform"].get(condition_id) or records["diffusive"].get(
            condition_id
        )
        values = source.get("axis_values", {}) if source else {}
        n = int(values.get("n_flagella", -1))
        if n not in cfg.allowed_n_flagella:
            continue
        result: dict[str, Any] = {
            "condition_id": condition_id,
            "n_flagella": n,
            "attach_seed": values.get("attach_seed"),
            "phase_seed": values.get("phase_seed"),
        }
        per: dict[str, dict[str, Any]] = {}
        for profile in profiles:
            record, row = (
                records[profile].get(condition_id),
                summaries[profile].get(condition_id),
            )
            per[profile] = (
                {"availability": "missing_condition"}
                if record is None or row is None
                else _metrics(profiles[profile], record, row)
            )
            for key, value in per[profile].items():
                result[f"{profile}_{key}"] = value
        result["pair_status"] = (
            "paired"
            if all(per[p].get("availability") == "available" for p in profiles)
            else "unpaired"
        )
        result["paired_strict_pass"] = result["pair_status"] == "paired" and all(
            bool(per[p].get("strict_pass")) for p in profiles
        )
        for metric in METRICS:
            result[f"delta_uniform_minus_diffusive_{metric}"] = _finite(
                per["uniform"].get(metric)
            ) - _finite(per["diffusive"].get(metric))
        paired.append(result)
    _write(cfg.output_dir / "paired_conditions.csv", paired)
    usable = [row for row in paired if row["paired_strict_pass"]]
    aggregate: list[dict[str, Any]] = []
    for n in cfg.allowed_n_flagella:
        subset = [row for row in usable if row["n_flagella"] == n]
        for metric in METRICS:
            values = np.asarray(
                [
                    _finite(row[f"delta_uniform_minus_diffusive_{metric}"])
                    for row in subset
                ]
            )
            aggregate.append(
                {
                    "n_flagella": n,
                    "metric": metric,
                    "paired_strict_run_count": len(subset),
                    "delta_mean": float(np.nanmean(values))
                    if np.isfinite(values).any()
                    else float("nan"),
                    "delta_std": float(np.nanstd(values))
                    if np.isfinite(values).any()
                    else float("nan"),
                }
            )
    _write(cfg.output_dir / "paired_aggregate.csv", aggregate)
    figure, axes = plt.subplots(2, 4, figsize=(15, 7))
    for axis, metric in zip(axes.flat, METRICS):
        for n in cfg.allowed_n_flagella:
            values = [
                _finite(row[f"delta_uniform_minus_diffusive_{metric}"])
                for row in usable
                if row["n_flagella"] == n
            ]
            axis.scatter([n] * len(values), values, alpha=0.7)
        axis.axhline(0, color="black", linewidth=0.8)
        axis.set(title=metric, xlabel="n_flagella", ylabel="uniform - diffusive")
    figure.tight_layout()
    figure.savefig(cfg.output_dir / "paired_delta.png", dpi=150)
    plt.close(figure)
    manifest = {
        "kind": "phase2_issue203_paired_comparison",
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "uniform_run_dir": str(cfg.uniform_run_dir),
        "diffusive_run_dir": str(cfg.diffusive_run_dir),
        "condition_count": len(paired),
        "paired_strict_pass_count": len(usable),
        "nonpass_or_unpaired_count": len(paired) - len(usable),
        "outputs": [
            "paired_conditions.csv",
            "paired_aggregate.csv",
            "paired_delta.png",
        ],
    }
    (cfg.output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return cfg.output_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--uniform-run-dir", type=Path)
    parser.add_argument("--diffusive-run-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    config = ComparisonConfig(
        uniform_run_dir=args.uniform_run_dir or config.uniform_run_dir,
        diffusive_run_dir=args.diffusive_run_dir or config.diffusive_run_dir,
        output_dir=args.output_dir or config.output_dir,
        allowed_n_flagella=config.allowed_n_flagella,
        overwrite=args.overwrite or config.overwrite,
    )
    print(analyze(config))
