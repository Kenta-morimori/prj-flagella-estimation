"""Compact contact/stability analysis for a profile-by-``dt_star`` matrix."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
import json
import math
from pathlib import Path
import re
from typing import Any
from zoneinfo import ZoneInfo

import matplotlib.pyplot as plt
import numpy as np
import yaml

from sim_swim.analysis.issue203_torque_profile_comparison import _strict


_DIAGNOSTIC_SEED_CASE = re.compile(
    r"^nf(?P<n>\d{2})__as(?P<attach>\d{3})__ps(?P<phase>\d{3})$"
)


def _reference_condition_id(seed_case: str) -> str:
    """Map a new n-first diagnostic ID to the historical reference ID."""

    match = _DIAGNOSTIC_SEED_CASE.fullmatch(seed_case)
    if match is None:
        raise ValueError(f"invalid n-first diagnostic seed case: {seed_case}")
    return "as{attach}__ps{phase}__nf{n}".format(**match.groupdict())


@dataclass(frozen=True)
class ContactConfig:
    uniform_reference_run_dir: Path
    diffusive_reference_run_dir: Path
    diagnostic_run_dir: Path
    output_dir: Path
    seed_cases: tuple[str, ...]
    profiles: tuple[str, ...]
    dt_stars: tuple[float, ...]
    overwrite: bool = False


def load_config(path: Path) -> ContactConfig:
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if type(raw.get("overwrite", False)) is not bool:
        raise ValueError("overwrite must be a boolean")
    return ContactConfig(
        uniform_reference_run_dir=Path(str(raw["uniform_reference_run_dir"])),
        diffusive_reference_run_dir=Path(str(raw["diffusive_reference_run_dir"])),
        diagnostic_run_dir=Path(str(raw["diagnostic_run_dir"])),
        output_dir=Path(str(raw["output_dir"])),
        seed_cases=tuple(str(value) for value in raw["seed_cases"]),
        profiles=tuple(str(value) for value in raw["profiles"]),
        dt_stars=tuple(float(value) for value in raw["dt_stars"]),
        overwrite=bool(raw.get("overwrite", False)),
    )


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _rows(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["condition_id"]: row for row in csv.DictReader(handle)}


def _condition_dir(root: Path, record: dict[str, Any]) -> Path:
    logical_id = str(record["condition_id"])
    value = Path(str(record.get("output_dir", "")))
    candidates = (
        root / "conditions" / logical_id,
        root / logical_id,
        root / value.name,
        value,
    )
    for candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    raise FileNotFoundError(f"condition artifact is missing: {root} / {logical_id}")


def _number(value: Any) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return result if math.isfinite(result) else float("nan")


def _metric(summary: dict[str, Any], name: str, *, statistic: str = "max") -> float:
    item = dict(summary.get("all_step_metrics", {}) or {}).get(name)
    if isinstance(item, dict):
        return _number(item.get(statistic))
    item = dict(summary.get("extrema", {}) or {}).get(name)
    return _number(item.get("value")) if isinstance(item, dict) else float("nan")


def _first_fail(summary: dict[str, Any]) -> tuple[str, float]:
    for gate_name in ("finite", "shape_nonbody", "shape_body"):
        gate = dict(summary.get("gates", {}) or {}).get(gate_name, {})
        when = _number(gate.get("first_observed_fail_t_s"))
        if math.isfinite(when):
            return str(gate.get("first_failure_category") or gate_name), when
    return "none", float("nan")


def _min_bead_distances(positions: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return per-frame body-flag and inter-flagellum bead minima.

    The bead layout is the fixed 2010 project layout (15 body + 3 * 11 flagella).
    These are a fast screening metric; exact segment distances are evaluated at
    the selected frame by the follow-up diagnostic.
    """

    if positions.ndim != 3 or positions.shape[1:] != (48, 3):
        raise ValueError(
            f"expected 2010 n=3 archive shape (N, 48, 3), got {positions.shape}"
        )
    body = positions[:, :15]
    flags = [
        positions[:, 15 + index * 11 : 15 + (index + 1) * 11] for index in range(3)
    ]
    body_flag = np.full(len(positions), np.inf)
    for flag in flags:
        distances = np.linalg.norm(body[:, :, None, :] - flag[:, None, :, :], axis=-1)
        body_flag = np.minimum(body_flag, np.min(distances, axis=(1, 2)))
    flag_flag = np.full(len(positions), np.inf)
    for left in range(3):
        for right in range(left + 1, 3):
            distances = np.linalg.norm(
                flags[left][:, :, None, :] - flags[right][:, None, :, :], axis=-1
            )
            flag_flag = np.minimum(flag_flag, np.min(distances, axis=(1, 2)))
    return body_flag, flag_flag


def _record(
    *, root: Path, record: dict[str, Any], source: str, profile: str, dt_star: float
) -> dict[str, Any]:
    condition_dir = _condition_dir(root, record)
    summary_path = condition_dir / "run_summary.json"
    archive_path = condition_dir / "state_archive.npz"
    if not summary_path.is_file() or not archive_path.is_file():
        raise FileNotFoundError(
            f"missing compact artifacts for {record['condition_id']}"
        )
    summary = _read_json(summary_path)
    strict_pass, strict_reason = _strict(summary)
    failure_category, failure_t_s = _first_fail(summary)
    with np.load(archive_path, allow_pickle=False) as archive:
        times = np.asarray(archive["t"], dtype=float)
        positions = np.asarray(archive["bead_positions_um"], dtype=float)
    body_flag, flag_flag = _min_bead_distances(positions)
    selected_index = (
        int(np.argmin(np.abs(times - failure_t_s)))
        if math.isfinite(failure_t_s)
        else int(np.argmin(body_flag))
    )
    values = dict(record.get("axis_values", {}) or {})
    return {
        "condition_id": record["condition_id"],
        "source": source,
        "profile": profile,
        "dt_star": dt_star,
        "seed_case": "nf{n_flagella:02d}__as{attach_seed:03d}__ps{phase_seed:03d}".format(
            n_flagella=int(values["n_flagella"]),
            attach_seed=int(values["attach_seed"]),
            phase_seed=int(values["phase_seed"]),
        ),
        "strict_pass": strict_pass,
        "strict_reason": strict_reason,
        "first_fail_category": failure_category,
        "first_fail_t_s": failure_t_s,
        "archive_selected_t_s": float(times[selected_index]),
        "body_flag_bead_min_um": float(np.min(body_flag)),
        "body_flag_bead_min_t_s": float(times[int(np.argmin(body_flag))]),
        "body_flag_bead_distance_at_selected_t_um": float(body_flag[selected_index]),
        "flag_flag_bead_min_um": float(np.min(flag_flag)),
        "flag_flag_bead_min_t_s": float(times[int(np.argmin(flag_flag))]),
        "flag_flag_bead_distance_at_selected_t_um": float(flag_flag[selected_index]),
        "hook_len_rel_err_max": _metric(summary, "hook_len_rel_err_max"),
        "flag_bond_rel_err_max": _metric(summary, "flag_bond_rel_err_max"),
        "body_spring_max_stretch_ratio": _metric(
            summary, "body_spring_max_stretch_ratio"
        ),
        "body_bend_max_error_deg": _metric(summary, "body_bend_max_error_deg"),
        "repulsion_mean_body_max_N": _metric(summary, "F_repulsion_mean_body"),
        "repulsion_mean_flag_max_N": _metric(summary, "F_repulsion_mean_flag"),
        "repulsion_basal_max_N": _metric(summary, "local_F_repulsion_basal_region"),
    }


def _campaign_records(root: Path) -> dict[str, dict[str, Any]]:
    manifest = _read_json(root / "run_manifest.json")
    return {str(row["condition_id"]): row for row in manifest.get("conditions", [])}


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _plot(rows: list[dict[str, Any]], output: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharex=True)
    metrics = (
        ("body_flag_bead_min_um", "min body–flag bead distance [µm]"),
        ("hook_len_rel_err_max", "max hook relative error"),
        ("first_fail_t_s", "first failure time [s]"),
    )
    for axis, (metric, label) in zip(axes, metrics):
        for profile, color in (("uniform", "tab:orange"), ("diffusive", "tab:blue")):
            subset = [row for row in rows if row["profile"] == profile]
            axis.scatter(
                [row["dt_star"] for row in subset],
                [row[metric] for row in subset],
                label=profile,
                color=color,
            )
        axis.set(xscale="log", xlabel="dt_star", ylabel=label)
        axis.grid(alpha=0.25)
    axes[0].legend()
    figure.tight_layout()
    figure.savefig(output, dpi=160)
    plt.close(figure)


def analyze(config: ContactConfig) -> Path:
    if config.output_dir.exists():
        if not config.overwrite:
            raise FileExistsError(config.output_dir)
        import shutil

        shutil.rmtree(config.output_dir)
    config.output_dir.mkdir(parents=True)
    references = {
        "uniform": _campaign_records(config.uniform_reference_run_dir),
        "diffusive": _campaign_records(config.diffusive_reference_run_dir),
    }
    diagnostic = _campaign_records(config.diagnostic_run_dir)
    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    for seed_case in config.seed_cases:
        for profile in config.profiles:
            for dt_star in config.dt_stars:
                if math.isclose(dt_star, 1.0e-3, rel_tol=0.0, abs_tol=1e-15):
                    record = references[profile].get(_reference_condition_id(seed_case))
                    root = (
                        config.uniform_reference_run_dir
                        if profile == "uniform"
                        else config.diffusive_reference_run_dir
                    )
                    source = "reused_reference"
                else:
                    condition_id = f"{seed_case}__{profile}__dt{dt_star:.0e}".replace(
                        "e-0", "e-"
                    )
                    record = diagnostic.get(condition_id)
                    root = config.diagnostic_run_dir
                    source = "new_diagnostic"
                if record is None:
                    missing.append(f"{seed_case} {profile} dt={dt_star:.0e}")
                    continue
                rows.append(
                    _record(
                        root=root,
                        record=record,
                        source=source,
                        profile=profile,
                        dt_star=dt_star,
                    )
                )
    if missing:
        raise RuntimeError(
            "required contact diagnostic conditions are missing: " + "; ".join(missing)
        )
    _write_csv(config.output_dir / "contact_stability_conditions.csv", rows)
    _plot(rows, config.output_dir / "contact_stability.png")
    (config.output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "kind": "phase2_torque_profile_dt_contact_diagnostic",
                "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
                "condition_count": len(rows),
                "reused_reference_count": sum(
                    row["source"] == "reused_reference" for row in rows
                ),
                "new_diagnostic_count": sum(
                    row["source"] == "new_diagnostic" for row in rows
                ),
                "outputs": [
                    "contact_stability_conditions.csv",
                    "contact_stability.png",
                ],
                "interpretation": "proximity is diagnostic only; it is not a penetration gate or a profile adoption decision.",
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return config.output_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--diagnostic-run-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    config = load_config(args.config)
    config = ContactConfig(
        **{
            **config.__dict__,
            "diagnostic_run_dir": args.diagnostic_run_dir or config.diagnostic_run_dir,
            "output_dir": args.output_dir or config.output_dir,
            "overwrite": args.overwrite or config.overwrite,
        }
    )
    print(analyze(config))
