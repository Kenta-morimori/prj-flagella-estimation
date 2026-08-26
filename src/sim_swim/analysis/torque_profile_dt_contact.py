"""Distributed contact/stability extraction and strict fragment integration."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any, Iterable
from zoneinfo import ZoneInfo

import matplotlib.pyplot as plt
import numpy as np
import yaml

from sim_swim.analysis.issue203_torque_profile_comparison import _strict

FRAGMENT_CSV = "contact_stability_fragment.csv"
FRAGMENT_MANIFEST = "fragment_manifest.json"
EXPECTED_SOURCES = {
    1.0e-3: "reused_reference",
    3.0e-4: "new_diagnostic",
    1.0e-4: "new_diagnostic",
}


@dataclass(frozen=True)
class CombineConfig:
    output_dir: Path
    seed_cases: tuple[str, ...]
    profiles: tuple[str, ...]
    dt_stars: tuple[float, ...]
    overwrite: bool = False


def load_config(path: Path) -> CombineConfig:
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if type(raw.get("overwrite", False)) is not bool:
        raise ValueError("overwrite must be a boolean")
    return CombineConfig(
        Path(str(raw["output_dir"])),
        tuple(map(str, raw["seed_cases"])),
        tuple(map(str, raw["profiles"])),
        tuple(map(float, raw["dt_stars"])),
        bool(raw.get("overwrite", False)),
    )


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_info() -> dict[str, str | None]:
    def value(*args: str) -> str | None:
        result = subprocess.run(
            ["git", *args], text=True, capture_output=True, check=False
        )
        return result.stdout.strip() if result.returncode == 0 else None

    return {
        "commit": value("rev-parse", "HEAD"),
        "branch": value("branch", "--show-current"),
    }


def _condition_dir(root: Path, record: dict[str, Any]) -> Path:
    logical_id, value = (
        str(record["condition_id"]),
        Path(str(record.get("output_dir", ""))),
    )
    for candidate in (
        root / "conditions" / logical_id,
        root / logical_id,
        root / value.name,
        value,
    ):
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
    if positions.ndim != 3 or positions.shape[1:] != (48, 3):
        raise ValueError(
            f"expected 2010 n=3 archive shape (N, 48, 3), got {positions.shape}"
        )
    body, flags = (
        positions[:, :15],
        [positions[:, 15 + i * 11 : 26 + i * 11] for i in range(3)],
    )
    body_flag, flag_flag = (
        np.full(len(positions), np.inf),
        np.full(len(positions), np.inf),
    )
    for flag in flags:
        body_flag = np.minimum(
            body_flag,
            np.min(
                np.linalg.norm(body[:, :, None, :] - flag[:, None, :, :], axis=-1),
                axis=(1, 2),
            ),
        )
    for left in range(3):
        for right in range(left + 1, 3):
            flag_flag = np.minimum(
                flag_flag,
                np.min(
                    np.linalg.norm(
                        flags[left][:, :, None, :] - flags[right][:, None, :, :],
                        axis=-1,
                    ),
                    axis=(1, 2),
                ),
            )
    return body_flag, flag_flag


def _seed_case(values: dict[str, Any]) -> str:
    return "nf{n_flagella:02d}__as{attach_seed:03d}__ps{phase_seed:03d}".format(
        n_flagella=int(values["n_flagella"]),
        attach_seed=int(values["attach_seed"]),
        phase_seed=int(values["phase_seed"]),
    )


def _record(
    root: Path,
    record: dict[str, Any],
    source: str,
    profile: str,
    dt_star: float,
    run_manifest_sha256: str,
) -> tuple[dict[str, Any], dict[str, str]]:
    condition_dir = _condition_dir(root, record)
    summary_path, archive_path = (
        condition_dir / "run_summary.json",
        condition_dir / "state_archive.npz",
    )
    if not summary_path.is_file() or not archive_path.is_file():
        raise FileNotFoundError(
            f"missing compact artifacts for {record['condition_id']}"
        )
    summary = _read_json(summary_path)
    strict_pass, strict_reason = _strict(summary)
    failure_category, failure_t_s = _first_fail(summary)
    with np.load(archive_path, allow_pickle=False) as archive:
        times, positions = (
            np.asarray(archive["t"], dtype=float),
            np.asarray(archive["bead_positions_um"], dtype=float),
        )
    body_flag, flag_flag = _min_bead_distances(positions)
    selected = (
        int(np.argmin(np.abs(times - failure_t_s)))
        if math.isfinite(failure_t_s)
        else int(np.argmin(body_flag))
    )
    values = dict(record.get("axis_values", {}) or {})
    provenance = {
        "run_manifest_sha256": run_manifest_sha256,
        "run_summary_sha256": _sha256(summary_path),
        "state_archive_sha256": _sha256(archive_path),
    }
    row: dict[str, Any] = {
        "condition_id": record["condition_id"],
        "source": source,
        "profile": profile,
        "dt_star": dt_star,
        "seed_case": _seed_case(values),
        "strict_pass": strict_pass,
        "strict_reason": strict_reason,
        "first_fail_category": failure_category,
        "first_fail_t_s": failure_t_s,
        "archive_selected_t_s": float(times[selected]),
        "body_flag_bead_min_um": float(np.min(body_flag)),
        "body_flag_bead_min_t_s": float(times[int(np.argmin(body_flag))]),
        "body_flag_bead_distance_at_selected_t_um": float(body_flag[selected]),
        "flag_flag_bead_min_um": float(np.min(flag_flag)),
        "flag_flag_bead_min_t_s": float(times[int(np.argmin(flag_flag))]),
        "flag_flag_bead_distance_at_selected_t_um": float(flag_flag[selected]),
        "hook_len_rel_err_max": _metric(summary, "hook_len_rel_err_max"),
        "flag_bond_rel_err_max": _metric(summary, "flag_bond_rel_err_max"),
        "body_spring_max_stretch_ratio": _metric(
            summary, "body_spring_max_stretch_ratio"
        ),
        "body_bend_max_error_deg": _metric(summary, "body_bend_max_error_deg"),
        "repulsion_mean_body_max_N": _metric(summary, "F_repulsion_mean_body"),
        "repulsion_mean_flag_max_N": _metric(summary, "F_repulsion_mean_flag"),
        "repulsion_basal_max_N": _metric(summary, "local_F_repulsion_basal_region"),
        **provenance,
    }
    return row, provenance


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=sorted({key for row in rows for key in row})
        )
        writer.writeheader()
        writer.writerows(rows)


def extract(
    *,
    run_dir: Path,
    output_dir: Path,
    source: str,
    seed_cases: Iterable[str],
    profile: str | None = None,
    dt_star: float | None = None,
    overwrite: bool = False,
) -> Path:
    if source not in {"new_diagnostic", "reused_reference"}:
        raise ValueError("source must be new_diagnostic or reused_reference")
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    seed_set, campaign, rows, artifacts = (
        set(seed_cases),
        _read_json(manifest_path),
        [],
        {},
    )
    run_manifest_sha256 = _sha256(manifest_path)
    for item in campaign.get("conditions", []):
        values = dict(item.get("axis_values", {}) or {})
        if _seed_case(values) not in seed_set:
            continue
        actual_profile = str(
            values.get(
                "torque_profile",
                values.get("torque_distribution_profile", profile or ""),
            )
        )
        actual_dt = _number(values.get("dt_star", dt_star))
        if not actual_profile or not math.isfinite(actual_dt):
            raise ValueError(
                f"profile/dt_star unavailable for {item['condition_id']}; provide --profile and --dt-star"
            )
        if profile is not None and actual_profile != profile:
            raise ValueError(
                f"profile mismatch for {item['condition_id']}: {actual_profile} != {profile}"
            )
        if dt_star is not None and not math.isclose(
            actual_dt, dt_star, rel_tol=0.0, abs_tol=1e-15
        ):
            raise ValueError(
                f"dt_star mismatch for {item['condition_id']}: {actual_dt} != {dt_star}"
            )
        row, hashes = _record(
            run_dir, item, source, actual_profile, actual_dt, run_manifest_sha256
        )
        rows.append(row)
        artifacts[str(item["condition_id"])] = hashes
    if not rows:
        raise RuntimeError("no requested conditions were found")
    csv_path = output_dir / FRAGMENT_CSV
    _write_csv(csv_path, rows)
    (output_dir / FRAGMENT_MANIFEST).write_text(
        json.dumps(
            {
                "kind": "phase2_torque_profile_dt_contact_fragment",
                "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
                "source": source,
                "source_run_dir": str(run_dir),
                "source_run_manifest_sha256": run_manifest_sha256,
                "fragment_csv": FRAGMENT_CSV,
                "fragment_csv_sha256": _sha256(csv_path),
                "condition_count": len(rows),
                "conditions": artifacts,
                "git": _git_info(),
                "archive_transfer": "not_required",
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return output_dir


def _plot(rows: list[dict[str, Any]], output: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharex=True)
    for axis, (metric, label) in zip(
        axes,
        (
            ("body_flag_bead_min_um", "min body–flag bead distance [µm]"),
            ("hook_len_rel_err_max", "max hook relative error"),
            ("first_fail_t_s", "first failure time [s]"),
        ),
    ):
        for profile, color in (("uniform", "tab:orange"), ("diffusive", "tab:blue")):
            subset = [row for row in rows if row["profile"] == profile]
            axis.scatter(
                [float(row["dt_star"]) for row in subset],
                [_number(row[metric]) for row in subset],
                label=profile,
                color=color,
            )
        axis.set(xscale="log", xlabel="dt_star", ylabel=label)
        axis.grid(alpha=0.25)
    axes[0].legend()
    figure.tight_layout()
    figure.savefig(output, dpi=160)
    plt.close(figure)


def combine(config: CombineConfig, fragment_dirs: Iterable[Path]) -> Path:
    expected = {
        (seed, profile, dt)
        for seed in config.seed_cases
        for profile in config.profiles
        for dt in config.dt_stars
    }
    rows: list[dict[str, Any]] = []
    for directory in fragment_dirs:
        manifest_path, csv_path = (
            directory / FRAGMENT_MANIFEST,
            directory / FRAGMENT_CSV,
        )
        if not manifest_path.is_file() or not csv_path.is_file():
            raise FileNotFoundError(f"fragment files missing in {directory}")
        manifest = _read_json(manifest_path)
        if manifest.get(
            "kind"
        ) != "phase2_torque_profile_dt_contact_fragment" or manifest.get(
            "fragment_csv_sha256"
        ) != _sha256(csv_path):
            raise ValueError(f"fragment provenance mismatch: {directory}")
        with csv_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                hashes = dict(
                    manifest.get("conditions", {}).get(row["condition_id"], {}) or {}
                )
                if any(
                    row.get(key) != value for key, value in hashes.items()
                ) or not all(hashes.values()):
                    raise ValueError(
                        f"condition provenance mismatch: {directory}/{row['condition_id']}"
                    )
                if row.get("source") != manifest.get("source"):
                    raise ValueError(
                        f"source provenance mismatch: {directory}/{row['condition_id']}"
                    )
                rows.append(row)
    found: dict[tuple[str, str, float], dict[str, Any]] = {}
    for row in rows:
        key = (row["seed_case"], row["profile"], _number(row["dt_star"]))
        if key not in expected:
            raise ValueError(f"unexpected condition: {key}")
        if row["source"] != EXPECTED_SOURCES.get(key[2]):
            raise ValueError(f"source mismatch for {key}")
        if key in found:
            raise ValueError(f"duplicate condition: {key}")
        found[key] = row
    if set(found) != expected:
        raise ValueError(
            f"missing required conditions: {sorted(expected - set(found))}"
        )
    if config.output_dir.exists():
        if not config.overwrite:
            raise FileExistsError(config.output_dir)
        import shutil

        shutil.rmtree(config.output_dir)
    config.output_dir.mkdir(parents=True)
    ordered = [found[key] for key in sorted(found)]
    _write_csv(config.output_dir / "contact_stability_conditions.csv", ordered)
    _plot(ordered, config.output_dir / "contact_stability.png")
    (config.output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "kind": "phase2_torque_profile_dt_contact_diagnostic",
                "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
                "condition_count": len(ordered),
                "reused_reference_count": sum(
                    row["source"] == "reused_reference" for row in ordered
                ),
                "new_diagnostic_count": sum(
                    row["source"] == "new_diagnostic" for row in ordered
                ),
                "input_fragments": [str(item) for item in fragment_dirs],
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
    commands = parser.add_subparsers(dest="command", required=True)
    extract_parser = commands.add_parser("extract")
    extract_parser.add_argument("--run-dir", type=Path, required=True)
    extract_parser.add_argument("--output-dir", type=Path, required=True)
    extract_parser.add_argument("--source", required=True)
    extract_parser.add_argument("--seed-case", action="append", required=True)
    extract_parser.add_argument("--profile")
    extract_parser.add_argument("--dt-star", type=float)
    extract_parser.add_argument("--overwrite", action="store_true")
    combine_parser = commands.add_parser("combine")
    combine_parser.add_argument("--config", type=Path, required=True)
    combine_parser.add_argument(
        "--fragment-dir", type=Path, action="append", required=True
    )
    combine_parser.add_argument("--output-dir", type=Path)
    combine_parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    if args.command == "extract":
        print(
            extract(
                run_dir=args.run_dir,
                output_dir=args.output_dir,
                source=args.source,
                seed_cases=args.seed_case,
                profile=args.profile,
                dt_star=args.dt_star,
                overwrite=args.overwrite,
            )
        )
        return
    config = load_config(args.config)
    config = CombineConfig(
        **{
            **config.__dict__,
            "output_dir": args.output_dir or config.output_dir,
            "overwrite": args.overwrite or config.overwrite,
        }
    )
    print(combine(config, args.fragment_dir))
