#!/usr/bin/env python3
"""Build a Phase 2.8 flagella-count behavior dataset from a run manifest."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.cli_profiles import (
    key_value_args_to_cli_args,
    split_config_key,
)
from sim_swim.analysis.multi_run_campaign import apply_campaign_cli_overrides, load_yaml


METADATA_FIELDS = [
    "sample_id",
    "dataset_id",
    "n_flagella",
    "seed",
    "attach_seed",
    "phase_seed",
    "duration_s",
    "dt_star",
    "torque_Nm",
    "force_distribution",
    "condition_tag",
]

QUALITY_FIELDS = [
    "quality_class",
    "shape_pass",
    "relaxed_pass",
    "use_for_analysis",
    "use_for_ml_candidate",
    "review_required",
    "valid_duration_s",
]

FEATURE_FIELDS = [
    "cell_displacement",
    "cell_path_length",
    "cell_mean_speed",
    "cell_speed_std",
    "cell_speed_cv",
    "cell_straightness",
    "cell_axis_angle_change",
    "cell_axis_angle_std",
    "cell_angular_velocity_mean",
    "cell_angular_velocity_std",
    "cell_angular_velocity_rms",
    "cell_wobble",
    "flagella_axis_alignment",
    "flagella_axis_spread",
    "flagella_axis_pair_angle_mean",
    "flagella_axis_pair_angle_max",
    "flagella_axis_rear_alignment",
    "cell_flagella_axis_angle",
    "cell_flagella_axis_angle_std",
    "cell_flagella_axis_stability",
    "hook_drift",
    "hook_wrapped",
    "flyaway",
    "abnormal_rotation",
    "first_fail_category",
]

BOOKKEEPING_FIELDS = [
    "run_status",
    "raw_dir",
    "step_summary_csv",
    "timeseries_csv",
    "step_count",
    "missing_value_count",
]

SUMMARY_FIELDS = METADATA_FIELDS + QUALITY_FIELDS + FEATURE_FIELDS + BOOKKEEPING_FIELDS

FLAGELLA_AXIS_FIELDS = [
    "flagella_axis_alignment",
    "flagella_axis_spread",
    "flagella_axis_pair_angle_mean",
    "flagella_axis_pair_angle_max",
    "flagella_axis_rear_alignment",
]

FLAGELLA_RELATION_FIELDS = [
    "cell_flagella_axis_angle",
    "cell_flagella_axis_angle_std",
    "cell_flagella_axis_stability",
]

FLAGELLA_TIMESERIES_FIELDS = [
    "flag_helix_axis_alignment_order",
    "flag_helix_axis_mean_deviation_deg_max",
    "flag_helix_axis_pair_angle_deg_mean",
    "flag_helix_axis_pair_angle_deg_max",
    "flag_helix_axis_rearward_projection_min",
    "bundle_axis_vs_body_axis_angle_deg",
    "bundle_axis_vs_rear_angle_deg",
]


def _load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _write_yaml(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _to_float(value: Any, default: float = float("nan")) -> float:
    try:
        if value in (None, ""):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _to_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _finite(values: list[float]) -> list[float]:
    return [float(v) for v in values if math.isfinite(float(v))]


def _mean(values: list[float]) -> float:
    finite = _finite(values)
    return sum(finite) / len(finite) if finite else float("nan")


def _std(values: list[float]) -> float:
    finite = _finite(values)
    if len(finite) < 2:
        return 0.0 if len(finite) == 1 else float("nan")
    avg = sum(finite) / len(finite)
    return math.sqrt(sum((v - avg) ** 2 for v in finite) / (len(finite) - 1))


def _path_length_um(rows: list[dict[str, str]]) -> float:
    if len(rows) < 2:
        return float("nan")
    length = 0.0
    prev_t = _to_float(rows[0].get("t_s"))
    for row in rows[1:]:
        t_s = _to_float(row.get("t_s"))
        speed = _to_float(row.get("body_speed_um_s"))
        if math.isfinite(t_s) and math.isfinite(prev_t) and math.isfinite(speed):
            length += max(t_s - prev_t, 0.0) * speed
        prev_t = t_s
    return length


def _angular_velocity_values(rows: list[dict[str, str]]) -> list[float]:
    out: list[float] = []
    for row in rows:
        dt_s = _to_float(row.get("dt_s"))
        step_angle_deg = _to_float(row.get("body_axis_step_angle_deg"))
        if math.isfinite(dt_s) and dt_s > 0.0 and math.isfinite(step_angle_deg):
            out.append(math.radians(step_angle_deg) / dt_s)
    return out


def _quality(last: dict[str, str] | None, rows: list[dict[str, str]]) -> dict[str, Any]:
    if last is None:
        return {
            "quality_class": "missing_raw",
            "shape_pass": False,
            "relaxed_pass": False,
            "use_for_analysis": False,
            "use_for_ml_candidate": False,
            "review_required": False,
            "valid_duration_s": float("nan"),
        }

    finite_pass = _to_bool(last.get("finite_pass"))
    shape_pass = _to_bool(last.get("shape_pass_nonbody_strict"))
    relaxed_pass = _to_bool(last.get("shape_pass_nonbody_hook_len_relaxed"))
    if not finite_pass:
        quality_class = "invalid_numeric"
    elif shape_pass:
        quality_class = "strict_pass"
    elif relaxed_pass:
        quality_class = "relaxed_pass"
    else:
        quality_class = "fail"
    use_for_analysis = bool(finite_pass and relaxed_pass)
    return {
        "quality_class": quality_class,
        "shape_pass": shape_pass,
        "relaxed_pass": relaxed_pass,
        "use_for_analysis": use_for_analysis,
        "use_for_ml_candidate": use_for_analysis,
        "review_required": bool(finite_pass and not shape_pass),
        "valid_duration_s": _to_float(last.get("t_s")),
    }


def summarize_sample(
    *,
    sample: dict[str, Any],
    dataset_id: str,
    timeseries_csv: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[str]]:
    raw_dir = Path(
        sample.get("raw_dir") or (sample.get("outputs", {}) or {}).get("raw_dir", "")
    )
    step_summary_csv = Path(
        (sample.get("outputs", {}) or {}).get("step_summary_csv")
        or raw_dir / "step_summary.csv"
    )
    rows = _read_csv(step_summary_csv) if step_summary_csv.is_file() else []
    last = rows[-1] if rows else None

    metadata = {
        "sample_id": sample.get("sample_id", ""),
        "dataset_id": dataset_id,
        "n_flagella": sample.get("n_flagella", ""),
        "seed": sample.get("seed", ""),
        "attach_seed": sample.get("attach_seed", sample.get("seed", "")),
        "phase_seed": sample.get("phase_seed", sample.get("seed", "")),
        "duration_s": sample.get("duration_s", ""),
        "dt_star": sample.get("dt_star", ""),
        "torque_Nm": sample.get("torque_Nm", ""),
        "force_distribution": sample.get("force_distribution", ""),
        "condition_tag": sample.get("condition_tag", ""),
    }
    quality = _quality(last, rows)

    body_speeds = [_to_float(row.get("body_speed_um_s")) for row in rows]
    body_speed_mean = _mean(body_speeds)
    body_speed_std = _std(body_speeds)
    angular_velocity = _angular_velocity_values(rows)
    relation_angles = [
        _to_float(row.get("bundle_axis_vs_body_axis_angle_deg")) for row in rows
    ]
    relation_std = _std(relation_angles)
    path_length = _path_length_um(rows)
    displacement = _to_float(last.get("body_displacement_um")) if last else float("nan")
    valid_duration = _to_float(quality["valid_duration_s"])

    features = {
        "cell_displacement": displacement,
        "cell_path_length": path_length,
        "cell_mean_speed": (
            displacement / valid_duration
            if math.isfinite(displacement) and valid_duration > 0.0
            else float("nan")
        ),
        "cell_speed_std": body_speed_std,
        "cell_speed_cv": (
            body_speed_std / body_speed_mean
            if math.isfinite(body_speed_std) and abs(body_speed_mean) > 1.0e-30
            else float("nan")
        ),
        "cell_straightness": (
            displacement / path_length
            if math.isfinite(displacement)
            and math.isfinite(path_length)
            and path_length > 0.0
            else float("nan")
        ),
        "cell_axis_angle_change": (
            _to_float(last.get("body_axis_cumulative_angle_deg"))
            if last
            else float("nan")
        ),
        "cell_axis_angle_std": _std(
            [_to_float(row.get("body_axis_cumulative_angle_deg")) for row in rows]
        ),
        "cell_angular_velocity_mean": _mean(angular_velocity),
        "cell_angular_velocity_std": _std(angular_velocity),
        "cell_angular_velocity_rms": (
            _to_float(last.get("body_angular_velocity_rms_rad_s"))
            if last
            else float("nan")
        ),
        "cell_wobble": (
            _to_float(last.get("body_axis_wobble_rms_deg")) if last else float("nan")
        ),
        "flagella_axis_alignment": (
            _to_float(last.get("flag_helix_axis_alignment_order"))
            if last
            else float("nan")
        ),
        "flagella_axis_spread": (
            _to_float(last.get("flag_helix_axis_mean_deviation_deg_max"))
            if last
            else float("nan")
        ),
        "flagella_axis_pair_angle_mean": (
            _to_float(last.get("flag_helix_axis_pair_angle_deg_mean"))
            if last
            else float("nan")
        ),
        "flagella_axis_pair_angle_max": (
            _to_float(last.get("flag_helix_axis_pair_angle_deg_max"))
            if last
            else float("nan")
        ),
        "flagella_axis_rear_alignment": (
            _to_float(last.get("flag_helix_axis_rearward_projection_min"))
            if last
            else float("nan")
        ),
        "cell_flagella_axis_angle": (
            _to_float(last.get("bundle_axis_vs_body_axis_angle_deg"))
            if last
            else float("nan")
        ),
        "cell_flagella_axis_angle_std": relation_std,
        "cell_flagella_axis_stability": (
            1.0 / (1.0 + relation_std) if math.isfinite(relation_std) else float("nan")
        ),
        "hook_drift": _to_float(last.get("hook_len_rel_err_max"))
        if last
        else float("nan"),
        "hook_wrapped": (
            str(last.get("first_fail_category_nonbody", "")).strip() == "hook"
            if last
            else False
        ),
        "flyaway": not bool(quality["shape_pass"])
        and str(last.get("first_fail_category_nonbody", "") if last else "").strip()
        in {"finite", "nan", "inf"},
        "abnormal_rotation": False,
        "first_fail_category": (
            last.get("first_fail_category_nonbody", "") if last else "missing_raw"
        ),
    }

    n_flagella = int(_to_float(metadata["n_flagella"], default=0.0))
    if n_flagella <= 1:
        for key in FLAGELLA_AXIS_FIELDS + FLAGELLA_RELATION_FIELDS:
            features[key] = float("nan")

    summary = {
        **metadata,
        **quality,
        **features,
        "run_status": sample.get("status", ""),
        "raw_dir": str(raw_dir),
        "step_summary_csv": str(step_summary_csv),
        "timeseries_csv": str(timeseries_csv),
        "step_count": len(rows),
    }
    summary["missing_value_count"] = sum(
        1
        for key in FEATURE_FIELDS
        if isinstance(summary.get(key), float) and not math.isfinite(summary[key])
    )

    ts_rows: list[dict[str, Any]] = []
    for row in rows:
        processed = dict(row)
        if n_flagella <= 1:
            for key in FLAGELLA_TIMESERIES_FIELDS:
                if key in processed:
                    processed[key] = float("nan")
        ts_rows.append({**metadata, **processed})
    fieldnames = list(metadata.keys()) + (list(rows[0].keys()) if rows else [])
    return summary, ts_rows, fieldnames


def _git_info() -> dict[str, Any]:
    def run(cmd: list[str]) -> str:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True).strip()

    try:
        commit = run(["git", "rev-parse", "HEAD"])
        return {
            "commit": commit,
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


def _now_jst() -> str:
    return datetime.now(ZoneInfo("Asia/Tokyo")).isoformat()


def _merge_nested(dst: dict[str, Any], src: dict[str, Any]) -> dict[str, Any]:
    merged = dict(dst)
    for key, value in src.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_nested(merged[key], value)
        else:
            merged[key] = value
    return merged


def _resolve_run_manifest_path(
    *,
    campaign: dict[str, Any],
    run_manifest_path: Path | None,
    run_dir: Path | None,
) -> Path:
    if run_manifest_path is not None and run_dir is not None:
        raise ValueError("Use either run_manifest=PATH or run_dir=PATH, not both")
    if run_manifest_path is not None:
        return run_manifest_path
    if run_dir is not None:
        return run_dir / "run_manifest.json"
    output_cfg = dict(campaign.get("output", {}) or {})
    if bool(output_cfg.get("timestamp_subdir", True)):
        raise ValueError(
            "run_manifest=... or run_dir=... is required when "
            "output.timestamp_subdir is true"
        )
    return Path(str(output_cfg["base_dir"])) / "run_manifest.json"


def _condition_to_sample(
    *,
    condition: dict[str, Any],
    campaign: dict[str, Any],
    sample_id_template: str,
) -> dict[str, Any]:
    axis_values = dict(condition.get("axis_values", {}) or {})
    effective_overrides = _merge_nested(
        dict(campaign.get("base_overrides", {}) or {}),
        dict(condition.get("config_overrides", {}) or {}),
    )
    time_cfg = dict(effective_overrides.get("time", {}) or {})
    motor_cfg = dict(effective_overrides.get("motor", {}) or {})
    seed_cfg = dict(effective_overrides.get("seed", {}) or {})
    n_flagella = int(axis_values.get("n_flagella", 0))
    attach_seed = int(axis_values.get("attach_seed", seed_cfg.get("attach_seed", 0)))
    phase_seed = int(axis_values.get("phase_seed", seed_cfg.get("phase_seed", 0)))
    sample_id = sample_id_template.format(
        n_flagella=n_flagella,
        attach_seed=attach_seed,
        phase_seed=phase_seed,
        condition_id=condition.get("condition_id", ""),
    )
    raw_dir = Path(str(condition["output_dir"]))
    step_summary_csv = raw_dir / "step_summary.csv"
    return {
        "sample_id": sample_id,
        "condition_id": condition.get("condition_id", ""),
        "condition_tag": condition.get("condition_label", ""),
        "n_flagella": n_flagella,
        "seed": attach_seed,
        "attach_seed": attach_seed,
        "phase_seed": phase_seed,
        "duration_s": time_cfg.get("duration_s", ""),
        "dt_star": time_cfg.get("dt_star", ""),
        "torque_Nm": motor_cfg.get("torque_Nm", ""),
        "force_distribution": motor_cfg.get("force_distribution", ""),
        "status": "completed" if step_summary_csv.is_file() else "missing_raw",
        "sample_dir": str(raw_dir),
        "raw_dir": str(raw_dir),
        "config_overrides": effective_overrides,
        "outputs": {
            "raw_dir": str(raw_dir),
            "step_summary_csv": str(step_summary_csv),
            "state_archive_npz": str(raw_dir / "state_archive.npz"),
            "trajectory_csv": str(raw_dir / "trajectory.csv"),
        },
    }


def _samples_from_multi_run_manifest(
    *,
    run_manifest: dict[str, Any],
    campaign: dict[str, Any],
    sample_id_template: str,
) -> list[dict[str, Any]]:
    if "conditions" not in run_manifest:
        raise ValueError("Multi-run run_manifest.json must contain conditions")
    return [
        _condition_to_sample(
            condition=condition,
            campaign=campaign,
            sample_id_template=sample_id_template,
        )
        for condition in run_manifest.get("conditions", [])
    ]


def build_dataset(
    *,
    campaign_config_path: Path,
    run_manifest_path: Path | None,
    run_dir: Path | None = None,
    overwrite: bool,
    cli_overrides: list[str] | None = None,
) -> Path:
    cli_overrides = cli_overrides or []
    raw_campaign = load_yaml(campaign_config_path)
    campaign = apply_campaign_cli_overrides(raw_campaign, cli_overrides)
    dataset_cfg = dict(campaign.get("dataset", {}) or {})
    if not dataset_cfg:
        raise ValueError("campaign config must contain a dataset section")
    dataset_id = str(dataset_cfg["dataset_id"])
    dataset_dir = Path(str(dataset_cfg["output_dir"]))
    if dataset_dir.exists() and overwrite:
        shutil.rmtree(dataset_dir)
    dataset_dir.mkdir(parents=True, exist_ok=True)
    effective_campaign_config_path = dataset_dir / "campaign_config_used.yaml"
    _write_yaml(effective_campaign_config_path, campaign)

    run_manifest_path = _resolve_run_manifest_path(
        campaign=campaign,
        run_manifest_path=run_manifest_path,
        run_dir=run_dir,
    )
    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    feature_schema_path = Path(
        str(
            dataset_cfg.get(
                "feature_schema",
                "conf/phase2_analysis/flagella_count_behavior_features.yaml",
            )
        )
    )
    shutil.copyfile(feature_schema_path, dataset_dir / "feature_schema_used.yaml")
    if str(dataset_cfg.get("timeseries_sampling", "all_steps")) != "all_steps":
        raise ValueError("Only dataset.timeseries_sampling=all_steps is supported")
    sample_id_template = str(
        dataset_cfg.get(
            "sample_id_template",
            "nf{n_flagella:02d}_as{attach_seed:03d}_ps{phase_seed:03d}",
        )
    )
    samples = _samples_from_multi_run_manifest(
        run_manifest=run_manifest,
        campaign=campaign,
        sample_id_template=sample_id_template,
    )
    base_cfg = _load_yaml(Path(str(campaign["base_config"])))
    config_dir = dataset_dir / "configs"
    for sample in samples:
        sample_id = str(sample["sample_id"])
        sample_config_path = config_dir / f"{sample_id}.yaml"
        sample_config = _merge_nested(
            base_cfg,
            dict(sample.get("config_overrides", {}) or {}),
        )
        _write_yaml(sample_config_path, sample_config)
        sample["config_path"] = str(sample_config_path)

    timeseries_dir = dataset_dir / "timeseries"
    summaries: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    for sample in samples:
        sample_id = str(sample.get("sample_id", ""))
        timeseries_csv = timeseries_dir / f"{sample_id}.csv"
        summary, ts_rows, ts_fieldnames = summarize_sample(
            sample=sample,
            dataset_id=dataset_id,
            timeseries_csv=timeseries_csv,
        )
        _write_csv(timeseries_csv, ts_rows, ts_fieldnames or METADATA_FIELDS)
        summaries.append(summary)
        qc_rows.append(
            {
                "sample_id": summary["sample_id"],
                "n_flagella": summary["n_flagella"],
                "seed": summary["seed"],
                "quality_class": summary["quality_class"],
                "shape_pass": summary["shape_pass"],
                "relaxed_pass": summary["relaxed_pass"],
                "use_for_analysis": summary["use_for_analysis"],
                "review_required": summary["review_required"],
                "first_fail_category": summary["first_fail_category"],
                "step_count": summary["step_count"],
                "missing_value_count": summary["missing_value_count"],
            }
        )

    summary_path = dataset_dir / "summary.csv"
    qc_path = dataset_dir / "qc_summary.csv"
    _write_csv(summary_path, summaries, SUMMARY_FIELDS)
    _write_csv(
        qc_path,
        qc_rows,
        [
            "sample_id",
            "n_flagella",
            "seed",
            "quality_class",
            "shape_pass",
            "relaxed_pass",
            "use_for_analysis",
            "review_required",
            "first_fail_category",
            "step_count",
            "missing_value_count",
        ],
    )

    manifest = {
        "dataset_id": dataset_id,
        "run_batch_id": run_manifest.get("output_root", ""),
        "created_at": _now_jst(),
        "campaign_config": str(campaign_config_path),
        "effective_campaign_config": campaign,
        "effective_campaign_config_yaml": str(effective_campaign_config_path),
        "cli_overrides": list(cli_overrides),
        "run_manifest": str(run_manifest_path),
        "feature_schema_source": str(feature_schema_path),
        "outputs": {
            "dataset_dir": str(dataset_dir),
            "summary_csv": str(summary_path),
            "qc_summary_csv": str(qc_path),
            "timeseries_dir": str(timeseries_dir),
            "feature_schema_used_yaml": str(dataset_dir / "feature_schema_used.yaml"),
        },
        "sample_count": len(summaries),
        "quality_counts": {
            key: sum(1 for row in summaries if row["quality_class"] == key)
            for key in sorted({str(row["quality_class"]) for row in summaries})
        },
        "timeseries_sampling": dataset_cfg.get("timeseries_sampling", "all_steps"),
        "samples": samples,
        "git": _git_info(),
    }
    manifest_path = dataset_dir / "dataset_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return dataset_dir


def main(argv: list[str] | None = None) -> None:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    config_from_key, parser_argv = split_config_key(raw_argv)
    parser_argv = key_value_args_to_cli_args(parser_argv)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
    )
    parser.add_argument("--run-manifest", type=Path, default=None)
    parser.add_argument("--run-dir", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "overrides",
        nargs="*",
        help=(
            "Optional campaign overrides as key=value "
            "(e.g. dataset.dataset_id=my_dataset)"
        ),
    )
    args = parser.parse_args(parser_argv)
    if config_from_key is not None and args.config is not None:
        parser.error("Use either config=PATH or --config PATH (not both)")
    config = config_from_key or args.config
    if config is None:
        parser.error("config=PATH or --config PATH is required")

    dataset_dir = build_dataset(
        campaign_config_path=config,
        run_manifest_path=args.run_manifest,
        run_dir=args.run_dir,
        overwrite=bool(args.overwrite),
        cli_overrides=args.overrides,
    )
    print(dataset_dir)


if __name__ == "__main__":
    main()
