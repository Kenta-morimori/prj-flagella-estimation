"""Aggregate reusable Phase 2 model-development evaluation screens.

This module deliberately separates numerical/physical QC from downstream
swimming-feature analysis.  It reads completed multi-run campaigns only; it
never starts simulations.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.analysis.flagella_count_behavior import load_state_archive
from sim_swim.sim.debug_summary import (
    NONBODY_FLAG_BEND_ERR_MAX_DEG_LIMIT,
    NONBODY_FLAG_BOND_REL_ERR_MAX_LIMIT,
    NONBODY_FLAG_TORSION_ERR_MAX_DEG_LIMIT,
    NONBODY_HOOK_REL_ERR_MAX_LIMIT,
)


SCREEN_METRICS: tuple[tuple[str, str], ...] = (
    ("hook_angle_err_max_deg", "max hook angle error [deg]"),
    ("hook_len_rel_err_max", "max hook length relative error"),
    ("flag_bond_rel_err_max", "max flag bond relative error"),
    ("body_spring_max_stretch_ratio", "max body spring stretch ratio"),
    ("motor_force_balance_residual_ratio", "max motor force residual ratio"),
    ("motor_torque_balance_residual_ratio", "max motor torque residual ratio"),
    ("wall_time_s", "wall time [s]"),
    ("steps_per_s", "steps/s"),
)
LONG_DURATION_REQUIRED_ARTIFACTS = (
    "run_summary.json",
    "performance.json",
    "state_archive.npz",
)
PARTIAL_ARCHIVE_NAME = "state_archive.partial.npz"
OPERATIONAL_LOG_NAMES = ("run.log", "render.log", "stdout.log", "stderr.log")


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _bool(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"1", "true", "yes"}


def _maximum(summary: dict[str, Any], name: str) -> float:
    try:
        return float(dict(summary.get("all_step_metrics", {}) or {})[name]["max"])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def _gate_failed(summary: dict[str, Any], name: str) -> bool:
    return bool(
        dict(summary.get("gates", {}).get(name, {}) or {}).get("any_fail", True)
    )


def _all_finite(values: list[float]) -> bool:
    return all(math.isfinite(value) for value in values)


def screen_status(summary: dict[str, Any], *, qc: dict[str, Any]) -> str:
    """Return PASS/FAIL using physical gates without hook-angle screening.

    Hook angle remains a diagnostic metric.  Its old 30-degree generic gate
    is intentionally excluded: it detects an initial geometric mismatch, not
    an irreversible shape failure.  Hook length, flag geometry, body and
    motor balance remain required.
    """

    if _gate_failed(summary, "finite") or _gate_failed(summary, "shape_body"):
        return "fail"
    hook = [
        _maximum(summary, "local_attach_first_rel_err"),
        _maximum(summary, "hook_len_rel_err_max"),
    ]
    flag = [
        _maximum(summary, "flag_bond_rel_err_max"),
        _maximum(summary, "flag_bend_err_max_deg"),
        _maximum(summary, "flag_torsion_err_max_deg"),
    ]
    motor = [
        _maximum(summary, "motor_force_balance_residual_ratio"),
        _maximum(summary, "motor_torque_balance_residual_ratio"),
    ]
    if not _all_finite([*hook, *flag, *motor]):
        return "fail"
    if motor[0] > float(qc["max_motor_force_balance_residual_ratio"]) or motor[
        1
    ] > float(qc["max_motor_torque_balance_residual_ratio"]):
        return "fail"
    if max(hook) > NONBODY_HOOK_REL_ERR_MAX_LIMIT:
        return "fail"
    if (
        flag[0] > NONBODY_FLAG_BOND_REL_ERR_MAX_LIMIT
        or flag[1] > NONBODY_FLAG_BEND_ERR_MAX_DEG_LIMIT
        or flag[2] > NONBODY_FLAG_TORSION_ERR_MAX_DEG_LIMIT
    ):
        return "fail"
    return "pass"


def _expected_conditions(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {
        str(condition["condition_id"]): condition
        for condition in build_campaign_conditions(normalize_campaign_config(config))
    }


def _record_key(record: dict[str, Any], axis_names: list[str]) -> tuple[Any, ...]:
    """Return the declared evaluation-axis identity for one condition."""

    axes = dict(record.get("axis_values", {}) or {})
    time = dict(record.get("time", {}) or {})
    values: list[Any] = []
    for name in axis_names:
        value = axes[name] if name in axes else time.get(name)
        if value is None:
            raise KeyError(f"Missing evaluation axis {name}")
        if name in {"motor_torque", "dt_star"}:
            values.append(float(value))
        elif name == "n_flagella":
            values.append(int(value))
        else:
            values.append(str(value))
    return tuple(values)


def _development_contract(config: dict[str, Any]) -> dict[str, Any]:
    contract = dict(config.get("development_evaluation", {}) or {})
    if contract.get("stage") not in {"short_screen", "long_duration"}:
        raise ValueError(
            "development_evaluation.stage must be short_screen or long_duration"
        )
    if int(contract.get("expected_condition_count", 0)) <= 0:
        raise ValueError("development_evaluation.expected_condition_count is required")
    if not isinstance(contract.get("axes"), list) or not contract["axes"]:
        raise ValueError("development_evaluation.axes is required")
    qc = dict(contract.get("qc", {}) or {})
    for name in (
        "max_motor_force_balance_residual_ratio",
        "max_motor_torque_balance_residual_ratio",
    ):
        if not math.isfinite(float(qc.get(name, float("nan")))):
            raise ValueError(f"development_evaluation.qc.{name} is required")
    if not isinstance(contract.get("accepted_source_campaigns"), list):
        raise ValueError("development_evaluation.accepted_source_campaigns is required")
    return contract


def _reused_source_spec(
    contract: dict[str, Any], campaign_config: str
) -> tuple[dict[str, Any], list[str]]:
    """Return declared defaults for a completed historical source campaign."""

    sources = dict(contract.get("reused_source_campaigns", {}) or {})
    source = dict(sources.get(campaign_config, {}) or {})
    defaults = dict(source.get("axis_defaults", {}) or {})
    paths = [str(path) for path in source.get("relaxed_config_override_paths", [])]
    return defaults, paths


def _normalize_reused_source_record(
    record: dict[str, Any], *, axis_defaults: dict[str, Any]
) -> dict[str, Any]:
    """Add only contract-declared axis defaults to a historical record."""

    normalized = copy.deepcopy(record)
    axes = dict(normalized.get("axis_values", {}) or {})
    axes.update(axis_defaults)
    normalized["axis_values"] = axes
    return normalized


def _without_override_paths(value: Any, paths: list[str]) -> Any:
    """Drop declared compatibility-only dotted paths before comparison."""

    normalized = copy.deepcopy(value)
    for path in paths:
        current = normalized
        parts = path.split(".")
        for part in parts[:-1]:
            if not isinstance(current, dict) or part not in current:
                current = None
                break
            current = current[part]
        if isinstance(current, dict):
            current.pop(parts[-1], None)
    return _drop_empty_mappings(normalized)


def _drop_empty_mappings(value: Any) -> Any:
    """Remove empty mappings left by a declared compatibility-only override."""

    if isinstance(value, dict):
        return {
            key: item
            for key, raw_item in value.items()
            if (item := _drop_empty_mappings(raw_item)) != {}
        }
    if isinstance(value, list):
        return [_drop_empty_mappings(item) for item in value]
    return value


def _source_rows(run_dir: Path) -> dict[str, dict[str, str]]:
    with (run_dir / "summary.csv").open(encoding="utf-8", newline="") as handle:
        return {str(row["condition_id"]): row for row in csv.DictReader(handle)}


def _campaign_git_provenance(
    *, run_dir: Path, manifest: dict[str, Any]
) -> dict[str, Any]:
    """Return verified Git provenance, including parallel aggregate campaigns.

    Aggregated campaign manifests deliberately contain only simulation metadata.
    Their fixed-commit provenance is stored in the parent parallel job manifest.
    Accept that layout only when the recorded porcelain status proves a clean
    checkout; otherwise never infer Git state from the current workspace.
    """

    git = dict(manifest.get("git", {}) or {})
    if str(git.get("commit") or "") and git.get("is_clean") is True:
        return git
    job_manifest = next(
        (
            path
            for path in (
                run_dir / "job_manifest.json",
                run_dir.parent / "job_manifest.json",
            )
            if path.is_file()
        ),
        None,
    )
    if job_manifest is None:
        raise ValueError(f"Invalid Git provenance in {run_dir}")
    job = _read_json(job_manifest)
    queued_git = dict(dict(job.get("provenance", {}) or {}).get("git", {}) or {})
    status = str(queued_git.get("status") or "")
    commit = str(queued_git.get("commit") or "")
    if not commit or not status.startswith("## ") or "\n" in status:
        raise ValueError(f"Invalid Git provenance in {run_dir}")
    return {
        "commit": commit,
        "is_clean": True,
        "source": "parallel_job_manifest",
        "status": status,
    }


def _resolve_condition_output_dir(
    *, run_dir: Path, record: dict[str, Any], source_condition_id: str
) -> Path:
    """Resolve a condition directory after an archive is moved between hosts."""

    configured = Path(str(record["output_dir"])).expanduser()
    candidates = [configured]
    for condition_id in (source_condition_id, str(record["condition_id"])):
        candidates.extend(
            (run_dir / condition_id, run_dir / "conditions" / condition_id)
        )
    for candidate in candidates:
        if (candidate / "run_summary.json").is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        "Could not resolve synchronized condition output for "
        f"{source_condition_id} under {run_dir}"
    )


def _profile_key(profile: dict[str, Any]) -> tuple[Any, ...]:
    return tuple(profile.get(key) for key in ("year", "variant", "resolution"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _first_failure(summary: dict[str, Any]) -> tuple[str, float | str]:
    for gate_name in ("finite", "shape_body", "shape_nonbody"):
        gate = dict(summary.get("gates", {}).get(gate_name, {}) or {})
        if gate.get("any_fail", False):
            first_failure_at = gate.get("first_observed_fail_t_s")
            if first_failure_at is None:
                first_failure_at = gate.get(
                    "first_failure_t_s", gate.get("first_failure_step", "")
                )
            return (
                str(gate.get("first_failure_category") or gate_name),
                first_failure_at,
            )
    return ("", "")


def _row(
    record: dict[str, Any],
    source_row: dict[str, str],
    *,
    qc: dict[str, Any],
    stage: str,
) -> dict[str, Any]:
    output_dir = Path(str(record["output_dir"])).resolve()
    summary = _read_json(output_dir / "run_summary.json")
    performance = _read_json(output_dir / "performance.json")
    axes = dict(record.get("axis_values", {}) or {})
    time = dict(record.get("time", {}) or {})
    first_failure_category, first_failure_at = _first_failure(summary)
    values: dict[str, Any] = {
        "condition_id": str(record["condition_id"]),
        "n_flagella": int(axes["n_flagella"]),
        "attachment_pattern": str(axes.get("attachment_pattern", "")),
        "attachment_pattern_label": str(
            dict(record.get("axis_labels", {}) or {}).get("attachment_pattern", "")
        ),
        "attachment_slots": json.dumps(axes.get("attachment_slots", [])),
        "screen_status": screen_status(summary, qc=qc),
        "raw_nonbody_any_fail": _gate_failed(summary, "shape_nonbody"),
        "raw_first_failure_category": str(
            dict(summary.get("gates", {}).get("shape_nonbody", {}) or {}).get(
                "first_failure_category"
            )
            or ""
        ),
        "wall_time_s": float(performance["wall_time_s"]),
        "steps_per_s": float(performance["steps_per_s"]),
        "source_output_dir": str(output_dir),
        "source_condition_id": str(
            record.get("source_condition_id", record["condition_id"])
        ),
        "source_campaign": str(record.get("source_campaign", "")),
        "source_git_commit": str(record.get("source_git_commit", "")),
        "source_reused": bool(record.get("source_reused", False)),
        "first_failure_category": first_failure_category,
        "first_failure_at": first_failure_at,
    }
    if stage == "short_screen":
        values.update(
            {
                "torque_Nm": float(axes["motor_torque"]),
                "dt_star": float(time["dt_star"]),
                "dt_internal_s": float(time["dt_internal_s"]),
            }
        )
    elif stage == "long_duration":
        values.update(
            {
                "duration_s": float(time.get("duration_s", 0.0)),
                "full_ring_rotation_equivalent": int(axes["n_flagella"]) == 6,
                "run_summary_sha256": _sha256(output_dir / "run_summary.json"),
                "performance_sha256": _sha256(output_dir / "performance.json"),
                "state_archive_sha256": _sha256(output_dir / "state_archive.npz"),
            }
        )
    for metric, _ in SCREEN_METRICS[:6]:
        values[metric] = _maximum(summary, metric)
    values["local_attach_first_rel_err"] = _maximum(
        summary, "local_attach_first_rel_err"
    )
    values["flag_bend_err_max_deg"] = _maximum(summary, "flag_bend_err_max_deg")
    values["flag_torsion_err_max_deg"] = _maximum(summary, "flag_torsion_err_max_deg")
    values["source_completed"] = source_row.get("completed", "")
    return values


def _equal_values(left: Any, right: Any) -> bool:
    """Compare manifest values while tolerating JSON float representation."""

    if isinstance(left, dict) and isinstance(right, dict):
        return set(left) == set(right) and all(
            _equal_values(left[key], right[key]) for key in left
        )
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(
            _equal_values(a, b) for a, b in zip(left, right, strict=True)
        )
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return math.isclose(float(left), float(right), rel_tol=1e-12, abs_tol=1e-30)
    return left == right


def _require_completed_source(
    *, summary: dict[str, Any], source_row: dict[str, str], condition_id: str
) -> None:
    if not _bool(source_row.get("completed", "")):
        raise ValueError(f"Source condition is not completed: {condition_id}")
    execution = dict(summary.get("execution", {}) or {})
    if execution.get("status") != "completed":
        raise ValueError(f"Source run_summary is not completed: {condition_id}")


def collect_rows(
    *, config: dict[str, Any], run_dirs: list[Path]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Load completed cells and validate profile, provenance and grid coverage."""

    expected = _expected_conditions(config)
    contract = _development_contract(config)
    stage = str(contract["stage"])
    contract_axes = [str(axis) for axis in contract["axes"]]
    expected_by_key = {
        _record_key(condition, contract_axes): condition_id
        for condition_id, condition in expected.items()
    }
    if len(expected) != int(contract["expected_condition_count"]):
        raise ValueError(
            "development_evaluation.expected_condition_count does not match sweep"
        )
    expected_profile = dict(
        load_yaml(Path(str(config["base_config"]))).get("model_profile") or {}
    )
    expected_base_config = str(config["base_config"])
    accepted_campaigns = {str(item) for item in contract["accepted_source_campaigns"]}
    qc = dict(contract["qc"])
    records_by_id: dict[str, tuple[dict[str, Any], dict[str, str]]] = {}
    provenance: list[dict[str, Any]] = []
    for run_dir in run_dirs:
        manifest = _read_json(run_dir / "run_manifest.json")
        profile = dict(manifest.get("model_profile", {}) or {})
        if profile != expected_profile:
            raise ValueError(f"Model profile mismatch in {run_dir}: {profile}")
        if str(manifest.get("base_config") or "") != expected_base_config:
            raise ValueError(f"Base config mismatch in {run_dir}")
        campaign_config = str(manifest.get("campaign_config") or "")
        if campaign_config not in accepted_campaigns:
            raise ValueError(
                f"Unaccepted source campaign in {run_dir}: {campaign_config}"
            )
        git = _campaign_git_provenance(run_dir=run_dir, manifest=manifest)
        source_rows = _source_rows(run_dir)
        axis_defaults, relaxed_paths = _reused_source_spec(contract, campaign_config)
        provenance.append(
            {
                "run_dir": str(run_dir.resolve()),
                "git": git,
                "model_profile": profile,
                "source_campaign": campaign_config,
                "reused_axis_defaults": axis_defaults,
            }
        )
        for raw_record in manifest.get("conditions", []) or []:
            record = _normalize_reused_source_record(
                dict(raw_record), axis_defaults=axis_defaults
            )
            source_condition_id = str(record["condition_id"])
            try:
                condition_id = expected_by_key[_record_key(record, contract_axes)]
            except KeyError as error:
                raise ValueError(
                    f"Unexpected condition {source_condition_id} in {run_dir}"
                ) from error
            if condition_id in records_by_id:
                raise ValueError(
                    f"Duplicate condition {condition_id} across input runs"
                )
            if source_condition_id not in source_rows:
                raise ValueError(
                    f"Missing summary row for {source_condition_id} in {run_dir}"
                )
            expected_record = expected[condition_id]
            if not _equal_values(
                _without_override_paths(record.get("config_overrides"), relaxed_paths),
                _without_override_paths(
                    expected_record["config_overrides"], relaxed_paths
                ),
            ):
                raise ValueError(
                    f"Config override mismatch for {source_condition_id} in {run_dir}"
                )
            output_dir = _resolve_condition_output_dir(
                run_dir=run_dir,
                record=record,
                source_condition_id=source_condition_id,
            )
            source_summary = _read_json(output_dir / "run_summary.json")
            _require_completed_source(
                summary=source_summary,
                source_row=source_rows[source_condition_id],
                condition_id=source_condition_id,
            )
            if stage == "long_duration":
                final_archive = output_dir / "state_archive.npz"
                if (
                    not final_archive.is_file()
                    and (output_dir / PARTIAL_ARCHIVE_NAME).is_file()
                ):
                    raise ValueError(
                        "Partial checkpoint cannot be used as a long-duration "
                        f"source: {source_condition_id}"
                    )
                for name in LONG_DURATION_REQUIRED_ARTIFACTS:
                    if not (output_dir / name).is_file():
                        raise FileNotFoundError(
                            f"Missing required long-duration artifact for {source_condition_id}: {name}"
                        )
                expected_hashes = dict(record.get("artifact_sha256", {}) or {})
                for name, expected_hash in expected_hashes.items():
                    if _sha256(output_dir / name) != str(expected_hash):
                        raise ValueError(
                            f"SHA-256 mismatch for {source_condition_id}: {name}"
                        )
            canonical_record = dict(record)
            canonical_record["output_dir"] = str(output_dir)
            canonical_record["condition_id"] = condition_id
            canonical_record["source_condition_id"] = source_condition_id
            canonical_record["axis_values"] = dict(expected_record["axis_values"])
            canonical_record["axis_labels"] = dict(expected_record["axis_labels"])
            canonical_record["source_campaign"] = campaign_config
            canonical_record["source_git_commit"] = str(git["commit"])
            canonical_record["source_reused"] = bool(axis_defaults)
            records_by_id[condition_id] = (
                canonical_record,
                source_rows[source_condition_id],
            )
    missing = sorted(set(expected) - set(records_by_id))
    if missing:
        raise ValueError("Missing expected conditions: " + ", ".join(missing))
    rows = [
        _row(*records_by_id[condition_id], qc=qc, stage=stage)
        for condition_id in sorted(expected)
    ]
    rows.sort(
        key=(
            (lambda row: (row["n_flagella"], row["torque_Nm"], row["dt_star"]))
            if stage == "short_screen"
            else (lambda row: (row["n_flagella"], row["attachment_pattern"]))
        )
    )
    return rows, provenance


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _plot_count(rows: list[dict[str, Any]], output_path: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    torques = sorted({float(row["torque_Nm"]) for row in rows})
    dts = sorted({float(row["dt_star"]) for row in rows})
    by_cell = {(float(row["torque_Nm"]), float(row["dt_star"])): row for row in rows}
    panels = (("screen_status", "screen status"),) + SCREEN_METRICS
    figure, axes = plt.subplots(3, 3, figsize=(15, 13), constrained_layout=True)
    for axis, (metric, label) in zip(axes.flat, panels, strict=False):
        matrix = np.full((len(torques), len(dts)), np.nan)
        for row_index, torque in enumerate(torques):
            for column_index, dt in enumerate(dts):
                value = by_cell[(torque, dt)][metric]
                matrix[row_index, column_index] = (
                    {"fail": 0.0, "pass": 1.0}[value]
                    if metric == "screen_status"
                    else float(value)
                )
        if metric == "screen_status":
            image = axis.imshow(
                matrix,
                aspect="auto",
                cmap=ListedColormap(["#c53d3d", "#278b6e"]),
                norm=BoundaryNorm([-0.5, 0.5, 1.5], 2),
            )
            colorbar = figure.colorbar(image, ax=axis, shrink=0.8)
            colorbar.set_ticks([0, 1], labels=["FAIL", "PASS"])
        else:
            image = axis.imshow(matrix, aspect="auto", cmap="viridis")
            figure.colorbar(image, ax=axis, shrink=0.8)
        axis.set_title(label, fontsize=10)
        axis.set_xticks(range(len(dts)), [f"{dt:.0e}" for dt in dts])
        axis.set_yticks(
            range(len(torques)), [f"{torque / 1e-20:g}" for torque in torques]
        )
        axis.set_xlabel("dt_star")
        axis.set_ylabel("torque [1e-20 N m / flagellum]")
        for (row_index, column_index), value in np.ndenumerate(matrix):
            text = (
                {0.0: "FAIL", 1.0: "PASS"}[value]
                if metric == "screen_status"
                else f"{value:.3g}"
            )
            axis.text(
                column_index, row_index, text, ha="center", va="center", fontsize=8
            )
    axis = axes.flat[-1]
    axis.axis("off")
    figure.suptitle(
        f"Model-development 1tau screen (n_flagella={rows[0]['n_flagella']})"
    )
    figure.savefig(output_path, dpi=220)
    plt.close(figure)


def _plot_attachment_patterns(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Render a sparse n-by-canonical-slot-pattern QC heatmap."""

    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    patterns = sorted(
        {str(row["attachment_pattern"]) for row in rows},
        key=lambda pattern: (
            next(
                int(row["n_flagella"])
                for row in rows
                if row["attachment_pattern"] == pattern
            ),
            pattern,
        ),
    )
    labels = {
        str(row["attachment_pattern"]): str(row["attachment_pattern_label"])
        for row in rows
    }
    counts = sorted({int(row["n_flagella"]) for row in rows})
    by_cell = {
        (int(row["n_flagella"]), str(row["attachment_pattern"])): row for row in rows
    }
    panels = (("screen_status", "screen status"),) + SCREEN_METRICS
    figure, axes = plt.subplots(3, 3, figsize=(22, 13), constrained_layout=True)
    for axis, (metric, label) in zip(axes.flat, panels, strict=False):
        matrix = np.full((len(counts), len(patterns)), np.nan)
        for row_index, count in enumerate(counts):
            for column_index, pattern in enumerate(patterns):
                row = by_cell.get((count, pattern))
                if row is None:
                    continue
                value = row[metric]
                matrix[row_index, column_index] = (
                    {"fail": 0.0, "pass": 1.0}[value]
                    if metric == "screen_status"
                    else float(value)
                )
        if metric == "screen_status":
            cmap = ListedColormap(["#c53d3d", "#278b6e"])
            cmap.set_bad("#e5e7eb")
            image = axis.imshow(
                matrix,
                aspect="auto",
                cmap=cmap,
                norm=BoundaryNorm([-0.5, 0.5, 1.5], 2),
            )
            colorbar = figure.colorbar(image, ax=axis, shrink=0.8)
            colorbar.set_ticks([0, 1], labels=["FAIL", "PASS"])
        else:
            cmap = plt.get_cmap("viridis").copy()
            cmap.set_bad("#e5e7eb")
            image = axis.imshow(matrix, aspect="auto", cmap=cmap)
            figure.colorbar(image, ax=axis, shrink=0.8)
        axis.set_title(label, fontsize=10)
        axis.set_xticks(range(len(patterns)), [labels[pattern] for pattern in patterns])
        axis.set_yticks(range(len(counts)), [f"n={count}" for count in counts])
        axis.tick_params(axis="x", labelrotation=60, labelsize=7)
        axis.set_xlabel("canonical attachment slots (C6 rotations identified)")
        axis.set_ylabel("n_flagella")
        for (row_index, column_index), value in np.ndenumerate(matrix):
            if not math.isfinite(float(value)):
                continue
            text = (
                {0.0: "FAIL", 1.0: "PASS"}[value]
                if metric == "screen_status"
                else f"{value:.3g}"
            )
            axis.text(
                column_index, row_index, text, ha="center", va="center", fontsize=7
            )
    axes.flat[-1].axis("off")
    figure.suptitle("Model-development attachment-pattern QC (masked = not applicable)")
    figure.savefig(output_path, dpi=220)
    plt.close(figure)


def _write_replay_input(
    *, output_dir: Path, run_dirs: list[Path], config: dict[str, Any]
) -> Path:
    """Build a replay-only manifest retaining absolute source archives."""

    replay_input = output_dir / "replay_input"
    replay_input.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    summary_rows: list[dict[str, str]] = []
    base_config: str | None = None
    contract = _development_contract(config)
    contract_axes = [str(axis) for axis in contract["axes"]]
    expected = _expected_conditions(config)
    expected_by_key = {
        _record_key(condition, contract_axes): condition_id
        for condition_id, condition in expected.items()
    }
    for run_dir in run_dirs:
        manifest = _read_json(run_dir / "run_manifest.json")
        campaign_config = str(manifest.get("campaign_config") or "")
        axis_defaults, _ = _reused_source_spec(contract, campaign_config)
        git = _campaign_git_provenance(run_dir=run_dir, manifest=manifest)
        if base_config is None:
            source_config = manifest.get("source_config_path")
            if source_config is not None and Path(str(source_config)).is_file():
                base_config = str(source_config)
            else:
                base_config = str(manifest.get("base_config") or "")
        source_rows = _source_rows(run_dir)
        for raw_record in manifest.get("conditions", []) or []:
            record = _normalize_reused_source_record(
                dict(raw_record), axis_defaults=axis_defaults
            )
            source_condition_id = str(record["condition_id"])
            condition_id = expected_by_key[_record_key(record, contract_axes)]
            record["condition_id"] = condition_id
            record["source_condition_id"] = source_condition_id
            record["axis_values"] = dict(expected[condition_id]["axis_values"])
            record["axis_labels"] = dict(expected[condition_id]["axis_labels"])
            record["source_campaign"] = campaign_config
            record["source_git_commit"] = str(git.get("commit") or "")
            record["source_reused"] = bool(axis_defaults)
            source_output_dir = _resolve_condition_output_dir(
                run_dir=run_dir,
                record=record,
                source_condition_id=source_condition_id,
            )
            record["output_dir"] = str(source_output_dir)
            geometry_path = source_output_dir / "initial_geometry_summary.json"
            if geometry_path.is_file():
                geometry = dict(_read_json(geometry_path).get("geometry", {}) or {})
                if geometry:
                    record["geometry"] = {"actual": geometry}
            records.append(record)
            summary_row = dict(source_rows[source_condition_id])
            summary_row["condition_id"] = str(record["condition_id"])
            summary_rows.append(summary_row)
    records.sort(key=lambda record: str(record["condition_id"]))
    summary_rows.sort(key=lambda row: str(row["condition_id"]))
    _write_csv(replay_input / "summary.csv", summary_rows)
    (replay_input / "run_manifest.json").write_text(
        json.dumps(
            {
                "kind": "model_development_replay_input",
                "base_config": base_config
                or str(records[0].get("source_config_path") or ""),
                "condition_order": [record["condition_id"] for record in records],
                "conditions": records,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return replay_input


def _render_replays(
    rows: list[dict[str, Any]], *, replay_input: Path, output_dir: Path, stage: str
) -> None:
    from sim_swim.analysis.phase2_replay import main as replay_main

    groups: list[tuple[str, list[str]]]
    if rows and rows[0].get("attachment_pattern"):
        groups = [
            (
                f"nf{n_flagella:02d}/attachment_patterns",
                [
                    row["condition_id"]
                    for row in rows
                    if int(row["n_flagella"]) == n_flagella
                ],
            )
            for n_flagella in sorted({int(row["n_flagella"]) for row in rows})
        ]
    elif stage == "short_screen":
        groups = [
            (
                f"dt{dt_star:.0e}/nf{n_flagella:02d}",
                [
                    row["condition_id"]
                    for row in rows
                    if float(row["dt_star"]) == dt_star
                    and int(row["n_flagella"]) == n_flagella
                ],
            )
            for dt_star in sorted({float(row["dt_star"]) for row in rows})
            for n_flagella in sorted({int(row["n_flagella"]) for row in rows})
        ]
    elif stage == "long_duration":
        groups = [
            (
                f"nf{int(row['n_flagella']):02d}/as{int(row['attach_seed']):03d}_ps{int(row['phase_seed']):03d}",
                [row["condition_id"]],
            )
            for row in rows
        ]
    for group_name, selected in groups:
        destination = output_dir / "replay" / group_name
        args = [
            "--input-dir",
            str(replay_input),
            "--camera-envelope-input-dir",
            str(replay_input),
            "--output-dir",
            str(destination),
            "--view",
            "3d+2d",
            "--mode",
            "render-only",
            "--camera-3d",
            "fixed",
            "--camera-2d",
            "fixed",
            "--view-range-mode",
            "campaign-envelope",
            "--target-frame-count",
            "41",
            "--max-panels-per-grid",
            "5",
            "--overwrite",
        ]
        for condition_id in selected:
            args.extend(["--condition-id", condition_id])
        replay_main(args)


def _body_axis_and_roll(
    beads_um: np.ndarray, *, n_prism: int
) -> tuple[np.ndarray, float]:
    """Return the body axis and a reproducible roll marker angle from archive beads."""

    layers = beads_um.shape[0] // n_prism
    first = beads_um[:n_prism]
    last = beads_um[(layers - 1) * n_prism : layers * n_prism]
    axis = np.mean(last, axis=0) - np.mean(first, axis=0)
    axis /= max(float(np.linalg.norm(axis)), 1e-18)
    reference = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(axis, reference))) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(axis, reference)
    e1 /= max(float(np.linalg.norm(e1)), 1e-18)
    e2 = np.cross(axis, e1)
    marker = first[0] - np.mean(first, axis=0)
    return axis, float(np.arctan2(np.dot(marker, e2), np.dot(marker, e1)))


def _stall_rows(
    *, rows: list[dict[str, Any]], config: dict[str, Any]
) -> list[dict[str, Any]]:
    """Derive configured stall candidates from completed long-duration archives."""

    diagnostic = dict(
        config.get("development_evaluation", {}).get("stall_diagnostic", {}) or {}
    )
    if not diagnostic:
        return []
    window_s = float(diagnostic["window_ms"]) / 1000.0
    threshold = float(diagnostic["relative_median_threshold"])
    base = load_yaml(Path(str(config["base_config"])))
    n_prism = int(base["body"]["prism"]["n_prism"])
    result: list[dict[str, Any]] = []
    for row in rows:
        states = load_state_archive(
            Path(str(row["source_output_dir"])) / "state_archive.npz"
        )
        times = np.asarray([state.t for state in states], dtype=float)
        speeds = np.asarray([np.linalg.norm(state.velocity_um_s) for state in states])
        axes, phases = zip(
            *[
                _body_axis_and_roll(
                    np.asarray(state.bead_positions_um)[
                        : int(base["model_profile"]["body_beads"])
                    ],
                    n_prism=n_prism,
                )
                for state in states
            ],
            strict=True,
        )
        phases_unwrapped = np.unwrap(np.asarray(phases, dtype=float))
        roll_hz = np.zeros_like(times)
        if len(times) > 1:
            roll_hz[1:] = np.abs(np.diff(phases_unwrapped) / np.diff(times)) / (
                2.0 * np.pi
            )
        speed_median = float(np.nanmedian(speeds[1:]))
        roll_median = float(np.nanmedian(roll_hz[1:]))
        start = float(times[0])
        while start < float(times[-1]):
            mask = (times >= start) & (times < start + window_s)
            if np.count_nonzero(mask) < 2:
                start += window_s
                continue
            indices = np.flatnonzero(mask)
            axis_angle = math.degrees(
                math.acos(
                    float(
                        np.clip(np.dot(axes[indices[0]], axes[indices[-1]]), -1.0, 1.0)
                    )
                )
            )
            mean_speed = float(np.mean(speeds[mask]))
            mean_roll = float(np.mean(roll_hz[mask]))
            is_stall = bool(
                speed_median > 0.0
                and roll_median > 0.0
                and mean_speed < threshold * speed_median
                and mean_roll < threshold * roll_median
            )
            result.append(
                {
                    "condition_id": row["condition_id"],
                    "window_start_s": start,
                    "window_end_s": min(start + window_s, float(times[-1])),
                    "mean_body_speed_um_s": mean_speed,
                    "mean_body_roll_hz": mean_roll,
                    "body_axis_angle_change_deg": axis_angle,
                    "speed_median_um_s": speed_median,
                    "roll_median_hz": roll_median,
                    "stall_candidate": is_stall,
                }
            )
            start += window_s
    return result


def _on_off_comparison_rows(
    *, rows: list[dict[str, Any]], stall_rows: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """Summarize a complete body--flagella-repulsion ON/OFF pairing."""

    arms: dict[str, dict[str, dict[str, Any]]] = {}
    windows_by_id: dict[str, list[dict[str, Any]]] = {}
    for window in stall_rows:
        windows_by_id.setdefault(str(window["condition_id"]), []).append(window)
    for row in rows:
        condition_id = str(row["condition_id"])
        for label in ("bfon", "bfoff"):
            suffix = f"__{label}"
            if condition_id.endswith(suffix):
                pair_id = condition_id.removesuffix(suffix)
                arms.setdefault(pair_id, {})[label] = row
                break
    comparisons: list[dict[str, Any]] = []
    for pair_id, pair in sorted(arms.items()):
        if set(pair) != {"bfon", "bfoff"}:
            raise ValueError(f"Incomplete ON/OFF pair: {pair_id}")
        values: dict[str, Any] = {
            "attachment_pattern": pair_id,
            "n_flagella": int(pair["bfon"]["n_flagella"]),
            "attachment_slots": pair["bfon"]["attachment_slots"],
        }
        for label, prefix in (("bfon", "on"), ("bfoff", "off")):
            row = pair[label]
            windows = windows_by_id.get(str(row["condition_id"]), [])
            if not windows:
                raise ValueError(f"Missing stall windows for {row['condition_id']}")
            stall_windows = [item for item in windows if _bool(item["stall_candidate"])]
            values.update(
                {
                    f"{prefix}_condition_id": row["condition_id"],
                    f"{prefix}_source_campaign": row["source_campaign"],
                    f"{prefix}_git_commit": row["source_git_commit"],
                    f"{prefix}_state_archive_sha256": row["state_archive_sha256"],
                    f"{prefix}_speed_median_um_s": float(
                        windows[0]["speed_median_um_s"]
                    ),
                    f"{prefix}_roll_median_hz": float(windows[0]["roll_median_hz"]),
                    f"{prefix}_axis_step_angle_median_deg": float(
                        np.median(
                            [
                                float(item["body_axis_angle_change_deg"])
                                for item in windows
                            ]
                        )
                    ),
                    f"{prefix}_stall_window_count": len(stall_windows),
                    f"{prefix}_stall_total_s": sum(
                        float(item["window_end_s"]) - float(item["window_start_s"])
                        for item in stall_windows
                    ),
                    f"{prefix}_stall_fraction": len(stall_windows) / len(windows),
                    f"{prefix}_strict_status": row["screen_status"],
                }
            )
        for metric in (
            "speed_median_um_s",
            "roll_median_hz",
            "axis_step_angle_median_deg",
            "stall_total_s",
            "stall_fraction",
        ):
            values[f"off_minus_on_{metric}"] = float(values[f"off_{metric}"]) - float(
                values[f"on_{metric}"]
            )
        comparisons.append(values)
    if len(comparisons) * 2 != len(rows):
        raise ValueError("Every evaluated condition must belong to one ON/OFF pair")
    return comparisons


def _plot_on_off_comparison(rows: list[dict[str, Any]], output_path: Path) -> None:
    """Plot attachment-level OFF-minus-ON diagnostics with zero-centered scales."""

    import matplotlib.pyplot as plt

    metrics = (
        ("off_minus_on_speed_median_um_s", "Δ body speed [µm/s]"),
        ("off_minus_on_roll_median_hz", "Δ body roll [Hz]"),
        ("off_minus_on_stall_total_s", "Δ stall time [s]"),
        ("off_minus_on_axis_step_angle_median_deg", "Δ axis step angle [deg]"),
    )
    values = np.asarray(
        [[float(row[key]) for key, _ in metrics] for row in rows], dtype=float
    )
    limits = np.maximum(np.max(np.abs(values), axis=0), 1.0e-12)
    figure, axes = plt.subplots(
        1, len(metrics), figsize=(14, max(4.5, len(rows) * 0.34))
    )
    labels = [f"n={row['n_flagella']} {row['attachment_pattern']}" for row in rows]
    for index, (axis, (_, title)) in enumerate(zip(axes, metrics, strict=True)):
        image = axis.imshow(
            values[:, [index]], cmap="coolwarm", vmin=-limits[index], vmax=limits[index]
        )
        axis.set_title(title, fontsize=9)
        axis.set_xticks([])
        axis.set_yticks(range(len(rows)))
        axis.set_yticklabels(labels if index == 0 else [], fontsize=7)
        figure.colorbar(image, ax=axis, fraction=0.08, pad=0.04)
    figure.suptitle(
        "Body--flagella repulsion OFF − ON (diagnostic-only; strict QC remains FAIL)",
        fontsize=11,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=220)
    plt.close(figure)


def _render_on_off_pair_replays(
    comparison_rows: list[dict[str, Any]], *, replay_input: Path, output_dir: Path
) -> None:
    """Render each ON/OFF pair using the same fixed cameras and sampling."""

    from sim_swim.analysis.phase2_replay import main as replay_main

    for row in comparison_rows:
        destination = output_dir / "pair_replay" / str(row["attachment_pattern"])
        replay_main(
            [
                "--input-dir",
                str(replay_input),
                "--camera-envelope-input-dir",
                str(replay_input),
                "--output-dir",
                str(destination),
                "--view",
                "3d+2d",
                "--mode",
                "render-only",
                "--camera-3d",
                "fixed",
                "--camera-2d",
                "fixed",
                "--view-range-mode",
                "campaign-envelope",
                "--target-frame-count",
                "41",
                "--max-panels-per-grid",
                "2",
                "--overwrite",
                "--condition-id",
                str(row["on_condition_id"]),
                "--condition-id",
                str(row["off_condition_id"]),
            ]
        )


def _output_hashes(output_dir: Path) -> dict[str, str]:
    return {
        str(path.relative_to(output_dir)): _sha256(path)
        for path in sorted(output_dir.rglob("*"))
        if path.is_file() and path.name not in OPERATIONAL_LOG_NAMES
    }


def build_evaluation(
    *,
    config_path: Path,
    run_dirs: list[Path],
    output_dir: Path,
    render_replay: bool = False,
    dry_run: bool = False,
) -> dict[str, Path]:
    config = load_yaml(config_path)
    contract = _development_contract(config)
    stage = str(contract["stage"])
    rows, provenance = collect_rows(config=config, run_dirs=run_dirs)
    if dry_run:
        return {"validated": config_path}
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "summary.csv"
    _write_csv(summary_path, rows)
    outputs: dict[str, Path] = {"summary_csv": summary_path}
    if rows and rows[0].get("attachment_pattern"):
        heatmap_dir = output_dir / "heatmaps"
        heatmap_dir.mkdir(exist_ok=True)
        path = heatmap_dir / "attachment_patterns_qc.png"
        _plot_attachment_patterns(rows, path)
        outputs["attachment_pattern_heatmap"] = path
    elif stage == "short_screen":
        heatmap_dir = output_dir / "heatmaps"
        heatmap_dir.mkdir(exist_ok=True)
        for n_flagella in sorted({int(row["n_flagella"]) for row in rows}):
            path = heatmap_dir / f"nf{n_flagella:02d}_screen.png"
            _plot_count(
                [row for row in rows if int(row["n_flagella"]) == n_flagella], path
            )
            outputs[f"heatmap_nf{n_flagella:02d}"] = path
    if stage == "long_duration":
        window_rows = [
            {
                "condition_id": row["condition_id"],
                "window_start_s": 0.0,
                "window_end_s": row["duration_s"],
                "strict_status": row["screen_status"],
                "first_failure_category": row["first_failure_category"],
                "first_failure_at": row["first_failure_at"],
            }
            for row in rows
        ]
        window_path = output_dir / "window_qc.csv"
        _write_csv(window_path, window_rows)
        outputs["window_qc_csv"] = window_path
        stall_rows = _stall_rows(rows=rows, config=config)
        if stall_rows:
            stall_path = output_dir / "stall_summary.csv"
            _write_csv(stall_path, stall_rows)
            outputs["stall_summary_csv"] = stall_path
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(2, 1, sharex=True, figsize=(9, 5))
            for condition_id in sorted(
                {str(item["condition_id"]) for item in stall_rows}
            ):
                selected = [
                    item for item in stall_rows if item["condition_id"] == condition_id
                ]
                x = [float(item["window_start_s"]) for item in selected]
                axes[0].plot(
                    x,
                    [item["mean_body_speed_um_s"] for item in selected],
                    label=condition_id,
                )
                axes[1].plot(
                    x,
                    [item["mean_body_roll_hz"] for item in selected],
                    label=condition_id,
                )
            axes[0].set_ylabel("body speed [µm/s]")
            axes[1].set_ylabel("body roll [Hz]")
            axes[1].set_xlabel("time [s]")
            axes[0].legend(fontsize=6, ncol=2)
            fig.tight_layout()
            stall_plot = output_dir / "stall_timeseries.png"
            fig.savefig(stall_plot, dpi=150)
            plt.close(fig)
            outputs["stall_timeseries_png"] = stall_plot
        if stall_rows and {"bfon", "bfoff"}.issubset(
            {str(row["condition_id"]).split("__")[-1] for row in rows}
        ):
            comparison_rows = _on_off_comparison_rows(rows=rows, stall_rows=stall_rows)
            comparison_path = output_dir / "comparison_summary.csv"
            _write_csv(comparison_path, comparison_rows)
            outputs["comparison_summary_csv"] = comparison_path
            comparison_heatmap = heatmap_dir / "on_off_comparison.png"
            _plot_on_off_comparison(comparison_rows, comparison_heatmap)
            outputs["on_off_comparison_heatmap"] = comparison_heatmap
    replay_input = _write_replay_input(
        output_dir=output_dir, run_dirs=run_dirs, config=config
    )
    outputs["replay_input"] = replay_input
    if rows and rows[0].get("attachment_pattern"):
        from sim_swim.analysis.attachment_slot_map import render_attachment_slot_map

        outputs.update(
            render_attachment_slot_map(replay_input, output_dir / "attachment_slots")
        )
    if render_replay:
        _render_replays(
            rows, replay_input=replay_input, output_dir=output_dir, stage=stage
        )
        outputs["replay"] = output_dir / "replay"
        if stage == "long_duration" and "comparison_rows" in locals():
            _render_on_off_pair_replays(
                comparison_rows, replay_input=replay_input, output_dir=output_dir
            )
            outputs["pair_replay"] = output_dir / "pair_replay"
    if "comparison_rows" in locals():
        visualization_path = output_dir / "visualization_manifest.json"
        visualization_path.write_text(
            json.dumps(
                {
                    "kind": "on_off_visualization_bundle",
                    "condition_pair_count": len(comparison_rows),
                    "comparison_rows": comparison_rows,
                    "strict_qc_note": "diagnostic-only; strict QC status is not changed",
                    "outputs": {key: str(value) for key, value in outputs.items()},
                    "sha256": _output_hashes(output_dir),
                },
                ensure_ascii=False,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        outputs["visualization_manifest"] = visualization_path
    manifest = {
        "kind": "model_development_evaluation",
        "config": str(config_path),
        "stage": stage,
        "condition_count": len(rows),
        "status_counts": {
            status: sum(row["screen_status"] == status for row in rows)
            for status in ("pass", "fail")
        },
        "qc_policy": "hook_angle_err_max_deg is diagnostic-only; all other required QC remains PASS/FAIL",
        "full_ring_rotation_equivalent_condition_ids": [
            row["condition_id"]
            for row in rows
            if row.get("full_ring_rotation_equivalent")
        ],
        "attachment_patterns": [
            {
                "condition_id": row["condition_id"],
                "n_flagella": row["n_flagella"],
                "attachment_pattern": row["attachment_pattern"],
                "attachment_slots": json.loads(row["attachment_slots"]),
            }
            for row in rows
            if row.get("attachment_pattern")
        ],
        "long_duration_artifact_policy": (
            {
                "required_completed_artifacts": list(LONG_DURATION_REQUIRED_ARTIFACTS),
                "partial_checkpoint_artifacts": [
                    "progress.json",
                    "diagnostic_samples.csv",
                    PARTIAL_ARCHIVE_NAME,
                    "trajectory.partial.csv",
                ],
                "partial_checkpoint_policy": "diagnostic-only; never a completed archive or acceptance input",
                "excluded_operational_logs": list(OPERATIONAL_LOG_NAMES),
            }
            if stage == "long_duration"
            else None
        ),
        "provenance": provenance,
        "outputs": {key: str(value) for key, value in outputs.items()},
    }
    manifest_text = json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    manifest_path = output_dir / "evaluation_manifest.json"
    manifest_path.write_text(manifest_text, encoding="utf-8")
    (output_dir / "manifest.json").write_text(manifest_text, encoding="utf-8")
    (output_dir / "run.log").write_text(
        "model-development-evaluation completed\n"
        f"config={config_path}\n"
        f"condition_count={len(rows)}\n"
        f"input_run_dirs=" + ",".join(str(path.resolve()) for path in run_dirs) + "\n",
        encoding="utf-8",
    )
    outputs["manifest"] = manifest_path
    return outputs


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--render-replay", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    outputs = build_evaluation(
        config_path=args.config,
        run_dirs=args.run_dir,
        output_dir=args.output_dir,
        render_replay=args.render_replay,
        dry_run=args.dry_run,
    )
    for path in outputs.values():
        print(path)


if __name__ == "__main__":
    main()
