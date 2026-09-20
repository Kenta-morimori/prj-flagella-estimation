"""Aggregate reusable Phase 2 model-development evaluation screens.

This module deliberately separates numerical/physical QC from downstream
swimming-feature analysis.  It reads completed multi-run campaigns only; it
never starts simulations.
"""

from __future__ import annotations

import argparse
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
        values.append(
            float(value) if name in {"motor_torque", "dt_star"} else int(value)
        )
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


def _source_rows(run_dir: Path) -> dict[str, dict[str, str]]:
    with (run_dir / "summary.csv").open(encoding="utf-8", newline="") as handle:
        return {str(row["condition_id"]): row for row in csv.DictReader(handle)}


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
            return (
                str(gate.get("first_failure_category") or gate_name),
                gate.get("first_failure_t_s", gate.get("first_failure_step", "")),
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
    else:
        values.update(
            {
                "attach_seed": int(axes["attach_seed"]),
                "phase_seed": int(axes["phase_seed"]),
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
        git = dict(manifest.get("git", {}) or {})
        if not str(git.get("commit") or "") or git.get("is_clean") is not True:
            raise ValueError(f"Invalid Git provenance in {run_dir}")
        source_rows = _source_rows(run_dir)
        provenance.append(
            {
                "run_dir": str(run_dir.resolve()),
                "git": dict(manifest.get("git", {}) or {}),
                "model_profile": profile,
            }
        )
        for raw_record in manifest.get("conditions", []) or []:
            record = dict(raw_record)
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
                record.get("config_overrides"), expected_record["config_overrides"]
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
                for name in (
                    "run_summary.json",
                    "performance.json",
                    "state_archive.npz",
                ):
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
            else (
                lambda row: (row["n_flagella"], row["attach_seed"], row["phase_seed"])
            )
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
    expected_by_key = {
        _record_key(condition, contract_axes): condition_id
        for condition_id, condition in _expected_conditions(config).items()
    }
    for run_dir in run_dirs:
        manifest = _read_json(run_dir / "run_manifest.json")
        if base_config is None:
            source_config = manifest.get("source_config_path") or manifest.get(
                "base_config"
            )
            if source_config is not None:
                base_config = str(source_config)
        source_rows = _source_rows(run_dir)
        for raw_record in manifest.get("conditions", []) or []:
            record = dict(raw_record)
            source_condition_id = str(record["condition_id"])
            record["condition_id"] = expected_by_key[_record_key(record, contract_axes)]
            record["source_condition_id"] = source_condition_id
            record["output_dir"] = str(
                _resolve_condition_output_dir(
                    run_dir=run_dir,
                    record=record,
                    source_condition_id=source_condition_id,
                )
            )
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
    if stage == "short_screen":
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
    else:
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
    if stage == "short_screen":
        heatmap_dir = output_dir / "heatmaps"
        heatmap_dir.mkdir(exist_ok=True)
        for n_flagella in sorted({int(row["n_flagella"]) for row in rows}):
            path = heatmap_dir / f"nf{n_flagella:02d}_screen.png"
            _plot_count(
                [row for row in rows if int(row["n_flagella"]) == n_flagella], path
            )
            outputs[f"heatmap_nf{n_flagella:02d}"] = path
    else:
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
    replay_input = _write_replay_input(
        output_dir=output_dir, run_dirs=run_dirs, config=config
    )
    outputs["replay_input"] = replay_input
    if render_replay:
        _render_replays(
            rows, replay_input=replay_input, output_dir=output_dir, stage=stage
        )
        outputs["replay"] = output_dir / "replay"
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
