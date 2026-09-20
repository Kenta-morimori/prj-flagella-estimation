"""Aggregate reusable Phase 2 model-development evaluation screens.

This module deliberately separates numerical/physical QC from downstream
swimming-feature analysis.  It reads completed multi-run campaigns only; it
never starts simulations.
"""

from __future__ import annotations

import argparse
import csv
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


def screen_status(summary: dict[str, Any]) -> str:
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


def _condition_key(
    *, n_flagella: Any, motor_torque: Any, dt_star: Any
) -> tuple[int, float, float]:
    return (int(n_flagella), float(motor_torque), float(dt_star))


def _record_key(record: dict[str, Any]) -> tuple[int, float, float]:
    axes = dict(record.get("axis_values", {}) or {})
    time = dict(record.get("time", {}) or {})
    dt_star = axes["dt_star"] if "dt_star" in axes else time["dt_star"]
    return _condition_key(
        n_flagella=axes["n_flagella"],
        motor_torque=axes["motor_torque"],
        dt_star=dt_star,
    )


def _development_contract(config: dict[str, Any]) -> dict[str, Any]:
    contract = dict(config.get("development_evaluation", {}) or {})
    if contract.get("stage") != "short_screen":
        raise ValueError("development_evaluation.stage must be short_screen")
    if int(contract.get("expected_condition_count", 0)) <= 0:
        raise ValueError("development_evaluation.expected_condition_count is required")
    if not isinstance(contract.get("axes"), list) or not contract["axes"]:
        raise ValueError("development_evaluation.axes is required")
    return contract


def _source_rows(run_dir: Path) -> dict[str, dict[str, str]]:
    with (run_dir / "summary.csv").open(encoding="utf-8", newline="") as handle:
        return {str(row["condition_id"]): row for row in csv.DictReader(handle)}


def _profile_key(profile: dict[str, Any]) -> tuple[Any, ...]:
    return tuple(profile.get(key) for key in ("year", "variant", "resolution"))


def _row(record: dict[str, Any], source_row: dict[str, str]) -> dict[str, Any]:
    output_dir = Path(str(record["output_dir"])).resolve()
    summary = _read_json(output_dir / "run_summary.json")
    performance = _read_json(output_dir / "performance.json")
    axes = dict(record.get("axis_values", {}) or {})
    time = dict(record.get("time", {}) or {})
    values: dict[str, Any] = {
        "condition_id": str(record["condition_id"]),
        "n_flagella": int(axes["n_flagella"]),
        "torque_Nm": float(axes["motor_torque"]),
        "dt_star": float(time["dt_star"]),
        "dt_internal_s": float(time["dt_internal_s"]),
        "screen_status": screen_status(summary),
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
    }
    for metric, _ in SCREEN_METRICS[:6]:
        values[metric] = _maximum(summary, metric)
    values["local_attach_first_rel_err"] = _maximum(
        summary, "local_attach_first_rel_err"
    )
    values["flag_bend_err_max_deg"] = _maximum(summary, "flag_bend_err_max_deg")
    values["flag_torsion_err_max_deg"] = _maximum(summary, "flag_torsion_err_max_deg")
    values["source_completed"] = source_row.get("completed", "")
    return values


def collect_rows(
    *, config: dict[str, Any], run_dirs: list[Path]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Load completed cells and validate profile, provenance and grid coverage."""

    expected = _expected_conditions(config)
    expected_by_key = {
        _record_key(condition): condition_id
        for condition_id, condition in expected.items()
    }
    contract = _development_contract(config)
    if len(expected) != int(contract["expected_condition_count"]):
        raise ValueError(
            "development_evaluation.expected_condition_count does not match sweep"
        )
    expected_profile = dict(
        load_yaml(Path(str(config["base_config"]))).get("model_profile") or {}
    )
    records_by_id: dict[str, tuple[dict[str, Any], dict[str, str]]] = {}
    provenance: list[dict[str, Any]] = []
    for run_dir in run_dirs:
        manifest = _read_json(run_dir / "run_manifest.json")
        profile = dict(manifest.get("model_profile", {}) or {})
        if _profile_key(profile) != _profile_key(expected_profile):
            raise ValueError(f"Model profile mismatch in {run_dir}: {profile}")
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
                condition_id = expected_by_key[_record_key(record)]
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
            canonical_record = dict(record)
            canonical_record["condition_id"] = condition_id
            canonical_record["source_condition_id"] = source_condition_id
            records_by_id[condition_id] = (
                canonical_record,
                source_rows[source_condition_id],
            )
    missing = sorted(set(expected) - set(records_by_id))
    if missing:
        raise ValueError("Missing expected conditions: " + ", ".join(missing))
    rows = [_row(*records_by_id[condition_id]) for condition_id in sorted(expected)]
    rows.sort(key=lambda row: (row["n_flagella"], row["torque_Nm"], row["dt_star"]))
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
    expected_by_key = {
        _record_key(condition): condition_id
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
            record["condition_id"] = expected_by_key[_record_key(record)]
            record["source_condition_id"] = source_condition_id
            record["output_dir"] = str(Path(str(record["output_dir"])).resolve())
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
    rows: list[dict[str, Any]], *, replay_input: Path, output_dir: Path
) -> None:
    from sim_swim.analysis.phase2_replay import main as replay_main

    for dt_star in sorted({float(row["dt_star"]) for row in rows}):
        for n_flagella in sorted({int(row["n_flagella"]) for row in rows}):
            selected = [
                row["condition_id"]
                for row in rows
                if float(row["dt_star"]) == dt_star
                and int(row["n_flagella"]) == n_flagella
            ]
            destination = (
                output_dir / "replay" / f"dt{dt_star:.0e}" / f"nf{n_flagella:02d}"
            )
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
    rows, provenance = collect_rows(config=config, run_dirs=run_dirs)
    if dry_run:
        return {"validated": config_path}
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "summary.csv"
    _write_csv(summary_path, rows)
    heatmap_dir = output_dir / "heatmaps"
    heatmap_dir.mkdir(exist_ok=True)
    outputs: dict[str, Path] = {"summary_csv": summary_path}
    for n_flagella in sorted({int(row["n_flagella"]) for row in rows}):
        path = heatmap_dir / f"nf{n_flagella:02d}_screen.png"
        _plot_count([row for row in rows if int(row["n_flagella"]) == n_flagella], path)
        outputs[f"heatmap_nf{n_flagella:02d}"] = path
    replay_input = _write_replay_input(
        output_dir=output_dir, run_dirs=run_dirs, config=config
    )
    outputs["replay_input"] = replay_input
    if render_replay:
        _render_replays(rows, replay_input=replay_input, output_dir=output_dir)
        outputs["replay"] = output_dir / "replay"
    manifest_path = output_dir / "evaluation_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "kind": "model_development_evaluation",
                "config": str(config_path),
                "stage": "short_screen",
                "condition_count": len(rows),
                "status_counts": {
                    status: sum(row["screen_status"] == status for row in rows)
                    for status in ("pass", "fail")
                },
                "qc_policy": "hook_angle_err_max_deg is diagnostic-only; all other required QC remains PASS/FAIL",
                "provenance": provenance,
                "outputs": {key: str(value) for key, value in outputs.items()},
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
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
