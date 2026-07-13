#!/usr/bin/env python3
"""Plot summary metrics for generic Phase 2 multi-run campaigns."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

from sim_swim.analysis.multi_run_campaign import (
    apply_campaign_cli_overrides,
    load_yaml,
)

matplotlib.use("Agg")

DEFAULT_METRICS = (
    "first_fail_t_s",
    "max_flag_bond_rel_err",
    "hook_len_rel_err_max",
    "axis_center_to_body_roll_ratio_mean",
)


def _get_plt():
    import matplotlib.pyplot as plt

    return plt


def _parse_args(argv: list[str] | None = None) -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-config", type=Path, required=True)
    parser.add_argument("--summary-csv", type=Path, default=None)
    parser.add_argument("--run-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--x-axis", type=str, default=None)
    parser.add_argument("--y-axis", type=str, default=None)
    parser.add_argument("--metrics", type=str, default=None)
    args, passthrough = parser.parse_known_args(argv)
    return args, passthrough


def _resolve_input_paths(args: argparse.Namespace, campaign: dict[str, Any]) -> None:
    if args.summary_csv is None:
        if args.run_dir is None:
            output_cfg = dict(campaign.get("output", {}) or {})
            if bool(output_cfg.get("timestamp_subdir", True)):
                raise ValueError(
                    "summary_csv=... or run_dir=... is required when "
                    "output.timestamp_subdir is true"
                )
            args.run_dir = Path(str(output_cfg["base_dir"]))
        args.summary_csv = args.run_dir / "summary.csv"
    if args.output_dir is None and args.run_dir is not None:
        args.output_dir = args.run_dir / "plots"


def _load_rows(summary_csv: Path) -> list[dict[str, str]]:
    if not summary_csv.is_file():
        raise FileNotFoundError(summary_csv)
    with summary_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    if not rows:
        raise ValueError(f"summary.csv is empty: {summary_csv}")
    return rows


def _axis_lookup(campaign: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {axis["name"]: axis for axis in campaign.get("sweep", {}).get("axes", [])}


def _metric_names(args: argparse.Namespace, campaign: dict[str, Any]) -> list[str]:
    if args.metrics:
        return [item.strip() for item in args.metrics.split(",") if item.strip()]
    metrics = list(campaign.get("plot", {}).get("metrics", []) or [])
    return metrics or list(DEFAULT_METRICS)


def _axis_names(
    args: argparse.Namespace, campaign: dict[str, Any]
) -> tuple[str, str | None]:
    plot_cfg = dict(campaign.get("plot", {}) or {})
    x_axis = args.x_axis or plot_cfg.get("default_x_axis")
    y_axis = args.y_axis if args.y_axis is not None else plot_cfg.get("default_y_axis")
    if not x_axis:
        raise ValueError("plot.default_x_axis or --x-axis is required")
    return str(x_axis), (None if y_axis in {"", None} else str(y_axis))


def _axis_label_column(axis_name: str) -> str:
    return f"axis_{axis_name}_label"


def _axis_index_column(axis_name: str) -> str:
    return f"axis_{axis_name}_index"


def _axis_labels(rows: list[dict[str, str]], axis: dict[str, Any]) -> list[str]:
    labels = list(axis["labels"])
    index_column = _axis_index_column(axis["name"])
    label_column = _axis_label_column(axis["name"])
    seen = {
        int(float(row[index_column]))
        for row in rows
        if row.get(index_column) not in {"", None}
    }
    if seen:
        return [labels[index] for index in sorted(seen)]
    present = [row.get(label_column, "") for row in rows if row.get(label_column, "")]
    return present or labels


def _parse_float(value: str | None) -> float:
    try:
        return float(value) if value not in {None, ""} else float("nan")
    except ValueError:
        return float("nan")


def _write_normalized_csv(rows: list[dict[str, str]], output_dir: Path) -> Path:
    out_path = output_dir / "plot_data.csv"
    fieldnames = list(rows[0].keys())
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return out_path


def _save_line_plot(
    *,
    rows: list[dict[str, str]],
    x_axis: dict[str, Any],
    metric: str,
    output_dir: Path,
) -> Path:
    plt = _get_plt()
    x_labels = _axis_labels(rows, x_axis)
    values = np.full(len(x_labels), np.nan, dtype=float)
    by_label = {label: index for index, label in enumerate(x_labels)}
    for row in rows:
        label = row.get(_axis_label_column(x_axis["name"]), "")
        if label not in by_label:
            continue
        values[by_label[label]] = _parse_float(row.get(metric))
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    ax.plot(np.arange(len(x_labels)), values, marker="o", color="#2b6cb0")
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=20, ha="right")
    ax.set_xlabel(x_axis["short_name"])
    ax.set_ylabel(metric)
    ax.set_title(f"{metric} vs {x_axis['short_name']}")
    out_path = output_dir / f"{metric}_vs_{x_axis['name']}.png"
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    return out_path


def _save_heatmap(
    *,
    rows: list[dict[str, str]],
    x_axis: dict[str, Any],
    y_axis: dict[str, Any],
    metric: str,
    output_dir: Path,
) -> Path:
    plt = _get_plt()
    x_labels = _axis_labels(rows, x_axis)
    y_labels = _axis_labels(rows, y_axis)
    x_index = {label: idx for idx, label in enumerate(x_labels)}
    y_index = {label: idx for idx, label in enumerate(y_labels)}
    matrix = np.full((len(y_labels), len(x_labels)), np.nan, dtype=float)
    for row in rows:
        x_label = row.get(_axis_label_column(x_axis["name"]), "")
        y_label = row.get(_axis_label_column(y_axis["name"]), "")
        if x_label not in x_index or y_label not in y_index:
            continue
        matrix[y_index[y_label], x_index[x_label]] = _parse_float(row.get(metric))
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_yticklabels(y_labels)
    ax.set_xlabel(x_axis["short_name"])
    ax.set_ylabel(y_axis["short_name"])
    ax.set_title(f"{metric} heatmap")
    fig.colorbar(image, ax=ax)
    out_path = output_dir / f"{metric}_heatmap.png"
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    return out_path


def main(argv: list[str] | None = None) -> None:
    args, passthrough = _parse_args(argv)
    campaign = apply_campaign_cli_overrides(
        load_yaml(args.campaign_config), passthrough
    )
    _resolve_input_paths(args, campaign)
    rows = _load_rows(args.summary_csv)
    output_dir = args.output_dir or (args.summary_csv.parent / "plots")
    output_dir.mkdir(parents=True, exist_ok=True)

    x_axis_name, y_axis_name = _axis_names(args, campaign)
    axes = _axis_lookup(campaign)
    if x_axis_name not in axes:
        raise ValueError(f"Unknown x-axis: {x_axis_name}")
    if y_axis_name is not None and y_axis_name not in axes:
        raise ValueError(f"Unknown y-axis: {y_axis_name}")

    normalized_csv = _write_normalized_csv(rows, output_dir)
    generated = [normalized_csv]
    metrics = _metric_names(args, campaign)
    for metric in metrics:
        if y_axis_name is None:
            generated.append(
                _save_line_plot(
                    rows=rows,
                    x_axis=axes[x_axis_name],
                    metric=metric,
                    output_dir=output_dir,
                )
            )
        else:
            generated.append(
                _save_heatmap(
                    rows=rows,
                    x_axis=axes[x_axis_name],
                    y_axis=axes[y_axis_name],
                    metric=metric,
                    output_dir=output_dir,
                )
            )
    for path in generated:
        print(path)


if __name__ == "__main__":
    main()
