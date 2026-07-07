#!/usr/bin/env python3
"""Plot generic Phase 2 shape-stability grid heatmaps."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np

matplotlib.use("Agg")

CATEGORY_ORDER = ["none", "hook", "flag", "finite"]
CATEGORY_LABELS = {
    "none": "stable",
    "hook": "hook",
    "flag": "flag",
    "finite": "numeric/other",
}
NUMERIC_METRICS = (
    "hook_len_rel_err_max",
    "local_attach_first_rel_err",
    "local_attach_first_vs_body_axis_err_deg",
    "local_first_second_rel_err",
    "local_attach_frame_position_rel_err",
    "local_attach_frame_position_angle_err_deg",
    "local_attach_frame_tangent_angle_err_deg",
)


def _get_plt():
    import matplotlib.pyplot as plt

    return plt


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-csv", type=Path, required=True)
    parser.add_argument(
        "--mode",
        choices=("body-first-grid", "first-second-grid", "attach-frame-grid"),
        required=True,
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args(argv)


def _parse_bool(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _parse_float(value: str | float | None) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")


def _legacy_category_rank(row: dict[str, str]) -> int:
    if _parse_bool(row.get("final_shape_pass_nonbody")):
        return CATEGORY_ORDER.index("none")
    category = str(row.get("final_first_fail_category_nonbody", "")).strip()
    if category == "none":
        return CATEGORY_ORDER.index("finite")
    if category in CATEGORY_ORDER:
        return CATEGORY_ORDER.index(category)
    return CATEGORY_ORDER.index("finite")


def _category_rank(row: dict[str, str]) -> int:
    first_fail_category = row.get("first_fail_category_nonbody")
    if first_fail_category is None:
        return _legacy_category_rank(row)

    category = str(first_fail_category).strip()
    first_fail_t = _parse_float(row.get("first_fail_t_s"))
    has_first_fail = np.isfinite(first_fail_t)
    if not has_first_fail and category in {"", "none"}:
        return CATEGORY_ORDER.index("none")
    if category in {"", "none"}:
        return CATEGORY_ORDER.index("finite")
    if category in CATEGORY_ORDER:
        return CATEGORY_ORDER.index(category)
    return CATEGORY_ORDER.index("finite")


def _load_rows(summary_csv: Path, mode: str) -> list[dict[str, str]]:
    if not summary_csv.is_file():
        raise FileNotFoundError(summary_csv)
    with summary_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]
    rows = [row for row in rows if row.get("mode", mode) == mode]
    if not rows:
        raise ValueError(f"No rows found for mode={mode!r} in {summary_csv}")
    return rows


def _summary_csv_candidates(summary_csv: Path) -> list[Path]:
    candidates: list[Path] = []
    for root in (summary_csv.parent, summary_csv.parent.parent):
        if not root.is_dir():
            continue
        candidates.extend(root.glob("summary.csv"))
        candidates.extend(root.glob("*/summary.csv"))
    return sorted(set(candidates))


def _missing_summary_message(summary_csv: Path) -> str:
    message = f"Summary CSV not found: {summary_csv}"
    candidates = _summary_csv_candidates(summary_csv)
    if candidates:
        candidate_lines = "\n".join(f"  - {path}" for path in candidates)
        message = f"{message}\nCandidate summary CSV files:\n{candidate_lines}"
    return message


def _axes_for_rows(
    rows: list[dict[str, str]], mode: str
) -> tuple[list[float], list[str], str, str]:
    if mode == "body-first-grid":
        x_values = sorted(
            {_parse_float(row["local_attach_first_spring_scale"]) for row in rows}
        )
        y_values_float = sorted(
            {
                _parse_float(row["local_attach_first_body_axis_angle_scale"])
                for row in rows
            }
        )
        return (
            x_values,
            [f"{value:g}" for value in y_values_float],
            "local_attach_first_spring_scale",
            "local_attach_first_body_axis_angle_scale",
        )

    if mode == "attach-frame-grid":
        x_values = sorted(
            {_parse_float(row["local_attach_frame_position_scale"]) for row in rows}
        )
        y_values_float = sorted(
            {_parse_float(row["local_attach_frame_tangent_scale"]) for row in rows}
        )
        return (
            x_values,
            [f"{value:g}" for value in y_values_float],
            "local_attach_frame_position_scale",
            "local_attach_frame_tangent_scale",
        )

    x_values = sorted(
        {_parse_float(row["local_first_second_spring_scale"]) for row in rows}
    )
    y_labels = sorted(
        {
            "af="
            f"{_parse_float(row['local_attach_first_spring_scale']):g}, "
            "axis="
            f"{_parse_float(row['local_attach_first_body_axis_angle_scale']):g}"
            for row in rows
        }
    )
    return (
        x_values,
        y_labels,
        "local_first_second_spring_scale",
        "fixed_attach_first_and_body_axis_angle_scales",
    )


def _row_y_label(row: dict[str, str], mode: str) -> str:
    if mode == "body-first-grid":
        return f"{_parse_float(row['local_attach_first_body_axis_angle_scale']):g}"
    if mode == "attach-frame-grid":
        return f"{_parse_float(row['local_attach_frame_tangent_scale']):g}"
    return (
        f"af={_parse_float(row['local_attach_first_spring_scale']):g}, "
        f"axis={_parse_float(row['local_attach_first_body_axis_angle_scale']):g}"
    )


def _build_matrix(
    rows: list[dict[str, str]],
    mode: str,
    *,
    metric: str | None,
) -> tuple[np.ndarray, list[float], list[str], str, str]:
    x_values, y_labels, x_name, y_name = _axes_for_rows(rows, mode)
    matrix = np.full((len(y_labels), len(x_values)), np.nan, dtype=float)
    for row in rows:
        x = _parse_float(
            row[
                "local_attach_first_spring_scale"
                if mode == "body-first-grid"
                else "local_attach_frame_position_scale"
                if mode == "attach-frame-grid"
                else "local_first_second_spring_scale"
            ]
        )
        y = _row_y_label(row, mode)
        x_idx = x_values.index(x)
        y_idx = y_labels.index(y)
        matrix[y_idx, x_idx] = (
            _parse_float(row.get(metric, ""))
            if metric is not None
            else _category_rank(row)
        )
    return matrix, x_values, y_labels, x_name, y_name


def _save_category_heatmap(
    matrix: np.ndarray,
    x_values: list[float],
    y_labels: list[str],
    x_name: str,
    y_name: str,
    out_path: Path,
) -> None:
    plt = _get_plt()
    cmap = ListedColormap(["#2f855a", "#dd6b20", "#c53030", "#4a5568"])
    norm = BoundaryNorm(np.arange(-0.5, len(CATEGORY_ORDER) + 0.5, 1.0), cmap.N)
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap=cmap, norm=norm)
    ax.set_xlabel(x_name)
    ax.set_ylabel(y_name)
    ax.set_title("Shape-stability first-fail category")
    ax.set_xticks(np.arange(len(x_values)))
    ax.set_xticklabels([f"{value:g}" for value in x_values], rotation=30, ha="right")
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_yticklabels(y_labels)
    cbar = fig.colorbar(image, ax=ax, ticks=list(range(len(CATEGORY_ORDER))))
    cbar.ax.set_yticklabels([CATEGORY_LABELS[key] for key in CATEGORY_ORDER])
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _save_pass_fail_heatmap(
    matrix: np.ndarray,
    x_values: list[float],
    y_labels: list[str],
    x_name: str,
    y_name: str,
    out_path: Path,
) -> None:
    plt = _get_plt()
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    image = ax.imshow(
        matrix,
        origin="lower",
        aspect="auto",
        cmap=ListedColormap(["#c53030", "#2f855a"]),
        norm=BoundaryNorm([-0.5, 0.5, 1.5], 2),
    )
    ax.set_xlabel(x_name)
    ax.set_ylabel(y_name)
    ax.set_title("Shape-stability pass/fail")
    ax.set_xticks(np.arange(len(x_values)))
    ax.set_xticklabels([f"{value:g}" for value in x_values], rotation=30, ha="right")
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_yticklabels(y_labels)
    cbar = fig.colorbar(image, ax=ax, ticks=[0, 1])
    cbar.ax.set_yticklabels(["fail", "pass"])
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _save_numeric_heatmap(
    matrix: np.ndarray,
    x_values: list[float],
    y_labels: list[str],
    x_name: str,
    y_name: str,
    metric: str,
    out_path: Path,
) -> None:
    plt = _get_plt()
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap="viridis")
    ax.set_xlabel(x_name)
    ax.set_ylabel(y_name)
    ax.set_title(f"Shape-stability {metric}")
    ax.set_xticks(np.arange(len(x_values)))
    ax.set_xticklabels([f"{value:g}" for value in x_values], rotation=30, ha="right")
    ax.set_yticks(np.arange(len(y_labels)))
    ax.set_yticklabels(y_labels)
    fig.colorbar(image, ax=ax)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _write_normalized_csv(rows: list[dict[str, str]], out_path: Path) -> None:
    fields = [
        "condition_id",
        "mode",
        "final_shape_pass_nonbody",
        "final_first_fail_category_nonbody",
        "first_fail_category_nonbody",
        "local_attach_first_spring_scale",
        "local_attach_first_body_axis_angle_scale",
        "local_first_second_spring_scale",
        "local_attach_frame_position_scale",
        "local_attach_frame_tangent_scale",
        *NUMERIC_METRICS,
        "hook_len_rel_err_max_flag_id",
        "hook_len_rel_err_max_attach_body_bead_index",
        "hook_len_rel_err_max_flag_first_bead_index",
        "hook_len_rel_err_max_len_over_b",
        "hook_len_rel_err_per_flag",
        "local_attach_first_rel_err_per_flag",
        "local_first_second_rel_err_per_flag",
        "local_attach_first_vs_body_axis_err_deg_per_flag",
        "local_attach_frame_position_rel_err_per_flag",
        "local_attach_frame_position_angle_err_deg_per_flag",
        "local_attach_frame_tangent_angle_err_deg_per_flag",
        "first_fail_t_s",
        "first_fail_hook_len_rel_err_max",
        "first_fail_hook_len_rel_err_max_flag_id",
        "first_fail_local_attach_first_rel_err",
        "first_fail_local_first_second_rel_err",
        "first_fail_local_attach_frame_position_rel_err",
        "first_fail_local_attach_frame_position_angle_err_deg",
        "first_fail_local_attach_frame_tangent_angle_err_deg",
        "max_hook_len_rel_err_t_s",
        "max_hook_len_rel_err",
        "max_hook_len_rel_err_flag_id",
        "max_hook_len_rel_err_attach_body_bead_index",
        "max_hook_len_rel_err_flag_first_bead_index",
        "max_hook_len_rel_err_len_over_b",
        "max_hook_local_attach_frame_position_rel_err",
        "max_hook_local_attach_frame_position_angle_err_deg",
        "max_hook_local_attach_frame_tangent_angle_err_deg",
        "net_abs_flag_helix_spin_revolutions",
        "flag_helix_spin_direction_consistency",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    try:
        rows = _load_rows(args.summary_csv, args.mode)
    except FileNotFoundError as exc:
        raise SystemExit(_missing_summary_message(args.summary_csv)) from exc

    args.output_dir.mkdir(parents=True, exist_ok=True)

    normalized_csv = args.output_dir / "heatmap_data.csv"
    _write_normalized_csv(rows, normalized_csv)

    category_matrix, x_values, y_labels, x_name, y_name = _build_matrix(
        rows, args.mode, metric=None
    )
    _save_category_heatmap(
        category_matrix,
        x_values,
        y_labels,
        x_name,
        y_name,
        args.output_dir / "first_fail_category_heatmap.png",
    )
    pass_fail_matrix = np.where(
        category_matrix == CATEGORY_ORDER.index("none"), 1.0, 0.0
    )
    _save_pass_fail_heatmap(
        pass_fail_matrix,
        x_values,
        y_labels,
        x_name,
        y_name,
        args.output_dir / "shape_pass_fail_heatmap.png",
    )

    for metric in NUMERIC_METRICS:
        matrix, x_values, y_labels, x_name, y_name = _build_matrix(
            rows, args.mode, metric=metric
        )
        _save_numeric_heatmap(
            matrix,
            x_values,
            y_labels,
            x_name,
            y_name,
            metric,
            args.output_dir / f"{metric}_heatmap.png",
        )

    print(f"Saved normalized CSV to {normalized_csv}")
    print(f"Saved heatmaps to {args.output_dir}")


if __name__ == "__main__":
    main()
