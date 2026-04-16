#!/usr/bin/env python3
"""Plot torque x scale collapse maps from sweep summary CSVs.

This script is intended for PhaseB1 and later only. It reads the long-form
summary emitted by `scripts.run_motor_scale_sweep`, pivots the table into a
torque-by-scale grid, and writes a category heatmap plus body / hook / flagella
pass-fail heatmaps, along with a normalized CSV.

Example:
  uv run python -m scripts.plot_motor_scale_collapse_heatmap \
        --summary-csv \
        outputs/phaseb1_torque_scale_smoke/local_hook_scale_sweep_summary.csv
"""

from __future__ import annotations

import csv
import argparse
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap

matplotlib.use("Agg")


def _get_plt():
    import matplotlib.pyplot as plt

    return plt


CATEGORY_ORDER = [
    "none",
    "hook",
    "body_spring",
    "body_bend",
    "body_centerline",
    "body_area",
    "body_nonfinite",
    "flag",
    "flag_nonfinite",
    "finite",
]

CATEGORY_LABELS = {
    "none": "stable",
    "hook": "hook",
    "body_spring": "body spring",
    "body_bend": "body bend",
    "body_centerline": "body centerline",
    "body_area": "body area",
    "body_nonfinite": "body nan/inf",
    "flag": "flag",
    "flag_nonfinite": "flag nan/inf",
    "finite": "numeric fail",
}

BODY_FAIL_CATEGORIES = {
    "body_spring",
    "body_bend",
    "body_centerline",
    "body_area",
    "body_nonfinite",
}
HOOK_FAIL_CATEGORIES = {"hook"}
FLAGELLA_FAIL_CATEGORIES = {"flag", "flag_nonfinite"}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=("Plot torque x scale collapse heatmaps from sweep summaries.")
    )
    parser.add_argument(
        "--summary-csv",
        required=True,
        type=Path,
        help="Long-form summary CSV emitted by run_motor_scale_sweep.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=("Directory for plots and normalized CSV. Defaults next to the input."),
    )
    parser.add_argument(
        "--title",
        type=str,
        default=None,
        help="Optional figure title prefix.",
    )
    return parser.parse_args()


def _collapse_rank(row: dict[str, str]) -> int:
    if str(row.get("shape_pass", "")).strip().lower() in {"true", "1"}:
        return CATEGORY_ORDER.index("none")

    category = str(row.get("first_fail_category", "none")).strip()
    if category in CATEGORY_ORDER:
        return CATEGORY_ORDER.index(category)
    return CATEGORY_ORDER.index("finite")


def _make_category_matrix(
    rows: list[dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    scale_col = "value"
    torque_col = "torque_Nm"
    if not rows:
        raise ValueError("summary CSV is empty")
    if torque_col not in rows[0]:
        raise ValueError(f"Missing required column: {torque_col}")

    torque_vals = np.array(sorted({float(row[torque_col]) for row in rows}))
    scale_vals = np.array(sorted({float(row[scale_col]) for row in rows}))

    mat = np.full((torque_vals.size, scale_vals.size), np.nan, dtype=float)
    for row in rows:
        t_idx = int(np.where(torque_vals == float(row[torque_col]))[0][0])
        s_idx = int(np.where(scale_vals == float(row[scale_col]))[0][0])
        mat[t_idx, s_idx] = _collapse_rank(row)

    return mat, torque_vals, scale_vals


def _make_component_matrix(
    rows: list[dict[str, str]],
    fail_categories: set[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    category_matrix, torque_vals, scale_vals = _make_category_matrix(rows)
    component_matrix = np.ones_like(category_matrix, dtype=float)

    for row in rows:
        first_fail_category = str(row.get("first_fail_category", "none")).strip()
        if first_fail_category in fail_categories:
            t_idx = int(np.where(torque_vals == float(row["torque_Nm"]))[0][0])
            s_idx = int(np.where(scale_vals == float(row["value"]))[0][0])
            component_matrix[t_idx, s_idx] = 0.0

    return component_matrix, torque_vals, scale_vals


def _save_heatmap(
    matrix: np.ndarray,
    torque_vals: np.ndarray,
    scale_vals: np.ndarray,
    out_path: Path,
    *,
    title: str,
    cmap: ListedColormap,
    norm: BoundaryNorm,
    cbar_ticks: list[int],
) -> None:
    plt = _get_plt()
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    image = ax.imshow(
        matrix,
        origin="lower",
        aspect="auto",
        cmap=cmap,
        norm=norm,
    )
    ax.set_xlabel("local scale value")
    ax.set_ylabel("torque_Nm")
    ax.set_title(title)
    ax.set_xticks(np.arange(scale_vals.size))
    ax.set_xticklabels([f"{value:g}" for value in scale_vals], rotation=45, ha="right")
    ax.set_yticks(np.arange(torque_vals.size))
    ax.set_yticklabels([f"{value:.2e}" for value in torque_vals])

    cbar = fig.colorbar(image, ax=ax, ticks=cbar_ticks)
    cbar.ax.set_yticklabels(
        [CATEGORY_LABELS[CATEGORY_ORDER[tick]] for tick in cbar_ticks]
    )

    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _save_pass_fail_heatmap(
    matrix: np.ndarray,
    torque_vals: np.ndarray,
    scale_vals: np.ndarray,
    out_path: Path,
    *,
    title: str,
) -> None:
    plt = _get_plt()
    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    image = ax.imshow(
        matrix,
        origin="lower",
        aspect="auto",
        cmap=ListedColormap(["#c53030", "#2f855a"]),
        norm=BoundaryNorm([-0.5, 0.5, 1.5], 2),
    )
    ax.set_xlabel("local scale value")
    ax.set_ylabel("torque_Nm")
    ax.set_title(title)
    ax.set_xticks(np.arange(scale_vals.size))
    ax.set_xticklabels([f"{value:g}" for value in scale_vals], rotation=45, ha="right")
    ax.set_yticks(np.arange(torque_vals.size))
    ax.set_yticklabels([f"{value:.2e}" for value in torque_vals])

    cbar = fig.colorbar(image, ax=ax, ticks=[0, 1])
    cbar.ax.set_yticklabels(["fail", "pass"])

    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _make_body_shape_pass_matrix(
    rows: list[dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract body_shape_pass column directly for body component pass/fail judgment.

    This function judges body pass/fail based on the body_shape_pass column,
    independent of the first_fail_category priority ordering.
    """
    scale_col = "value"
    torque_col = "torque_Nm"
    if not rows:
        raise ValueError("summary CSV is empty")
    if torque_col not in rows[0]:
        raise ValueError(f"Missing required column: {torque_col}")
    if "body_shape_pass" not in rows[0]:
        raise ValueError("Missing required column: body_shape_pass")

    torque_vals = np.array(sorted({float(row[torque_col]) for row in rows}))
    scale_vals = np.array(sorted({float(row[scale_col]) for row in rows}))

    mat = np.full((torque_vals.size, scale_vals.size), np.nan, dtype=float)
    for row in rows:
        t_idx = int(np.where(torque_vals == float(row[torque_col]))[0][0])
        s_idx = int(np.where(scale_vals == float(row[scale_col]))[0][0])
        # 1.0 = pass, 0.0 = fail
        body_pass = str(row.get("body_shape_pass", "")).strip().lower() in {"true", "1"}
        mat[t_idx, s_idx] = 1.0 if body_pass else 0.0

    return mat, torque_vals, scale_vals


def _make_hook_pass_matrix(
    rows: list[dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract hook pass/fail from first_fail_category."""
    category_matrix, torque_vals, scale_vals = _make_category_matrix(rows)
    component_matrix = np.ones_like(category_matrix, dtype=float)

    for row in rows:
        first_fail_category = str(row.get("first_fail_category", "none")).strip()
        if first_fail_category in HOOK_FAIL_CATEGORIES:
            t_idx = int(np.where(torque_vals == float(row["torque_Nm"]))[0][0])
            s_idx = int(np.where(scale_vals == float(row["value"]))[0][0])
            component_matrix[t_idx, s_idx] = 0.0

    return component_matrix, torque_vals, scale_vals


def _make_flagella_pass_matrix(
    rows: list[dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract flagella pass/fail from first_fail_category."""
    category_matrix, torque_vals, scale_vals = _make_category_matrix(rows)
    component_matrix = np.ones_like(category_matrix, dtype=float)

    for row in rows:
        first_fail_category = str(row.get("first_fail_category", "none")).strip()
        if first_fail_category in FLAGELLA_FAIL_CATEGORIES:
            t_idx = int(np.where(torque_vals == float(row["torque_Nm"]))[0][0])
            s_idx = int(np.where(scale_vals == float(row["value"]))[0][0])
            component_matrix[t_idx, s_idx] = 0.0

    return component_matrix, torque_vals, scale_vals


def _save_combined_pass_fail_heatmap(
    body_matrix: np.ndarray,
    hook_matrix: np.ndarray,
    flagella_matrix: np.ndarray,
    torque_vals: np.ndarray,
    scale_vals: np.ndarray,
    out_path: Path,
    *,
    title: str,
) -> None:
    """Save a 1x3 subplot heatmap showing body, hook, and flagella pass/fail."""
    plt = _get_plt()
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

    cmap = ListedColormap(["#c53030", "#2f855a"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5], 2)

    components = [
        ("body", body_matrix),
        ("hook", hook_matrix),
        ("flagella", flagella_matrix),
    ]

    for ax, (label, matrix) in zip(axes, components):
        image = ax.imshow(
            matrix,
            origin="lower",
            aspect="auto",
            cmap=cmap,
            norm=norm,
        )
        ax.set_xlabel("local scale value")
        ax.set_ylabel("torque_Nm")
        ax.set_title(f"{label} pass/fail")
        ax.set_xticks(np.arange(scale_vals.size))
        ax.set_xticklabels(
            [f"{value:g}" for value in scale_vals], rotation=45, ha="right"
        )
        ax.set_yticks(np.arange(torque_vals.size))
        ax.set_yticklabels([f"{value:.2e}" for value in torque_vals])

        cbar = fig.colorbar(image, ax=ax, ticks=[0, 1])
        cbar.ax.set_yticklabels(["fail", "pass"])

    fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _has_flagella_failure(rows: list[dict[str, str]]) -> bool:
    """Check if any row has flagella failure.

    Flagella heatmap is output only if at least one run shows flagella failure
    (first_fail_category in FLAGELLA_FAIL_CATEGORIES). This correctly skips
    output for minimal_basal_stub mode where flagella is not deployed.
    """
    return any(
        str(row.get("first_fail_category", "none")).strip() in FLAGELLA_FAIL_CATEGORIES
        for row in rows
    )


def main() -> None:
    args = _parse_args()
    with args.summary_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    if not rows:
        raise ValueError("summary CSV is empty")

    if "torque_Nm" not in rows[0]:
        raise ValueError("summary CSV must include torque_Nm")
    if "shape_pass" not in rows[0] or "first_fail_category" not in rows[0]:
        raise ValueError("summary CSV must include shape_pass and first_fail_category")
    if "body_shape_pass" not in rows[0]:
        raise ValueError("summary CSV must include body_shape_pass")

    output_dir = args.output_dir or args.summary_csv.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    matrix, torque_vals, scale_vals = _make_category_matrix(rows)

    normalized_csv = output_dir / f"{args.summary_csv.stem}_collapse_matrix.csv"
    with normalized_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "torque_Nm",
                "value",
                "shape_pass",
                "first_fail_category",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "torque_Nm": float(row["torque_Nm"]),
                    "value": float(row["value"]),
                    "shape_pass": str(row["shape_pass"]).strip().lower()
                    in {"true", "1"},
                    "first_fail_category": str(row["first_fail_category"]),
                }
            )

    cmap = ListedColormap(
        [
            "#2f855a",  # stable
            "#dd6b20",  # hook
            "#b7791f",  # body spring
            "#d53f8c",  # body bend
            "#805ad5",  # body centerline
            "#3182ce",  # body area
            "#4a5568",  # body nan/inf
            "#c53030",  # flag
            "#2c5282",  # flag nan/inf
            "#718096",  # numeric fail
        ]
    )
    norm = BoundaryNorm(
        np.arange(-0.5, len(CATEGORY_ORDER) + 0.5, 1.0),
        cmap.N,
    )

    title_prefix = args.title or args.summary_csv.stem
    category_plot = output_dir / f"{args.summary_csv.stem}_category_heatmap.png"
    combined_pass_fail_plot = (
        output_dir / f"{args.summary_csv.stem}_combined_pass_fail_heatmap.png"
    )
    flagella_plot = (
        output_dir / f"{args.summary_csv.stem}_flagella_pass_fail_heatmap.png"
    )

    _save_heatmap(
        matrix,
        torque_vals,
        scale_vals,
        category_plot,
        title=f"{title_prefix}: first-fail category",
        cmap=cmap,
        norm=norm,
        cbar_ticks=list(range(len(CATEGORY_ORDER))),
    )

    # Generate component pass/fail matrices using independent judgment axes
    body_matrix, body_torque_vals, body_scale_vals = _make_body_shape_pass_matrix(rows)
    hook_matrix, hook_torque_vals, hook_scale_vals = _make_hook_pass_matrix(rows)
    flagella_matrix, flag_torque_vals, flag_scale_vals = _make_flagella_pass_matrix(
        rows
    )

    # Save combined 1x3 subplot heatmap (primary visualization)
    _save_combined_pass_fail_heatmap(
        body_matrix,
        hook_matrix,
        flagella_matrix,
        torque_vals,
        scale_vals,
        combined_pass_fail_plot,
        title=f"{title_prefix}: component pass/fail (independent judgment axes)",
    )

    # Conditionally save flagella heatmap only if flagella failures exist
    if _has_flagella_failure(rows):
        _save_pass_fail_heatmap(
            flagella_matrix,
            flag_torque_vals,
            flag_scale_vals,
            flagella_plot,
            title=f"{title_prefix}: flagella pass/fail",
        )
        print(f"Saved flagella pass/fail heatmap to {flagella_plot}")
    else:
        print("Skipped flagella pass/fail heatmap because no flag failures were found")

    print(f"Saved normalized CSV to {normalized_csv}")
    print(f"Saved category heatmap to {category_plot}")
    print(
        f"Saved combined pass/fail heatmap (1x3 subplots) to {combined_pass_fail_plot}"
    )


if __name__ == "__main__":
    main()
