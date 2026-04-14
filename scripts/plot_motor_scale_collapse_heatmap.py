#!/usr/bin/env python3
"""Plot torque x scale collapse maps from sweep summary CSVs.

This script is intended for PhaseB1 and later only. It reads the long-form
summary emitted by `scripts.run_motor_scale_sweep`, pivots the table into a
torque-by-scale grid, and writes compact heatmap figures plus a normalized CSV.

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
    pass_plot = output_dir / f"{args.summary_csv.stem}_shape_pass_heatmap.png"

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

    pass_matrix = np.where(
        np.isfinite(matrix),
        (matrix == 0).astype(float),
        np.nan,
    )
    pass_cmap = ListedColormap(["#c53030", "#2f855a"])
    pass_norm = BoundaryNorm([-0.5, 0.5, 1.5], pass_cmap.N)

    _save_heatmap(
        pass_matrix,
        torque_vals,
        scale_vals,
        pass_plot,
        title=f"{title_prefix}: shape pass",
        cmap=pass_cmap,
        norm=pass_norm,
        cbar_ticks=[0, 1],
    )

    print(f"Saved normalized CSV to {normalized_csv}")
    print(f"Saved category heatmap to {category_plot}")
    print(f"Saved shape-pass heatmap to {pass_plot}")


if __name__ == "__main__":
    main()
