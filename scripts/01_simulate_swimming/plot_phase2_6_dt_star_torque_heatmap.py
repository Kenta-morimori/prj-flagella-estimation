#!/usr/bin/env python3
"""Plot Phase 2.6 torque x dt_star heatmaps for baseline local scales."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np

matplotlib.use("Agg")


CATEGORY_ORDER = [
    "none",
    "motor_no_rotation",
    "hook",
    "flag",
    "finite",
]
CATEGORY_LABELS = {
    "none": "stable",
    "motor_no_rotation": "no rotation",
    "hook": "hook",
    "flag": "flag",
    "finite": "numeric/other",
}


def _get_plt():
    import matplotlib.pyplot as plt

    return plt


def _parse_float_list(text: str | None) -> set[float] | None:
    if text is None:
        return None
    values = {float(item.strip()) for item in text.split(",") if item.strip()}
    return values or None


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot Phase 2.6 torque x dt_star heatmaps from torque model "
            "evaluation summaries."
        )
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        action="append",
        required=True,
        help="Summary CSV emitted by run_phase2_6_torque_model_evaluation.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for heatmap outputs.",
    )
    parser.add_argument(
        "--torques",
        type=str,
        default=None,
        help="Optional comma-separated torque_Nm filter.",
    )
    parser.add_argument(
        "--dt-stars",
        type=str,
        default=None,
        help="Optional comma-separated dt_star filter.",
    )
    return parser.parse_args()


def _parse_bool(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _category_rank(row: dict[str, str]) -> int:
    if _parse_bool(row.get("shape_pass")):
        return CATEGORY_ORDER.index("none")
    category = str(row.get("first_fail_category", "finite")).strip()
    if category in CATEGORY_ORDER:
        return CATEGORY_ORDER.index(category)
    return CATEGORY_ORDER.index("finite")


def _is_baseline_row(row: dict[str, str]) -> bool:
    target = str(row.get("target", "")).strip()
    if target != "all_local_scales":
        return False
    try:
        return np.isclose(float(row.get("value", "nan")), 1.0)
    except ValueError:
        return False


def _collect_rows(
    summary_paths: list[Path],
    *,
    torques_filter: set[float] | None,
    dt_stars_filter: set[float] | None,
) -> list[dict[str, str]]:
    by_key: dict[tuple[float, float], dict[str, str]] = {}
    for path in summary_paths:
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if not _is_baseline_row(row):
                    continue
                torque = float(row["torque_Nm"])
                dt_star = float(row["dt_star"])
                if torques_filter is not None and torque not in torques_filter:
                    continue
                if dt_stars_filter is not None and dt_star not in dt_stars_filter:
                    continue
                by_key[(torque, dt_star)] = dict(row)

    rows = list(by_key.values())
    rows.sort(key=lambda r: (float(r["torque_Nm"]), float(r["dt_star"])))
    return rows


def _build_matrices(
    rows: list[dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if not rows:
        raise ValueError("No baseline rows matched the requested filters.")
    torques = np.array(sorted({float(row["torque_Nm"]) for row in rows}))
    dt_stars = np.array(sorted({float(row["dt_star"]) for row in rows}))
    category = np.full((torques.size, dt_stars.size), np.nan, dtype=float)
    pass_fail = np.full((torques.size, dt_stars.size), np.nan, dtype=float)
    for row in rows:
        t_idx = int(np.where(torques == float(row["torque_Nm"]))[0][0])
        d_idx = int(np.where(dt_stars == float(row["dt_star"]))[0][0])
        category[t_idx, d_idx] = _category_rank(row)
        pass_fail[t_idx, d_idx] = 1.0 if _parse_bool(row.get("shape_pass")) else 0.0
    return category, torques, dt_stars, pass_fail


def _save_category_heatmap(
    matrix: np.ndarray,
    torques: np.ndarray,
    dt_stars: np.ndarray,
    out_path: Path,
) -> None:
    plt = _get_plt()
    cmap = ListedColormap(["#2f855a", "#718096", "#dd6b20", "#c53030", "#4a5568"])
    cmap.set_bad("#edf2f7")
    norm = BoundaryNorm(np.arange(-0.5, len(CATEGORY_ORDER) + 0.5, 1.0), cmap.N)
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap=cmap, norm=norm)
    ax.set_xlabel("dt_star")
    ax.set_ylabel("torque_Nm")
    ax.set_title("Phase 2.6 dt_star x torque: first-fail category")
    ax.set_xticks(np.arange(dt_stars.size))
    ax.set_xticklabels([f"{value:.0e}" for value in dt_stars])
    ax.set_yticks(np.arange(torques.size))
    ax.set_yticklabels([f"{torque:.2e}" for torque in torques])
    cbar = fig.colorbar(image, ax=ax, ticks=list(range(len(CATEGORY_ORDER))))
    cbar.ax.set_yticklabels([CATEGORY_LABELS[key] for key in CATEGORY_ORDER])
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _save_pass_fail_heatmap(
    matrix: np.ndarray,
    torques: np.ndarray,
    dt_stars: np.ndarray,
    out_path: Path,
) -> None:
    plt = _get_plt()
    cmap = ListedColormap(["#c53030", "#2f855a"])
    cmap.set_bad("#edf2f7")
    norm = BoundaryNorm([-0.5, 0.5, 1.5], 2)
    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap=cmap, norm=norm)
    ax.set_xlabel("dt_star")
    ax.set_ylabel("torque_Nm")
    ax.set_title("Phase 2.6 dt_star x torque: shape pass/fail")
    ax.set_xticks(np.arange(dt_stars.size))
    ax.set_xticklabels([f"{value:.0e}" for value in dt_stars])
    ax.set_yticks(np.arange(torques.size))
    ax.set_yticklabels([f"{torque:.2e}" for torque in torques])
    cbar = fig.colorbar(image, ax=ax, ticks=[0, 1])
    cbar.ax.set_yticklabels(["fail", "pass"])
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _write_normalized_csv(rows: list[dict[str, str]], out_path: Path) -> None:
    fields = [
        "torque_Nm",
        "dt_star",
        "shape_pass",
        "first_fail_category",
        "first_fail_reason",
        "net_abs_flag_helix_spin_revolutions",
        "flag_helix_spin_direction_consistency",
        "max_flag_bond_rel_err",
        "max_flag_bend_err_deg",
        "max_flag_torsion_err_deg",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> None:
    args = _parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = _collect_rows(
        args.summary_csv,
        torques_filter=_parse_float_list(args.torques),
        dt_stars_filter=_parse_float_list(args.dt_stars),
    )
    category, torques, dt_stars, pass_fail = _build_matrices(rows)
    normalized_csv = args.output_dir / "phase2_6_dt_star_torque_heatmap.csv"
    category_png = args.output_dir / "phase2_6_dt_star_torque_category_heatmap.png"
    pass_fail_png = args.output_dir / "phase2_6_dt_star_torque_pass_fail_heatmap.png"
    _write_normalized_csv(rows, normalized_csv)
    _save_category_heatmap(category, torques, dt_stars, category_png)
    _save_pass_fail_heatmap(pass_fail, torques, dt_stars, pass_fail_png)
    print(f"Saved normalized CSV to {normalized_csv}")
    print(f"Saved category heatmap to {category_png}")
    print(f"Saved pass/fail heatmap to {pass_fail_png}")


if __name__ == "__main__":
    main()
