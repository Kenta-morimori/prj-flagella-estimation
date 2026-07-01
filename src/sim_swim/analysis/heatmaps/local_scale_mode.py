#!/usr/bin/env python3
"""Plot Phase 2.6 torque x local-scale-mode heatmaps."""

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
TARGET_LABELS = {
    "local_spring_scale": "spring",
    "local_bend_scale": "bend",
    "local_torsion_scale": "torsion",
    "local_hook_scale": "hook",
}
MODE_ORDER = ["all=1", "spring=2", "bend=2", "torsion=2", "hook=2"]


def _get_plt():
    import matplotlib.pyplot as plt

    return plt


def _parse_float_list(text: str | None) -> set[float] | None:
    if text is None:
        return None
    values = {float(item.strip()) for item in text.split(",") if item.strip()}
    return values or None


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Phase 2.6 torque x local-scale-mode heatmaps."
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
        "--mode-value",
        type=float,
        default=2.0,
        help="Local scale value used for one-factor modes.",
    )
    return parser.parse_args(argv)


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


def _mode_for_row(row: dict[str, str], *, mode_value: float) -> str | None:
    target = str(row.get("target", ""))
    value = float(row.get("value", "nan"))
    if target == "all_local_scales" and np.isclose(value, 1.0):
        return "all=1"
    if target in TARGET_LABELS and np.isclose(value, mode_value):
        return f"{TARGET_LABELS[target]}={mode_value:g}"
    return None


def _collect_rows(
    summary_paths: list[Path],
    *,
    torques_filter: set[float] | None,
    mode_value: float,
) -> list[dict[str, str]]:
    by_key: dict[tuple[float, str], dict[str, str]] = {}
    for path in summary_paths:
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                torque = float(row["torque_Nm"])
                if torques_filter is not None and torque not in torques_filter:
                    continue
                mode = _mode_for_row(row, mode_value=mode_value)
                if mode is None:
                    continue
                row = dict(row)
                row["mode"] = mode
                by_key[(torque, mode)] = row

    rows = list(by_key.values())
    rows.sort(key=lambda r: (float(r["torque_Nm"]), MODE_ORDER.index(r["mode"])))
    return rows


def _build_matrices(
    rows: list[dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray]:
    if not rows:
        raise ValueError("No rows matched the requested filters.")
    torques = np.array(sorted({float(row["torque_Nm"]) for row in rows}))
    modes = [mode for mode in MODE_ORDER if any(row["mode"] == mode for row in rows)]
    category = np.full((torques.size, len(modes)), np.nan, dtype=float)
    pass_fail = np.full((torques.size, len(modes)), np.nan, dtype=float)
    for row in rows:
        t_idx = int(np.where(torques == float(row["torque_Nm"]))[0][0])
        m_idx = modes.index(row["mode"])
        category[t_idx, m_idx] = _category_rank(row)
        pass_fail[t_idx, m_idx] = 1.0 if _parse_bool(row.get("shape_pass")) else 0.0
    return category, torques, modes, pass_fail


def _save_category_heatmap(
    matrix: np.ndarray,
    torques: np.ndarray,
    modes: list[str],
    out_path: Path,
) -> None:
    plt = _get_plt()
    cmap = ListedColormap(["#2f855a", "#718096", "#dd6b20", "#c53030", "#4a5568"])
    norm = BoundaryNorm(np.arange(-0.5, len(CATEGORY_ORDER) + 0.5, 1.0), cmap.N)
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap=cmap, norm=norm)
    ax.set_xlabel("local scale mode")
    ax.set_ylabel("torque_Nm")
    ax.set_title("Phase 2.6 local-scale mode: first-fail category")
    ax.set_xticks(np.arange(len(modes)))
    ax.set_xticklabels(modes, rotation=30, ha="right")
    ax.set_yticks(np.arange(torques.size))
    ax.set_yticklabels([f"{torque:.2e}" for torque in torques])
    cbar = fig.colorbar(image, ax=ax, ticks=list(range(len(CATEGORY_ORDER))))
    cbar.ax.set_yticklabels([CATEGORY_LABELS[key] for key in CATEGORY_ORDER])
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _save_pass_fail_heatmap(
    matrix: np.ndarray,
    torques: np.ndarray,
    modes: list[str],
    out_path: Path,
) -> None:
    plt = _get_plt()
    cmap = ListedColormap(["#c53030", "#2f855a"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5], 2)
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap=cmap, norm=norm)
    ax.set_xlabel("local scale mode")
    ax.set_ylabel("torque_Nm")
    ax.set_title("Phase 2.6 local-scale mode: shape pass/fail")
    ax.set_xticks(np.arange(len(modes)))
    ax.set_xticklabels(modes, rotation=30, ha="right")
    ax.set_yticks(np.arange(torques.size))
    ax.set_yticklabels([f"{torque:.2e}" for torque in torques])
    cbar = fig.colorbar(image, ax=ax, ticks=[0, 1])
    cbar.ax.set_yticklabels(["fail", "pass"])
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _write_normalized_csv(rows: list[dict[str, str]], out_path: Path) -> None:
    fields = [
        "torque_Nm",
        "mode",
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


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    mode_value = float(args.mode_value)
    if not np.isclose(mode_value, 2.0):
        raise SystemExit(
            f"--mode-value={mode_value:g} is not supported yet (expected 2.0). "
            "Update MODE_ORDER/sorting logic to enable other values."
        )

    rows = _collect_rows(
        args.summary_csv,
        torques_filter=_parse_float_list(args.torques),
        mode_value=mode_value,
    )
    category, torques, modes, pass_fail = _build_matrices(rows)
    normalized_csv = args.output_dir / "heatmap_data.csv"
    category_png = args.output_dir / "local_scale_mode_category_heatmap.png"
    pass_fail_png = args.output_dir / "local_scale_mode_pass_fail_heatmap.png"
    _write_normalized_csv(rows, normalized_csv)
    _save_category_heatmap(category, torques, modes, category_png)
    _save_pass_fail_heatmap(pass_fail, torques, modes, pass_fail_png)
    print(f"Saved normalized CSV to {normalized_csv}")
    print(f"Saved category heatmap to {category_png}")
    print(f"Saved pass/fail heatmap to {pass_fail_png}")


if __name__ == "__main__":
    main()
