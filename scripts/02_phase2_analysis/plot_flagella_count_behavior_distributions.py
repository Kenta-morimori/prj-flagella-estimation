#!/usr/bin/env python3
"""Plot Phase 2.8 flagella-count behavior feature distributions."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
from typing import Any

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "flagella_estimation_matplotlib"),
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import yaml


DEFAULT_DATASET_ROOT = Path("outputs/phase2_analysis/flagella_count_behavior/datasets")

EXCLUDED_FEATURE_CATEGORIES = {"metadata", "quality"}
QC_COLUMNS = [
    "quality_class",
    "use_for_analysis",
    "shape_pass",
    "relaxed_pass",
    "review_required",
    "first_fail_category",
]
PATH_AND_STATUS_COLUMNS = {
    "run_status",
    "raw_dir",
    "step_summary_csv",
    "timeseries_csv",
}
BOOL_TRUE = {"true", "1", "yes"}
BOOL_FALSE = {"false", "0", "no"}
QUALITY_COLORS = {
    "strict_pass": "#2563eb",
    "relaxed_pass": "#f59e0b",
    "fail": "#dc2626",
    "invalid_numeric": "#7c3aed",
    "missing_raw": "#6b7280",
}


def resolve_dataset_dir(
    *,
    dataset_dir: Path | None,
    dataset_id: str | None,
    dataset_root: Path = DEFAULT_DATASET_ROOT,
) -> Path:
    if dataset_dir is not None:
        return dataset_dir
    if dataset_id:
        return dataset_root / dataset_id
    raise ValueError("Either --dataset-dir or --dataset-id is required")


def _load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _read_inputs(dataset_dir: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    summary_csv = dataset_dir / "summary.csv"
    manifest_json = dataset_dir / "dataset_manifest.json"
    feature_schema_yaml = dataset_dir / "feature_schema_used.yaml"
    missing = [
        str(path)
        for path in (summary_csv, manifest_json, feature_schema_yaml)
        if not path.is_file()
    ]
    if missing:
        raise FileNotFoundError(f"Missing required dataset input(s): {missing}")

    summary = pd.read_csv(summary_csv)
    if summary.empty:
        raise ValueError(f"summary.csv is empty: {summary_csv}")
    for column in ("n_flagella", "quality_class", "use_for_analysis"):
        if column not in summary.columns:
            raise ValueError(f"summary.csv must include {column}")

    json.loads(manifest_json.read_text(encoding="utf-8"))
    schema = _load_yaml(feature_schema_yaml)
    return summary, schema


def _output_paths(dataset_dir: Path) -> dict[str, Path]:
    return {
        "distribution_dir": dataset_dir / "plots/distributions",
        "qc_dir": dataset_dir / "plots/qc",
        "analysis_dir": dataset_dir / "analysis",
    }


def _prepare_outputs(dataset_dir: Path, *, overwrite: bool) -> dict[str, Path]:
    paths = _output_paths(dataset_dir)
    existing = [path for path in paths.values() if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Analysis outputs already exist; pass --overwrite to regenerate: "
            + ", ".join(str(path) for path in existing)
        )
    for path in existing:
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def _schema_feature_categories(schema: dict[str, Any]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    categories = schema.get("feature_categories", {}) or {}
    for category, raw in categories.items():
        if category in EXCLUDED_FEATURE_CATEGORIES:
            continue
        variables = raw.get("variables", []) if isinstance(raw, dict) else []
        out[str(category)] = [str(value) for value in variables]
    return out


def _is_bool_like(series: pd.Series) -> bool:
    values = series.dropna().astype(str).str.strip().str.lower()
    values = values[values != ""]
    if values.empty:
        return False
    return bool(values.isin(BOOL_TRUE | BOOL_FALSE).all())


def _coerce_numeric(series: pd.Series) -> pd.Series:
    if _is_bool_like(series):
        return (
            series.astype(str)
            .str.strip()
            .str.lower()
            .map({**dict.fromkeys(BOOL_TRUE, 1.0), **dict.fromkeys(BOOL_FALSE, 0.0)})
        )
    return pd.to_numeric(series, errors="coerce")


def _analysis_feature_categories(
    summary: pd.DataFrame,
    schema: dict[str, Any],
) -> tuple[dict[str, list[str]], pd.DataFrame]:
    schema_categories = _schema_feature_categories(schema)
    numeric_columns: dict[str, pd.Series] = {}
    categories: dict[str, list[str]] = {}

    for category, variables in schema_categories.items():
        category_features: list[str] = []
        for feature in variables:
            if feature not in summary.columns or feature in PATH_AND_STATUS_COLUMNS:
                continue
            if feature in QC_COLUMNS:
                continue
            numeric = _coerce_numeric(summary[feature])
            if numeric.notna().any() or summary[feature].notna().any():
                numeric_columns[feature] = numeric
                category_features.append(feature)
        if category_features:
            categories[category] = category_features

    if not categories:
        fallback_features: list[str] = []
        for feature in summary.columns:
            if feature in PATH_AND_STATUS_COLUMNS or feature in QC_COLUMNS:
                continue
            if feature in {"sample_id", "dataset_id", "condition_tag"}:
                continue
            numeric = _coerce_numeric(summary[feature])
            if numeric.notna().any():
                numeric_columns[feature] = numeric
                fallback_features.append(feature)
        categories["features"] = fallback_features

    numeric_df = pd.DataFrame(numeric_columns, index=summary.index)
    return categories, numeric_df


def _parse_bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(BOOL_TRUE)


def _summary_base(summary: pd.DataFrame) -> pd.DataFrame:
    base = summary.copy()
    base["n_flagella"] = pd.to_numeric(base["n_flagella"], errors="coerce")
    base["quality_class"] = base["quality_class"].astype(str)
    base["use_for_analysis_bool"] = _parse_bool_series(base["use_for_analysis"])
    return base


def _stats(values: pd.Series) -> dict[str, float]:
    clean = values.dropna()
    if clean.empty:
        return {
            "mean": math.nan,
            "std": math.nan,
            "median": math.nan,
            "q25": math.nan,
            "q75": math.nan,
            "min": math.nan,
            "max": math.nan,
        }
    return {
        "mean": float(clean.mean()),
        "std": float(clean.std(ddof=1)) if len(clean) > 1 else 0.0,
        "median": float(clean.median()),
        "q25": float(clean.quantile(0.25)),
        "q75": float(clean.quantile(0.75)),
        "min": float(clean.min()),
        "max": float(clean.max()),
    }


def build_feature_summary(
    summary: pd.DataFrame,
    numeric_df: pd.DataFrame,
    features: list[str],
) -> pd.DataFrame:
    base = _summary_base(summary)
    rows: list[dict[str, Any]] = []
    for feature in features:
        for (n_flagella, quality_class), group in base.groupby(
            ["n_flagella", "quality_class"],
            dropna=False,
            sort=True,
        ):
            values = numeric_df.loc[group.index, feature]
            row = {
                "feature": feature,
                "n_flagella": n_flagella,
                "quality_class": quality_class,
                "count": int(values.notna().sum()),
                "nan_count": int(values.isna().sum()),
            }
            row.update(_stats(values))
            rows.append(row)
    return pd.DataFrame(rows)


def build_nan_summary(
    summary: pd.DataFrame,
    numeric_df: pd.DataFrame,
    features: list[str],
) -> pd.DataFrame:
    base = _summary_base(summary)
    rows: list[dict[str, Any]] = []
    for feature in features:
        for n_flagella, group in base.groupby("n_flagella", dropna=False, sort=True):
            values = numeric_df.loc[group.index, feature]
            total = int(len(values))
            nan_count = int(values.isna().sum())
            rows.append(
                {
                    "feature": feature,
                    "n_flagella": n_flagella,
                    "total_count": total,
                    "nan_count": nan_count,
                    "nan_fraction": nan_count / total if total else math.nan,
                }
            )
        values = numeric_df[feature]
        total = int(len(values))
        nan_count = int(values.isna().sum())
        rows.append(
            {
                "feature": feature,
                "n_flagella": "all",
                "total_count": total,
                "nan_count": nan_count,
                "nan_fraction": nan_count / total if total else math.nan,
            }
        )
    return pd.DataFrame(rows)


def build_quality_summary(summary: pd.DataFrame) -> pd.DataFrame:
    base = _summary_base(summary)
    rows: list[dict[str, Any]] = []
    for column in QC_COLUMNS + ["missing_value_count"]:
        if column not in base.columns:
            continue
        values = base[column].astype(str)
        for n_flagella, group in base.assign(_value=values).groupby(
            "n_flagella",
            dropna=False,
            sort=True,
        ):
            total = int(len(group))
            counts = group["_value"].value_counts(dropna=False).sort_index()
            for category, count in counts.items():
                rows.append(
                    {
                        "metric": column,
                        "n_flagella": n_flagella,
                        "category": category,
                        "count": int(count),
                        "total_count": total,
                        "fraction": int(count) / total if total else math.nan,
                    }
                )
    return pd.DataFrame(rows)


def build_feature_screening(
    summary: pd.DataFrame,
    numeric_df: pd.DataFrame,
    features: list[str],
) -> pd.DataFrame:
    base = _summary_base(summary)
    screen_base = base[base["use_for_analysis_bool"]]
    if screen_base.empty:
        screen_base = base

    rows: list[dict[str, Any]] = []
    for feature in features:
        values = numeric_df.loc[screen_base.index, feature]
        by_group = values.groupby(screen_base["n_flagella"])
        group_means = by_group.mean().dropna()
        group_stds = by_group.std(ddof=0).dropna()
        group_counts = by_group.count()
        between_range = (
            float(group_means.max() - group_means.min())
            if len(group_means) >= 2
            else math.nan
        )
        within_std = float(group_stds.mean()) if not group_stds.empty else math.nan
        if math.isfinite(within_std) and within_std > 0.0:
            ratio = between_range / within_std
        elif math.isfinite(between_range):
            ratio = math.inf
        else:
            ratio = math.nan

        notes: list[str] = []
        if len(group_means) < 2:
            notes.append("insufficient_groups")
        if not group_counts.empty and int(group_counts.max()) < 2:
            notes.append("single_sample_per_group")
        nan_fraction = float(numeric_df[feature].isna().mean())
        if nan_fraction >= 0.5:
            notes.append("high_nan_fraction")
        if math.isinf(ratio):
            notes.append("within_group_std_zero")

        rows.append(
            {
                "feature": feature,
                "between_group_range": between_range,
                "within_group_std_mean": within_std,
                "range_over_within_std": ratio,
                "nan_count": int(numeric_df[feature].isna().sum()),
                "notes": ";".join(notes) if notes else "",
            }
        )
    return pd.DataFrame(rows).sort_values(
        by=["range_over_within_std", "between_group_range"],
        ascending=[False, False],
        na_position="last",
    )


def _safe_name(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return safe.strip("_") or "plot"


def _feature_label(feature: str) -> str:
    return feature.replace("_", "\n")


def _plot_category(
    *,
    summary: pd.DataFrame,
    numeric_df: pd.DataFrame,
    category: str,
    features: list[str],
    output_path: Path,
    analysis_only: bool,
) -> None:
    base = _summary_base(summary)
    if analysis_only:
        base = base[base["use_for_analysis_bool"]]
    if base.empty:
        return

    ncols = 2 if len(features) > 1 else 1
    nrows = math.ceil(len(features) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.2 * ncols, max(3.2 * nrows, 3.2)),
        squeeze=False,
        constrained_layout=True,
    )
    handles: dict[str, Any] = {}
    for ax, feature in zip(axes.flat, features, strict=False):
        values = numeric_df.loc[base.index, feature]
        for offset_idx, (row_idx, row) in enumerate(base.iterrows()):
            value = values.loc[row_idx]
            if pd.isna(value):
                continue
            quality = str(row["quality_class"])
            color = QUALITY_COLORS.get(quality, "#111827")
            marker = "o" if bool(row["use_for_analysis_bool"]) else "x"
            x_value = float(row["n_flagella"]) + ((offset_idx % 7) - 3) * 0.025
            plotted = ax.scatter(
                x_value,
                float(value),
                color=color,
                marker=marker,
                s=42,
                alpha=0.9,
                label=quality,
            )
            handles.setdefault(quality, plotted)
        ax.set_title(_feature_label(feature), fontsize=9)
        ax.set_xlabel("n_flagella")
        ax.grid(True, axis="y", alpha=0.25)
        ax.tick_params(axis="x", labelsize=8)
        ax.tick_params(axis="y", labelsize=8)
    for ax in axes.flat[len(features) :]:
        ax.axis("off")

    title_suffix = "use_for_analysis only" if analysis_only else "all samples"
    fig.suptitle(f"{category}: {title_suffix}", fontsize=12)
    if handles:
        fig.legend(
            handles.values(),
            handles.keys(),
            loc="outside lower center",
            ncols=min(len(handles), 4),
            fontsize=8,
        )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_distribution_plots(
    *,
    summary: pd.DataFrame,
    numeric_df: pd.DataFrame,
    categories: dict[str, list[str]],
    output_dir: Path,
) -> list[Path]:
    paths: list[Path] = []
    for category, features in categories.items():
        existing_features = [feature for feature in features if feature in numeric_df]
        if not existing_features:
            continue
        for analysis_only, suffix in (
            (False, "all_samples"),
            (True, "use_for_analysis_true"),
        ):
            output_path = output_dir / f"{_safe_name(category)}_{suffix}.png"
            _plot_category(
                summary=summary,
                numeric_df=numeric_df,
                category=category,
                features=existing_features,
                output_path=output_path,
                analysis_only=analysis_only,
            )
            if output_path.is_file():
                paths.append(output_path)
    return paths


def _plot_categorical_qc(summary: pd.DataFrame, column: str, output_path: Path) -> None:
    base = _summary_base(summary)
    table = pd.crosstab(base["n_flagella"], base[column].astype(str))
    if table.empty:
        return
    ax = table.plot(kind="bar", stacked=True, figsize=(7.0, 4.2), width=0.75)
    ax.set_xlabel("n_flagella")
    ax.set_ylabel("sample count")
    ax.set_title(column)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(title=column, fontsize=8, title_fontsize=8)
    ax.figure.tight_layout()
    ax.figure.savefig(output_path, dpi=180)
    plt.close(ax.figure)


def _plot_missing_values(summary: pd.DataFrame, output_path: Path) -> None:
    if "missing_value_count" not in summary.columns:
        return
    base = _summary_base(summary)
    values = pd.to_numeric(base["missing_value_count"], errors="coerce")
    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    for idx, row in base.iterrows():
        value = values.loc[idx]
        if pd.isna(value):
            continue
        color = QUALITY_COLORS.get(str(row["quality_class"]), "#111827")
        marker = "o" if bool(row["use_for_analysis_bool"]) else "x"
        ax.scatter(row["n_flagella"], value, color=color, marker=marker, s=48)
    ax.set_xlabel("n_flagella")
    ax.set_ylabel("missing_value_count")
    ax.set_title("missing_value_count by n_flagella")
    ax.grid(True, axis="y", alpha=0.25)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_qc_plots(summary: pd.DataFrame, output_dir: Path) -> list[Path]:
    paths: list[Path] = []
    for column in QC_COLUMNS:
        if column not in summary.columns:
            continue
        output_path = output_dir / f"{_safe_name(column)}_by_n_flagella.png"
        _plot_categorical_qc(summary, column, output_path)
        if output_path.is_file():
            paths.append(output_path)
    missing_path = output_dir / "missing_value_count_by_n_flagella.png"
    _plot_missing_values(summary, missing_path)
    if missing_path.is_file():
        paths.append(missing_path)
    return paths


def analyze_dataset(
    *,
    dataset_dir: Path,
    overwrite: bool = False,
) -> dict[str, Any]:
    summary, schema = _read_inputs(dataset_dir)
    paths = _prepare_outputs(dataset_dir, overwrite=overwrite)
    categories, numeric_df = _analysis_feature_categories(summary, schema)
    features = [
        feature
        for category_features in categories.values()
        for feature in category_features
        if feature in numeric_df
    ]
    features = list(dict.fromkeys(features))
    if not features:
        raise ValueError("No numeric feature columns are available for analysis")

    feature_summary = build_feature_summary(summary, numeric_df, features)
    nan_summary = build_nan_summary(summary, numeric_df, features)
    quality_summary = build_quality_summary(summary)
    feature_screening = build_feature_screening(summary, numeric_df, features)

    analysis_dir = paths["analysis_dir"]
    feature_summary_path = analysis_dir / "feature_summary_by_n_flagella.csv"
    nan_summary_path = analysis_dir / "nan_summary.csv"
    quality_summary_path = analysis_dir / "quality_summary.csv"
    feature_screening_path = analysis_dir / "feature_screening_summary.csv"
    feature_summary.to_csv(feature_summary_path, index=False)
    nan_summary.to_csv(nan_summary_path, index=False)
    quality_summary.to_csv(quality_summary_path, index=False)
    feature_screening.to_csv(feature_screening_path, index=False)

    distribution_plots = write_distribution_plots(
        summary=summary,
        numeric_df=numeric_df,
        categories=categories,
        output_dir=paths["distribution_dir"],
    )
    qc_plots = write_qc_plots(summary, paths["qc_dir"])

    return {
        "dataset_dir": dataset_dir,
        "analysis_csvs": [
            feature_summary_path,
            nan_summary_path,
            quality_summary_path,
            feature_screening_path,
        ],
        "distribution_plots": distribution_plots,
        "qc_plots": qc_plots,
        "features": features,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", type=Path, default=None)
    parser.add_argument("--dataset-id", type=str, default=None)
    parser.add_argument(
        "--dataset-root",
        type=Path,
        default=DEFAULT_DATASET_ROOT,
        help="Root used with --dataset-id.",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    try:
        dataset_dir = resolve_dataset_dir(
            dataset_dir=args.dataset_dir,
            dataset_id=args.dataset_id,
            dataset_root=args.dataset_root,
        )
        result = analyze_dataset(dataset_dir=dataset_dir, overwrite=args.overwrite)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc

    print(result["dataset_dir"])
    for path in result["analysis_csvs"]:
        print(path)
    print(f"distribution_plots={len(result['distribution_plots'])}")
    print(f"qc_plots={len(result['qc_plots'])}")


if __name__ == "__main__":
    main()
