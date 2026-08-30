"""Compare Issue #204 and #215 body--flagella axis-angle artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

import matplotlib.pyplot as plt
import numpy as np
import yaml


DOMAINS = ("3d", "2d")
ANGLE = "body_flagella_axis_angle_deg"


@dataclass(frozen=True)
class AxisAngleComparisonConfig:
    reference_2s_dir: Path
    reference_5s_dir: Path
    output_dir: Path
    compare_until_s: float = 2.0
    post_start_s: float = 2.0
    overwrite: bool = False


def load_config(path: Path) -> AxisAngleComparisonConfig:
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if type(raw.get("overwrite", False)) is not bool:
        raise ValueError("overwrite must be a boolean")
    return AxisAngleComparisonConfig(
        reference_2s_dir=Path(str(raw["reference_2s_dir"])),
        reference_5s_dir=Path(str(raw["reference_5s_dir"])),
        output_dir=Path(str(raw["output_dir"])),
        compare_until_s=float(raw.get("compare_until_s", 2.0)),
        post_start_s=float(raw.get("post_start_s", 2.0)),
        overwrite=bool(raw.get("overwrite", False)),
    )


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    keys = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def _finite(value: Any) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return parsed if math.isfinite(parsed) else float("nan")


def _truth(value: Any) -> bool:
    return str(value).lower() in {"true", "1"}


def _analysis_dir(reference_dir: Path) -> Path:
    path = reference_dir / "analysis" / "motion_features"
    if not path.is_dir():
        raise FileNotFoundError(f"missing motion-feature analysis: {path}")
    return path


def _series(reference_dir: Path, domain: str) -> list[dict[str, str]]:
    return _read_rows(_analysis_dir(reference_dir) / f"time_series_{domain}.csv")


def _windows(reference_dir: Path, domain: str) -> list[dict[str, str]]:
    return _read_rows(_analysis_dir(reference_dir) / f"window_features_{domain}.csv")


def _comparison_rows(
    rows_2s: list[dict[str, str]],
    rows_5s: list[dict[str, str]],
    domain: str,
    until_s: float,
) -> list[dict[str, Any]]:
    def indexed(rows: list[dict[str, str]]) -> dict[tuple[str, float], dict[str, str]]:
        return {
            (row["sample_id"], round(_finite(row["t_s"]), 9)): row
            for row in rows
            if _finite(row["t_s"]) <= until_s + 1e-12
        }

    left, right = indexed(rows_2s), indexed(rows_5s)
    results: list[dict[str, Any]] = []
    for sample_id in sorted({key[0] for key in left} | {key[0] for key in right}):
        common = sorted(
            key for key in left.keys() & right.keys() if key[0] == sample_id
        )
        source = next((row for row in rows_5s if row["sample_id"] == sample_id), None)
        values = np.asarray(
            [_finite(left[key][ANGLE]) - _finite(right[key][ANGLE]) for key in common]
        )
        valid = values[np.isfinite(values)]
        results.append(
            {
                "domain": domain.upper(),
                "sample_id": sample_id,
                "n_flagella": source.get("n_flagella", "") if source else "",
                "attach_seed": source.get("attach_seed", "") if source else "",
                "phase_seed": source.get("phase_seed", "") if source else "",
                "two_s_run_strict_pass": _truth(
                    next(
                        (
                            row.get("run_strict_pass")
                            for row in rows_2s
                            if row["sample_id"] == sample_id
                        ),
                        False,
                    )
                ),
                "five_s_run_strict_pass": _truth(
                    next(
                        (
                            row.get("run_strict_pass")
                            for row in rows_5s
                            if row["sample_id"] == sample_id
                        ),
                        False,
                    )
                ),
                "two_s_available": any(key[0] == sample_id for key in left),
                "five_s_available": any(key[0] == sample_id for key in right),
                "common_timepoint_count": len(common),
                "valid_pair_count": len(valid),
                "mean_abs_difference_deg": float(np.mean(np.abs(valid)))
                if len(valid)
                else float("nan"),
                "max_abs_difference_deg": float(np.max(np.abs(valid)))
                if len(valid)
                else float("nan"),
                "rmse_difference_deg": float(np.sqrt(np.mean(valid**2)))
                if len(valid)
                else float("nan"),
            }
        )
        results[-1]["comparison_strict_pass"] = bool(
            results[-1]["two_s_run_strict_pass"]
            and results[-1]["five_s_run_strict_pass"]
        )
    return results


def _post_windows(
    rows: list[dict[str, str]], domain: str, start_s: float
) -> list[dict[str, Any]]:
    output = []
    for row in rows:
        if _finite(row["t_start_s"]) >= start_s - 1e-12:
            output.append({"domain": domain.upper(), **row})
    return output


def _usable(row: dict[str, str]) -> bool:
    return _truth(row.get("run_strict_pass"))


def _mean_by_time(rows: list[dict[str, str]]) -> tuple[list[float], list[float]]:
    """Aggregate seed trajectories once per time point, not once per lookup."""
    grouped: dict[float, list[float]] = {}
    for row in rows:
        grouped.setdefault(round(_finite(row["t_s"]), 9), []).append(
            _finite(row[ANGLE])
        )
    times = sorted(grouped)
    means = [
        float(np.nanmean(grouped[time]))
        if np.isfinite(grouped[time]).any()
        else float("nan")
        for time in times
    ]
    return times, means


def _plot_trajectories(rows: list[dict[str, str]], domain: str, output: Path) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(12, 7), sharex=True, sharey=True)
    for n, axis in zip((1, 2, 3, 4), axes.flat, strict=True):
        subset = [row for row in rows if int(row["n_flagella"]) == n and _usable(row)]
        by_sample: dict[str, list[dict[str, str]]] = {}
        for row in subset:
            by_sample.setdefault(row["sample_id"], []).append(row)
        for values in by_sample.values():
            values.sort(key=lambda row: _finite(row["t_s"]))
            axis.plot(
                [_finite(row["t_s"]) for row in values],
                [_finite(row[ANGLE]) for row in values],
                color="C0",
                alpha=0.22,
            )
        times, means = _mean_by_time(subset)
        axis.plot(times, means, color="C0", linewidth=2.2, label="5 s mean")
        axis.axvline(2.0, color="black", linestyle="--", linewidth=0.8)
        axis.set(
            title=f"n={n}", xlabel="time (s)", ylabel="body--flagella axis angle (deg)"
        )
        axis.grid(alpha=0.25)
    axes[0, 0].legend()
    figure.suptitle(f"{domain.upper()} body--flagella axis angle (5 s campaign)")
    figure.tight_layout()
    figure.savefig(output, dpi=150)
    plt.close(figure)


def _plot_overlay(
    rows_2s: list[dict[str, str]],
    rows_5s: list[dict[str, str]],
    domain: str,
    until_s: float,
    output: Path,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(12, 7), sharex=True, sharey=True)
    for n, axis in zip((1, 2, 3, 4), axes.flat, strict=True):
        for label, rows, color in (("2 s", rows_2s, "C1"), ("5 s", rows_5s, "C0")):
            subset = [
                row
                for row in rows
                if int(row["n_flagella"]) == n
                and _usable(row)
                and _finite(row["t_s"]) <= until_s + 1e-12
            ]
            times, means = _mean_by_time(subset)
            axis.plot(times, means, color=color, linewidth=2.0, label=label)
        axis.set(
            title=f"n={n}", xlabel="time (s)", ylabel="body--flagella axis angle (deg)"
        )
        axis.grid(alpha=0.25)
    axes[0, 0].legend()
    figure.suptitle(f"{domain.upper()} 0--{until_s:g} s reference overlay")
    figure.tight_layout()
    figure.savefig(output, dpi=150)
    plt.close(figure)


def analyze(cfg: AxisAngleComparisonConfig) -> Path:
    if cfg.output_dir.exists():
        if not cfg.overwrite:
            raise FileExistsError(cfg.output_dir)
        shutil.rmtree(cfg.output_dir)
    cfg.output_dir.mkdir(parents=True)
    agreement: list[dict[str, Any]] = []
    post: list[dict[str, Any]] = []
    plots: list[str] = []
    for domain in DOMAINS:
        rows_2s, rows_5s = (
            _series(cfg.reference_2s_dir, domain),
            _series(cfg.reference_5s_dir, domain),
        )
        agreement.extend(
            _comparison_rows(rows_2s, rows_5s, domain, cfg.compare_until_s)
        )
        post.extend(
            _post_windows(
                _windows(cfg.reference_5s_dir, domain), domain, cfg.post_start_s
            )
        )
        trajectory = f"{domain.upper()}_axis_angle_0_5s.png"
        overlay = f"{domain.upper()}_axis_angle_0_{cfg.compare_until_s:g}s_overlay.png"
        _plot_trajectories(rows_5s, domain, cfg.output_dir / trajectory)
        _plot_overlay(
            rows_2s, rows_5s, domain, cfg.compare_until_s, cfg.output_dir / overlay
        )
        plots.extend((trajectory, overlay))
    _write_rows(cfg.output_dir / "overlap_0_2s_agreement.csv", agreement)
    _write_rows(cfg.output_dir / "post_2s_window_axis_angle.csv", post)
    valid = [row for row in agreement if int(row["valid_pair_count"]) > 0]
    summary = {
        "kind": "phase2_issue215_axis_angle_comparison",
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "reference_2s_dir": str(cfg.reference_2s_dir),
        "reference_5s_dir": str(cfg.reference_5s_dir),
        "compare_interval_s": [0.0, cfg.compare_until_s],
        "post_window_interval_s": [cfg.post_start_s, 5.0],
        "agreement_row_count": len(agreement),
        "agreement_with_pairs_count": len(valid),
        "post_window_row_count": len(post),
        "strict_nonpass_post_window_count": sum(
            not _truth(row.get("run_strict_pass")) for row in post
        ),
        "outputs": [
            "overlap_0_2s_agreement.csv",
            "post_2s_window_axis_angle.csv",
            *plots,
            "review.md",
        ],
    }
    (cfg.output_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    review = [
        "# Issue #215 軸角度比較レビュー",
        "",
        "本資料は診断結果の記録であり、model・dataset・ML policyの採択判断を含まない。",
        "",
        f"- 0–{cfg.compare_until_s:g} s 対応条件: {len(valid)}/{len(agreement)}",
        f"- {cfg.post_start_s:g}–5 s window行数: {len(post)}",
        f"- strict non-PASS window行数: {summary['strict_nonpass_post_window_count']}",
        "",
        "詳細値はCSV、時系列の確認はPNGを参照する。",
    ]
    (cfg.output_dir / "review.md").write_text(
        "\n".join(review) + "\n", encoding="utf-8"
    )
    return cfg.output_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    cfg = load_config(args.config)
    if args.overwrite:
        cfg = AxisAngleComparisonConfig(**{**cfg.__dict__, "overwrite": True})
    print(analyze(cfg))


if __name__ == "__main__":
    main()
