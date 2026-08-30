from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from sim_swim.analysis.axis_angle_comparison import (
    AxisAngleComparisonConfig,
    _usable,
    analyze,
)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _reference(root: Path, *, duration_s: int) -> None:
    series: list[dict[str, object]] = []
    windows: list[dict[str, object]] = []
    for n in (1, 2, 4):
        sample_id = f"as000__ps000__nf{n:02d}"
        strict = n == 1
        for t_s in range(duration_s + 1):
            series.append(
                {
                    "sample_id": sample_id,
                    "n_flagella": n,
                    "attach_seed": 0,
                    "phase_seed": 0,
                    "run_strict_pass": strict,
                    "t_s": t_s,
                    "body_flagella_axis_angle_deg": 10.0 + n + t_s,
                }
            )
        for t_start_s in range(duration_s):
            windows.append(
                {
                    "sample_id": sample_id,
                    "n_flagella": n,
                    "attach_seed": 0,
                    "phase_seed": 0,
                    "run_strict_pass": strict,
                    "requested_duration_s": 1.0,
                    "window_index": t_start_s,
                    "t_start_s": t_start_s,
                    "t_end_s": t_start_s + 1,
                    "window_shape_pass": strict,
                    "window_first_fail_t_s": "" if strict else t_start_s,
                    "window_first_fail_category": "none" if strict else "shape",
                    "mean_body_flagella_axis_angle_deg": 10.5 + n + t_start_s,
                }
            )
    analysis = root / "analysis" / "motion_features"
    for domain in ("3d", "2d"):
        _write_csv(analysis / f"time_series_{domain}.csv", series)
        _write_csv(analysis / f"window_features_{domain}.csv", windows)


@pytest.mark.light
def test_axis_angle_comparison_writes_overlap_post_window_and_review(
    tmp_path: Path,
) -> None:
    reference_2s, reference_5s, output = (
        tmp_path / "2s",
        tmp_path / "5s",
        tmp_path / "comparison",
    )
    _reference(reference_2s, duration_s=2)
    _reference(reference_5s, duration_s=5)

    analyze(
        AxisAngleComparisonConfig(
            reference_2s_dir=reference_2s,
            reference_5s_dir=reference_5s,
            output_dir=output,
        )
    )

    agreement = list(csv.DictReader((output / "overlap_0_2s_agreement.csv").open()))
    assert len(agreement) == 6
    assert {row["domain"] for row in agreement} == {"2D", "3D"}
    assert all(int(row["common_timepoint_count"]) == 3 for row in agreement)
    assert all(float(row["rmse_difference_deg"]) == 0.0 for row in agreement)
    assert {
        "two_s_run_strict_pass",
        "five_s_run_strict_pass",
        "comparison_strict_pass",
    } <= set(agreement[0])
    assert {row["comparison_strict_pass"] for row in agreement} == {"False", "True"}
    post = list(csv.DictReader((output / "post_2s_window_axis_angle.csv").open()))
    assert len(post) == 18
    assert all(float(row["t_start_s"]) >= 2.0 for row in post)
    summary = json.loads((output / "summary.json").read_text())
    assert summary["agreement_with_pairs_count"] == 6
    assert summary["strict_nonpass_post_window_count"] == 12
    assert (output / "3D_axis_angle_0_5s.png").is_file()
    assert (output / "2D_axis_angle_0_2s_overlay.png").is_file()
    assert "採択判断を含まない" in (output / "review.md").read_text()


@pytest.mark.light
def test_axis_angle_plot_eligibility_excludes_any_strict_nonpass() -> None:
    assert _usable({"n_flagella": "4", "run_strict_pass": "true"})
    assert not _usable({"n_flagella": "2", "run_strict_pass": "false"})
