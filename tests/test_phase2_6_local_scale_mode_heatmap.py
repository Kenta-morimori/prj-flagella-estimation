from __future__ import annotations

import csv
from pathlib import Path

from sim_swim.analysis.heatmaps import local_scale_mode


def test_phase2_6_local_scale_mode_heatmap_outputs_files(tmp_path: Path) -> None:
    summary_csv = tmp_path / "summary.csv"
    fields = [
        "stage",
        "target",
        "value",
        "torque_Nm",
        "shape_pass",
        "first_fail_category",
        "first_fail_reason",
        "net_abs_flag_helix_spin_revolutions",
        "flag_helix_spin_direction_consistency",
        "max_flag_bond_rel_err",
        "max_flag_bend_err_deg",
        "max_flag_torsion_err_deg",
    ]
    rows = [
        {
            "stage": "baseline_all_local_scales_1",
            "target": "all_local_scales",
            "value": "1.0",
            "torque_Nm": "3.0e-20",
            "shape_pass": "True",
            "first_fail_category": "none",
            "first_fail_reason": "none",
        },
        {
            "stage": "one_factor_scale",
            "target": "local_spring_scale",
            "value": "2.0",
            "torque_Nm": "3.0e-20",
            "shape_pass": "True",
            "first_fail_category": "none",
            "first_fail_reason": "none",
        },
        {
            "stage": "baseline_all_local_scales_1",
            "target": "all_local_scales",
            "value": "1.0",
            "torque_Nm": "3.5e-20",
            "shape_pass": "False",
            "first_fail_category": "flag",
            "first_fail_reason": "shape_pass_nonbody=False",
        },
        {
            "stage": "one_factor_scale",
            "target": "local_spring_scale",
            "value": "2.0",
            "torque_Nm": "3.5e-20",
            "shape_pass": "False",
            "first_fail_category": "flag",
            "first_fail_reason": "shape_pass_nonbody=False",
        },
    ]
    with summary_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    out_dir = tmp_path / "plots"
    local_scale_mode.main(
        [
            "--summary-csv",
            str(summary_csv),
            "--output-dir",
            str(out_dir),
        ],
    )

    normalized_csv = out_dir / "heatmap_data.csv"
    assert normalized_csv.is_file()
    assert (out_dir / "local_scale_mode_category_heatmap.png").is_file()
    assert (out_dir / "local_scale_mode_pass_fail_heatmap.png").is_file()

    with normalized_csv.open("r", encoding="utf-8", newline="") as handle:
        normalized_rows = list(csv.DictReader(handle))
    assert {(row["torque_Nm"], row["mode"]) for row in normalized_rows} == {
        ("3.0e-20", "all=1"),
        ("3.0e-20", "spring=2"),
        ("3.5e-20", "all=1"),
        ("3.5e-20", "spring=2"),
    }
