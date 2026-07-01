from __future__ import annotations

import csv
from pathlib import Path

from sim_swim.analysis.heatmaps import dt_star_torque


def test_phase2_6_dt_star_torque_heatmap_outputs_files(tmp_path: Path) -> None:
    summary_csv = tmp_path / "summary.csv"
    fields = [
        "target",
        "value",
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
    rows = [
        {
            "target": "all_local_scales",
            "value": "1.0",
            "torque_Nm": "2.0e-20",
            "dt_star": "1.0e-4",
            "shape_pass": "True",
            "first_fail_category": "none",
            "first_fail_reason": "none",
        },
        {
            "target": "all_local_scales",
            "value": "1.0",
            "torque_Nm": "2.0e-20",
            "dt_star": "1.0e-3",
            "shape_pass": "False",
            "first_fail_category": "flag",
            "first_fail_reason": "shape_pass_nonbody=False",
        },
        {
            "target": "local_spring_scale",
            "value": "2.0",
            "torque_Nm": "2.0e-20",
            "dt_star": "1.0e-4",
            "shape_pass": "True",
            "first_fail_category": "none",
            "first_fail_reason": "diagnostic non-baseline row",
        },
    ]
    with summary_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    out_dir = tmp_path / "plots"
    dt_star_torque.main(
        [
            "--summary-csv",
            str(summary_csv),
            "--output-dir",
            str(out_dir),
        ],
    )

    normalized_csv = out_dir / "heatmap_data.csv"
    assert normalized_csv.is_file()
    assert (out_dir / "dt_star_torque_category_heatmap.png").is_file()
    assert (out_dir / "dt_star_torque_pass_fail_heatmap.png").is_file()

    with normalized_csv.open("r", encoding="utf-8", newline="") as handle:
        normalized_rows = list(csv.DictReader(handle))
    assert {(row["torque_Nm"], row["dt_star"]) for row in normalized_rows} == {
        ("2.0e-20", "1.0e-4"),
        ("2.0e-20", "1.0e-3"),
    }
