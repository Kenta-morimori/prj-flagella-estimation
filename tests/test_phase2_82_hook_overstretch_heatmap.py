from __future__ import annotations

import csv
from pathlib import Path

import pytest

from sim_swim.analysis.heatmaps import hook_overstretch as script


def test_phase2_82_category_rank_uses_first_fail_category() -> None:
    assert script._category_rank(
        {
            "final_shape_pass_nonbody": "True",
            "final_first_fail_category_nonbody": "none",
            "first_fail_category_nonbody": "hook",
            "first_fail_t_s": "0.12",
        }
    ) == script.CATEGORY_ORDER.index("hook")
    assert script._category_rank(
        {
            "final_shape_pass_nonbody": "False",
            "final_first_fail_category_nonbody": "hook",
            "first_fail_category_nonbody": "flag",
            "first_fail_t_s": "0.23",
        }
    ) == script.CATEGORY_ORDER.index("flag")


def test_phase2_82_category_rank_keeps_legacy_summary_fallback() -> None:
    assert script._category_rank(
        {
            "final_shape_pass_nonbody": "True",
            "final_first_fail_category_nonbody": "none",
        }
    ) == script.CATEGORY_ORDER.index("none")
    assert script._category_rank(
        {
            "final_shape_pass_nonbody": "False",
            "final_first_fail_category_nonbody": "hook",
        }
    ) == script.CATEGORY_ORDER.index("hook")
    assert script._category_rank(
        {
            "final_shape_pass_nonbody": "True",
            "final_first_fail_category_nonbody": "none",
            "first_fail_category_nonbody": "none",
            "first_fail_t_s": "",
        }
    ) == script.CATEGORY_ORDER.index("none")
    assert script._category_rank(
        {
            "final_shape_pass_nonbody": "False",
            "final_first_fail_category_nonbody": "none",
            "first_fail_category_nonbody": "none",
            "first_fail_t_s": "0.34",
        }
    ) == script.CATEGORY_ORDER.index("finite")


def test_phase2_82_body_first_heatmap_outputs_files(tmp_path: Path) -> None:
    summary_csv = tmp_path / "summary.csv"
    fields = [
        "condition_id",
        "mode",
        "final_shape_pass_nonbody",
        "final_first_fail_category_nonbody",
        "local_attach_first_spring_scale",
        "local_attach_first_body_axis_angle_scale",
        "local_first_second_spring_scale",
        "hook_len_rel_err_max",
        "local_attach_first_rel_err",
        "local_attach_first_vs_body_axis_err_deg",
        "local_first_second_rel_err",
        "net_abs_flag_helix_spin_revolutions",
        "flag_helix_spin_direction_consistency",
    ]
    rows = [
        {
            "condition_id": "af1_axis1_fs1",
            "mode": "body-first-grid",
            "final_shape_pass_nonbody": "True",
            "final_first_fail_category_nonbody": "none",
            "local_attach_first_spring_scale": "1.0",
            "local_attach_first_body_axis_angle_scale": "1.0",
            "local_first_second_spring_scale": "1.0",
            "hook_len_rel_err_max": "0.4",
            "local_attach_first_rel_err": "0.4",
            "local_attach_first_vs_body_axis_err_deg": "20.0",
            "local_first_second_rel_err": "0.1",
            "net_abs_flag_helix_spin_revolutions": "0.05",
            "flag_helix_spin_direction_consistency": "0.5",
        },
        {
            "condition_id": "af2_axis1_fs1",
            "mode": "body-first-grid",
            "final_shape_pass_nonbody": "False",
            "final_first_fail_category_nonbody": "hook",
            "local_attach_first_spring_scale": "2.0",
            "local_attach_first_body_axis_angle_scale": "1.0",
            "local_first_second_spring_scale": "1.0",
            "hook_len_rel_err_max": "1.2",
            "local_attach_first_rel_err": "1.2",
            "local_attach_first_vs_body_axis_err_deg": "25.0",
            "local_first_second_rel_err": "0.2",
            "net_abs_flag_helix_spin_revolutions": "0.06",
            "flag_helix_spin_direction_consistency": "0.6",
        },
        {
            "condition_id": "af1_axis2_fs1",
            "mode": "body-first-grid",
            "final_shape_pass_nonbody": "True",
            "final_first_fail_category_nonbody": "none",
            "local_attach_first_spring_scale": "1.0",
            "local_attach_first_body_axis_angle_scale": "2.0",
            "local_first_second_spring_scale": "1.0",
            "hook_len_rel_err_max": "0.2",
            "local_attach_first_rel_err": "0.2",
            "local_attach_first_vs_body_axis_err_deg": "2.0",
            "local_first_second_rel_err": "0.1",
            "net_abs_flag_helix_spin_revolutions": "0.04",
            "flag_helix_spin_direction_consistency": "0.4",
        },
    ]
    with summary_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    output_dir = tmp_path / "plots"
    script.main(
        [
            "--summary-csv",
            str(summary_csv),
            "--mode",
            "body-first-grid",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert (output_dir / "heatmap_data.csv").is_file()
    assert (output_dir / "first_fail_category_heatmap.png").is_file()
    assert (output_dir / "shape_pass_fail_heatmap.png").is_file()
    assert (output_dir / "hook_len_rel_err_max_heatmap.png").is_file()
    assert (
        output_dir / "local_attach_first_vs_body_axis_err_deg_heatmap.png"
    ).is_file()


def test_phase2_82_first_second_heatmap_outputs_files(tmp_path: Path) -> None:
    summary_csv = tmp_path / "summary.csv"
    fields = [
        "condition_id",
        "mode",
        "final_shape_pass_nonbody",
        "final_first_fail_category_nonbody",
        "local_attach_first_spring_scale",
        "local_attach_first_body_axis_angle_scale",
        "local_first_second_spring_scale",
        "hook_len_rel_err_max",
        "local_attach_first_rel_err",
        "local_attach_first_vs_body_axis_err_deg",
        "local_first_second_rel_err",
    ]
    rows = [
        {
            "condition_id": "af2_axis2_fs1",
            "mode": "first-second-grid",
            "final_shape_pass_nonbody": "True",
            "final_first_fail_category_nonbody": "none",
            "local_attach_first_spring_scale": "2.0",
            "local_attach_first_body_axis_angle_scale": "2.0",
            "local_first_second_spring_scale": "1.0",
            "hook_len_rel_err_max": "0.3",
            "local_attach_first_rel_err": "0.3",
            "local_attach_first_vs_body_axis_err_deg": "3.0",
            "local_first_second_rel_err": "0.2",
        },
        {
            "condition_id": "af2_axis2_fs2",
            "mode": "first-second-grid",
            "final_shape_pass_nonbody": "True",
            "final_first_fail_category_nonbody": "none",
            "local_attach_first_spring_scale": "2.0",
            "local_attach_first_body_axis_angle_scale": "2.0",
            "local_first_second_spring_scale": "2.0",
            "hook_len_rel_err_max": "0.25",
            "local_attach_first_rel_err": "0.25",
            "local_attach_first_vs_body_axis_err_deg": "2.0",
            "local_first_second_rel_err": "0.1",
        },
    ]
    with summary_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    output_dir = tmp_path / "plots"
    script.main(
        [
            "--summary-csv",
            str(summary_csv),
            "--mode",
            "first-second-grid",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert (output_dir / "heatmap_data.csv").is_file()
    assert (output_dir / "first_fail_category_heatmap.png").is_file()
    assert (output_dir / "local_first_second_rel_err_heatmap.png").is_file()


def test_phase2_82_attach_frame_heatmap_outputs_files(tmp_path: Path) -> None:
    summary_csv = tmp_path / "summary.csv"
    fields = [
        "condition_id",
        "mode",
        "final_shape_pass_nonbody",
        "final_first_fail_category_nonbody",
        "local_attach_first_spring_scale",
        "local_attach_first_body_axis_angle_scale",
        "local_first_second_spring_scale",
        "local_attach_frame_position_scale",
        "local_attach_frame_tangent_scale",
        "hook_len_rel_err_max",
        "local_attach_first_rel_err",
        "local_attach_first_vs_body_axis_err_deg",
        "local_first_second_rel_err",
        "local_attach_frame_position_rel_err",
        "local_attach_frame_position_angle_err_deg",
        "local_attach_frame_tangent_angle_err_deg",
    ]
    rows = [
        {
            "condition_id": "af3_axis1p25_fs1p25_fp1_ft1",
            "mode": "attach-frame-grid",
            "final_shape_pass_nonbody": "True",
            "final_first_fail_category_nonbody": "none",
            "local_attach_first_spring_scale": "3.0",
            "local_attach_first_body_axis_angle_scale": "1.25",
            "local_first_second_spring_scale": "1.25",
            "local_attach_frame_position_scale": "1.0",
            "local_attach_frame_tangent_scale": "1.0",
            "hook_len_rel_err_max": "0.3",
            "local_attach_first_rel_err": "0.3",
            "local_attach_first_vs_body_axis_err_deg": "2.0",
            "local_first_second_rel_err": "0.2",
            "local_attach_frame_position_rel_err": "0.1",
            "local_attach_frame_position_angle_err_deg": "3.0",
            "local_attach_frame_tangent_angle_err_deg": "4.0",
        },
        {
            "condition_id": "af3_axis1p25_fs1p25_fp2_ft1",
            "mode": "attach-frame-grid",
            "final_shape_pass_nonbody": "False",
            "final_first_fail_category_nonbody": "hook",
            "local_attach_first_spring_scale": "3.0",
            "local_attach_first_body_axis_angle_scale": "1.25",
            "local_first_second_spring_scale": "1.25",
            "local_attach_frame_position_scale": "2.0",
            "local_attach_frame_tangent_scale": "1.0",
            "hook_len_rel_err_max": "1.1",
            "local_attach_first_rel_err": "0.9",
            "local_attach_first_vs_body_axis_err_deg": "3.0",
            "local_first_second_rel_err": "0.25",
            "local_attach_frame_position_rel_err": "0.4",
            "local_attach_frame_position_angle_err_deg": "5.0",
            "local_attach_frame_tangent_angle_err_deg": "6.0",
        },
    ]
    with summary_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    output_dir = tmp_path / "plots"
    script.main(
        [
            "--summary-csv",
            str(summary_csv),
            "--mode",
            "attach-frame-grid",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert (output_dir / "heatmap_data.csv").is_file()
    assert (output_dir / "first_fail_category_heatmap.png").is_file()
    assert (output_dir / "local_attach_frame_position_rel_err_heatmap.png").is_file()


def test_phase2_82_position_only_heatmap_outputs_files(tmp_path: Path) -> None:
    summary_csv = tmp_path / "summary.csv"
    fields = [
        "condition_id",
        "mode",
        "final_shape_pass_nonbody",
        "final_first_fail_category_nonbody",
        "first_fail_category_nonbody",
        "first_fail_t_s",
        "local_attach_first_spring_scale",
        "local_attach_first_body_axis_angle_scale",
        "local_first_second_spring_scale",
        "local_attach_frame_position_scale",
        "local_attach_frame_tangent_scale",
        "hook_len_rel_err_max",
        "local_attach_first_rel_err",
        "local_attach_first_vs_body_axis_err_deg",
        "local_first_second_rel_err",
        "local_attach_frame_position_rel_err",
        "local_attach_frame_position_angle_err_deg",
        "local_attach_frame_tangent_angle_err_deg",
        "max_flag_bond_rel_err",
        "body_roll_net_abs_revolutions",
        "axis_center_to_body_roll_ratio_mean",
    ]
    rows = [
        {
            "condition_id": "no_frame",
            "mode": "position-only-grid",
            "final_shape_pass_nonbody": "False",
            "final_first_fail_category_nonbody": "hook",
            "first_fail_category_nonbody": "hook",
            "first_fail_t_s": "0.0449",
            "local_attach_first_spring_scale": "1.0",
            "local_attach_first_body_axis_angle_scale": "1.0",
            "local_first_second_spring_scale": "1.0",
            "local_attach_frame_position_scale": "1.0",
            "local_attach_frame_tangent_scale": "1.0",
            "hook_len_rel_err_max": "1.3",
            "local_attach_first_rel_err": "1.3",
            "local_attach_first_vs_body_axis_err_deg": "12.0",
            "local_first_second_rel_err": "0.2",
            "local_attach_frame_position_rel_err": "0.8",
            "local_attach_frame_position_angle_err_deg": "20.0",
            "local_attach_frame_tangent_angle_err_deg": "9.0",
            "max_flag_bond_rel_err": "0.5",
            "body_roll_net_abs_revolutions": "0.01",
            "axis_center_to_body_roll_ratio_mean": "100.0",
        },
        {
            "condition_id": "fp1p25",
            "mode": "position-only-grid",
            "final_shape_pass_nonbody": "False",
            "final_first_fail_category_nonbody": "flag",
            "first_fail_category_nonbody": "flag",
            "first_fail_t_s": "0.3295",
            "local_attach_first_spring_scale": "1.0",
            "local_attach_first_body_axis_angle_scale": "1.0",
            "local_first_second_spring_scale": "1.0",
            "local_attach_frame_position_scale": "1.25",
            "local_attach_frame_tangent_scale": "1.0",
            "hook_len_rel_err_max": "0.12",
            "local_attach_first_rel_err": "0.12",
            "local_attach_first_vs_body_axis_err_deg": "4.0",
            "local_first_second_rel_err": "0.3",
            "local_attach_frame_position_rel_err": "0.1",
            "local_attach_frame_position_angle_err_deg": "5.0",
            "local_attach_frame_tangent_angle_err_deg": "7.0",
            "max_flag_bond_rel_err": "1.1",
            "body_roll_net_abs_revolutions": "0.0135",
            "axis_center_to_body_roll_ratio_mean": "151.2",
        },
    ]
    with summary_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    output_dir = tmp_path / "plots"
    script.main(
        [
            "--summary-csv",
            str(summary_csv),
            "--mode",
            "position-only-grid",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert (output_dir / "heatmap_data.csv").is_file()
    assert (output_dir / "first_fail_category_heatmap.png").is_file()
    assert (output_dir / "shape_pass_fail_heatmap.png").is_file()
    assert (output_dir / "max_flag_bond_rel_err_heatmap.png").is_file()
    assert (output_dir / "body_roll_net_abs_revolutions_heatmap.png").is_file()
    assert (output_dir / "axis_center_to_body_roll_ratio_mean_heatmap.png").is_file()


def test_phase2_82_heatmap_missing_summary_lists_candidates(tmp_path: Path) -> None:
    phase2_dir = tmp_path / "outputs" / "hook_overstretch"
    candidate_dir = phase2_dir / "first_second_grid_af3_axis1p25"
    candidate_dir.mkdir(parents=True)
    candidate_csv = candidate_dir / "summary.csv"
    candidate_csv.write_text(
        "condition_id,mode,final_shape_pass_nonbody\n",
        encoding="utf-8",
    )
    output_dir = phase2_dir / "first_second_grid" / "plots"
    missing_csv = phase2_dir / "first_second_grid" / "summary.csv"

    with pytest.raises(SystemExit) as excinfo:
        script.main(
            [
                "--summary-csv",
                str(missing_csv),
                "--mode",
                "first-second-grid",
                "--output-dir",
                str(output_dir),
            ]
        )

    message = str(excinfo.value)
    assert f"Summary CSV not found: {missing_csv}" in message
    assert str(candidate_csv) in message
    assert not output_dir.exists()
