from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest


def _load_plot_script():
    script_path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "01_simulate_swimming"
        / "plot_phase2_82_hook_overstretch_heatmap.py"
    )
    spec = importlib.util.spec_from_file_location("phase2_82_heatmap", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_phase2_82_body_first_heatmap_outputs_files(
    tmp_path: Path, monkeypatch
) -> None:
    script = _load_plot_script()
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
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_phase2_82_hook_overstretch_heatmap",
            "--summary-csv",
            str(summary_csv),
            "--mode",
            "body-first-grid",
            "--output-dir",
            str(output_dir),
        ],
    )

    script.main()

    assert (output_dir / "phase2_82_hook_overstretch_heatmap.csv").is_file()
    assert (output_dir / "phase2_82_first_fail_category_heatmap.png").is_file()
    assert (output_dir / "phase2_82_shape_pass_fail_heatmap.png").is_file()
    assert (output_dir / "phase2_82_hook_len_rel_err_max_heatmap.png").is_file()
    assert (
        output_dir / "phase2_82_local_attach_first_vs_body_axis_err_deg_heatmap.png"
    ).is_file()


def test_phase2_82_first_second_heatmap_outputs_files(
    tmp_path: Path, monkeypatch
) -> None:
    script = _load_plot_script()
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
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_phase2_82_hook_overstretch_heatmap",
            "--summary-csv",
            str(summary_csv),
            "--mode",
            "first-second-grid",
            "--output-dir",
            str(output_dir),
        ],
    )

    script.main()

    assert (output_dir / "phase2_82_hook_overstretch_heatmap.csv").is_file()
    assert (output_dir / "phase2_82_first_fail_category_heatmap.png").is_file()
    assert (output_dir / "phase2_82_local_first_second_rel_err_heatmap.png").is_file()


def test_phase2_82_attach_frame_heatmap_outputs_files(
    tmp_path: Path, monkeypatch
) -> None:
    script = _load_plot_script()
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
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_phase2_82_hook_overstretch_heatmap",
            "--summary-csv",
            str(summary_csv),
            "--mode",
            "attach-frame-grid",
            "--output-dir",
            str(output_dir),
        ],
    )

    script.main()

    assert (output_dir / "phase2_82_hook_overstretch_heatmap.csv").is_file()
    assert (output_dir / "phase2_82_first_fail_category_heatmap.png").is_file()
    assert (
        output_dir / "phase2_82_local_attach_frame_position_rel_err_heatmap.png"
    ).is_file()


def test_phase2_82_heatmap_missing_summary_lists_candidates(
    tmp_path: Path, monkeypatch
) -> None:
    script = _load_plot_script()
    phase2_dir = tmp_path / "outputs" / "phase2_82"
    candidate_dir = phase2_dir / "first_second_grid_af3_axis1p25"
    candidate_dir.mkdir(parents=True)
    candidate_csv = candidate_dir / "phase2_82_hook_scale_sweep_summary.csv"
    candidate_csv.write_text(
        "condition_id,mode,final_shape_pass_nonbody\n",
        encoding="utf-8",
    )
    output_dir = phase2_dir / "first_second_grid" / "plots"
    missing_csv = (
        phase2_dir / "first_second_grid" / "phase2_82_hook_scale_sweep_summary.csv"
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_phase2_82_hook_overstretch_heatmap",
            "--summary-csv",
            str(missing_csv),
            "--mode",
            "first-second-grid",
            "--output-dir",
            str(output_dir),
        ],
    )

    with pytest.raises(SystemExit) as excinfo:
        script.main()

    message = str(excinfo.value)
    assert f"Summary CSV not found: {missing_csv}" in message
    assert str(candidate_csv) in message
    assert not output_dir.exists()
