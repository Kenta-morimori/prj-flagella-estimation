from __future__ import annotations

import importlib.util
import csv
import sys
from pathlib import Path


def _load_plot_script():
    script_path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "plot_motor_scale_collapse_heatmap.py"
    )
    spec = importlib.util.spec_from_file_location("plot_script", script_path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_plot_motor_scale_collapse_heatmap(tmp_path: Path, monkeypatch) -> None:
    plot_script = _load_plot_script()
    summary_csv = tmp_path / "local_hook_scale_sweep_summary.csv"
    summary_csv.write_text(
        "\n".join(
            [
                (
                    "target,torque_Nm,value,output_dir,pos_all_finite,any_nan,"
                    "any_inf,finite_pass,shape_pass_nonbody,body_shape_pass,"
                    "shape_pass,first_fail_category"
                ),
                (
                    "local_hook_scale,4e-21,8.0,run_a,True,False,False,True,"
                    "False,False,False,hook"
                ),
                (
                    "local_hook_scale,4e-21,4.0,run_b,True,False,False,True,"
                    "False,False,False,body_area"
                ),
                (
                    "local_hook_scale,4e-21,0.0,run_c,True,False,False,True,"
                    "True,True,True,none"
                ),
            ]
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_motor_scale_collapse_heatmap",
            "--summary-csv",
            str(summary_csv),
            "--output-dir",
            str(tmp_path / "plots"),
        ],
    )

    plot_script.main()

    output_dir = tmp_path / "plots"
    normalized_csv = output_dir / "local_hook_scale_sweep_summary_collapse_matrix.csv"
    category_png = output_dir / "local_hook_scale_sweep_summary_category_heatmap.png"
    combined_pass_fail_png = (
        output_dir / "local_hook_scale_sweep_summary_combined_pass_fail_heatmap.png"
    )
    flagella_png = (
        output_dir / "local_hook_scale_sweep_summary_flagella_pass_fail_heatmap.png"
    )

    assert normalized_csv.is_file()
    assert category_png.is_file()
    assert combined_pass_fail_png.is_file()
    assert not flagella_png.exists()

    with normalized_csv.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 3
    assert rows[0]["torque_Nm"] == "4e-21"
    assert rows[0]["first_fail_category"] == "hook"
    assert rows[1]["first_fail_category"] == "body_area"
    assert rows[2]["first_fail_category"] == "none"


def test_plot_motor_scale_collapse_heatmap_emits_flagella_plot(
    tmp_path: Path, monkeypatch
) -> None:
    plot_script = _load_plot_script()
    summary_csv = tmp_path / "local_hook_scale_sweep_summary.csv"
    summary_csv.write_text(
        "\n".join(
            [
                (
                    "target,torque_Nm,value,output_dir,pos_all_finite,any_nan,"
                    "any_inf,finite_pass,shape_pass_nonbody,body_shape_pass,"
                    "shape_pass,first_fail_category"
                ),
                (
                    "local_hook_scale,4e-21,8.0,run_a,True,False,False,True,"
                    "False,False,False,flag"
                ),
                (
                    "local_hook_scale,4e-21,0.0,run_b,True,False,False,True,"
                    "True,True,True,none"
                ),
            ]
        ),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_motor_scale_collapse_heatmap",
            "--summary-csv",
            str(summary_csv),
            "--output-dir",
            str(tmp_path / "plots"),
        ],
    )

    plot_script.main()

    output_dir = tmp_path / "plots"
    flagella_png = (
        output_dir / "local_hook_scale_sweep_summary_flagella_pass_fail_heatmap.png"
    )

    assert flagella_png.is_file()
