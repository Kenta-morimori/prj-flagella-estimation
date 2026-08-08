from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

from sim_swim.analysis.run_summary import (
    MAX_RUN_SUMMARY_BYTES,
    build_run_summary,
    write_run_summary,
)


def _write_csv(
    path: Path, fieldnames: list[str], rows: list[dict[str, object]]
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _step_rows(passes: list[bool]) -> list[dict[str, object]]:
    return [
        {
            "step": index,
            "t_s": index * 0.1,
            "finite_pass": True,
            "shape_pass_nonbody": value,
            "first_fail_category_nonbody": "none" if value else "hook",
            "hook_len_rel_err_max": index / 10,
            "flag_bond_rel_err_max": index / 100,
            "flag_bend_err_max_deg": index,
            "flag_torsion_err_max_deg": index * 2,
        }
        for index, value in enumerate(passes)
    ]


def _run_dir(tmp_path: Path, passes: list[bool]) -> Path:
    run_dir = tmp_path / "condition_001"
    run_dir.mkdir()
    rows = _step_rows(passes)
    _write_csv(run_dir / "step_summary.csv", list(rows[0]), rows)
    return run_dir


def test_run_summary_distinguishes_short_fail_recovery_and_missing_body(
    tmp_path: Path,
) -> None:
    run_dir = _run_dir(tmp_path, [True, False, True, False, False, False, True])

    summary = build_run_summary(
        run_dir,
        time_manifest={"total_steps": 7, "final_step_summary_t_s": 0.6},
    )

    assert summary["execution"]["status"] == "completed"
    gate = summary["gates"]["shape_nonbody"]
    assert gate["final_pass"] is True
    assert gate["recovered_after_fail"] is True
    assert gate["episode_count"] == 2
    assert gate["episodes"][0]["persistent_observed"] is False
    assert gate["episodes"][0]["end_t_s"] == pytest.approx(0.1)
    assert gate["episodes"][0]["next_observed_pass_t_s"] == pytest.approx(0.2)
    assert gate["episodes"][1]["persistent_observed"] is True
    assert summary["gates"]["shape_body"] == {
        "status": "unavailable",
        "reason": "body_constraint_diagnostics.csv is missing or empty",
    }


def test_run_summary_marks_partial_and_handles_body_gate(tmp_path: Path) -> None:
    run_dir = _run_dir(tmp_path, [True, True])
    body_rows = [
        {
            "step": 0,
            "t_s": 0.0,
            "body_spring_max_stretch_ratio": 0.9,
            "body_bend_max_error_deg": 10.0,
            "body_centerline_max_deviation_um": 0.1,
            "body_triangle_area_min": 2.0,
        },
        {
            "step": 1,
            "t_s": 0.1,
            "body_spring_max_stretch_ratio": 1.1,
            "body_bend_max_error_deg": 10.0,
            "body_centerline_max_deviation_um": 0.1,
            "body_triangle_area_min": 2.0,
        },
    ]
    _write_csv(
        run_dir / "body_constraint_diagnostics.csv", list(body_rows[0]), body_rows
    )

    summary = build_run_summary(
        run_dir,
        time_manifest={"total_steps": 3, "final_step_summary_t_s": 0.2},
    )

    assert summary["execution"]["status"] == "partial"
    body = summary["gates"]["shape_body"]
    assert body["final_pass"] is False
    assert body["episodes"][0]["category_counts"] == {"body_spring": 1}


def test_run_summary_caps_episode_storage_and_preserves_file(tmp_path: Path) -> None:
    passes = [value for _ in range(40) for value in (True, False)]
    run_dir = _run_dir(tmp_path, passes)
    source = (run_dir / "step_summary.csv").read_bytes()

    output = write_run_summary(run_dir)
    summary = json.loads(output.read_text(encoding="utf-8"))

    gate = summary["gates"]["shape_nonbody"]
    assert gate["episode_count"] == 40
    assert gate["stored_episode_count"] == 32
    assert gate["omitted_episode_count"] == 8
    assert (run_dir / "step_summary.csv").read_bytes() == source
    assert output.stat().st_size <= MAX_RUN_SUMMARY_BYTES
    with pytest.raises(FileExistsError):
        write_run_summary(run_dir)


def _load_inspector() -> object:
    path = Path("scripts/02_phase2_analysis/inspect_step_summary.py")
    spec = importlib.util.spec_from_file_location("inspect_step_summary", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_inspector_returns_only_requested_columns_and_rows(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    run_dir = _run_dir(tmp_path, [True, False, False, True])
    write_run_summary(run_dir)
    inspector = _load_inspector()

    inspector.main(
        [
            "--input-dir",
            str(run_dir),
            "--gate",
            "shape_nonbody",
            "--episode",
            "1",
            "--columns",
            "t_s,shape_pass_nonbody",
            "--max-rows",
            "1",
        ]
    )

    result = json.loads(capsys.readouterr().out)
    assert result["returned_row_count"] == 1
    assert result["truncated"] is True
    assert list(result["rows"][0]) == ["t_s", "shape_pass_nonbody"]
