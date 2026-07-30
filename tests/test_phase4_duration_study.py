from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import subprocess

import numpy as np
import pytest

from flagella_estimation.phase4.duration_study import (
    DurationStudyConfig,
    analyze_duration_seed_study,
    load_duration_study_config,
)
from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.sim.core import SimulationState


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_fixture(dataset_dir: Path) -> None:
    summary_rows: list[dict[str, object]] = []
    for attach_seed in (0, 1):
        for phase_seed in (0, 1):
            sample_id = f"nf01_as{attach_seed:03d}_ps{phase_seed:03d}"
            raw_dir = dataset_dir / "raw" / sample_id
            states = [
                SimulationState(
                    t=index * 0.1,
                    position_um=(index * (0.01 + attach_seed * 0.002), 0.0, 0.0),
                    quaternion=(0.0, 0.0, 0.0, 1.0),
                    velocity_um_s=(0.1 + phase_seed * 0.02, 0.0, 0.0),
                    omega_rad_s=(0.0, 0.0, 0.2 + attach_seed * 0.1),
                    bead_positions_um=np.asarray(
                        [[0.0, 0.0, 0.0], [0.5, 0.1 + phase_seed * 0.05, 0.0]]
                    ),
                )
                for index in range(21)
            ]
            save_state_archive(raw_dir / "state_archive.npz", states)
            step_rows = []
            for index, state in enumerate(states):
                strict_pass = not (
                    attach_seed == 1 and phase_seed == 1 and state.t >= 1.0
                )
                step_rows.append(
                    {
                        "t_s": state.t,
                        "finite_pass": True,
                        "shape_pass_nonbody_strict": strict_pass,
                        "shape_pass_nonbody_hook_len_relaxed": True,
                        "first_fail_category_nonbody_strict": (
                            "flag" if not strict_pass else "none"
                        ),
                        "flag_helix_axis_alignment_order": 0.7 + attach_seed * 0.1,
                        "bundle_axis_vs_body_axis_angle_deg": 10.0 + phase_seed,
                        "hook_len_rel_err_max": index * 0.001,
                    }
                )
            _write_csv(raw_dir / "step_summary.csv", step_rows)
            run_pass = not (attach_seed == 1 and phase_seed == 1)
            summary_rows.append(
                {
                    "sample_id": sample_id,
                    "dataset_id": "v1_r1_duration_3s",
                    "n_flagella": 1,
                    "attach_seed": attach_seed,
                    "phase_seed": phase_seed,
                    "quality_class": "strict_pass" if run_pass else "relaxed_pass",
                    "shape_pass": run_pass,
                    "relaxed_pass": True,
                    "raw_dir": raw_dir,
                }
            )
    _write_csv(dataset_dir / "summary.csv", summary_rows)
    (dataset_dir / "dataset_manifest.json").write_text(
        json.dumps(
            {
                "effective_campaign_config": {
                    "metadata": {
                        "dataset_version": "v1",
                        "dataset_revision": "r1",
                    },
                    "dataset": {
                        "dataset_version": "v1",
                        "dataset_revision": "r1",
                    },
                }
            }
        ),
        encoding="utf-8",
    )


@pytest.mark.light
def test_duration_study_writes_window_qc_and_seed_artifacts(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "dataset"
    output_dir = tmp_path / "study"
    _write_fixture(dataset_dir)

    result = analyze_duration_seed_study(
        DurationStudyConfig(
            dataset_dir=dataset_dir,
            output_dir=output_dir,
            durations_s=(0.5, 1.0),
            frame_rate_hz=2.0,
            crop_size_px=32,
            pixel_size_um=0.1,
        )
    )

    assert result == output_dir
    rows_3d = list(
        csv.DictReader(
            (output_dir / "window_features_3d.csv").open(
                "r", encoding="utf-8", newline=""
            )
        )
    )
    rows_2d = list(
        csv.DictReader(
            (output_dir / "window_features_2d.csv").open(
                "r", encoding="utf-8", newline=""
            )
        )
    )
    assert rows_3d and len(rows_3d) == len(rows_2d)
    assert {row["dataset_revision"] for row in rows_3d} == {"r1"}
    assert {row["group_key"] for row in rows_3d} == {
        f"phase2:v1:nf01_as{attach:03d}_ps{phase:03d}"
        for attach in (0, 1)
        for phase in (0, 1)
    }
    failed = [
        row
        for row in rows_3d
        if row["run_id"] == "nf01_as001_ps001" and float(row["t_start_s"]) >= 1.0
    ]
    assert failed
    assert {row["window_shape_pass"] for row in failed} == {"False"}
    assert {row["window_first_fail_category"] for row in failed} == {"flag"}

    seed_rows = list(
        csv.DictReader(
            (output_dir / "seed_effects.csv").open("r", encoding="utf-8", newline="")
        )
    )
    assert seed_rows
    assert {row["qc_scope"] for row in seed_rows} == {
        "all_labeled",
        "strict_run_and_windows",
    }
    assert {
        row["run_count"] for row in seed_rows if row["qc_scope"] == "all_labeled"
    } == {"4"}
    assert all(0.0 <= float(row["attach_eta_squared"]) <= 1.0 for row in seed_rows)
    run_means = list(
        csv.DictReader(
            (output_dir / "run_feature_means.csv").open(
                "r", encoding="utf-8", newline=""
            )
        )
    )
    assert run_means
    assert {
        row["all_windows_shape_pass"]
        for row in run_means
        if row["run_id"] == "nf01_as001_ps001"
    } == {"False"}
    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["dataset_version"] == "v1"
    assert manifest["dataset_revision"] == "r1"
    assert manifest["within_run_windows_are_independent"] is False
    assert len(manifest["outputs"]["plots"]) >= 6
    plot_names = {Path(path).name for path in manifest["outputs"]["plots"]}
    assert "2d_centroid_step_mean_seed_heatmap.png" in plot_names
    assert not any(
        name.startswith("2d_centroid_step_mean_seed_heatmap_") for name in plot_names
    )


@pytest.mark.light
def test_duration_study_config_accepts_key_value_overrides(tmp_path: Path) -> None:
    config_path = tmp_path / "study.yaml"
    config_path.write_text(
        "\n".join(
            [
                "dataset_dir: input",
                "output_dir: output",
                "study:",
                "  durations_s: [0.5]",
                "clip:",
                "  frame_rate_hz: 25.0",
            ]
        ),
        encoding="utf-8",
    )

    cfg = load_duration_study_config(
        config_path,
        [
            "study.durations_s=[0.25,0.5,1.0]",
            "clip.frame_rate_hz=20.0",
            "overwrite=true",
        ],
    )

    assert cfg.durations_s == (0.25, 0.5, 1.0)
    assert cfg.frame_rate_hz == pytest.approx(20.0)
    assert cfg.overwrite is True


@pytest.mark.light
def test_duration_study_shell_full_path_supports_empty_overwrite(
    tmp_path: Path,
) -> None:
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    command_log = tmp_path / "commands.log"
    fake_uv = fake_bin / "uv"
    fake_uv.write_text(
        '#!/usr/bin/env bash\nprintf \'%s\\n\' "$*" >> "${COMMAND_LOG}"\n',
        encoding="utf-8",
    )
    fake_uv.chmod(0o755)
    env = {
        **os.environ,
        "PATH": f"{fake_bin}:{os.environ['PATH']}",
        "COMMAND_LOG": str(command_log),
    }
    script = (
        Path(__file__).parents[1]
        / "scripts"
        / "04_phase4"
        / "run_duration_seed_study.sh"
    )

    completed = subprocess.run(
        ["bash", str(script), "full"],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    commands = command_log.read_text(encoding="utf-8").splitlines()
    assert len(commands) == 3
    assert all("overwrite=true" not in command for command in commands)
    assert "output_dir=" not in commands[2]
