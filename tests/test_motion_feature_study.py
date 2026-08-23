from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.analysis.motion_feature_study import (
    MotionFeatureStudyConfig,
    _axis_velocity,
    _basis,
    _observation_frame_rate,
    _plot_rows,
    _plot_series,
    _time_windows,
    analyze_motion_feature_study,
    load_config,
)
from sim_swim.sim.core import SimulationState


def _csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _fixture(
    root: Path,
    *,
    broken_n4: bool = False,
    quaternion: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 1.0),
    fps_out_2d: float | None = 2.0,
) -> None:
    conditions = []
    for n in (1, 4):
        sample = f"as000__ps000__nf{n:02d}"
        raw = root / sample
        raw.mkdir(parents=True)
        states = [
            SimulationState(
                t=index * 0.2,
                position_um=(index * 0.2, 0.0, 0.0),
                quaternion=quaternion,
                velocity_um_s=(1.0, 0.0, 0.0),
                omega_rad_s=(0.0, 0.0, 0.0),
                bead_positions_um=np.asarray([[0.0, 0.0, 0.0]]),
            )
            for index in range(11)
        ]
        save_state_archive(raw / "state_archive.npz", states)
        _csv(
            raw / "step_summary.csv",
            [
                {
                    "t_s": state.t,
                    "finite_pass": True,
                    "shape_pass_nonbody_strict": not (
                        broken_n4 and n == 4 and round(state.t, 1) == 0.6
                    ),
                    "first_fail_category_nonbody_strict": "none",
                }
                for state in states
            ],
        )
        axis_rows = []
        for state in states:
            for flag_id in range(n):
                axis_rows.append(
                    {
                        "t_s": state.t,
                        "flag_id": flag_id,
                        "axis_dir_x": -1.0,
                        "axis_dir_y": 0.0,
                        "axis_dir_z": 0.0,
                    }
                )
        _csv(raw / "flag_helix_axis_diagnostics.csv", axis_rows)
        conditions.append(
            {
                "condition_id": sample,
                "output_dir": str(raw),
                "config_overrides": (
                    {"output_sampling": {"fps_out_2d": fps_out_2d}}
                    if fps_out_2d is not None
                    else {}
                ),
                "axis_values": {"n_flagella": n, "attach_seed": 0, "phase_seed": 0},
            }
        )
    (root / "run_manifest.json").write_text(
        json.dumps({"conditions": conditions}), encoding="utf-8"
    )


@pytest.mark.light
def test_motion_feature_study_raw_campaign_writes_feature_contract(
    tmp_path: Path,
) -> None:
    run_dir, output_dir = tmp_path / "run", tmp_path / "analysis"
    _fixture(run_dir)
    analyze_motion_feature_study(
        MotionFeatureStudyConfig(
            run_dir=run_dir,
            output_dir=output_dir,
            durations_s=(1.0,),
        )
    )
    rows = list(csv.DictReader((output_dir / "window_features_2d.csv").open()))
    assert rows and {row["n_flagella"] for row in rows} == {"1", "4"}
    assert all(float(row["mean_speed_um_s"]) == pytest.approx(1.0) for row in rows)
    series3d = list(csv.DictReader((output_dir / "time_series_3d.csv").open()))
    series2d = list(csv.DictReader((output_dir / "time_series_2d.csv").open()))
    assert len(series3d) == 22
    assert len(series2d) == 10
    assert {row["sampling_source"] for row in series3d} == {"simulation_step"}
    assert {row["sampling_source"] for row in series2d} == {"2d_output_frame"}
    assert all(float(row["observation_frame_rate_hz"]) == 2.0 for row in series2d)
    assert (
        json.loads((output_dir / "manifest.json").read_text())["observability"][
            "2d_body_axis_angular_velocity"
        ]
        == "pixel_observable"
    )
    provenance = json.loads((output_dir / "manifest.json").read_text())[
        "analysis_provenance"
    ]
    assert provenance["analysis_git_commit"]
    assert provenance["python_version"]
    assert provenance["dependencies"]["numpy"]
    assert (output_dir / "time_series" / "3D_speed.png").is_file()
    assert (output_dir / "windows" / "2D_speed.png").is_file()


@pytest.mark.light
def test_motion_feature_study_excludes_strict_failing_n4_from_plots(
    tmp_path: Path,
) -> None:
    run_dir, output_dir = tmp_path / "run", tmp_path / "analysis"
    _fixture(run_dir, broken_n4=True)
    analyze_motion_feature_study(
        MotionFeatureStudyConfig(
            run_dir=run_dir,
            output_dir=output_dir,
            durations_s=(1.0,),
            observation_frame_rate_hz=2.0,
        )
    )
    manifest = json.loads((output_dir / "manifest.json").read_text())
    assert manifest["plot_exclusions"]["n4_strict_failure"] == ["as000__ps000__nf04"]
    rows = list(csv.DictReader((output_dir / "window_features_3d.csv").open()))
    assert {row["n_flagella"] for row in rows} == {"1", "4"}


@pytest.mark.light
def test_motion_feature_axes_are_axial_and_projection_basis_is_orthonormal() -> None:
    vectors = np.asarray([[1.0, 0.0], [-1.0, 0.0]])
    assert _axis_velocity(vectors, np.asarray([0.0, 1.0]))[1] == pytest.approx(0.0)
    basis = _basis(((2.0, 0.0, 0.0), (1.0, 3.0, 0.0)))
    assert basis @ basis.T == pytest.approx(np.eye(2))


@pytest.mark.light
def test_motion_feature_windows_use_physical_time_boundaries() -> None:
    t = np.arange(0.0, 2.00004 + 1e-12, 4e-5)
    windows = _time_windows(t, 0.5)
    assert [(start, end) for start, end, _ in windows] == pytest.approx(
        [(0.0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0)]
    )
    assert [sl.stop - sl.start for _, _, sl in windows] == [12500] * 4


@pytest.mark.light
def test_motion_feature_study_rejects_ambiguous_frame_rate_config(
    tmp_path: Path,
) -> None:
    path = tmp_path / "analysis.yaml"
    path.write_text("study:\n  frame_rate_hz: 25\n", encoding="utf-8")
    with pytest.raises(ValueError, match="study.frame_rate_hz is ambiguous"):
        load_config(path)
    path.write_text("plot:\n  flagella_axis_time_bin_s: 0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="flagella_axis_time_bin_s must be > 0"):
        load_config(path)
    path.write_text("overwrite: 'false'\n", encoding="utf-8")
    with pytest.raises(ValueError, match="overwrite must be a boolean"):
        load_config(path)
    with pytest.raises(ValueError, match="overwrite must be a boolean"):
        load_config(None, ["overwrite=flase"])


@pytest.mark.light
def test_motion_feature_explicit_observation_rate_must_match_manifest(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    _fixture(run_dir)
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    cfg = MotionFeatureStudyConfig(
        run_dir=run_dir, output_dir=tmp_path / "analysis", observation_frame_rate_hz=5.0
    )
    with pytest.raises(ValueError, match="must match"):
        _observation_frame_rate(manifest, cfg)
    _fixture(tmp_path / "run_no_rate", fps_out_2d=None)
    no_rate_manifest = json.loads(
        (tmp_path / "run_no_rate" / "run_manifest.json").read_text()
    )
    assert _observation_frame_rate(cfg=cfg, manifest=no_rate_manifest) == 5.0


@pytest.mark.light
def test_motion_feature_plot_distinguishes_individual_and_mean(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from matplotlib.axes import Axes

    original = Axes.plot
    styles: list[tuple[str | None, float | None]] = []

    def capture(self: Axes, *args: object, **kwargs: object) -> object:
        styles.append((kwargs.get("linestyle"), kwargs.get("linewidth")))
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "plot", capture)
    rows = [
        {
            "sample_id": sample,
            "n_flagella": 1,
            "t_s": t,
            "speed_um_s": 1.0,
            "body_axis_angular_velocity_rad_s": 1.0,
            "mean_flagella_axis_angular_velocity_rad_s": 1.0,
            "body_flagella_axis_angle_deg": 1.0,
        }
        for sample in ("a", "b")
        for t in (0.0, 1.0)
    ]
    _plot_series(rows, tmp_path, "3D", 0.02)
    assert ("--", 0.9) in styles
    assert (None, 2.4) in styles


@pytest.mark.light
def test_flagella_axis_plot_is_binned_without_changing_other_features() -> None:
    rows = [
        {
            "sample_id": "sample",
            "n_flagella": 1,
            "t_s": t,
            "mean_flagella_axis_angular_velocity_rad_s": value,
            "speed_um_s": value,
        }
        for t, value in ((0.0, 1.0), (0.01, 3.0), (0.02, 5.0))
    ]
    reduced = _plot_rows(rows, "mean_flagella_axis_angular_velocity_rad_s", "t_s", 0.02)
    assert [row["mean_flagella_axis_angular_velocity_rad_s"] for row in reduced] == [
        pytest.approx(2.0),
        pytest.approx(5.0),
    ]
    assert _plot_rows(rows, "speed_um_s", "t_s", 0.02) is rows
    assert (
        _plot_rows(rows, "mean_flagella_axis_angular_velocity_rad_s", "t_s", None)
        is rows
    )


@pytest.mark.light
def test_flagella_axis_plot_bin_tolerates_decimal_boundaries() -> None:
    rows = [
        {
            "sample_id": "sample",
            "n_flagella": 1,
            "t_s": t,
            "mean_flagella_axis_angular_velocity_rad_s": value,
        }
        for t, value in ((0.56, 1.0), (0.58, 3.0), (0.6, 5.0))
    ]
    reduced = _plot_rows(rows, "mean_flagella_axis_angular_velocity_rad_s", "t_s", 0.02)
    assert [row["t_s"] for row in reduced] == pytest.approx([0.56, 0.58, 0.6])


@pytest.mark.light
def test_motion_feature_2d_axis_is_nan_when_silhouette_is_circular(
    tmp_path: Path,
) -> None:
    run_dir, output_dir = tmp_path / "run", tmp_path / "analysis"
    # Rotate the body long axis from world X onto world Z; it is circular in XY.
    _fixture(
        run_dir,
        quaternion=(0.0, math.sqrt(0.5), 0.0, math.sqrt(0.5)),
    )
    analyze_motion_feature_study(
        MotionFeatureStudyConfig(
            run_dir=run_dir, output_dir=output_dir, durations_s=(1.0,)
        )
    )
    rows = list(csv.DictReader((output_dir / "time_series_2d.csv").open()))
    assert all(
        math.isnan(float(row["body_axis_angular_velocity_rad_s"])) for row in rows
    )
    assert all(math.isnan(float(row["body_flagella_axis_angle_deg"])) for row in rows)


@pytest.mark.light
def test_motion_feature_camera_basis_controls_2d_silhouette(
    tmp_path: Path,
) -> None:
    run_dir, output_dir = tmp_path / "run", tmp_path / "analysis"
    # The world-Z body axis is circular in XY but observable in an XZ camera.
    _fixture(
        run_dir,
        quaternion=(0.0, math.sqrt(0.5), 0.0, math.sqrt(0.5)),
    )
    analyze_motion_feature_study(
        MotionFeatureStudyConfig(
            run_dir=run_dir,
            output_dir=output_dir,
            durations_s=(1.0,),
            projection_basis=((1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
        )
    )
    rows = list(csv.DictReader((output_dir / "time_series_2d.csv").open()))
    assert all(
        math.isfinite(float(row["body_flagella_axis_angle_deg"])) for row in rows
    )


@pytest.mark.light
def test_motion_feature_defaults_to_unbinned_plot_and_protects_existing_output(
    tmp_path: Path,
) -> None:
    run_dir, output_dir = tmp_path / "run", tmp_path / "analysis"
    _fixture(run_dir)
    config = MotionFeatureStudyConfig(run_dir=run_dir, output_dir=output_dir)
    assert config.flagella_axis_plot_bin_s is None
    analyze_motion_feature_study(config)
    with pytest.raises(FileExistsError, match="--overwrite"):
        analyze_motion_feature_study(config)
    analyze_motion_feature_study(
        MotionFeatureStudyConfig(run_dir=run_dir, output_dir=output_dir, overwrite=True)
    )


@pytest.mark.light
@pytest.mark.parametrize("target", ("run", "condition"))
def test_motion_feature_overwrite_rejects_raw_source_overlap(
    tmp_path: Path, target: str
) -> None:
    run_dir = tmp_path / "run"
    _fixture(run_dir)
    output_dir = run_dir if target == "run" else run_dir / "as000__ps000__nf01"
    with pytest.raises(ValueError, match="output_dir must not"):
        analyze_motion_feature_study(
            MotionFeatureStudyConfig(
                run_dir=run_dir, output_dir=output_dir, overwrite=True
            )
        )
    assert (run_dir / "as000__ps000__nf01" / "state_archive.npz").is_file()
