import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from sim_swim.analysis.issue203_composite_replay import (
    _write_frames_with_fallback,
    composite_manifest_path,
    nominal_segment_weights,
    reconstructed_segment_weights,
)
from sim_swim.analysis.issue203_torque_profile_comparison import _flag_axes, load_config
from sim_swim.analysis.torque_profile_dt_contact import (
    ContactConfig,
    _min_bead_distances,
    analyze as analyze_contact,
)
from sim_swim.analysis.motion_feature_study import (
    load_config as load_motion_feature_config,
)
from sim_swim.analysis.phase2_replay import _archive_path
from sim_swim.analysis.multi_run_campaign import (
    apply_campaign_cli_overrides,
    build_campaign_conditions,
    load_yaml,
)


def test_issue203_campaign_is_the_fixed_27_condition_uniform_grid() -> None:
    path = Path(
        "conf/phase2_multi_run/2010_project_uniform_torque_profile_2s_issue203.yaml"
    )
    config = load_yaml(path)
    conditions = build_campaign_conditions(apply_campaign_cli_overrides(config, []))
    assert len(conditions) == 27
    assert config["base_config"] == "conf/sim_swim_2010.yaml"
    assert config["base_overrides"]["motor.torque_distribution_profile"] == "uniform"
    assert (
        config["base_overrides"]["motor.force_distribution"]
        == "root_torque_segment_couples"
    )
    assert config["base_overrides"]["time.integration.dt_star"] == 1.0e-3


def test_linked_axis_preserves_correlated_seed_cases() -> None:
    config = {
        "base_overrides": {},
        "sweep": {
            "axes": {
                "seed_case": {
                    "keys": {
                        "attach_seed": "seed.attach_seed",
                        "phase_seed": "seed.phase_seed",
                        "n_flagella": "flagella.n_flagella",
                    },
                    "ids": ["as000__ps000__nf03", "as001__ps000__nf03"],
                    "values": [
                        {"attach_seed": 0, "phase_seed": 0, "n_flagella": 3},
                        {"attach_seed": 1, "phase_seed": 0, "n_flagella": 3},
                    ],
                },
                "dt_star": {
                    "key": "time.integration.dt_star",
                    "values": [1.0e-3, 1.0e-4],
                    "ids": ["dt1e3", "dt1e4"],
                },
            }
        },
    }

    conditions = build_campaign_conditions(apply_campaign_cli_overrides(config, []))

    assert len(conditions) == 4
    assert conditions[0]["axis_values"] == {
        "attach_seed": 0,
        "phase_seed": 0,
        "n_flagella": 3,
        "seed_case": "as000__ps000__nf03",
        "dt_star": 1.0e-3,
    }
    assert conditions[-1]["config_overrides"]["seed"]["attach_seed"] == 1


def test_issue203_uniform_weight_panel_uses_segments_and_sums_to_one() -> None:
    weights = nominal_segment_weights("uniform", 10)
    assert weights.shape == (10,)
    assert np.allclose(weights, 0.1)
    assert weights.sum() == 1.0


def test_issue203_diffusive_weight_panel_reconstructs_dynamic_segment_weights() -> None:
    weights = reconstructed_segment_weights(
        "diffusive",
        4,
        times_s=np.asarray([0.0, 0.01, 0.1]),
        dt_s=1.0e-3,
        torque_Nm=2.5e-20,
    )

    assert all(np.isclose(weight.sum(), 1.0) for weight in weights)
    assert np.allclose(weights[0], 0.25)
    assert weights[-1][0] > weights[-1][-1]


def test_issue203_flag_axis_metrics_accept_large_archive_layout() -> None:
    positions = np.zeros((3, 5, 3), dtype=float)
    positions[:, 2, 0] = 1.0
    positions[:, 4, 1] = 1.0

    mean, alignment, spread = _flag_axes(
        positions, [np.asarray([0, 1, 2]), np.asarray([0, 3, 4])]
    )

    assert mean.shape == (3, 3)
    assert np.all(np.isfinite(alignment))
    assert np.all(np.isfinite(spread))


def test_torque_profile_dt_contact_distance_screen_detects_body_flag_proximity() -> (
    None
):
    positions = np.full((2, 48, 3), 10.0)
    positions[:, :15] = 0.0
    positions[:, 15:48] = 2.0
    positions[1, 15] = np.asarray([0.01, 0.0, 0.0])

    body_flag, flag_flag = _min_bead_distances(positions)

    assert body_flag[0] == 2.0 * np.sqrt(3.0)
    assert body_flag[1] == 0.01
    assert np.all(np.isfinite(flag_flag))


def test_torque_profile_dt_contact_analysis_keeps_reference_provenance(
    tmp_path: Path,
) -> None:
    root = tmp_path / "uniform"
    condition_id = "as000__ps000__nf03"
    run_dir = root / condition_id
    run_dir.mkdir(parents=True)
    (root / "run_manifest.json").write_text(
        json.dumps(
            {
                "conditions": [
                    {
                        "condition_id": condition_id,
                        "output_dir": str(run_dir),
                        "axis_values": {
                            "attach_seed": 0,
                            "phase_seed": 0,
                            "n_flagella": 3,
                        },
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "run_summary.json").write_text(
        json.dumps(
            {
                "execution": {"status": "completed"},
                "gates": {
                    name: {
                        "status": "available",
                        "any_fail": False,
                        "final_pass": True,
                    }
                    for name in ("finite", "shape_nonbody", "shape_body")
                },
                "all_step_metrics": {},
            }
        ),
        encoding="utf-8",
    )
    positions = np.full((2, 48, 3), 2.0)
    positions[:, :15] = 0.0
    np.savez(
        run_dir / "state_archive.npz",
        t=np.asarray([0.0, 2.0]),
        bead_positions_um=positions,
    )
    output = analyze_contact(
        ContactConfig(
            uniform_reference_run_dir=root,
            diffusive_reference_run_dir=root,
            diagnostic_run_dir=root,
            output_dir=tmp_path / "analysis",
            seed_cases=(condition_id,),
            profiles=("uniform",),
            dt_stars=(1.0e-3,),
        )
    )

    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["condition_count"] == 1
    assert manifest["reused_reference_count"] == 1


def test_issue203_composite_manifest_is_condition_scoped(tmp_path: Path) -> None:
    condition_id = "as000__ps000__nf01"
    assert composite_manifest_path(tmp_path, condition_id) == (
        tmp_path / "as000__ps000__nf01_composite_manifest.json"
    )


def test_issue203_composite_retries_writer_when_opened_codec_is_not_decodable(
    monkeypatch, tmp_path: Path
) -> None:
    attempted: list[str] = []

    class Writer:
        def write(self, image: np.ndarray) -> None:
            del image

        def release(self) -> None:
            return None

    def fake_open(path: Path, **kwargs: object) -> SimpleNamespace:
        del path
        attempted.append(kwargs["codec_candidates"][0])  # type: ignore[index]
        return SimpleNamespace(writer=Writer())

    monkeypatch.setattr(
        "sim_swim.analysis.issue203_composite_replay.open_mp4_writer", fake_open
    )
    monkeypatch.setattr(
        "sim_swim.analysis.issue203_composite_replay.shutil.which", lambda _: None
    )
    monkeypatch.setattr(
        "sim_swim.analysis.issue203_composite_replay._written_frame_count",
        lambda _: 1 if attempted[-1] == "mp4v" else 0,
    )

    codec, codecs = _write_frames_with_fallback(
        tmp_path / "movie.mp4",
        [np.zeros((2, 2, 3), dtype=np.uint8)],
        fps=10.0,
        frame_size=(2, 2),
    )

    assert codec == "mp4v"
    assert codecs == ("avc1", "H264", "mp4v")


def test_issue203_comparison_config_keeps_the_existing_diffusive_reference() -> None:
    config = load_config(
        Path("conf/phase2_analysis/issue203_uniform_paired_comparison.yaml")
    )
    assert config.allowed_n_flagella == (1, 2, 3)
    assert "2010_project_tau_linked_2s" in str(config.diffusive_run_dir)


def test_issue203_motion_feature_config_matches_reference_windows() -> None:
    config = load_motion_feature_config(
        Path("conf/phase2_analysis/issue203_uniform_motion_feature_study.yaml")
    )
    assert config.allowed_n_flagella == (1, 2, 3)
    assert config.durations_s == (0.25, 0.5, 1.0)


def test_replay_resolves_parallel_aggregate_condition_symlink(tmp_path: Path) -> None:
    condition_id = "as000__ps000__nf01"
    archive = tmp_path / "conditions" / condition_id / "state_archive.npz"
    archive.parent.mkdir(parents=True)
    archive.touch()

    resolved = _archive_path(
        tmp_path,
        {"condition_id": condition_id, "output_dir": "/stale/campaign/nf01"},
    )

    assert resolved == archive
