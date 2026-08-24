from pathlib import Path

import numpy as np

from sim_swim.analysis.issue203_composite_replay import (
    composite_manifest_path,
    nominal_segment_weights,
)
from sim_swim.analysis.issue203_torque_profile_comparison import load_config
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


def test_issue203_uniform_weight_panel_uses_segments_and_sums_to_one() -> None:
    weights = nominal_segment_weights("uniform", 10)
    assert weights.shape == (10,)
    assert np.allclose(weights, 0.1)
    assert weights.sum() == 1.0


def test_issue203_composite_manifest_is_condition_scoped(tmp_path: Path) -> None:
    condition_id = "as000__ps000__nf01"
    assert composite_manifest_path(tmp_path, condition_id) == (
        tmp_path / "as000__ps000__nf01_composite_manifest.json"
    )


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
