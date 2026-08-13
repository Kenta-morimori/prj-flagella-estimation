from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_phase2_and_dataset_cli_responsibilities_are_separated() -> None:
    simulation_dir = ROOT / "scripts/01_simulate_swimming"
    analysis_dir = ROOT / "scripts/02_phase2_analysis"
    dataset_dir = ROOT / "scripts/03_dataset_building"

    assert (simulation_dir / "run_sweep.py").is_file()
    assert (simulation_dir / "run_multi_run.py").is_file()
    assert not (simulation_dir / "plan_2010_torque_dt_stability.py").exists()
    assert not (simulation_dir / "render_shape_stability_grid_replay.py").exists()

    assert (analysis_dir / "plan_2010_torque_dt_stability.py").is_file()
    assert (analysis_dir / "plot_heatmap.py").is_file()
    assert (analysis_dir / "render_phase2_replay.py").is_file()
    assert (dataset_dir / "build_dataset.py").is_file()
    assert (dataset_dir / "build_clip_dataset.py").is_file()
    assert not (ROOT / "scripts/03_phase3").exists()


def test_moved_phase2_implementations_live_under_src() -> None:
    analysis_dir = ROOT / "src/sim_swim/analysis"

    assert (analysis_dir / "behavior_dataset.py").is_file()
    assert (analysis_dir / "behavior_dataset_replay.py").is_file()
    assert (analysis_dir / "behavior_dataset_distributions.py").is_file()
    assert (analysis_dir / "behavior_dataset_separability.py").is_file()
    assert (analysis_dir / "phase2_replay.py").is_file()
