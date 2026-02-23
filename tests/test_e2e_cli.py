from __future__ import annotations

import importlib
import importlib.util
from pathlib import Path

import yaml


def test_script_generates_outputs(tmp_path: Path, monkeypatch) -> None:
    run_context = importlib.import_module("sim_swim.core.run_context")
    monkeypatch.setattr(run_context, "_require_clean_git", lambda: None)

    script_path = (
        Path(__file__).resolve().parents[1] / "scripts" / "01_simulate_swimming.py"
    )
    spec = importlib.util.spec_from_file_location("phase2_script", script_path)
    assert spec is not None
    assert spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    cfg = {
        "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
        "body": {
            "prism": {
                "n_prism": 3,
                "dz_over_b": 0.5,
                "radius_over_b": 0.5,
                "axis": "x",
            },
            "length_total_um": 2.0,
        },
        "flagella": {
            "n_flagella": 3,
            "placement_mode": "uniform",
            "discretization": {"ds_over_b": 0.58},
            "bond_L_over_b": 0.58,
            "length_over_b": 2.32,
            "helix_init": {"radius_over_b": 0.2, "pitch_over_b": 1.0},
        },
        "fluid": {"viscosity_Pa_s": 0.001},
        "motor": {"torque_Nm": 4.0e-18, "reverse_n_flagella": 1},
        "potentials": {
            "spring": {"H_over_T_over_b": 10.0, "s": 0.1},
            "bend": {
                "kb_over_T": 20.0,
                "theta0_deg": {"normal": 142.0, "semicoiled": 90.0, "curly1": 105.0},
            },
            "torsion": {
                "kt_over_T": 10.0,
                "phi0_deg": {"normal": -60.0, "semicoiled": 65.0, "curly1": 120.0},
            },
            "spring_spring_repulsion": {
                "A_ss_over_T": 1.0,
                "a_ss_over_b": 0.2,
                "cutoff_over_b": 0.2,
            },
        },
        "hook": {"enabled": True, "threshold_deg": 90.0, "kb_over_T": 20.0},
        "run_tumble": {
            "run_tau": 20.0,
            "tumble_tau": 8.0,
            "semicoiled_tau": 4.0,
            "curly1_tau": 4.0,
        },
        "time": {"duration_s": 5.0e-5, "dt_s": 2.5e-7},
        "output_sampling": {"out_all_steps_3d": True, "fps_out_2d": 25.0},
        "brownian": {
            "enabled": False,
            "temperature_K": 298.0,
            "method": "cholesky",
            "jitter": 1.0e-20,
        },
        "render": {
            "image_size_px": 96,
            "pixel_size_um": 0.203,
            "flagella_linewidth_px": 2.0,
            "render_flagella": True,
            "save_frames_3d": True,
            "follow_camera_3d": True,
            "view_range_um": 8.0,
            "timestamp_3d": True,
            "timestamp_fmt": "t = {t:.3f} s",
            "label_flagella": True,
            "follow_camera_2d": False,
            "save_frames_2d": True,
        },
        "seed": {"global_seed": 0},
        "output": {"base_dir": str(tmp_path / "outputs")},
    }

    cfg_path = tmp_path / "sim_swim.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

    mod.main(
        config=cfg_path,
        duration_s=5.0e-5,
        fps_out=25.0,
        render_flagella=True,
        overrides=[],
    )

    run_logs = list((tmp_path / "outputs").rglob("run.log"))
    manifests = list((tmp_path / "outputs").rglob("manifest.json"))
    assert run_logs
    assert manifests

    latest = sorted((tmp_path / "outputs").rglob("manifest.json"))[-1].parent
    assert (latest / "render" / "movie_3d.mp4").is_file()
    assert (latest / "render" / "swim3d.mp4").is_file()
    assert (latest / "render" / "swim3d_final.png").is_file()
    assert any((latest / "render" / "frames_3d").glob("frame_*.png"))
    assert (latest / "render2d" / "projection.mp4").is_file()
    assert any((latest / "render2d" / "frames").glob("frame_*.png"))
