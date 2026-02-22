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
        "discretization": {"ds_um": 0.8},
        "body": {"length_total_um": 3.0, "diameter_um": 0.8, "bond_L_um": None},
        "flagella": {
            "n_flagella": 2,
            "length_um": 4.0,
            "bond_L_um": None,
            "pitch_um": 2.3,
            "radius_um": 0.2,
            "filament_diameter_um": 0.02,
            "placement_mode": "uniform",
            "helix_step_um": 0.2,
        },
        "scale": {"bead_radius_a_um": 0.15},
        "fluid": {"viscosity_Pa_s": 0.001},
        "motor": {"torque_Nm": 4.0e-18, "reverse_n_flagella": 1},
        "potentials": {
            "spring": {"H": 1.0e-4, "s_um": 0.15},
            "bend": {
                "kb": 2.0e-19,
                "theta0_deg": {"normal": 25.0, "semicoiled": 55.0, "curly1": 75.0},
            },
            "torsion": {
                "kt": 2.0e-19,
                "phi0_deg": {"normal": 15.0, "semicoiled": 95.0, "curly1": 145.0},
            },
            "spring_spring_repulsion": {
                "A_ss": 2.0e-19,
                "a_ss_um": 0.1,
                "cutoff_um": 0.6,
            },
        },
        "hook": {"enabled": True, "kb": 8.0e-20, "threshold_deg": 90.0},
        "run_tumble": {
            "run_tau": 0.2,
            "tumble_tau": 0.08,
            "semicoiled_tau": 0.03,
            "curly1_tau": 0.03,
        },
        "brownian": {
            "enabled": False,
            "temperature_K": 298.0,
            "method": "cholesky",
            "jitter": 1.0e-20,
        },
        "time": {"dt": 0.002, "fps_out": 10.0, "duration_s": 0.02},
        "render": {
            "image_size_px": 96,
            "pixel_size_um": 0.203,
            "render_flagella": False,
            "flagella_linewidth_px": 2.0,
        },
        "seed": {"global_seed": 0},
        "output": {"base_dir": str(tmp_path / "outputs")},
    }

    cfg_path = tmp_path / "sim_swim.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")

    mod.main(
        config=cfg_path,
        duration_s=0.02,
        fps_out=10.0,
        render_flagella=False,
        overrides=[],
    )

    run_logs = list((tmp_path / "outputs").rglob("run.log"))
    manifests = list((tmp_path / "outputs").rglob("manifest.json"))
    assert run_logs
    assert manifests
