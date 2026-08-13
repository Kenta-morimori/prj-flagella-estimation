"""Issue #190: plan a torque-linked 2010 project dt stability campaign."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.sim.params import SimulationConfig


def build_plan(config_path: Path) -> dict[str, Any]:
    raw = load_yaml(config_path)
    base_path = Path(str(raw["base_config"]))
    duration_tau = float(raw["duration_tau"])
    samples = int(raw.get("comparison_sample_count", 201))
    if duration_tau <= 0 or samples < 2:
        raise ValueError(
            "duration_tau must be positive and comparison_sample_count >= 2"
        )
    torques = [float(v) for v in raw["torques_Nm"]]
    dt_stars = [float(v) for v in raw["dt_stars"]]
    contract = dict(raw.get("campaign_contract", {}) or {})
    formal_reference = float(contract.get("formal_reference_dt_star", 0.0))
    if not torques or not dt_stars or any(v <= 0 for v in torques + dt_stars):
        raise ValueError("torques_Nm and dt_stars must contain positive values")
    base = load_yaml(base_path)
    conditions: list[dict[str, Any]] = []
    for torque in torques:
        for dt_star in dt_stars:
            overrides = {
                "time": {
                    "scale_policy": "reference_torque",
                    "duration": {"value": duration_tau, "unit": "tau"},
                    "integration": {"dt_star": dt_star},
                },
                "motor": {
                    "torque_Nm": torque,
                    "reference_torque_Nm": torque,
                    "torque_for_forces_override_Nm": 0.0,
                },
                "brownian": {"enabled": False},
                "output": {"policy": "compact"},
            }
            cfg = SimulationConfig.from_dict(base).with_overrides(overrides)
            if not (
                abs(cfg.motor_torque_Nm)
                == cfg.reference_torque_Nm
                == cfg.torque_for_forces_Nm
            ):
                raise ValueError(
                    "torque-linked condition is not equal across motor/reference/forces"
                )
            archive_interval_s = cfg.time.duration_s / float(samples - 1)
            if archive_interval_s <= 0.0:
                raise ValueError("derived comparison archive interval must be positive")
            overrides["output"]["archive_interval_s"] = archive_interval_s
            conditions.append(
                {
                    "condition_id": f"T{torque:.0e}_dt{dt_star:.0e}",
                    "torque_Nm_per_flagellum": torque,
                    "dt_star": dt_star,
                    "comparison_role": (
                        "formal_reference"
                        if math.isclose(dt_star, formal_reference, rel_tol=0.0)
                        else (
                            "screen_comparator"
                            if math.isclose(dt_star, min(dt_stars), rel_tol=0.0)
                            else "candidate"
                        )
                    ),
                    "comparison_sample_count": samples,
                    "comparison_archive_interval_s": archive_interval_s,
                    "config_overrides": overrides,
                    "execution_command": (
                        "uv run python -m scripts.01_simulate_swimming "
                        f"config={base_path} "
                        "time.scale_policy=reference_torque "
                        f"time.duration.value={duration_tau:g} time.duration.unit=tau "
                        f"time.integration.dt_star={dt_star:.12g} "
                        f"motor.torque_Nm={torque:.12g} "
                        f"motor.reference_torque_Nm={torque:.12g} "
                        "motor.torque_for_forces_override_Nm=0 "
                        "brownian.enabled=false output.policy=compact"
                        f" output.archive_interval_s={archive_interval_s:.12g}"
                    ),
                    "time": cfg.time_manifest(),
                }
            )
    return {
        "kind": "phase2_2010_torque_linked_dt_stability_plan",
        "contract_version": 1,
        "source_config": str(base_path),
        "duration_tau": duration_tau,
        "conditions": conditions,
        "interpretation_prohibitions": [
            "This plan does not change the supported 2010 project baseline.",
            "This plan does not establish 2015 adoption or long-time swimming stability.",
        ],
    }


def main(argv: list[str] | None = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    payload = json.dumps(build_plan(args.config), ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
