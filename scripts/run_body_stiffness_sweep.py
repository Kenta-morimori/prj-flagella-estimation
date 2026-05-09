#!/usr/bin/env python3
"""Sweep body_stiffness_scale and summarize observables.

Usage example:
  uv run python -m scripts.run_body_stiffness_sweep \
    --values 50,25,10,5 --torques 1e-21,1.1e-21,1.2e-21 --duration 0.05
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Any

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _parse_values(text: str) -> list[float]:
    return [float(v.strip()) for v in text.split(",") if v.strip()]


def _parse_torques(text: str) -> list[float]:
    vals = _parse_values(text)
    if any(v < 0.0 for v in vals):
        raise ValueError("Torque list must be non-negative.")
    return vals


def _base_cfg() -> dict[str, Any]:
    # Reuse base from run_motor_scale_sweep minimal config
    from scripts.run_motor_scale_sweep import _base_cfg as _b

    return _b()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run a sweep over body stiffness scale."
    )
    parser.add_argument(
        "--values",
        type=_parse_values,
        required=True,
        help="Comma-separated stiffness values",
    )
    parser.add_argument(
        "--torques",
        type=_parse_torques,
        default=None,
        help="Comma-separated torque_Nm values",
    )
    parser.add_argument(
        "--duration", type=float, default=0.05, help="Simulation duration in seconds"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="outputs/body_stiffness_sweep",
        help="Directory for results",
    )
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    summary_path = out_dir / "body_stiffness_sweep_summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "body_stiffness_scale",
            "torque_Nm",
            "output_dir",
            "finite_pass",
            "shape_pass_nonbody",
            "body_shape_pass",
            "shape_pass",
            "first_fail_category",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()

        base = _base_cfg()
        torques = (
            args.torques if args.torques is not None else [base["motor"]["torque_Nm"]]
        )

        for torque in torques:
            for val in args.values:
                cfg = base.copy()
                cfg["motor"]["torque_Nm"] = float(torque)
                cfg.setdefault("stiffness_scales", {})["body_stiffness_scale"] = float(
                    val
                )
                cfg["time"]["duration_s"] = float(args.duration)
                sim_cfg = SimulationConfig.from_dict(cfg)

                run_dir = out_dir / f"torque_{torque:.2e}" / f"body_stiffness_{val:g}"
                run_dir.mkdir(parents=True, exist_ok=True)
                sim = Simulator(sim_cfg)
                sim.run(sim_cfg.time.duration_s, step_summary_dir=run_dir)

                # read last step summary
                import csv as _csv

                last = {}
                step_path = run_dir / "step_summary.csv"
                if step_path.is_file():
                    with step_path.open("r", encoding="utf-8", newline="") as h:
                        rows = list(_csv.DictReader(h))
                        if rows:
                            last = rows[-1]

                body_shape = {}
                body_path = run_dir / "body_constraint_diagnostics.csv"
                if body_path.is_file():
                    with body_path.open("r", encoding="utf-8", newline="") as h:
                        rows = list(_csv.DictReader(h))
                        if rows:
                            # evaluate pass based on existing helper fields
                            spring = rows[-1].get("body_spring_max_stretch_ratio", "")
                            body_shape = {"body_shape_pass": bool(spring)}

                finite_pass = last.get("finite_pass", "")
                shape_pass_nonbody = last.get("shape_pass_nonbody", "")
                body_shape_pass = body_shape.get("body_shape_pass", "")
                shape_pass = bool(
                    finite_pass and shape_pass_nonbody and body_shape_pass
                )

                writer.writerow(
                    {
                        "body_stiffness_scale": float(val),
                        "torque_Nm": float(torque),
                        "output_dir": str(run_dir),
                        "finite_pass": finite_pass,
                        "shape_pass_nonbody": shape_pass_nonbody,
                        "body_shape_pass": body_shape_pass,
                        "shape_pass": shape_pass,
                        "first_fail_category": last.get("first_fail_category", ""),
                    }
                )

    print(f"Sweep summary saved to {summary_path}")


if __name__ == "__main__":
    main()
