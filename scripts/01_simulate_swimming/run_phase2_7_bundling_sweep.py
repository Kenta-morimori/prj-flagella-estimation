"""Phase 2.7 用の多べん毛 torque sweep。"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig, merge_overrides


SUMMARY_FIELDS = [
    "orientation_mode",
    "n_flagella",
    "torque_Nm",
    "duration_s",
    "dt_star",
    "output_dir",
    "phase27_class",
    "shape_pass_nonbody",
    "first_fail_category_nonbody",
    "net_abs_flag_helix_spin_revolutions",
    "bundle_axis_vs_rear_angle_deg",
    "bundle_axis_vs_body_axis_angle_deg",
    "bundle_rearward_projection",
    "bundle_participation_ratio",
    "bundle_independent_flagella_count",
    "bundle_tip_axis_dist_mean_um",
    "flag_tip_pair_dist_mean_um",
    "body_displacement_um",
    "body_speed_um_s",
    "body_axis_cumulative_angle_deg",
    "body_axis_wobble_rms_deg",
    "body_angular_velocity_rms_rad_s",
]

STEP_SUMMARY_METRIC_FIELDS = [
    key
    for key in SUMMARY_FIELDS
    if key
    not in {
        "orientation_mode",
        "n_flagella",
        "torque_Nm",
        "duration_s",
        "dt_star",
        "output_dir",
        "phase27_class",
    }
]


def _parse_csv_floats(text: str) -> list[float]:
    values = [float(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("at least one value is required")
    return values


def _parse_csv_ints(text: str) -> list[int]:
    values = [int(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("at least one value is required")
    return values


def _parse_csv_strs(text: str) -> list[str]:
    values = [part.strip() for part in text.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("at least one value is required")
    return values


def _load_yaml(path: Path) -> dict[str, Any]:
    raw_text = path.read_text(encoding="utf-8") if path.exists() else ""
    return yaml.safe_load(raw_text) or {}


def _step_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _last_step_row(path: Path) -> dict[str, str]:
    rows = _step_rows(path)
    if not rows:
        return {}
    return rows[-1]


def _as_float(row: dict[str, str], key: str, default: float = float("nan")) -> float:
    try:
        return float(row.get(key, default))
    except (TypeError, ValueError):
        return default


def _as_bool(row: dict[str, str], key: str) -> bool:
    return str(row.get(key, "")).strip().lower() in {"true", "1", "yes"}


def classify_phase27_condition(row: dict[str, str], n_flagella: int) -> str:
    """最終stepの診断値から Phase 2.7 用の粗分類を返す。"""

    if not _as_bool(row, "shape_pass_nonbody"):
        return "collapse"

    rear_angle = _as_float(row, "bundle_axis_vs_rear_angle_deg")
    rear_projection = _as_float(row, "bundle_rearward_projection")
    participation = _as_float(row, "bundle_participation_ratio")
    independent = int(_as_float(row, "bundle_independent_flagella_count", n_flagella))

    posterior_like = rear_angle <= 45.0 and rear_projection >= 0.5
    if posterior_like and participation >= 0.95 and independent == 0:
        return "posterior_bundle"
    if posterior_like and participation >= 0.5 and independent < n_flagella:
        return "partial_bundle"
    return "no_bundle"


def _net_abs_helix_spin_revolutions(rows: list[dict[str, str]]) -> float:
    if len(rows) < 2:
        return float("nan")
    first = _as_float(rows[0], "flag_helix_spin_phase_deg")
    last = _as_float(rows[-1], "flag_helix_spin_phase_deg")
    if not np.isfinite(first) or not np.isfinite(last):
        return float("nan")
    return abs(last - first) / 360.0


def _build_config(
    raw_cfg: dict[str, Any],
    *,
    orientation_mode: str,
    n_flagella: int,
    torque_Nm: float,
    duration_s: float,
    dt_star: float,
    overrides: list[str] | None,
) -> SimulationConfig:
    override_dict = merge_overrides({}, overrides)
    override_dict.setdefault("flagella", {})["initial_orientation_mode"] = (
        orientation_mode
    )
    override_dict.setdefault("flagella", {})["n_flagella"] = int(n_flagella)
    override_dict.setdefault("motor", {})["torque_Nm"] = float(torque_Nm)
    override_dict.setdefault("time", {})["duration_s"] = float(duration_s)
    override_dict.setdefault("time", {})["dt_star"] = float(dt_star)
    return SimulationConfig.from_dict(raw_cfg).with_overrides(override_dict)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("conf/sim_swim.yaml"))
    parser.add_argument(
        "--orientation-modes",
        type=_parse_csv_strs,
        default=["side_attach", "posterior_aligned"],
        help="Comma-separated initial orientation modes.",
    )
    parser.add_argument(
        "--n-flagella",
        type=_parse_csv_ints,
        default=[3],
        help="Comma-separated flagella counts.",
    )
    parser.add_argument(
        "--torques",
        type=_parse_csv_floats,
        default=[0.5e-20, 1.0e-20, 1.5e-20, 2.0e-20, 2.5e-20, 3.0e-20],
        help="Comma-separated motor.torque_Nm values.",
    )
    parser.add_argument("--duration-s", type=float, default=0.5)
    parser.add_argument("--dt-star", type=float, default=1.0e-4)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/phase2_7_bundling_sweep"),
    )
    parser.add_argument(
        "overrides",
        nargs="*",
        help="Optional config overrides such as output_sampling.out_all_steps_3d=false.",
    )
    args = parser.parse_args()

    raw_cfg = _load_yaml(args.config)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "phase2_7_bundling_sweep_summary.csv"

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()

        for orientation_mode in args.orientation_modes:
            for n_flagella in args.n_flagella:
                for torque in args.torques:
                    cfg = _build_config(
                        raw_cfg,
                        orientation_mode=orientation_mode,
                        n_flagella=int(n_flagella),
                        torque_Nm=float(torque),
                        duration_s=float(args.duration_s),
                        dt_star=float(args.dt_star),
                        overrides=args.overrides,
                    )
                    run_dir = (
                        out_dir
                        / f"orientation_{orientation_mode}"
                        / f"n_{int(n_flagella)}"
                        / f"torque_{float(torque):.2e}"
                    )
                    run_dir.mkdir(parents=True, exist_ok=True)
                    sim = Simulator(cfg)
                    sim.run(cfg.time.duration_s, step_summary_dir=run_dir)
                    step_rows = _step_rows(run_dir / "step_summary.csv")
                    last = step_rows[-1] if step_rows else {}
                    net_abs_spin = _net_abs_helix_spin_revolutions(step_rows)
                    phase27_class = classify_phase27_condition(last, int(n_flagella))

                    writer.writerow(
                        {
                            "orientation_mode": orientation_mode,
                            "n_flagella": int(n_flagella),
                            "torque_Nm": float(torque),
                            "duration_s": float(cfg.time.duration_s),
                            "dt_star": float(cfg.dt_star),
                            "output_dir": str(run_dir),
                            "phase27_class": phase27_class,
                            "net_abs_flag_helix_spin_revolutions": net_abs_spin,
                            **{
                                key: last.get(key, "")
                                for key in STEP_SUMMARY_METRIC_FIELDS
                                if key != "net_abs_flag_helix_spin_revolutions"
                            },
                        }
                    )

    print(summary_path)


if __name__ == "__main__":
    main()
