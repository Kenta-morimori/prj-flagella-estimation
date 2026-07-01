"""Phase 2.7 用の多べん毛 torque sweep。"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import yaml

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig, merge_overrides

matplotlib.use("Agg")

AXIS_ALIGNMENT_THRESHOLD_DEG = 15.0
STABLE_WINDOW_FRACTION = 0.8

SUMMARY_FIELDS = [
    "initial_helix_axis_from_rear_deg",
    "n_flagella",
    "torque_Nm",
    "duration_s",
    "final_t_s",
    "dt_star",
    "output_dir",
    "flag_helix_axis_timeseries_plot",
    "phase27_class",
    "phase27_class_hook_len_relaxed",
    "phase27_axis_alignment_threshold_deg",
    "phase27_axis_alignment_stable",
    "phase27_axis_alignment_stable_fraction",
    "shape_pass_nonbody",
    "first_fail_category_nonbody",
    "shape_pass_nonbody_strict",
    "first_fail_category_nonbody_strict",
    "shape_pass_nonbody_hook_len_relaxed",
    "first_fail_category_nonbody_hook_len_relaxed",
    "hook_len_strict_limit",
    "hook_len_relaxed_limit",
    "net_abs_flag_helix_spin_revolutions",
    "hook_len_rel_err_max",
    "hook_angle_err_max_deg",
    "local_attach_first_vs_body_axis_angle_deg",
    "local_attach_first_vs_body_axis_err_deg",
    "flag_helix_axis_vs_rear_angle_deg_min",
    "flag_helix_axis_vs_rear_angle_deg_mean",
    "flag_helix_axis_vs_rear_angle_deg_max",
    "flag_helix_axis_rearward_projection_min",
    "flag_helix_axis_fit_r2_min",
    "flag_helix_axis_degenerate_count",
    "flag_helix_axis_pair_angle_deg_mean",
    "flag_helix_axis_pair_angle_deg_max",
    "flag_helix_axis_mean_deviation_deg_max",
    "flag_helix_axis_alignment_order",
    "flag_flag_helix_bead_dist_min_um",
    "flag_flag_helix_close_pair_count",
    "flag_helix_bundle_radius_mean_um",
    "flag_helix_bundle_radius_max_um",
    "flag_bond_rel_err_max",
    "flag_bend_err_max_deg",
    "flag_torsion_err_max_deg",
    "bundle_axis_vs_rear_angle_deg",
    "bundle_axis_vs_body_axis_angle_deg",
    "bundle_rearward_projection",
    "bundle_participation_ratio",
    "bundle_independent_flagella_count",
    "bundle_tip_axis_dist_mean_um",
    "flag_tip_pair_dist_mean_um",
    "flag_flag_bead_pair_dist_min_um",
    "flag_flag_bead_pair_dist_mean_um",
    "flag_flag_close_pair_count",
    "flag_flag_repulsion_force_mean_N",
    "flag_flag_repulsion_force_max_N",
    "flag_flag_basal_repulsion_force_mean_N",
    "flag_flag_basal_repulsion_force_max_N",
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
        "initial_helix_axis_from_rear_deg",
        "n_flagella",
        "torque_Nm",
        "duration_s",
        "final_t_s",
        "dt_star",
        "output_dir",
        "flag_helix_axis_timeseries_plot",
        "phase27_class",
        "phase27_class_hook_len_relaxed",
        "phase27_axis_alignment_threshold_deg",
        "phase27_axis_alignment_stable",
        "phase27_axis_alignment_stable_fraction",
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


def _stable_window_rows(
    rows: list[dict[str, str]],
    stable_window_fraction: float = STABLE_WINDOW_FRACTION,
) -> list[dict[str, str]]:
    if not rows:
        return []
    times = np.asarray([_as_float(row, "t_s") for row in rows], dtype=float)
    finite_times = times[np.isfinite(times)]
    if finite_times.size == 0:
        return rows
    t0 = float(np.min(finite_times))
    t1 = float(np.max(finite_times))
    start = t0 + (1.0 - float(stable_window_fraction)) * max(t1 - t0, 0.0)
    stable_rows = [
        row
        for row in rows
        if np.isfinite(_as_float(row, "t_s")) and _as_float(row, "t_s") >= start
    ]
    return stable_rows or rows[-1:]


def _axis_alignment_stable_summary(
    rows: list[dict[str, str]],
    *,
    threshold_deg: float = AXIS_ALIGNMENT_THRESHOLD_DEG,
    stable_window_fraction: float = STABLE_WINDOW_FRACTION,
) -> tuple[bool, float]:
    stable_rows = _stable_window_rows(rows, stable_window_fraction)
    if not stable_rows:
        return False, float("nan")

    pass_values: list[bool] = []
    for row in stable_rows:
        pair_max = _as_float(row, "flag_helix_axis_pair_angle_deg_max")
        mean_dev_max = _as_float(row, "flag_helix_axis_mean_deviation_deg_max")
        if np.isfinite(mean_dev_max):
            pass_values.append(mean_dev_max <= float(threshold_deg))
            continue
        if np.isfinite(pair_max):
            pass_values.append(pair_max <= float(threshold_deg))
            continue
        else:
            pass_values.append(False)

    if not pass_values:
        return False, float("nan")
    stable_fraction = float(np.mean(np.asarray(pass_values, dtype=float)))
    return bool(all(pass_values)), stable_fraction


def _has_nonhook_shape_failure(rows: list[dict[str, str]], *, relaxed: bool) -> bool:
    shape_key = (
        "shape_pass_nonbody_hook_len_relaxed" if relaxed else "shape_pass_nonbody"
    )
    fail_key = (
        "first_fail_category_nonbody_hook_len_relaxed"
        if relaxed
        else "first_fail_category_nonbody"
    )
    for row in rows:
        if _as_bool(row, shape_key):
            continue
        fail_category = str(row.get(fail_key, "")).strip()
        if fail_category not in {"", "none", "hook"}:
            return True
    return False


def _has_hook_shape_failure(rows: list[dict[str, str]], *, relaxed: bool) -> bool:
    shape_key = (
        "shape_pass_nonbody_hook_len_relaxed" if relaxed else "shape_pass_nonbody"
    )
    fail_key = (
        "first_fail_category_nonbody_hook_len_relaxed"
        if relaxed
        else "first_fail_category_nonbody"
    )
    return any(
        (not _as_bool(row, shape_key)) and str(row.get(fail_key, "")).strip() == "hook"
        for row in rows
    )


def classify_phase27_condition(
    rows_or_row: list[dict[str, str]] | dict[str, str],
    n_flagella: int,
    *,
    relaxed: bool = False,
    threshold_deg: float = AXIS_ALIGNMENT_THRESHOLD_DEG,
    stable_window_fraction: float = STABLE_WINDOW_FRACTION,
) -> str:
    """Classify Phase 2.7 runs by stable flagellar helix-axis alignment."""

    rows = [rows_or_row] if isinstance(rows_or_row, dict) else list(rows_or_row)
    if not rows:
        return "collapse"
    if int(n_flagella) <= 1:
        return "axis_not_aligned"

    stable_rows = _stable_window_rows(rows, stable_window_fraction)
    if _has_nonhook_shape_failure(stable_rows, relaxed=relaxed):
        return "collapse"

    axis_stable, _stable_fraction = _axis_alignment_stable_summary(
        rows,
        threshold_deg=threshold_deg,
        stable_window_fraction=stable_window_fraction,
    )
    if not axis_stable:
        return "axis_not_aligned"

    if _has_hook_shape_failure(stable_rows, relaxed=relaxed):
        return "hook_wrapped_axis_aligned"
    return "axis_aligned_stable"


def classify_phase27_condition_hook_len_relaxed(
    rows_or_row: list[dict[str, str]] | dict[str, str],
    n_flagella: int,
) -> str:
    return classify_phase27_condition(rows_or_row, n_flagella, relaxed=True)


def _net_abs_helix_spin_revolutions(rows: list[dict[str, str]]) -> float:
    if len(rows) < 2:
        return float("nan")
    first = _as_float(rows[0], "flag_helix_spin_phase_deg")
    last = _as_float(rows[-1], "flag_helix_spin_phase_deg")
    if not np.isfinite(first) or not np.isfinite(last):
        return float("nan")
    return abs(last - first) / 360.0


def _read_axis_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _axis_vector_from_row(row: dict[str, str]) -> np.ndarray:
    return np.asarray(
        [
            _as_float(row, "axis_dir_x"),
            _as_float(row, "axis_dir_y"),
            _as_float(row, "axis_dir_z"),
        ],
        dtype=float,
    )


def _angle_deg_between(a: np.ndarray, b: np.ndarray) -> float:
    norm_a = float(np.linalg.norm(a))
    norm_b = float(np.linalg.norm(b))
    if norm_a <= 1.0e-18 or norm_b <= 1.0e-18:
        return float("nan")
    cos_angle = float(np.clip(np.dot(a, b) / (norm_a * norm_b), -1.0, 1.0))
    return float(np.rad2deg(np.arccos(cos_angle)))


def _axis_deviation_rows(axis_rows: list[dict[str, str]]) -> list[dict[str, float]]:
    by_time: dict[float, list[dict[str, str]]] = {}
    for row in axis_rows:
        t_s = _as_float(row, "t_s")
        if np.isfinite(t_s):
            by_time.setdefault(t_s, []).append(row)

    out: list[dict[str, float]] = []
    for t_s, rows in sorted(by_time.items()):
        vectors: dict[int, np.ndarray] = {}
        for row in rows:
            flag_id = int(_as_float(row, "flag_id", -1))
            vec = _axis_vector_from_row(row)
            norm = float(np.linalg.norm(vec))
            if flag_id >= 0 and norm > 1.0e-18 and np.isfinite(vec).all():
                vectors[flag_id] = vec / norm
        if not vectors:
            continue
        mean_axis = np.mean(np.asarray(list(vectors.values()), dtype=float), axis=0)
        mean_norm = float(np.linalg.norm(mean_axis))
        if mean_norm <= 1.0e-18:
            continue
        mean_axis = mean_axis / mean_norm
        for flag_id, vec in vectors.items():
            out.append(
                {
                    "t_s": float(t_s),
                    "flag_id": float(flag_id),
                    "mean_axis_deviation_deg": _angle_deg_between(vec, mean_axis),
                }
            )
    return out


def plot_flag_helix_axis_timeseries(
    *,
    axis_csv_path: Path,
    output_path: Path,
    threshold_deg: float = AXIS_ALIGNMENT_THRESHOLD_DEG,
    stable_window_fraction: float = STABLE_WINDOW_FRACTION,
) -> Path:
    """Plot per-flagellum rear angle and mean-axis deviation over time."""

    axis_rows = _read_axis_rows(axis_csv_path)
    if not axis_rows:
        raise ValueError(f"No axis diagnostics rows found: {axis_csv_path}")

    import matplotlib.pyplot as plt

    times = np.asarray([_as_float(row, "t_s") for row in axis_rows], dtype=float)
    finite_times = times[np.isfinite(times)]
    if finite_times.size == 0:
        raise ValueError("Axis diagnostics rows do not contain finite t_s values.")
    t0 = float(np.min(finite_times))
    t1 = float(np.max(finite_times))
    stable_start = t0 + (1.0 - float(stable_window_fraction)) * max(t1 - t0, 0.0)

    fig, (ax_rear, ax_dev) = plt.subplots(
        2,
        1,
        figsize=(8.0, 6.0),
        sharex=True,
        constrained_layout=True,
    )

    flag_ids = sorted({int(_as_float(row, "flag_id", -1)) for row in axis_rows})
    flag_ids = [flag_id for flag_id in flag_ids if flag_id >= 0]
    cmap = plt.get_cmap("tab10")
    for idx, flag_id in enumerate(flag_ids):
        color = cmap(idx % 10)
        rows = [
            row for row in axis_rows if int(_as_float(row, "flag_id", -1)) == flag_id
        ]
        rows.sort(key=lambda row: _as_float(row, "t_s"))
        t = [_as_float(row, "t_s") for row in rows]
        rear_angle = [
            _as_float(row, "flag_helix_axis_vs_rear_angle_deg") for row in rows
        ]
        ax_rear.plot(t, rear_angle, label=f"F{flag_id}", color=color, linewidth=1.6)

    deviation_rows = _axis_deviation_rows(axis_rows)
    for idx, flag_id in enumerate(flag_ids):
        color = cmap(idx % 10)
        rows = [row for row in deviation_rows if int(row["flag_id"]) == flag_id]
        rows.sort(key=lambda row: row["t_s"])
        if not rows:
            continue
        ax_dev.plot(
            [row["t_s"] for row in rows],
            [row["mean_axis_deviation_deg"] for row in rows],
            label=f"F{flag_id}",
            color=color,
            linewidth=1.6,
        )

    for ax in (ax_rear, ax_dev):
        ax.axvspan(stable_start, t1, color="0.9", alpha=0.5, zorder=-10)
        ax.grid(True, color="0.85", linewidth=0.8)

    ax_rear.axhline(
        0.0,
        color="0.35",
        linestyle="--",
        linewidth=1.0,
        label="rear axis (0 deg)",
    )
    ax_rear.axhline(
        90.0,
        color="0.35",
        linestyle=":",
        linewidth=1.0,
        label="side axis (90 deg)",
    )
    ax_rear.set_ylim(0.0, 90.0)
    ax_rear.set_ylabel("angle vs rear [deg]")
    ax_rear.set_title("Flagellar helix-axis direction")
    ax_rear.legend(loc="best", fontsize=8)

    ax_dev.axhline(
        float(threshold_deg),
        color="crimson",
        linestyle="--",
        linewidth=1.2,
        label=f"alignment threshold ({threshold_deg:g} deg)",
    )
    ax_dev.set_ylabel("deviation from mean axis [deg]")
    ax_dev.set_xlabel("time [s]")
    ax_dev.legend(loc="best", fontsize=8)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)
    return output_path


def _build_config(
    raw_cfg: dict[str, Any],
    *,
    initial_helix_axis_from_rear_deg: float,
    n_flagella: int,
    torque_Nm: float,
    duration_s: float,
    dt_star: float,
    overrides: list[str] | None,
) -> SimulationConfig:
    override_dict = merge_overrides({}, overrides)
    flagella_overrides = override_dict.setdefault("flagella", {})
    flagella_overrides["initial_helix_axis_from_rear_deg"] = float(
        initial_helix_axis_from_rear_deg
    )
    flagella_overrides["n_flagella"] = int(n_flagella)
    override_dict.setdefault("motor", {})["torque_Nm"] = float(torque_Nm)
    override_dict.setdefault("time", {})["duration_s"] = float(duration_s)
    override_dict.setdefault("time", {})["dt_star"] = float(dt_star)
    return SimulationConfig.from_dict(raw_cfg).with_overrides(override_dict)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("conf/sim_swim.yaml"))
    parser.add_argument(
        "--helix-axis-angles-deg",
        type=_parse_csv_floats,
        default=[0.0],
        help=(
            "Comma-separated flagella.initial_helix_axis_from_rear_deg values. "
            "0 aligns the helix center axis estimated from bead 2 onward to the "
            "body rear direction."
        ),
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
        "--no-stop-on-shape-fail",
        action="store_true",
        help=(
            "Deprecated no-op. Phase 2.7 axis-alignment sweeps run to full "
            "duration by default so temporal stability can be evaluated."
        ),
    )
    parser.add_argument(
        "--stop-on-shape-fail",
        action="store_true",
        help="Stop each condition when shape_pass_nonbody first fails.",
    )
    parser.add_argument(
        "overrides",
        nargs="*",
        help="Optional config overrides such as output_sampling.out_all_steps_3d=false.",
    )
    args = parser.parse_args(argv)

    raw_cfg = _load_yaml(args.config)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "summary.csv"

    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()

        conditions = [
            (helix_axis_angle_deg, n_flagella, torque)
            for helix_axis_angle_deg in args.helix_axis_angles_deg
            for n_flagella in args.n_flagella
            for torque in args.torques
        ]
        total = len(conditions)
        for index, (helix_axis_angle_deg, n_flagella, torque) in enumerate(
            conditions, start=1
        ):
            print(
                "[{}/{}] bundling_alignment angle={} n_flagella={} "
                "torque={:.3e}".format(
                    index,
                    total,
                    float(helix_axis_angle_deg),
                    int(n_flagella),
                    float(torque),
                ),
                flush=True,
            )
            cfg = _build_config(
                raw_cfg,
                initial_helix_axis_from_rear_deg=float(helix_axis_angle_deg),
                n_flagella=int(n_flagella),
                torque_Nm=float(torque),
                duration_s=float(args.duration_s),
                dt_star=float(args.dt_star),
                overrides=args.overrides,
            )
            angle_label = f"{float(helix_axis_angle_deg):g}deg"
            run_dir = (
                out_dir
                / f"helix_axis_angle_{angle_label}"
                / f"n_{int(n_flagella)}"
                / f"torque_{float(torque):.2e}"
            )
            run_dir.mkdir(parents=True, exist_ok=True)
            sim = Simulator(cfg)
            sim.run(
                cfg.time.duration_s,
                step_summary_dir=run_dir,
                stop_on_shape_fail=bool(args.stop_on_shape_fail),
            )
            step_rows = _step_rows(run_dir / "step_summary.csv")
            last = step_rows[-1] if step_rows else {}
            net_abs_spin = _net_abs_helix_spin_revolutions(step_rows)
            phase27_class = classify_phase27_condition(
                step_rows,
                int(n_flagella),
            )
            phase27_class_hook_len_relaxed = (
                classify_phase27_condition_hook_len_relaxed(
                    step_rows,
                    int(n_flagella),
                )
            )
            axis_alignment_stable, axis_alignment_stable_fraction = (
                _axis_alignment_stable_summary(step_rows)
            )
            axis_plot_path = run_dir / "flag_helix_axis_angles_timeseries.png"
            try:
                plot_flag_helix_axis_timeseries(
                    axis_csv_path=run_dir / "flag_helix_axis_diagnostics.csv",
                    output_path=axis_plot_path,
                )
                axis_plot_value = str(axis_plot_path)
            except ValueError:
                axis_plot_value = ""

            writer.writerow(
                {
                    "initial_helix_axis_from_rear_deg": float(helix_axis_angle_deg),
                    "n_flagella": int(n_flagella),
                    "torque_Nm": float(torque),
                    "duration_s": float(cfg.time.duration_s),
                    "final_t_s": last.get("t_s", ""),
                    "dt_star": float(cfg.dt_star),
                    "output_dir": str(run_dir),
                    "flag_helix_axis_timeseries_plot": axis_plot_value,
                    "phase27_class": phase27_class,
                    "phase27_class_hook_len_relaxed": (phase27_class_hook_len_relaxed),
                    "phase27_axis_alignment_threshold_deg": (
                        AXIS_ALIGNMENT_THRESHOLD_DEG
                    ),
                    "phase27_axis_alignment_stable": axis_alignment_stable,
                    "phase27_axis_alignment_stable_fraction": (
                        axis_alignment_stable_fraction
                    ),
                    "net_abs_flag_helix_spin_revolutions": net_abs_spin,
                    **{
                        key: last.get(key, "")
                        for key in STEP_SUMMARY_METRIC_FIELDS
                        if key != "net_abs_flag_helix_spin_revolutions"
                    },
                }
            )
            handle.flush()

    print(summary_path)


if __name__ == "__main__":
    main()
