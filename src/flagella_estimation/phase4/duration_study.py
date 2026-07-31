"""Window-duration and seed-effect analysis for Phase 4 pseudo clips."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
import json
import math
import platform
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from flagella_estimation.phase3.render import render_clip_array, select_frames
from flagella_estimation.phase3.windows import generate_windows
from flagella_estimation.phase4.baseline import FEATURE_NAMES, extract_clip_features
from sim_swim.analysis.flagella_count_behavior import load_state_archive


THREE_D_FEATURE_NAMES = (
    "cell_displacement_um",
    "cell_path_length_um",
    "cell_mean_speed_um_s",
    "cell_speed_std_um_s",
    "cell_speed_cv",
    "cell_straightness",
    "cell_angular_velocity_mean_rad_s",
    "cell_angular_velocity_std_rad_s",
    "cell_angular_velocity_rms_rad_s",
    "body_axis_angle_change_deg",
    "flagella_axis_alignment_mean",
    "cell_flagella_axis_angle_mean_deg",
    "cell_flagella_axis_angle_std_deg",
    "hook_len_rel_err_max",
)

PLOT_FEATURES = {
    "3d": (
        "cell_mean_speed_um_s",
        "cell_angular_velocity_rms_rad_s",
        "flagella_axis_alignment_mean",
    ),
    "2d": (
        "temporal_diff_mean",
        "centroid_step_mean",
        "radial_spread_mean",
    ),
}

COMMON_FIELDS = (
    "sample_id",
    "run_id",
    "group_key",
    "dataset_id",
    "dataset_version",
    "dataset_revision",
    "n_flagella",
    "attach_seed",
    "phase_seed",
    "requested_duration_s",
    "actual_duration_s",
    "window_index",
    "frame_count",
    "frame_rate_hz",
    "t_start_s",
    "t_end_s",
    "run_quality_class",
    "run_shape_pass",
    "run_relaxed_pass",
    "window_quality_class",
    "window_shape_pass",
    "window_relaxed_pass",
    "window_first_fail_t_s",
    "window_first_fail_category",
)


@dataclass(frozen=True)
class DurationStudyConfig:
    dataset_dir: Path
    output_dir: Path
    config_path: Path | None = None
    cli_overrides: tuple[str, ...] = ()
    durations_s: tuple[float, ...] = (0.25, 0.5, 1.0)
    frame_rate_hz: float = 25.0
    crop_size_px: int = 96
    pixel_size_um: float = 0.1
    allowed_n_flagella: tuple[int, ...] = (1, 2, 3)
    overwrite: bool = False


def default_output_dir() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / "phase4_duration_seed_study"
    )


def _coerce_override(raw: str) -> tuple[list[str], Any]:
    if "=" not in raw:
        raise ValueError(f"Invalid override; expected KEY=VALUE: {raw}")
    key, value = raw.split("=", 1)
    parts = key.strip().split(".")
    if not all(parts):
        raise ValueError(f"Invalid override path: {raw}")
    return parts, yaml.safe_load(value)


def load_duration_study_config(
    path: Path | None, overrides: list[str] | None = None
) -> DurationStudyConfig:
    data = (
        yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        if path is not None
        else {}
    )
    for raw in overrides or []:
        parts, value = _coerce_override(raw)
        node = data
        for part in parts[:-1]:
            node = node.setdefault(part, {})
            if not isinstance(node, dict):
                raise ValueError(f"Override path conflicts with scalar: {raw}")
        node[parts[-1]] = value
    clip = dict(data.get("clip", {}) or {})
    study = dict(data.get("study", {}) or {})
    output_dir = data.get("output_dir")
    return DurationStudyConfig(
        dataset_dir=Path(str(data.get("dataset_dir", ""))),
        output_dir=Path(str(output_dir)) if output_dir else default_output_dir(),
        config_path=path,
        cli_overrides=tuple(overrides or ()),
        durations_s=tuple(
            float(value) for value in study.get("durations_s", [0.25, 0.5, 1.0])
        ),
        frame_rate_hz=float(clip.get("frame_rate_hz", 25.0)),
        crop_size_px=int(clip.get("crop_size_px", 96)),
        pixel_size_um=float(clip.get("pixel_size_um", 0.1)),
        allowed_n_flagella=tuple(
            int(value) for value in study.get("allowed_n_flagella", [1, 2, 3])
        ),
        overwrite=bool(data.get("overwrite", False)),
    )


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(
    path: Path, rows: list[dict[str, Any]], fieldnames: tuple[str, ...] | list[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _to_float(value: Any, default: float = float("nan")) -> float:
    try:
        if value in (None, ""):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _to_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _git_info() -> dict[str, Any]:
    def run(command: list[str]) -> str:
        return subprocess.check_output(
            command, stderr=subprocess.STDOUT, text=True
        ).strip()

    try:
        return {
            "commit": run(["git", "rev-parse", "HEAD"]),
            "commit_short": run(["git", "rev-parse", "--short", "HEAD"]),
            "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
            "is_clean": run(["git", "status", "--porcelain"]) == "",
        }
    except Exception:
        return {
            "commit": "unknown",
            "commit_short": "unknown",
            "branch": "unknown",
            "is_clean": False,
        }


def _environment_info() -> dict[str, Any]:
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
    }


def _mean(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return float(np.mean(finite)) if finite else float("nan")


def _std(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return float("nan")
    return float(np.std(finite, ddof=1)) if len(finite) > 1 else 0.0


def _body_axis(quaternion: tuple[float, float, float, float]) -> np.ndarray:
    x, y, z, w = np.asarray(quaternion, dtype=float)
    norm = math.sqrt(x * x + y * y + z * z + w * w)
    if not math.isfinite(norm) or norm <= 0.0:
        return np.full(3, np.nan)
    x, y, z, w = (x / norm, y / norm, z / norm, w / norm)
    return np.asarray(
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y + z * w),
            2.0 * (x * z - y * w),
        ],
        dtype=float,
    )


def _angle_deg(first: np.ndarray, last: np.ndarray) -> float:
    if not np.isfinite(first).all() or not np.isfinite(last).all():
        return float("nan")
    denominator = float(np.linalg.norm(first) * np.linalg.norm(last))
    if denominator <= 0.0:
        return float("nan")
    cosine = float(np.clip(np.dot(first, last) / denominator, -1.0, 1.0))
    return float(np.rad2deg(np.arccos(cosine)))


def summarize_states_3d(states: list[Any]) -> dict[str, float]:
    if len(states) < 2:
        return {name: float("nan") for name in THREE_D_FEATURE_NAMES}
    positions = np.asarray([state.position_um for state in states], dtype=float)
    velocities = np.asarray([state.velocity_um_s for state in states], dtype=float)
    omega = np.linalg.norm(
        np.asarray([state.omega_rad_s for state in states], dtype=float), axis=1
    )
    segments = np.linalg.norm(np.diff(positions, axis=0), axis=1)
    displacement = float(np.linalg.norm(positions[-1] - positions[0]))
    path_length = float(np.sum(segments))
    speed = np.linalg.norm(velocities, axis=1)
    speed_mean = float(np.mean(speed))
    speed_std = float(np.std(speed, ddof=1)) if len(speed) > 1 else 0.0
    result = {
        "cell_displacement_um": displacement,
        "cell_path_length_um": path_length,
        "cell_mean_speed_um_s": speed_mean,
        "cell_speed_std_um_s": speed_std,
        "cell_speed_cv": speed_std / speed_mean
        if abs(speed_mean) > 1.0e-30
        else float("nan"),
        "cell_straightness": displacement / path_length
        if path_length > 0.0
        else float("nan"),
        "cell_angular_velocity_mean_rad_s": float(np.mean(omega)),
        "cell_angular_velocity_std_rad_s": float(np.std(omega, ddof=1))
        if len(omega) > 1
        else 0.0,
        "cell_angular_velocity_rms_rad_s": float(np.sqrt(np.mean(omega**2))),
        "body_axis_angle_change_deg": _angle_deg(
            _body_axis(states[0].quaternion), _body_axis(states[-1].quaternion)
        ),
    }
    return result


def _step_rows_for_window(
    rows: list[dict[str, str]], *, start_s: float, end_s: float
) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if start_s - 1.0e-12 <= _to_float(row.get("t_s")) < end_s - 1.0e-12
    ]


def _window_quality(rows: list[dict[str, str]]) -> dict[str, Any]:
    if not rows:
        return {
            "window_quality_class": "missing",
            "window_shape_pass": False,
            "window_relaxed_pass": False,
            "window_first_fail_t_s": float("nan"),
            "window_first_fail_category": "missing",
        }
    finite_pass = all(_to_bool(row.get("finite_pass")) for row in rows)
    shape_pass = finite_pass and all(
        _to_bool(row.get("shape_pass_nonbody_strict")) for row in rows
    )
    relaxed_pass = finite_pass and all(
        _to_bool(row.get("shape_pass_nonbody_hook_len_relaxed")) for row in rows
    )
    failed = next(
        (
            row
            for row in rows
            if not _to_bool(row.get("finite_pass"))
            or not _to_bool(row.get("shape_pass_nonbody_strict"))
        ),
        None,
    )
    if not finite_pass:
        quality_class = "invalid_numeric"
    elif shape_pass:
        quality_class = "strict_pass"
    elif relaxed_pass:
        quality_class = "relaxed_pass"
    else:
        quality_class = "fail"
    return {
        "window_quality_class": quality_class,
        "window_shape_pass": shape_pass,
        "window_relaxed_pass": relaxed_pass,
        "window_first_fail_t_s": _to_float(failed.get("t_s"))
        if failed
        else float("nan"),
        "window_first_fail_category": str(
            failed.get("first_fail_category_nonbody_strict")
            or failed.get("first_fail_category_nonbody")
            or "unknown"
        )
        if failed
        else "none",
    }


def _step_features(rows: list[dict[str, str]]) -> dict[str, float]:
    relation = [
        _to_float(row.get("bundle_axis_vs_body_axis_angle_deg")) for row in rows
    ]
    hook_error = [_to_float(row.get("hook_len_rel_err_max")) for row in rows]
    finite_hook = [value for value in hook_error if math.isfinite(value)]
    return {
        "flagella_axis_alignment_mean": _mean(
            [_to_float(row.get("flag_helix_axis_alignment_order")) for row in rows]
        ),
        "cell_flagella_axis_angle_mean_deg": _mean(relation),
        "cell_flagella_axis_angle_std_deg": _std(relation),
        "hook_len_rel_err_max": max(finite_hook) if finite_hook else float("nan"),
    }


def _common_row(
    sample: dict[str, str],
    *,
    requested_duration_s: float,
    actual_duration_s: float,
    window_index: int,
    frame_count: int,
    frame_rate_hz: float,
    t_start_s: float,
    t_end_s: float,
    quality: dict[str, Any],
    dataset_version: str,
    dataset_revision: str,
) -> dict[str, Any]:
    sample_id = str(sample["sample_id"])
    return {
        "sample_id": sample_id,
        "run_id": sample_id,
        "group_key": f"phase2:v1:{sample_id}",
        "dataset_id": sample.get("dataset_id", ""),
        "dataset_version": dataset_version,
        "dataset_revision": dataset_revision,
        "n_flagella": int(_to_float(sample.get("n_flagella"))),
        "attach_seed": int(_to_float(sample.get("attach_seed"))),
        "phase_seed": int(_to_float(sample.get("phase_seed"))),
        "requested_duration_s": requested_duration_s,
        "actual_duration_s": actual_duration_s,
        "window_index": window_index,
        "frame_count": frame_count,
        "frame_rate_hz": frame_rate_hz,
        "t_start_s": t_start_s,
        "t_end_s": t_end_s,
        "run_quality_class": sample.get("quality_class", ""),
        "run_shape_pass": _to_bool(sample.get("shape_pass")),
        "run_relaxed_pass": _to_bool(sample.get("relaxed_pass")),
        **quality,
    }


def _run_feature_mean_rows(
    rows: list[dict[str, Any]], feature_names: tuple[str, ...], domain: str
) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[tuple[float, bool]]] = {}
    for row in rows:
        for feature in feature_names:
            value = _to_float(row.get(feature))
            if not math.isfinite(value):
                continue
            key = (
                row["requested_duration_s"],
                row["actual_duration_s"],
                row["n_flagella"],
                row["attach_seed"],
                row["phase_seed"],
                feature,
            )
            grouped.setdefault(key, []).append((value, bool(row["window_shape_pass"])))
    output: list[dict[str, Any]] = []
    for key, values in sorted(grouped.items()):
        requested, actual, n_flagella, attach_seed, phase_seed, feature = key
        sample = next(
            row
            for row in rows
            if row["requested_duration_s"] == requested
            and row["actual_duration_s"] == actual
            and row["n_flagella"] == n_flagella
            and row["attach_seed"] == attach_seed
            and row["phase_seed"] == phase_seed
        )
        output.append(
            {
                "domain": domain,
                "feature": feature,
                "requested_duration_s": requested,
                "actual_duration_s": actual,
                "n_flagella": n_flagella,
                "attach_seed": attach_seed,
                "phase_seed": phase_seed,
                "run_id": sample["run_id"],
                "run_shape_pass": sample["run_shape_pass"],
                "all_windows_shape_pass": all(passed for _, passed in values),
                "window_count": len(values),
                "run_mean": float(np.mean([value for value, _ in values])),
            }
        )
    return output


def _seed_effect_rows(
    run_means: list[dict[str, Any]], *, qc_scope: str
) -> list[dict[str, Any]]:
    if qc_scope not in {"all_labeled", "strict_run_and_windows"}:
        raise ValueError(f"Unsupported seed-effect QC scope: {qc_scope}")
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in run_means:
        key = (
            row["domain"],
            row["requested_duration_s"],
            row["actual_duration_s"],
            row["n_flagella"],
            row["feature"],
        )
        groups.setdefault(key, []).append(row)

    output: list[dict[str, Any]] = []
    for (
        domain,
        requested,
        actual,
        n_flagella,
        feature,
    ), all_rows in sorted(groups.items()):
        expected_attach_levels = sorted({int(row["attach_seed"]) for row in all_rows})
        expected_phase_levels = sorted({int(row["phase_seed"]) for row in all_rows})
        selected_rows = (
            [
                row
                for row in all_rows
                if bool(row["run_shape_pass"]) and bool(row["all_windows_shape_pass"])
            ]
            if qc_scope == "strict_run_and_windows"
            else all_rows
        )
        values = [
            (
                int(row["attach_seed"]),
                int(row["phase_seed"]),
                float(row["run_mean"]),
            )
            for row in selected_rows
        ]
        if not values:
            continue
        observed = np.asarray([value for _, _, value in values], dtype=float)
        grand = float(np.mean(observed))
        total_ss = float(np.sum((observed - grand) ** 2))
        observed_cells = {(attach, phase) for attach, phase, _ in values}
        expected_cells = {
            (attach, phase)
            for attach in expected_attach_levels
            for phase in expected_phase_levels
        }
        factorial_complete = (
            len(values) == len(expected_cells) and observed_cells == expected_cells
        )
        attach_eta: float | None = None
        phase_eta: float | None = None
        residual_eta: float | None = None
        if factorial_complete:
            attach_ss = sum(
                sum(1 for a, _, _ in values if a == attach)
                * (np.mean([value for a, _, value in values if a == attach]) - grand)
                ** 2
                for attach in expected_attach_levels
            )
            phase_ss = sum(
                sum(1 for _, p, _ in values if p == phase)
                * (np.mean([value for _, p, value in values if p == phase]) - grand)
                ** 2
                for phase in expected_phase_levels
            )
            residual_ss = max(total_ss - float(attach_ss) - float(phase_ss), 0.0)
            if total_ss > 0.0:
                attach_eta = float(attach_ss / total_ss)
                phase_eta = float(phase_ss / total_ss)
                residual_eta = residual_ss / total_ss
            else:
                attach_eta = phase_eta = residual_eta = 0.0
        output.append(
            {
                "domain": domain,
                "qc_scope": qc_scope,
                "feature": feature,
                "requested_duration_s": requested,
                "actual_duration_s": actual,
                "n_flagella": n_flagella,
                "run_count": len(values),
                "attach_level_count": len(expected_attach_levels),
                "phase_level_count": len(expected_phase_levels),
                "factorial_complete": factorial_complete,
                "mean": grand,
                "between_run_std": float(np.std(observed, ddof=1))
                if len(observed) > 1
                else 0.0,
                "attach_eta_squared": attach_eta,
                "phase_eta_squared": phase_eta,
                "residual_eta_squared": residual_eta,
            }
        )
    return output


def _plot_seed_heatmaps(
    run_means: list[dict[str, Any]],
    *,
    output_dir: Path,
) -> list[str]:
    paths: list[str] = []
    for domain, features in PLOT_FEATURES.items():
        for feature in features:
            matching_feature = [
                row
                for row in run_means
                if row["domain"] == domain and row["feature"] == feature
            ]
            durations = sorted(
                {float(row["requested_duration_s"]) for row in matching_feature}
            )
            classes = sorted({int(row["n_flagella"]) for row in matching_feature})
            figure, axes = plt.subplots(
                len(durations),
                len(classes),
                figsize=(4.2 * len(classes), 3.6 * len(durations)),
                squeeze=False,
            )
            for duration_index, duration in enumerate(durations):
                duration_rows = [
                    row
                    for row in matching_feature
                    if float(row["requested_duration_s"]) == duration
                ]
                observed = np.asarray(
                    [float(row["run_mean"]) for row in duration_rows], dtype=float
                )
                vmin = float(np.min(observed))
                vmax = float(np.max(observed))
                if math.isclose(vmin, vmax):
                    vmax = vmin + 1.0
                image = None
                for class_index, n_flagella in enumerate(classes):
                    axis = axes[duration_index, class_index]
                    class_rows = [
                        row
                        for row in duration_rows
                        if int(row["n_flagella"]) == n_flagella
                    ]
                    if not class_rows:
                        axis.set_visible(False)
                        continue
                    attach_levels = sorted(
                        {int(row["attach_seed"]) for row in class_rows}
                    )
                    phase_levels = sorted(
                        {int(row["phase_seed"]) for row in class_rows}
                    )
                    matrix = np.full(
                        (len(attach_levels), len(phase_levels)), np.nan, dtype=float
                    )
                    strict = np.zeros_like(matrix, dtype=bool)
                    for row in class_rows:
                        attach_index = attach_levels.index(int(row["attach_seed"]))
                        phase_index = phase_levels.index(int(row["phase_seed"]))
                        matrix[attach_index, phase_index] = float(row["run_mean"])
                        strict[attach_index, phase_index] = bool(
                            row["run_shape_pass"]
                        ) and bool(row["all_windows_shape_pass"])
                    image = axis.imshow(matrix, vmin=vmin, vmax=vmax, cmap="viridis")
                    for row_index in range(matrix.shape[0]):
                        for column_index in range(matrix.shape[1]):
                            if not math.isfinite(matrix[row_index, column_index]):
                                continue
                            suffix = "" if strict[row_index, column_index] else " !"
                            axis.text(
                                column_index,
                                row_index,
                                f"{matrix[row_index, column_index]:.3g}{suffix}",
                                ha="center",
                                va="center",
                                color="white",
                                fontsize=8,
                            )
                    axis.set_xticks(range(len(phase_levels)), phase_levels)
                    axis.set_yticks(range(len(attach_levels)), attach_levels)
                    axis.set_xlabel("phase seed")
                    axis.set_ylabel("attach seed")
                    axis.set_title(f"requested {duration:g} s, n={n_flagella}")
                if image is not None:
                    figure.colorbar(
                        image,
                        ax=list(axes[duration_index]),
                        shrink=0.8,
                        label=f"{duration:g} s scale",
                    )
            figure.suptitle(f"{feature} by duration and seed (! = QC fail)")
            path = output_dir / f"{domain}_{feature}_seed_heatmap.png"
            figure.savefig(path, dpi=150, bbox_inches="tight")
            plt.close(figure)
            paths.append(str(path))
    return paths


def _plot_time_features(
    rows: list[dict[str, Any]],
    *,
    domain: str,
    features: tuple[str, ...],
    output_dir: Path,
) -> list[str]:
    paths: list[str] = []
    durations = sorted({float(row["requested_duration_s"]) for row in rows})
    colors = {1: "#0072B2", 2: "#D55E00", 3: "#009E73"}
    for feature in features:
        figure, axes = plt.subplots(
            1,
            len(durations),
            figsize=(5.0 * len(durations), 4.0),
            squeeze=False,
            sharey=True,
        )
        for axis, duration in zip(axes[0], durations, strict=True):
            duration_rows = [
                row for row in rows if float(row["requested_duration_s"]) == duration
            ]
            for n_flagella in sorted({int(row["n_flagella"]) for row in duration_rows}):
                class_rows = [
                    row for row in duration_rows if int(row["n_flagella"]) == n_flagella
                ]
                by_run: dict[str, list[dict[str, Any]]] = {}
                for row in class_rows:
                    by_run.setdefault(str(row["run_id"]), []).append(row)
                for run_rows in by_run.values():
                    run_rows.sort(key=lambda row: float(row["t_start_s"]))
                    axis.plot(
                        [float(row["t_start_s"]) for row in run_rows],
                        [_to_float(row.get(feature)) for row in run_rows],
                        color=colors.get(n_flagella, "#666666"),
                        alpha=0.22,
                        linewidth=0.8,
                    )
                time_values = sorted({float(row["t_start_s"]) for row in class_rows})
                means = [
                    _mean(
                        [
                            _to_float(row.get(feature))
                            for row in class_rows
                            if math.isclose(
                                float(row["t_start_s"]), time_value, abs_tol=1.0e-9
                            )
                        ]
                    )
                    for time_value in time_values
                ]
                axis.plot(
                    time_values,
                    means,
                    color=colors.get(n_flagella, "#666666"),
                    linewidth=2.2,
                    marker="o",
                    markersize=3,
                    label=f"n={n_flagella}",
                )
            failed_rows = [
                row
                for row in duration_rows
                if not bool(row["window_shape_pass"])
                and math.isfinite(_to_float(row.get(feature)))
            ]
            if failed_rows:
                axis.scatter(
                    [float(row["t_start_s"]) for row in failed_rows],
                    [_to_float(row.get(feature)) for row in failed_rows],
                    marker="x",
                    color="#111111",
                    s=18,
                    linewidths=0.8,
                    label="QC fail" if axis is axes[0, 0] else None,
                )
            axis.set_title(f"requested {duration:g} s")
            axis.set_xlabel("window start (s)")
            axis.set_ylabel(feature)
            axis.grid(alpha=0.25)
        axes[0, 0].legend(frameon=False)
        figure.tight_layout()
        path = output_dir / f"{domain}_{feature}_by_time.png"
        figure.savefig(path, dpi=150)
        plt.close(figure)
        paths.append(str(path))
    return paths


def analyze_duration_seed_study(cfg: DurationStudyConfig) -> Path:
    summary_path = cfg.dataset_dir / "summary.csv"
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    if cfg.output_dir.exists():
        if not cfg.overwrite:
            raise FileExistsError(cfg.output_dir)
        shutil.rmtree(cfg.output_dir)
    cfg.output_dir.mkdir(parents=True)
    plots_dir = cfg.output_dir / "plots"
    plots_dir.mkdir()

    dataset_manifest_path = cfg.dataset_dir / "dataset_manifest.json"
    dataset_manifest = (
        json.loads(dataset_manifest_path.read_text(encoding="utf-8"))
        if dataset_manifest_path.is_file()
        else {}
    )
    effective = dataset_manifest.get("effective_campaign_config", {})
    metadata = effective.get("metadata", {})
    dataset_cfg = effective.get("dataset", {})
    dataset_version = str(
        dataset_cfg.get("dataset_version", metadata.get("dataset_version", "v1"))
    )
    dataset_revision = str(
        dataset_cfg.get("dataset_revision", metadata.get("dataset_revision", ""))
    )

    rows_3d: list[dict[str, Any]] = []
    rows_2d: list[dict[str, Any]] = []
    eligible_samples = [
        row
        for row in _read_csv(summary_path)
        if int(_to_float(row.get("n_flagella"), default=-1)) in cfg.allowed_n_flagella
    ]
    if not eligible_samples:
        raise ValueError("No eligible samples found")

    for sample in eligible_samples:
        raw_dir = Path(str(sample.get("raw_dir", "")))
        archive_path = raw_dir / "state_archive.npz"
        step_path = raw_dir / "step_summary.csv"
        if not archive_path.is_file():
            raise FileNotFoundError(archive_path)
        if not step_path.is_file():
            raise FileNotFoundError(step_path)
        states_all = load_state_archive(archive_path)
        states_2d = select_frames(states_all, cfg.frame_rate_hz)
        step_rows = _read_csv(step_path)
        for requested_duration in cfg.durations_s:
            windows = generate_windows(
                source_frame_count=len(states_2d),
                frame_rate_hz=cfg.frame_rate_hz,
                duration_s=requested_duration,
                policy="non_overlap",
            )
            for window_index, window in enumerate(windows):
                frame_states = states_2d[window.start : window.end]
                if not frame_states:
                    continue
                actual_duration = window.frame_count / cfg.frame_rate_hz
                t_start = float(frame_states[0].t)
                t_end = t_start + actual_duration
                raw_states = [
                    state
                    for state in states_all
                    if t_start - 1.0e-12 <= float(state.t) < t_end - 1.0e-12
                ]
                step_window = _step_rows_for_window(
                    step_rows, start_s=t_start, end_s=t_end
                )
                quality = _window_quality(step_window)
                common = _common_row(
                    sample,
                    requested_duration_s=requested_duration,
                    actual_duration_s=actual_duration,
                    window_index=window_index,
                    frame_count=window.frame_count,
                    frame_rate_hz=cfg.frame_rate_hz,
                    t_start_s=t_start,
                    t_end_s=t_end,
                    quality=quality,
                    dataset_version=dataset_version,
                    dataset_revision=dataset_revision,
                )
                rows_3d.append(
                    {
                        **common,
                        **summarize_states_3d(raw_states),
                        **_step_features(step_window),
                    }
                )
                clip, _ = render_clip_array(
                    frame_states,
                    image_size_px=cfg.crop_size_px,
                    pixel_size_um=cfg.pixel_size_um,
                )
                rows_2d.append(
                    {
                        **common,
                        **dict(
                            zip(
                                FEATURE_NAMES,
                                extract_clip_features(clip).tolist(),
                                strict=True,
                            )
                        ),
                    }
                )

    features_3d_path = cfg.output_dir / "window_features_3d.csv"
    features_2d_path = cfg.output_dir / "window_features_2d.csv"
    run_means_path = cfg.output_dir / "run_feature_means.csv"
    seed_effects_path = cfg.output_dir / "seed_effects.csv"
    _write_csv(
        features_3d_path,
        rows_3d,
        [*COMMON_FIELDS, *THREE_D_FEATURE_NAMES],
    )
    _write_csv(
        features_2d_path,
        rows_2d,
        [*COMMON_FIELDS, *FEATURE_NAMES],
    )
    run_mean_rows = [
        *_run_feature_mean_rows(rows_3d, THREE_D_FEATURE_NAMES, "3d"),
        *_run_feature_mean_rows(rows_2d, FEATURE_NAMES, "2d"),
    ]
    _write_csv(
        run_means_path,
        run_mean_rows,
        [
            "domain",
            "feature",
            "requested_duration_s",
            "actual_duration_s",
            "n_flagella",
            "attach_seed",
            "phase_seed",
            "run_id",
            "run_shape_pass",
            "all_windows_shape_pass",
            "window_count",
            "run_mean",
        ],
    )
    seed_rows = [
        *(_seed_effect_rows(run_mean_rows, qc_scope="all_labeled")),
        *(_seed_effect_rows(run_mean_rows, qc_scope="strict_run_and_windows")),
    ]
    _write_csv(
        seed_effects_path,
        seed_rows,
        [
            "domain",
            "qc_scope",
            "feature",
            "requested_duration_s",
            "actual_duration_s",
            "n_flagella",
            "run_count",
            "attach_level_count",
            "phase_level_count",
            "factorial_complete",
            "mean",
            "between_run_std",
            "attach_eta_squared",
            "phase_eta_squared",
            "residual_eta_squared",
        ],
    )
    plot_paths = [
        *_plot_time_features(
            rows_3d,
            domain="3d",
            features=PLOT_FEATURES["3d"],
            output_dir=plots_dir,
        ),
        *_plot_seed_heatmaps(run_mean_rows, output_dir=plots_dir),
        *_plot_time_features(
            rows_2d,
            domain="2d",
            features=PLOT_FEATURES["2d"],
            output_dir=plots_dir,
        ),
    ]
    window_counts: dict[str, int] = {}
    for row in rows_2d:
        key = f"{float(row['requested_duration_s']):g}"
        window_counts[key] = window_counts.get(key, 0) + 1
    created_at = datetime.now(ZoneInfo("Asia/Tokyo")).isoformat()
    attach_seeds = sorted({int(row["attach_seed"]) for row in eligible_samples})
    phase_seeds = sorted({int(row["phase_seed"]) for row in eligible_samples})
    manifest = {
        "pipeline_name": "phase4_duration_seed_study",
        "created_at": created_at,
        "dataset_dir": str(cfg.dataset_dir),
        "dataset_version": dataset_version,
        "dataset_revision": dataset_revision,
        "durations_s_requested": list(cfg.durations_s),
        "frame_rate_hz": cfg.frame_rate_hz,
        "body_centered_2d": True,
        "study": {
            "allowed_n_flagella": list(cfg.allowed_n_flagella),
            "durations_s": list(cfg.durations_s),
        },
        "clip": {
            "frame_rate_hz": cfg.frame_rate_hz,
            "crop_size_px": cfg.crop_size_px,
            "pixel_size_um": cfg.pixel_size_um,
        },
        "seeds": {
            "attach": attach_seeds,
            "phase": phase_seeds,
        },
        "invocation": {
            "config_path": str(cfg.config_path) if cfg.config_path else None,
            "cli_overrides": list(cfg.cli_overrides),
        },
        "source_dataset_manifest": {
            "path": str(dataset_manifest_path),
            "git": dataset_manifest.get("git"),
        },
        "independent_unit": "run_id / group_key",
        "within_run_windows_are_independent": False,
        "sample_count": len(eligible_samples),
        "window_counts_2d": window_counts,
        "outputs": {
            "window_features_3d_csv": str(features_3d_path),
            "window_features_2d_csv": str(features_2d_path),
            "run_feature_means_csv": str(run_means_path),
            "seed_effects_csv": str(seed_effects_path),
            "plots": plot_paths,
        },
        "git": _git_info(),
        "environment": _environment_info(),
    }
    (cfg.output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (cfg.output_dir / "run.log").write_text(
        "\n".join(
            [
                f"created_at={created_at}",
                f"config_path={manifest['invocation']['config_path']}",
                f"cli_overrides={json.dumps(list(cfg.cli_overrides))}",
                f"dataset_dir={cfg.dataset_dir}",
                f"output_dir={cfg.output_dir}",
                f"attach_seeds={attach_seeds}",
                f"phase_seeds={phase_seeds}",
                f"sample_count={len(eligible_samples)}",
                f"git_commit={manifest['git']['commit']}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return cfg.output_dir
