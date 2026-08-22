"""Reusable 3D/2D motion-feature analysis for generic multi-run campaigns."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
import json
import math
from pathlib import Path
import shutil
from typing import Any
from zoneinfo import ZoneInfo

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yaml

from flagella_estimation.phase3.feature_comparison import body_axis_angles_rad
from flagella_estimation.phase3.render import render_clip_array
from flagella_estimation.phase3.windows import generate_windows


FEATURES = (
    "mean_speed_um_s",
    "mean_body_axis_angular_velocity_rad_s",
    "mean_flagella_axis_angular_velocity_rad_s",
    "mean_body_flagella_axis_angle_deg",
)

PLOT_NAMES = {
    "speed_um_s": "speed",
    "body_axis_angular_velocity_rad_s": "body_axis_angular_velocity",
    "mean_flagella_axis_angular_velocity_rad_s": "flagella_axis_angular_velocity",
    "body_flagella_axis_angle_deg": "body_flagella_axis_angle",
}

WINDOW_FEATURES = dict(zip(FEATURES, PLOT_NAMES.values(), strict=True))


@dataclass(frozen=True)
class MotionFeatureStudyConfig:
    run_dir: Path
    output_dir: Path
    config_path: Path | None = None
    durations_s: tuple[float, ...] = (0.25, 0.5, 1.0)
    frame_rate_hz: float = 25.0
    allowed_n_flagella: tuple[int, ...] = (1, 2, 3, 4)
    projection_basis: tuple[tuple[float, float, float], tuple[float, float, float]] = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
    )
    overwrite: bool = False


def _default_output_dir() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / "phase2_motion_feature_study"
    )


def _override(data: dict[str, Any], raw: str) -> None:
    if "=" not in raw:
        raise ValueError(f"Expected KEY=VALUE: {raw}")
    key, value = raw.split("=", 1)
    node = data
    parts = key.split(".")
    for part in parts[:-1]:
        node = node.setdefault(part, {})
        if not isinstance(node, dict):
            raise ValueError(f"Override conflicts with scalar: {raw}")
    node[parts[-1]] = yaml.safe_load(value)


def load_config(
    path: Path | None, overrides: list[str] | None = None
) -> MotionFeatureStudyConfig:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) if path else {}
    data = data or {}
    for raw in overrides or []:
        _override(data, raw)
    study, projection = data.get("study", {}) or {}, data.get("projection", {}) or {}
    basis = projection.get("basis", [[1, 0, 0], [0, 1, 0]])
    if len(basis) != 2 or any(len(row) != 3 for row in basis):
        raise ValueError("projection.basis must contain two 3D vectors")
    return MotionFeatureStudyConfig(
        run_dir=Path(str(data.get("run_dir", ""))),
        output_dir=Path(str(data.get("output_dir") or _default_output_dir())),
        config_path=path,
        durations_s=tuple(float(x) for x in study.get("durations_s", [0.25, 0.5, 1.0])),
        frame_rate_hz=float(study.get("frame_rate_hz", 25.0)),
        allowed_n_flagella=tuple(
            int(x) for x in study.get("allowed_n_flagella", [1, 2, 3, 4])
        ),
        projection_basis=tuple(tuple(float(v) for v in row) for row in basis),  # type: ignore[arg-type]
        overwrite=bool(data.get("overwrite", False)),
    )


def _basis(
    value: tuple[tuple[float, float, float], tuple[float, float, float]],
) -> np.ndarray:
    basis = np.asarray(value, dtype=float)
    if (
        not np.isfinite(basis).all()
        or np.linalg.norm(basis[0]) < 1e-12
        or np.linalg.norm(basis[1]) < 1e-12
    ):
        raise ValueError("projection basis must be finite and non-zero")
    basis[0] /= np.linalg.norm(basis[0])
    basis[1] -= np.dot(basis[1], basis[0]) * basis[0]
    if np.linalg.norm(basis[1]) < 1e-12:
        raise ValueError("projection basis vectors must be independent")
    basis[1] /= np.linalg.norm(basis[1])
    return basis


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _condition_dir(run_dir: Path, condition: dict[str, Any]) -> Path:
    recorded = Path(str(condition["output_dir"]))
    return recorded if recorded.is_dir() else run_dir / recorded.name


def _select_indices(t: np.ndarray, fps: float) -> np.ndarray:
    if fps <= 0 or not len(t):
        raise ValueError("frame_rate_hz must be > 0 and archive must not be empty")
    wanted = np.arange(float(t[0]), float(t[-1]) + 1e-12, 1.0 / fps)
    result = np.searchsorted(t, wanted, side="left")
    result = np.clip(result, 0, len(t) - 1)
    if result[-1] != len(t) - 1:
        result = np.append(result, len(t) - 1)
    return np.unique(result)


def _body_axes(quaternions: np.ndarray) -> np.ndarray:
    x, y, z, w = quaternions.T
    axes = np.column_stack(
        (1 - 2 * (y * y + z * z), 2 * (x * y + z * w), 2 * (x * z - y * y))
    )
    # Correct the third rotation-matrix component (kept separate for readable vectorisation).
    axes[:, 2] = 2 * (x * z - y * w)
    norms = np.linalg.norm(axes, axis=1)
    axes[norms > 0] /= norms[norms > 0, None]
    axes[norms <= 0] = np.nan
    return _continuous(axes)


def _continuous(vectors: np.ndarray) -> np.ndarray:
    result = np.asarray(vectors, dtype=float).copy()
    previous: np.ndarray | None = None
    for index, vector in enumerate(result):
        if not np.isfinite(vector).all() or np.linalg.norm(vector) < 1e-12:
            result[index] = np.nan
            continue
        result[index] = vector / np.linalg.norm(vector)
        if previous is not None and np.dot(previous, result[index]) < 0:
            result[index] *= -1
        previous = result[index]
    return result


def _axis_velocity(vectors: np.ndarray, t: np.ndarray) -> np.ndarray:
    output = np.full(len(t), np.nan)
    for index in range(1, len(t)):
        if not (
            np.isfinite(vectors[index - 1]).all() and np.isfinite(vectors[index]).all()
        ):
            continue
        dt = t[index] - t[index - 1]
        if dt <= 0:
            continue
        output[index] = (
            math.acos(
                float(np.clip(abs(np.dot(vectors[index - 1], vectors[index])), -1, 1))
            )
            / dt
        )
    return output


def _project(vectors: np.ndarray, basis: np.ndarray) -> np.ndarray:
    projected = vectors @ basis.T
    norms = np.linalg.norm(projected, axis=1)
    projected[norms > 1e-12] /= norms[norms > 1e-12, None]
    projected[norms <= 1e-12] = np.nan
    return projected


def _axis_from_angles(angles: np.ndarray) -> np.ndarray:
    return np.column_stack((np.cos(angles), np.sin(angles)))


def _pixel_axis_velocity(
    states: list[Any], t: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    clip, _ = render_clip_array(states, image_size_px=96, pixel_size_um=0.1)
    angles = body_axis_angles_rad(clip)
    return _axis_from_angles(angles), _axis_velocity(_axis_from_angles(angles), t)


def _flag_axes(axis_path: Path, times: np.ndarray, n_flagella: int) -> np.ndarray:
    keys = {round(float(t), 9): index for index, t in enumerate(times)}
    grouped: dict[int, list[np.ndarray]] = {index: [] for index in range(len(times))}
    with axis_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            index = keys.get(round(float(row["t_s"]), 9))
            if index is None:
                continue
            vector = np.asarray(
                [row["axis_dir_x"], row["axis_dir_y"], row["axis_dir_z"]], dtype=float
            )
            if np.isfinite(vector).all() and np.linalg.norm(vector) > 1e-12:
                grouped[index].append(vector / np.linalg.norm(vector))
    mean = np.full((len(times), 3), np.nan)
    for index, vectors in grouped.items():
        if len(vectors) != n_flagella:
            continue
        ref = vectors[0]
        aligned = np.asarray([v if np.dot(v, ref) >= 0 else -v for v in vectors])
        candidate = aligned.mean(axis=0)
        if np.linalg.norm(candidate) > 1e-12:
            mean[index] = candidate / np.linalg.norm(candidate)
    return _continuous(mean)


def _window_qc(path: Path, windows: list[tuple[float, float]]) -> list[dict[str, Any]]:
    output = [
        {
            "window_shape_pass": True,
            "window_first_fail_t_s": float("nan"),
            "window_first_fail_category": "none",
        }
        for _ in windows
    ]
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            t = float(row["t_s"])
            for index, (start, end) in enumerate(windows):
                if not start - 1e-12 <= t < end - 1e-12:
                    continue
                valid = str(row.get("finite_pass", "")).lower() in {
                    "true",
                    "1",
                } and str(row.get("shape_pass_nonbody_strict", "")).lower() in {
                    "true",
                    "1",
                }
                if not valid and output[index]["window_shape_pass"]:
                    output[index].update(
                        window_shape_pass=False,
                        window_first_fail_t_s=t,
                        window_first_fail_category=row.get(
                            "first_fail_category_nonbody_strict"
                        )
                        or "unknown",
                    )
    return output


def _strict_run_pass(path: Path) -> bool:
    """Return whether every recorded source step passes the strict QC."""
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            finite = str(row.get("finite_pass", "")).lower() in {"true", "1"}
            shape = str(row.get("shape_pass_nonbody_strict", "")).lower() in {
                "true",
                "1",
            }
            if not finite or not shape:
                return False
    return True


def _mean(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    return float(np.mean(values)) if len(values) else float("nan")


def _write(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = list(rows[0]) if rows else ["sample_id"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _plot_series(
    rows: list[dict[str, Any]], output_dir: Path, domain: str
) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []
    colors = {1: "#0072B2", 2: "#D55E00", 3: "#009E73", 4: "#CC79A7"}
    for feature, name in PLOT_NAMES.items():
        figure, axis = plt.subplots(figsize=(7, 4))
        for n in sorted({int(row["n_flagella"]) for row in rows}):
            subset = [row for row in rows if int(row["n_flagella"]) == n]
            for sample in sorted({str(row["sample_id"]) for row in subset}):
                part = sorted(
                    (row for row in subset if row["sample_id"] == sample),
                    key=lambda x: float(x["t_s"]),
                )
                axis.plot(
                    [float(x["t_s"]) for x in part],
                    [float(x[feature]) for x in part],
                    color=colors[n],
                    alpha=0.3,
                )
            times = sorted({float(row["t_s"]) for row in subset})
            axis.plot(
                times,
                [
                    _mean(
                        np.asarray(
                            [float(x[feature]) for x in subset if float(x["t_s"]) == t]
                        )
                    )
                    for t in times
                ],
                color=colors[n],
                label=f"n={n}" + (" diagnostic" if n == 4 else ""),
                linewidth=2.4,
            )
        axis.set(title=f"{domain}: {name}", xlabel="time (s)", ylabel=name)
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
        path = output_dir / f"{domain}_{name}.png"
        figure.tight_layout()
        figure.savefig(path, dpi=150)
        plt.close(figure)
        paths.append(str(path))
    return paths


def _plot_windows(
    rows: list[dict[str, Any]], output_dir: Path, domain: str
) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []
    colors = {1: "#0072B2", 2: "#D55E00", 3: "#009E73", 4: "#CC79A7"}
    durations = sorted({float(row["requested_duration_s"]) for row in rows})
    for feature, name in WINDOW_FEATURES.items():
        figure, axes = plt.subplots(
            1,
            len(durations),
            figsize=(5 * len(durations), 4),
            sharey=True,
            squeeze=False,
        )
        for axis, duration in zip(axes[0], durations, strict=True):
            duration_rows = [
                row for row in rows if float(row["requested_duration_s"]) == duration
            ]
            for n in sorted({int(row["n_flagella"]) for row in duration_rows}):
                subset = [row for row in duration_rows if int(row["n_flagella"]) == n]
                for sample in sorted({str(row["sample_id"]) for row in subset}):
                    part = sorted(
                        (row for row in subset if row["sample_id"] == sample),
                        key=lambda x: float(x["t_start_s"]),
                    )
                    axis.plot(
                        [float(x["t_start_s"]) for x in part],
                        [float(x[feature]) for x in part],
                        color=colors[n],
                        alpha=0.3,
                    )
                times = sorted({float(row["t_start_s"]) for row in subset})
                axis.plot(
                    times,
                    [
                        _mean(
                            np.asarray(
                                [
                                    float(x[feature])
                                    for x in subset
                                    if float(x["t_start_s"]) == t
                                ]
                            )
                        )
                        for t in times
                    ],
                    color=colors[n],
                    linewidth=2.4,
                    label=f"n={n}" + (" diagnostic" if n == 4 else ""),
                )
            axis.set(title=f"{duration:g} s", xlabel="window start (s)", ylabel=name)
            axis.grid(alpha=0.25)
        axes[0, 0].legend(frameon=False)
        figure.tight_layout()
        path = output_dir / f"{domain}_{name}.png"
        figure.savefig(path, dpi=150)
        plt.close(figure)
        paths.append(str(path))
    return paths


def analyze_motion_feature_study(cfg: MotionFeatureStudyConfig) -> Path:
    if not cfg.run_dir.is_dir():
        raise FileNotFoundError(cfg.run_dir)
    basis = _basis(cfg.projection_basis)
    manifest = _read_json(cfg.run_dir / "run_manifest.json")
    if cfg.output_dir.exists():
        if not cfg.overwrite:
            raise FileExistsError(cfg.output_dir)
        shutil.rmtree(cfg.output_dir)
    cfg.output_dir.mkdir(parents=True)
    series3d: list[dict[str, Any]] = []
    series2d: list[dict[str, Any]] = []
    windows3d: list[dict[str, Any]] = []
    windows2d: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []
    for condition in manifest.get("conditions", []):
        values = condition.get("axis_values", {})
        n = int(values.get("n_flagella", -1))
        if n not in cfg.allowed_n_flagella:
            continue
        sample_id, raw = (
            str(condition["condition_id"]),
            _condition_dir(cfg.run_dir, condition),
        )
        required = [
            raw / name
            for name in (
                "state_archive.npz",
                "step_summary.csv",
                "flag_helix_axis_diagnostics.csv",
            )
        ]
        if not all(path.is_file() for path in required):
            raise FileNotFoundError(f"{sample_id}: missing required raw artifact")
        run_strict_pass = _strict_run_pass(required[1])
        with np.load(required[0], allow_pickle=False) as archive:
            indices = _select_indices(
                np.asarray(archive["t"], dtype=float), cfg.frame_rate_hz
            )
            t = np.asarray(archive["t"], dtype=float)[indices]
            positions = np.asarray(archive["position_um"], dtype=float)[indices]
            quaternions = np.asarray(archive["quaternion"], dtype=float)[indices]
            from sim_swim.sim.core import SimulationState

            states = [
                SimulationState(
                    t=float(t[i]),
                    position_um=tuple(positions[i]),
                    quaternion=tuple(quaternions[i]),
                    velocity_um_s=(0.0, 0.0, 0.0),
                    omega_rad_s=(0.0, 0.0, 0.0),
                    bead_positions_um=np.empty((0, 3)),
                )
                for i in range(len(t))
            ]
        body3d = _body_axes(quaternions)
        flag3d = _flag_axes(required[2], t, n)
        body2d, body2d_omega = _pixel_axis_velocity(states, t)
        flag2d = _project(flag3d, basis)
        centroid2d = positions @ basis.T
        speed3d = np.r_[
            np.nan, np.linalg.norm(np.diff(positions, axis=0), axis=1) / np.diff(t)
        ]
        speed2d = np.r_[
            np.nan, np.linalg.norm(np.diff(centroid2d, axis=0), axis=1) / np.diff(t)
        ]
        body3d_omega, flag3d_omega, flag2d_omega = (
            _axis_velocity(body3d, t),
            _axis_velocity(flag3d, t),
            _axis_velocity(flag2d, t),
        )
        relation3d = np.rad2deg(
            np.arccos(np.clip(np.abs(np.sum(body3d * flag3d, axis=1)), -1, 1))
        )
        relation2d = np.rad2deg(
            np.arccos(np.clip(np.abs(np.sum(body2d * flag2d, axis=1)), -1, 1))
        )
        base = {
            "sample_id": sample_id,
            "n_flagella": n,
            "attach_seed": int(values.get("attach_seed", -1)),
            "phase_seed": int(values.get("phase_seed", -1)),
            "diagnostic_only": n == 4,
            "run_strict_pass": run_strict_pass,
        }
        for i in range(len(t)):
            series3d.append(
                {
                    **base,
                    "t_s": t[i],
                    "body_centroid_x_um": positions[i, 0],
                    "body_centroid_y_um": positions[i, 1],
                    "body_centroid_z_um": positions[i, 2],
                    "speed_um_s": speed3d[i],
                    "body_axis_angular_velocity_rad_s": body3d_omega[i],
                    "mean_flagella_axis_angular_velocity_rad_s": flag3d_omega[i],
                    "body_flagella_axis_angle_deg": relation3d[i],
                }
            )
            series2d.append(
                {
                    **base,
                    "t_s": t[i],
                    "projected_centroid_u_um": centroid2d[i, 0],
                    "projected_centroid_v_um": centroid2d[i, 1],
                    "speed_um_s": speed2d[i],
                    "body_axis_angular_velocity_rad_s": body2d_omega[i],
                    "mean_flagella_axis_angular_velocity_rad_s": flag2d_omega[i],
                    "body_flagella_axis_angle_deg": relation2d[i],
                }
            )
        for duration in cfg.durations_s:
            frame_windows = generate_windows(
                source_frame_count=len(t),
                frame_rate_hz=cfg.frame_rate_hz,
                duration_s=duration,
                policy="non_overlap",
            )
            ranges = [
                (
                    float(t[w.start]),
                    float(t[w.start] + w.frame_count / cfg.frame_rate_hz),
                )
                for w in frame_windows
            ]
            qcs = _window_qc(required[1], ranges)
            for index, window in enumerate(frame_windows):
                sl = slice(window.start, window.end)
                common = {
                    **base,
                    "requested_duration_s": duration,
                    "window_index": index,
                    "t_start_s": ranges[index][0],
                    "t_end_s": ranges[index][1],
                    **qcs[index],
                }
                for target, speed, body_omega, flag_omega, relation in (
                    (windows3d, speed3d, body3d_omega, flag3d_omega, relation3d),
                    (windows2d, speed2d, body2d_omega, flag2d_omega, relation2d),
                ):
                    target.append(
                        {
                            **common,
                            "mean_speed_um_s": _mean(speed[sl]),
                            "mean_body_axis_angular_velocity_rad_s": _mean(
                                body_omega[sl]
                            ),
                            "mean_flagella_axis_angular_velocity_rad_s": _mean(
                                flag_omega[sl]
                            ),
                            "mean_body_flagella_axis_angle_deg": _mean(relation[sl]),
                        }
                    )
    excluded_n4 = sorted(
        {
            str(row["sample_id"])
            for row in windows3d
            if int(row["n_flagella"]) == 4 and not bool(row["run_strict_pass"])
        }
    )
    plot_windows3d = [
        row for row in windows3d if str(row["sample_id"]) not in excluded_n4
    ]
    plot_windows2d = [
        row for row in windows2d if str(row["sample_id"]) not in excluded_n4
    ]
    plot_series3d = [
        row for row in series3d if str(row["sample_id"]) not in excluded_n4
    ]
    plot_series2d = [
        row for row in series2d if str(row["sample_id"]) not in excluded_n4
    ]
    _write(cfg.output_dir / "time_series_3d.csv", series3d)
    _write(cfg.output_dir / "time_series_2d.csv", series2d)
    _write(cfg.output_dir / "window_features_3d.csv", windows3d)
    _write(cfg.output_dir / "window_features_2d.csv", windows2d)
    time_dir, window_dir = cfg.output_dir / "time_series", cfg.output_dir / "windows"
    plots = {
        "time_series": _plot_series(plot_series3d, time_dir, "3D")
        + _plot_series(plot_series2d, time_dir, "2D"),
        "windows": _plot_windows(plot_windows3d, window_dir, "3D")
        + _plot_windows(plot_windows2d, window_dir, "2D"),
    }
    output = {
        "pipeline_name": "phase2_motion_feature_study",
        "run_dir": str(cfg.run_dir),
        "durations_s": list(cfg.durations_s),
        "frame_rate_hz": cfg.frame_rate_hz,
        "projection_basis": basis.tolist(),
        "n_flagella": list(cfg.allowed_n_flagella),
        "independent_unit": "condition/run",
        "observability": {
            "2d_body_axis_angular_velocity": "pixel_observable",
            "2d_speed_and_flagella_features": "projected_latent_simulation_gt",
        },
        "outputs": {
            "time_series_3d": str(cfg.output_dir / "time_series_3d.csv"),
            "time_series_2d": str(cfg.output_dir / "time_series_2d.csv"),
            "window_features_3d": str(cfg.output_dir / "window_features_3d.csv"),
            "window_features_2d": str(cfg.output_dir / "window_features_2d.csv"),
            "plots": plots,
        },
        "plot_exclusions": {"n4_strict_failure": excluded_n4},
        "skipped_conditions": skipped,
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
    }
    (cfg.output_dir / "manifest.json").write_text(
        json.dumps(output, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (cfg.output_dir / "run.log").write_text(
        f"run_dir={cfg.run_dir}\noutput_dir={cfg.output_dir}\n", encoding="utf-8"
    )
    return cfg.output_dir
