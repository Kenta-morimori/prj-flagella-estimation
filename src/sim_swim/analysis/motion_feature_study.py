"""Reusable 3D/2D motion-feature analysis for generic multi-run campaigns."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import yaml

from flagella_estimation.phase3.feature_comparison import body_axis_angles_rad
from flagella_estimation.phase3.render import render_clip_array


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
    observation_frame_rate_hz: float | None = None
    flagella_axis_plot_bin_s: float | None = None
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
    observation, plot = data.get("observation", {}) or {}, data.get("plot", {}) or {}
    if "frame_rate_hz" in study:
        raise ValueError(
            "study.frame_rate_hz is ambiguous; use observation.frame_rate_hz "
            "for 2D output-frame sampling"
        )
    if (
        plot.get("flagella_axis_time_bin_s") is not None
        and float(plot["flagella_axis_time_bin_s"]) <= 0
    ):
        raise ValueError("plot.flagella_axis_time_bin_s must be > 0 when set")
    basis = projection.get("basis", [[1, 0, 0], [0, 1, 0]])
    if len(basis) != 2 or any(len(row) != 3 for row in basis):
        raise ValueError("projection.basis must contain two 3D vectors")
    return MotionFeatureStudyConfig(
        run_dir=Path(str(data.get("run_dir", ""))),
        output_dir=Path(str(data.get("output_dir") or _default_output_dir())),
        config_path=path,
        durations_s=tuple(float(x) for x in study.get("durations_s", [0.25, 0.5, 1.0])),
        observation_frame_rate_hz=(
            float(observation["frame_rate_hz"])
            if observation.get("frame_rate_hz") is not None
            else None
        ),
        flagella_axis_plot_bin_s=(
            float(plot["flagella_axis_time_bin_s"])
            if plot.get("flagella_axis_time_bin_s") is not None
            else None
        ),
        allowed_n_flagella=tuple(
            int(x) for x in study.get("allowed_n_flagella", [1, 2, 3, 4])
        ),
        projection_basis=tuple(tuple(float(v) for v in row) for row in basis),  # type: ignore[arg-type]
        overwrite=_strict_bool(data.get("overwrite", False), "overwrite"),
    )


def _strict_bool(value: Any, key: str) -> bool:
    """Accept YAML booleans only; never silently coerce a typo to True."""
    if type(value) is not bool:
        raise ValueError(f"{key} must be a boolean (true or false)")
    return value


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


def _camera_basis(basis: np.ndarray) -> np.ndarray:
    """Return the world-to-camera rotation for the requested image-plane basis."""
    normal = np.cross(basis[0], basis[1])
    normal /= np.linalg.norm(normal)
    return np.vstack((basis, normal))


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _condition_dir(run_dir: Path, condition: dict[str, Any]) -> Path:
    recorded = Path(str(condition["output_dir"]))
    candidates = (
        recorded,
        run_dir / recorded.name,
        run_dir / "conditions" / recorded.name,
    )
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    # Preserve the canonical aggregate location in diagnostics.  This also
    # supports compact transfers that retain only ``<root>/<condition_id>``.
    return run_dir / "conditions" / recorded.name


def _select_indices(t: np.ndarray, fps: float) -> np.ndarray:
    if fps <= 0 or not len(t):
        raise ValueError(
            "observation frame rate must be > 0 and archive must not be empty"
        )
    wanted = np.arange(float(t[0]), float(t[-1]) + 1e-12, 1.0 / fps)
    result = np.searchsorted(t, wanted, side="left")
    result = np.clip(result, 0, len(t) - 1)
    if result[-1] != len(t) - 1:
        result = np.append(result, len(t) - 1)
    return np.unique(result)


def _time_windows(t: np.ndarray, duration_s: float) -> list[tuple[float, float, slice]]:
    """Return complete non-overlapping windows bounded in physical time."""
    if duration_s <= 0 or not len(t):
        raise ValueError("duration_s must be > 0 and source times must not be empty")
    start, final = float(t[0]), float(t[-1])
    starts = np.arange(start, final - duration_s + 1e-12, duration_s)
    windows: list[tuple[float, float, slice]] = []
    for window_start in starts:
        window_end = float(window_start + duration_s)
        left = int(np.searchsorted(t, window_start, side="left"))
        right = int(np.searchsorted(t, window_end, side="left"))
        if right > left:
            windows.append((float(window_start), window_end, slice(left, right)))
    return windows


def _observation_frame_rate(
    manifest: dict[str, Any], cfg: MotionFeatureStudyConfig
) -> float:
    values = {
        float(rate)
        for condition in manifest.get("conditions", [])
        if (
            rate := condition.get("config_overrides", {})
            .get("output_sampling", {})
            .get("fps_out_2d")
        )
        is not None
    }
    if len(values) > 1:
        raise ValueError(
            "run_manifest contains multiple output_sampling.fps_out_2d values; "
            "an analysis requires one observation frame rate"
        )
    if cfg.observation_frame_rate_hz is not None:
        if cfg.observation_frame_rate_hz <= 0:
            raise ValueError("observation.frame_rate_hz must be > 0")
        if values and not math.isclose(
            cfg.observation_frame_rate_hz,
            next(iter(values)),
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            raise ValueError(
                "observation.frame_rate_hz must match the run_manifest "
                "output_sampling.fps_out_2d"
            )
        return cfg.observation_frame_rate_hz
    if len(values) != 1:
        raise ValueError(
            "Could not determine one output_sampling.fps_out_2d from run_manifest; "
            "set observation.frame_rate_hz explicitly"
        )
    return values.pop()


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


def _quaternion_from_x_axis(axis: np.ndarray) -> tuple[float, float, float, float]:
    """Return an xyzw quaternion rotating the body x-axis onto ``axis``."""
    axis = np.asarray(axis, dtype=float)
    norm = float(np.linalg.norm(axis))
    if not np.isfinite(axis).all() or norm <= 1e-12:
        return (0.0, 0.0, 0.0, 1.0)
    target = axis / norm
    dot = float(np.clip(target[0], -1.0, 1.0))
    if dot <= -1.0 + 1e-12:
        return (0.0, 1.0, 0.0, 0.0)
    cross = np.cross(np.asarray([1.0, 0.0, 0.0]), target)
    quaternion = np.r_[cross, 1.0 + dot]
    quaternion /= np.linalg.norm(quaternion)
    return tuple(float(value) for value in quaternion)


def _camera_states(states: list[Any], basis: np.ndarray) -> list[Any]:
    """Express states in a camera whose image axes are ``basis``.

    The capsule renderer is intentionally XY-only.  Rotating the state first
    makes its silhouette and every projected latent feature share one camera.
    """
    from sim_swim.sim.core import SimulationState

    camera = _camera_basis(basis)
    output: list[SimulationState] = []
    for state in states:
        body_axis = _body_axes(np.asarray([state.quaternion], dtype=float))[0]
        output.append(
            SimulationState(
                t=float(state.t),
                position_um=tuple(camera @ np.asarray(state.position_um, dtype=float)),
                quaternion=_quaternion_from_x_axis(camera @ body_axis),
                velocity_um_s=tuple(
                    camera @ np.asarray(state.velocity_um_s, dtype=float)
                ),
                omega_rad_s=tuple(camera @ np.asarray(state.omega_rad_s, dtype=float)),
                bead_positions_um=np.empty((0, 3)),
            )
        )
    return output


def _pixel_axis_velocity(
    states: list[Any], t: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    clip, geometries = render_clip_array(states, image_size_px=96, pixel_size_um=0.1)
    angles = body_axis_angles_rad(clip)
    axis = _axis_from_angles(angles)
    axis[[geometry.body_axis_angle_rad is None for geometry in geometries]] = np.nan
    return axis, _axis_velocity(axis, t)


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


def _mean_trace(
    rows: list[dict[str, Any]], feature: str, time_key: str
) -> tuple[list[float], list[float]]:
    grouped: dict[float, list[float]] = {}
    for row in rows:
        grouped.setdefault(float(row[time_key]), []).append(float(row[feature]))
    times = sorted(grouped)
    return times, [_mean(np.asarray(grouped[t])) for t in times]


def _plot_rows(
    rows: list[dict[str, Any]], feature: str, time_key: str, bin_s: float | None
) -> list[dict[str, Any]]:
    """Reduce only dense flagella-axis traces for legible rendering."""
    if feature != "mean_flagella_axis_angular_velocity_rad_s" or bin_s is None:
        return rows
    if bin_s <= 0:
        raise ValueError("plot.flagella_axis_time_bin_s must be > 0")
    origin = min(float(row[time_key]) for row in rows)
    buckets: dict[tuple[str, int, int], list[float]] = {}
    for row in rows:
        ratio = (float(row[time_key]) - origin) / bin_s
        nearest = round(ratio)
        bin_index = (
            int(nearest)
            if math.isclose(ratio, nearest, rel_tol=0.0, abs_tol=1e-9)
            else int(math.floor(ratio))
        )
        key = (
            str(row["sample_id"]),
            int(row["n_flagella"]),
            bin_index,
        )
        buckets.setdefault(key, []).append(float(row[feature]))
    return [
        {
            "sample_id": sample_id,
            "n_flagella": n_flagella,
            time_key: origin + bin_index * bin_s,
            feature: _mean(np.asarray(values)),
        }
        for (sample_id, n_flagella, bin_index), values in sorted(buckets.items())
    ]


def _add_style_legend(axis: Any) -> None:
    n_legend = axis.legend(frameon=False, title="mean by n")
    axis.add_artist(n_legend)
    axis.legend(
        handles=[
            Line2D(
                [],
                [],
                color="#444444",
                linestyle="--",
                linewidth=0.9,
                label="individual run",
            ),
            Line2D([], [], color="#444444", linewidth=2.4, label="mean"),
        ],
        frameon=False,
        loc="upper right",
    )


def _plot_series(
    rows: list[dict[str, Any]],
    output_dir: Path,
    domain: str,
    flagella_axis_plot_bin_s: float | None,
) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[str] = []
    colors = {1: "#0072B2", 2: "#D55E00", 3: "#009E73", 4: "#CC79A7"}
    for feature, name in PLOT_NAMES.items():
        feature_rows = _plot_rows(rows, feature, "t_s", flagella_axis_plot_bin_s)
        figure, axis = plt.subplots(figsize=(7, 4))
        for n in sorted({int(row["n_flagella"]) for row in feature_rows}):
            subset = [row for row in feature_rows if int(row["n_flagella"]) == n]
            for sample in sorted({str(row["sample_id"]) for row in subset}):
                part = sorted(
                    (row for row in subset if row["sample_id"] == sample),
                    key=lambda x: float(x["t_s"]),
                )
                axis.plot(
                    [float(x["t_s"]) for x in part],
                    [float(x[feature]) for x in part],
                    color=colors[n],
                    alpha=0.35,
                    linestyle="--",
                    linewidth=0.9,
                )
            times, values = _mean_trace(subset, feature, "t_s")
            axis.plot(
                times,
                values,
                color=colors[n],
                label=f"n={n}" + (" diagnostic" if n == 4 else ""),
                linewidth=2.4,
            )
        axis.set(title=f"{domain}: {name}", xlabel="time (s)", ylabel=name)
        if (
            feature == "mean_flagella_axis_angular_velocity_rad_s"
            and flagella_axis_plot_bin_s is not None
        ):
            axis.set_title(
                f"{domain}: {name} ({flagella_axis_plot_bin_s:g} s bin mean)"
            )
        axis.grid(alpha=0.25)
        _add_style_legend(axis)
        path = output_dir / f"{domain}_{name}.png"
        figure.tight_layout()
        figure.savefig(path, dpi=150)
        plt.close(figure)
        paths.append(str(path))
    return paths


def _is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _validate_output_dir(
    output_dir: Path, run_dir: Path, manifest: dict[str, Any]
) -> None:
    """Reject overwrite targets that could erase raw source artifacts."""
    output = output_dir.resolve()
    run = run_dir.resolve()
    # An analysis child such as ``run/analysis/motion_features`` is valid, but
    # the run root itself and its ancestors would delete the source archive.
    if _is_within(run, output):
        raise ValueError(
            "output_dir must not equal or contain the raw run directory: "
            f"output_dir={output}, run_dir={run}"
        )
    for condition in manifest.get("conditions", []):
        raw = _condition_dir(run_dir, condition).resolve()
        if _is_within(raw, output) or _is_within(output, raw):
            raise ValueError(
                "output_dir must not overlap a raw condition directory: "
                f"output_dir={output}, condition_dir={raw}"
            )


def _analysis_provenance(config_path: Path | None) -> dict[str, Any]:
    config = config_path.resolve() if config_path is not None else None
    config_sha256 = (
        hashlib.sha256(config.read_bytes()).hexdigest()
        if config is not None and config.is_file()
        else None
    )
    try:
        git_commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        git_commit = None
    return {
        "config_path": str(config) if config is not None else None,
        "config_sha256": config_sha256,
        "analysis_git_commit": git_commit,
        "python_version": sys.version,
        "platform": platform.platform(),
        "dependencies": {
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
            "pyyaml": yaml.__version__,
        },
    }


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
                        alpha=0.35,
                        linestyle="--",
                        linewidth=0.9,
                    )
                times, values = _mean_trace(subset, feature, "t_start_s")
                axis.plot(
                    times,
                    values,
                    color=colors[n],
                    linewidth=2.4,
                    label=f"n={n}" + (" diagnostic" if n == 4 else ""),
                )
            axis.set(title=f"{duration:g} s", xlabel="window start (s)", ylabel=name)
            axis.grid(alpha=0.25)
        _add_style_legend(axes[0, 0])
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
    observation_frame_rate_hz = _observation_frame_rate(manifest, cfg)
    _validate_output_dir(cfg.output_dir, cfg.run_dir, manifest)
    if cfg.output_dir.exists():
        if not cfg.overwrite:
            raise FileExistsError(
                f"{cfg.output_dir} already exists; rerun with --overwrite to replace "
                "it, or set output_dir to a new location"
            )
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
            t3d = np.asarray(archive["t"], dtype=float)
            positions3d = np.asarray(archive["position_um"], dtype=float)
            quaternions3d = np.asarray(archive["quaternion"], dtype=float)
        observation_indices = _select_indices(t3d, observation_frame_rate_hz)
        t2d = t3d[observation_indices]
        positions2d = positions3d[observation_indices]
        quaternions2d = quaternions3d[observation_indices]
        from sim_swim.sim.core import SimulationState

        states = [
            SimulationState(
                t=float(t2d[i]),
                position_um=tuple(positions2d[i]),
                quaternion=tuple(quaternions2d[i]),
                velocity_um_s=(0.0, 0.0, 0.0),
                omega_rad_s=(0.0, 0.0, 0.0),
                bead_positions_um=np.empty((0, 3)),
            )
            for i in range(len(t2d))
        ]
        body3d = _body_axes(quaternions3d)
        flag3d = _flag_axes(required[2], t3d, n)
        body2d, body2d_omega = _pixel_axis_velocity(_camera_states(states, basis), t2d)
        flag2d = _project(flag3d[observation_indices], basis)
        centroid2d = positions2d @ basis.T
        speed3d = np.r_[
            np.nan,
            np.linalg.norm(np.diff(positions3d, axis=0), axis=1) / np.diff(t3d),
        ]
        speed2d = np.r_[
            np.nan, np.linalg.norm(np.diff(centroid2d, axis=0), axis=1) / np.diff(t2d)
        ]
        body3d_omega, flag3d_omega, flag2d_omega = (
            _axis_velocity(body3d, t3d),
            _axis_velocity(flag3d, t3d),
            _axis_velocity(flag2d, t2d),
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
        base3d = {**base, "sampling_source": "simulation_step"}
        base2d = {
            **base,
            "sampling_source": "2d_output_frame",
            "observation_frame_rate_hz": observation_frame_rate_hz,
        }
        for i in range(len(t3d)):
            series3d.append(
                {
                    **base3d,
                    "t_s": t3d[i],
                    "body_centroid_x_um": positions3d[i, 0],
                    "body_centroid_y_um": positions3d[i, 1],
                    "body_centroid_z_um": positions3d[i, 2],
                    "speed_um_s": speed3d[i],
                    "body_axis_angular_velocity_rad_s": body3d_omega[i],
                    "mean_flagella_axis_angular_velocity_rad_s": flag3d_omega[i],
                    "body_flagella_axis_angle_deg": relation3d[i],
                }
            )
        for i in range(len(t2d)):
            series2d.append(
                {
                    **base2d,
                    "t_s": t2d[i],
                    "projected_centroid_u_um": centroid2d[i, 0],
                    "projected_centroid_v_um": centroid2d[i, 1],
                    "speed_um_s": speed2d[i],
                    "body_axis_angular_velocity_rad_s": body2d_omega[i],
                    "mean_flagella_axis_angular_velocity_rad_s": flag2d_omega[i],
                    "body_flagella_axis_angle_deg": relation2d[i],
                }
            )
        for duration in cfg.durations_s:
            frame_windows = _time_windows(t3d, duration)
            ranges = [(start, end) for start, end, _ in frame_windows]
            qcs = _window_qc(required[1], ranges)
            for index, (start, end, slice3d) in enumerate(frame_windows):
                slice2d = slice(
                    int(np.searchsorted(t2d, start, side="left")),
                    int(np.searchsorted(t2d, end, side="left")),
                )
                common = {
                    "requested_duration_s": duration,
                    "window_index": index,
                    "t_start_s": start,
                    "t_end_s": end,
                    **qcs[index],
                }
                for (
                    target,
                    domain_base,
                    sample_slice,
                    speed,
                    body_omega,
                    flag_omega,
                    relation,
                ) in (
                    (
                        windows3d,
                        base3d,
                        slice3d,
                        speed3d,
                        body3d_omega,
                        flag3d_omega,
                        relation3d,
                    ),
                    (
                        windows2d,
                        base2d,
                        slice2d,
                        speed2d,
                        body2d_omega,
                        flag2d_omega,
                        relation2d,
                    ),
                ):
                    target.append(
                        {
                            **domain_base,
                            **common,
                            "sample_count": sample_slice.stop - sample_slice.start,
                            "mean_speed_um_s": _mean(speed[sample_slice]),
                            "mean_body_axis_angular_velocity_rad_s": _mean(
                                body_omega[sample_slice]
                            ),
                            "mean_flagella_axis_angular_velocity_rad_s": _mean(
                                flag_omega[sample_slice]
                            ),
                            "mean_body_flagella_axis_angle_deg": _mean(
                                relation[sample_slice]
                            ),
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
        "time_series": _plot_series(
            plot_series3d, time_dir, "3D", cfg.flagella_axis_plot_bin_s
        )
        + _plot_series(plot_series2d, time_dir, "2D", cfg.flagella_axis_plot_bin_s),
        "windows": _plot_windows(plot_windows3d, window_dir, "3D")
        + _plot_windows(plot_windows2d, window_dir, "2D"),
    }
    output = {
        "pipeline_name": "phase2_motion_feature_study",
        "run_dir": str(cfg.run_dir),
        "durations_s": list(cfg.durations_s),
        "sampling_contract": {
            "3d": "all saved simulation steps",
            "2d": {
                "source": "output_sampling.fps_out_2d",
                "frame_rate_hz": observation_frame_rate_hz,
            },
            "windows": "non-overlapping physical-time intervals [start_s, end_s)",
        },
        "plot_contract": {
            "flagella_axis_time_bin_s": cfg.flagella_axis_plot_bin_s,
            "raw_csv_is_unbinned": True,
        },
        "projection_basis": basis.tolist(),
        "analysis_provenance": _analysis_provenance(cfg.config_path),
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
