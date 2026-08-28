"""Campaign-level quantitative analysis for compact free-space RPY archives."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from sim_swim.analysis.flagella_count_behavior import load_state_archive
from sim_swim.analysis.hydrodynamics import (
    HydroArchive,
    load_hydro_archive,
    rpy_flow_velocity,
    velocity_contributions,
)


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _norm(values: np.ndarray) -> np.ndarray:
    return np.linalg.norm(values, axis=-1)


def _rotation_matrix_xyzw(quaternion: np.ndarray) -> np.ndarray:
    x, y, z, w = quaternion
    return np.asarray(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ],
        dtype=float,
    )


def _qc_fields(run_summary: dict[str, Any]) -> dict[str, Any]:
    gates = dict(run_summary.get("gates", {}) or {})
    nonbody = dict(gates.get("shape_nonbody", {}) or {})
    body = dict(gates.get("shape_body", {}) or {})
    return {
        "qc_execution_status": dict(run_summary.get("execution", {}) or {}).get(
            "status"
        ),
        "qc_nonbody_final_pass": nonbody.get("final_pass"),
        "qc_body_any_fail": body.get("any_fail"),
        "qc_nonbody_first_fail_t_s": nonbody.get("first_observed_fail_t_s"),
        "qc_body_first_fail_t_s": body.get("first_observed_fail_t_s"),
    }


def _source_body_velocity(archive: HydroArchive, sample_index: int) -> dict[str, float]:
    pos = archive.positions_m[sample_index]
    forces = archive.total_forces_N[sample_index]
    body = archive.bead_is_body
    total = velocity_contributions(
        pos,
        forces,
        bead_radius_m=archive.bead_radius_m,
        viscosity_Pa_s=archive.viscosity_Pa_s,
    )
    result = {
        "body_bead_velocity_total_um_s": float(
            np.mean(_norm(total["source"][body])) * 1e6
        ),
        "body_bead_velocity_self_um_s": float(
            np.mean(_norm(total["self"][body])) * 1e6
        ),
        "body_bead_velocity_other_um_s": float(
            np.mean(_norm(total["other"][body])) * 1e6
        ),
    }
    for source_name, mask in [("body", body), ("flagella", ~body)]:
        source = velocity_contributions(
            pos,
            forces,
            bead_radius_m=archive.bead_radius_m,
            viscosity_Pa_s=archive.viscosity_Pa_s,
            source_mask=mask,
        )
        result[f"body_bead_velocity_from_{source_name}_um_s"] = float(
            np.mean(_norm(source["source"][body])) * 1e6
        )
    return result


def analyze_condition(
    condition_dir: Path, condition: dict[str, Any]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Write an archive-aligned timeseries and return one comparison summary row."""
    archive = load_hydro_archive(condition_dir / "hydro_archive.npz")
    states = load_state_archive(condition_dir / "state_archive.npz")
    state_t = np.asarray([state.t for state in states], dtype=float)
    rows: list[dict[str, Any]] = []
    for index, sample_t in enumerate(archive.t_s):
        state = states[int(np.argmin(np.abs(state_t - sample_t)))]
        force = archive.total_forces_N[index]
        center = np.mean(archive.positions_m[index][archive.bead_is_body], axis=0)
        torque = np.sum(np.cross(archive.positions_m[index] - center, force), axis=0)
        row: dict[str, Any] = {
            "condition_id": condition["condition_id"],
            "t_s": float(sample_t),
            "body_translation_um_s": float(np.linalg.norm(state.velocity_um_s)),
            "body_rotation_rad_s": float(np.linalg.norm(state.omega_rad_s)),
            "net_force_N": float(np.linalg.norm(np.sum(force, axis=0))),
            "net_torque_Nm": float(np.linalg.norm(torque)),
        }
        row.update(_source_body_velocity(archive, index))
        rows.append(row)
    fields = list(rows[0])
    path = condition_dir / "hydrodynamics_timeseries.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    summary: dict[str, Any] = {
        "condition_id": condition["condition_id"],
        **condition.get("axis_values", {}),
    }
    for name in fields:
        if name in {"condition_id", "t_s"}:
            continue
        values = np.asarray([float(row[name]) for row in rows], dtype=float)
        summary[f"mean_{name}"] = float(np.mean(values))
        summary[f"max_{name}"] = float(np.max(values))
    summary.update(_qc_fields(_json(condition_dir / "run_summary.json")))
    summary["hydro_archive_format"] = archive.provenance.get("hydrodynamics", {}).get(
        "model", "free_space_rpy"
    )
    summary["hydro_sample_count"] = int(archive.t_s.size)
    return rows, summary


def _axial_slice(
    archive: HydroArchive, states: list[Any]
) -> tuple[np.ndarray, np.ndarray]:
    """Average full-run RPY velocity over a body-fixed x-y (long-axis) section."""
    extent_m = 5e-6
    grid_1d = np.linspace(-extent_m, extent_m, 21)
    x, y = np.meshgrid(grid_1d, grid_1d)
    local_points = np.column_stack((x.ravel(), y.ravel(), np.zeros(x.size)))
    state_t = np.asarray([state.t for state in states], dtype=float)
    mean_velocity = np.zeros_like(local_points)
    for index, sample_t in enumerate(archive.t_s):
        state = states[int(np.argmin(np.abs(state_t - sample_t)))]
        rotation = _rotation_matrix_xyzw(np.asarray(state.quaternion, dtype=float))
        center = np.mean(archive.positions_m[index][archive.bead_is_body], axis=0)
        world_points = center + local_points @ rotation.T
        world_velocity = rpy_flow_velocity(
            world_points,
            archive.positions_m[index],
            archive.total_forces_N[index],
            bead_radius_m=archive.bead_radius_m,
            viscosity_Pa_s=archive.viscosity_Pa_s,
        )
        mean_velocity += world_velocity @ rotation
    return local_points * 1e6, mean_velocity / max(archive.t_s.size, 1) * 1e6


def _plot_summary(rows: list[dict[str, Any]], path: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    for phase in sorted({int(row.get("phase_seed", 0)) for row in rows}):
        selected = [row for row in rows if int(row.get("phase_seed", 0)) == phase]
        selected.sort(key=lambda row: int(row.get("n_flagella", 0)))
        x = [int(row["n_flagella"]) for row in selected]
        axes[0].plot(
            x,
            [row["mean_body_translation_um_s"] for row in selected],
            marker="o",
            label=f"phase {phase}",
        )
        axes[1].plot(
            x,
            [row["max_net_force_N"] for row in selected],
            marker="o",
            label=f"phase {phase}",
        )
    axes[0].set(xlabel="n_flagella", ylabel="mean body translation [um/s]")
    axes[1].set(xlabel="n_flagella", ylabel="max net force residual [N]")
    axes[0].legend()
    axes[1].legend()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def _plot_slice(points_um: np.ndarray, velocity_um_s: np.ndarray, path: Path) -> None:
    figure, axis = plt.subplots(figsize=(5, 4), constrained_layout=True)
    axis.quiver(
        points_um[:, 0],
        points_um[:, 1],
        velocity_um_s[:, 0],
        velocity_um_s[:, 1],
        scale=None,
    )
    axis.set(
        aspect="equal",
        xlabel="body-fixed long axis x [um]",
        ylabel="body-fixed y [um]",
        title="full 0.5 s mean RPY velocity",
    )
    figure.savefig(path, dpi=160)
    plt.close(figure)


def analyze_campaign(run_dir: Path, output_dir: Path | None = None) -> Path:
    """Generate all quantitative #225 artifacts without resimulating."""
    manifest = _json(run_dir / "run_manifest.json")
    output = output_dir or run_dir / "analysis" / "hydrodynamics"
    output.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, Any]] = []
    representative: tuple[HydroArchive, list[Any]] | None = None
    for condition in manifest["conditions"]:
        condition_dir = Path(condition["output_dir"])
        _, summary = analyze_condition(condition_dir, condition)
        summaries.append(summary)
        if (
            representative is None
            and int(condition.get("axis_values", {}).get("n_flagella", 0)) == 3
        ):
            representative = (
                load_hydro_archive(condition_dir / "hydro_archive.npz"),
                load_state_archive(condition_dir / "state_archive.npz"),
            )
    fieldnames = list(summaries[0])
    with (output / "hydrodynamics_comparison.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summaries)
    _plot_summary(summaries, output / "flagella_count_comparison.png")
    if representative is not None:
        _plot_slice(
            *_axial_slice(*representative), output / "body_fixed_axial_flow.png"
        )
    (output / "analysis_manifest.json").write_text(
        json.dumps(
            {
                "kind": "free_space_rpy_hydrodynamics_analysis",
                "run_dir": str(run_dir),
                "comparison_csv": str(output / "hydrodynamics_comparison.csv"),
                "full_duration_average": True,
                "condition_count": len(summaries),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return output


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    print(analyze_campaign(args.run_dir, args.output_dir))
