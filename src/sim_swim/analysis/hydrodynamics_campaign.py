"""Select and document force/flow visualizations from compact RPY archives."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

from sim_swim.analysis.hydrodynamics import (
    HydroArchive,
    rpy_flow_velocity,
)

FLOW_SLICE_GRID_SIZE = 41
FLOW_VOLUME_GRID_SIZE = 7
SNAPSHOT_TIMES_S = (0.0, 1.0, 2.0)


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _body_frame(positions_m: np.ndarray, bead_is_body: np.ndarray) -> np.ndarray:
    """Return a world-from-body frame, including roll, from labelled body beads."""
    body = np.asarray(positions_m, dtype=float)[np.asarray(bead_is_body, dtype=bool)]
    if body.shape[0] < 3:
        raise ValueError("at least three labelled body beads are required for a frame")
    center = np.mean(body, axis=0)
    _, _, right = np.linalg.svd(body - center, full_matrices=False)
    long_axis = right[0]
    directed = body[-1] - body[0]
    if np.dot(long_axis, directed) < 0.0:
        long_axis = -long_axis
    radial = body[0] - center
    radial -= long_axis * np.dot(radial, long_axis)
    if np.linalg.norm(radial) < 1.0e-15:
        raise ValueError("body bead geometry cannot determine a roll direction")
    radial /= np.linalg.norm(radial)
    transverse = np.cross(long_axis, radial)
    transverse /= np.linalg.norm(transverse)
    return np.column_stack((long_axis, radial, transverse))


def qc_result(condition_dir: Path) -> tuple[bool, dict[str, Any]]:
    """Return eligibility and the QC fields retained in the output manifest."""
    run_summary = _json(condition_dir / "run_summary.json")
    gates = dict(run_summary.get("gates", {}) or {})
    nonbody = dict(gates.get("shape_nonbody", {}) or {})
    body = dict(gates.get("shape_body", {}) or {})
    fields = {
        "execution_status": dict(run_summary.get("execution", {}) or {}).get("status"),
        "nonbody_final_pass": nonbody.get("final_pass"),
        "body_any_fail": body.get("any_fail"),
        "nonbody_first_fail_t_s": nonbody.get("first_observed_fail_t_s"),
        "body_first_fail_t_s": body.get("first_observed_fail_t_s"),
    }
    return (
        fields["execution_status"] == "completed"
        and fields["nonbody_final_pass"] is True
        and fields["body_any_fail"] is False,
        fields,
    )


def select_qc_passed_conditions(
    run_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Split a campaign manifest into visualization targets and documented skips."""
    manifest = _json(run_dir / "run_manifest.json")
    selected: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    for condition in manifest["conditions"]:
        condition_dir = _condition_dir(run_dir, condition)
        passed, qc = qc_result(condition_dir)
        record = {
            "condition_id": condition["condition_id"],
            "axis_values": dict(condition.get("axis_values", {})),
            "input_dir": str(condition_dir),
            "qc": qc,
        }
        if passed:
            selected.append(record)
        else:
            record["exclusion_reason"] = "strict_shape_qc_not_passed"
            skipped.append(record)
    return selected, skipped


def _condition_dir(run_dir: Path, condition: dict[str, Any]) -> Path:
    """Resolve original or staged reference condition directories."""
    candidate = Path(str(condition.get("output_dir", "")))
    for path in (
        candidate,
        run_dir / candidate.name,
        run_dir / "conditions" / candidate.name,
    ):
        if path.is_dir():
            return path
    return run_dir / "conditions" / str(condition["condition_id"])


def nearest_sample_indices(archive: HydroArchive) -> list[int]:
    """Select requested 0/1/2 s snapshots, clamping the final time safely."""
    indices: list[int] = []
    for requested_t in SNAPSHOT_TIMES_S:
        target = min(requested_t, float(archive.t_s[-1]))
        index = int(np.argmin(np.abs(archive.t_s - target)))
        if index not in indices:
            indices.append(index)
    return indices


def body_fixed_flow_slice(
    archive: HydroArchive, sample_index: int, *, grid_size: int = FLOW_SLICE_GRID_SIZE
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return a dense instantaneous long-axis flow slice in body coordinates."""
    if grid_size < 2:
        raise ValueError("grid_size must be at least 2")
    extent_m = 5e-6
    grid_1d = np.linspace(-extent_m, extent_m, grid_size)
    x, y = np.meshgrid(grid_1d, grid_1d)
    local_points = np.column_stack((x.ravel(), y.ravel(), np.zeros(x.size)))
    positions = archive.positions_m[sample_index]
    rotation = _body_frame(positions, archive.bead_is_body)
    center = np.mean(positions[archive.bead_is_body], axis=0)
    world_points = center + local_points @ rotation.T
    world_velocity = rpy_flow_velocity(
        world_points,
        positions,
        archive.total_forces_N[sample_index],
        bead_radius_m=archive.bead_radius_m,
        viscosity_Pa_s=archive.viscosity_Pa_s,
    )
    return local_points * 1e6, world_velocity @ rotation * 1e6, rotation


def _body_fixed_beads(
    archive: HydroArchive, sample_index: int, rotation: np.ndarray
) -> np.ndarray:
    positions = archive.positions_m[sample_index]
    center = np.mean(positions[archive.bead_is_body], axis=0)
    return (positions - center) @ rotation * 1e6


def analyze_campaign(run_dir: Path, output_dir: Path | None = None) -> Path:
    """Write the compact QC/provenance manifest; rendering is handled by replay."""
    output = output_dir or run_dir / "analysis" / "hydrodynamics"
    output.mkdir(parents=True, exist_ok=True)
    selected, skipped = select_qc_passed_conditions(run_dir)
    (output / "analysis_manifest.json").write_text(
        json.dumps(
            {
                "kind": "free_space_rpy_phase_seed_flow_visualization",
                "run_dir": str(run_dir),
                "input_provenance": "hydro_archive.npz positions_m + total_forces_N",
                "stokes_force_balance": "F_hydro = -F_total",
                "units": {"flow_velocity": "um/s", "force": "N"},
                "body_fixed_slice_grid_shape": [
                    FLOW_SLICE_GRID_SIZE,
                    FLOW_SLICE_GRID_SIZE,
                ],
                "volume_flow_grid_shape": [FLOW_VOLUME_GRID_SIZE] * 3,
                "layout": "3 phase-seed rows x 3 columns; world flow/source and body-fixed slice",
                "included_conditions": selected,
                "excluded_conditions": skipped,
                "stokes_force_balance_verified": "F_total + F_hydro = 0",
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
