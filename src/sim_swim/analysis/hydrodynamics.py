"""Compact RPY hydrodynamics archives and post-processing primitives."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np

from sim_swim.dynamics.hydro_rpy import compute_rpy_mobility, compute_rpy_pair_mobility

HYDRO_ARCHIVE_FORMAT = "sim_swim.hydro_archive"
HYDRO_ARCHIVE_VERSION = 1


@dataclass(frozen=True)
class HydroSample:
    """One physically consistent force/position sample from an RPY step."""

    t_s: float
    positions_m: np.ndarray
    total_forces_N: np.ndarray


@dataclass(frozen=True)
class HydroArchive:
    """Portable input for velocity-contribution and flow-field reconstruction."""

    t_s: np.ndarray
    positions_m: np.ndarray
    total_forces_N: np.ndarray
    bead_is_body: np.ndarray
    bead_flagella_id: np.ndarray
    bead_radius_m: float
    viscosity_Pa_s: float
    provenance: dict[str, Any]


def save_hydro_archive(
    path: Path,
    samples: list[HydroSample],
    *,
    bead_is_body: np.ndarray,
    bead_flagella_id: np.ndarray,
    bead_radius_m: float,
    viscosity_Pa_s: float,
    provenance: dict[str, Any],
) -> None:
    """Save only fields needed to reconstruct RPY velocities after a run."""

    if not samples:
        raise ValueError("Cannot save an empty hydrodynamics archive")
    positions = np.asarray([sample.positions_m for sample in samples], dtype=float)
    forces = np.asarray([sample.total_forces_N for sample in samples], dtype=float)
    if (
        positions.shape != forces.shape
        or positions.ndim != 3
        or positions.shape[2] != 3
    ):
        raise ValueError("Hydrodynamics samples must have matching shape (S, N, 3)")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        archive_format=np.asarray(HYDRO_ARCHIVE_FORMAT),
        archive_version=np.asarray(HYDRO_ARCHIVE_VERSION, dtype=np.int64),
        t_s=np.asarray([sample.t_s for sample in samples], dtype=float),
        positions_m=positions,
        total_forces_N=forces,
        bead_is_body=np.asarray(bead_is_body, dtype=bool),
        bead_flagella_id=np.asarray(bead_flagella_id, dtype=np.int64),
        bead_radius_m=np.asarray(bead_radius_m, dtype=float),
        viscosity_Pa_s=np.asarray(viscosity_Pa_s, dtype=float),
        provenance_json=np.asarray(json.dumps(provenance, sort_keys=True)),
    )


def load_hydro_archive(path: Path) -> HydroArchive:
    """Load a versioned archive without pickle-backed metadata."""

    with np.load(path, allow_pickle=False) as data:
        if str(data["archive_format"].item()) != HYDRO_ARCHIVE_FORMAT:
            raise ValueError("Unsupported hydrodynamics archive format")
        if int(np.asarray(data["archive_version"]).item()) != HYDRO_ARCHIVE_VERSION:
            raise ValueError("Unsupported hydrodynamics archive version")
        return HydroArchive(
            t_s=np.asarray(data["t_s"], dtype=float),
            positions_m=np.asarray(data["positions_m"], dtype=float),
            total_forces_N=np.asarray(data["total_forces_N"], dtype=float),
            bead_is_body=np.asarray(data["bead_is_body"], dtype=bool),
            bead_flagella_id=np.asarray(data["bead_flagella_id"], dtype=np.int64),
            bead_radius_m=float(np.asarray(data["bead_radius_m"]).item()),
            viscosity_Pa_s=float(np.asarray(data["viscosity_Pa_s"]).item()),
            provenance=json.loads(str(data["provenance_json"].item())),
        )


def velocity_contributions(
    positions_m: np.ndarray,
    total_forces_N: np.ndarray,
    *,
    bead_radius_m: float,
    viscosity_Pa_s: float,
    source_mask: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Decompose bead velocities into total, self, and other-source RPY terms."""

    positions = np.asarray(positions_m, dtype=float)
    forces = np.asarray(total_forces_N, dtype=float)
    if (
        positions.shape != forces.shape
        or positions.ndim != 2
        or positions.shape[1] != 3
    ):
        raise ValueError("positions_m and total_forces_N must both have shape (N, 3)")
    mobility = compute_rpy_mobility(positions, bead_radius_m, viscosity_Pa_s)
    source_forces = (
        forces
        if source_mask is None
        else forces * np.asarray(source_mask, dtype=bool)[:, None]
    )
    source_velocity = (mobility @ source_forces.reshape(-1)).reshape(positions.shape)
    self_velocity = source_forces / (6.0 * np.pi * viscosity_Pa_s * bead_radius_m)
    return {
        "source": source_velocity,
        "self": self_velocity,
        "other": source_velocity - self_velocity,
    }


def rpy_flow_velocity(
    observation_points_m: np.ndarray,
    positions_m: np.ndarray,
    total_forces_N: np.ndarray,
    *,
    bead_radius_m: float,
    viscosity_Pa_s: float,
    source_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Reconstruct the RPY flow velocity at arbitrary fixed spatial points."""

    points = np.asarray(observation_points_m, dtype=float)
    positions = np.asarray(positions_m, dtype=float)
    forces = np.asarray(total_forces_N, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("observation_points_m must have shape (M, 3)")
    if (
        positions.shape != forces.shape
        or positions.ndim != 2
        or positions.shape[1] != 3
    ):
        raise ValueError("positions_m and total_forces_N must have shape (N, 3)")
    if source_mask is not None:
        forces = forces * np.asarray(source_mask, dtype=bool)[:, None]
    velocity = np.zeros_like(points)
    for source_position, source_force in zip(positions, forces, strict=True):
        for point_index, point in enumerate(points):
            velocity[point_index] += (
                compute_rpy_pair_mobility(
                    point - source_position, bead_radius_m, viscosity_Pa_s
                )
                @ source_force
            )
    return velocity
