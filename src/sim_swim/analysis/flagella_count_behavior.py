"""Archive helpers for Phase 2.8 flagella-count behavior runs."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np

from sim_swim.sim.core import SimulationState
from sim_swim.sim.params import merge_overrides

ARCHIVE_FORMAT = "sim_swim.flagella_count_behavior_state_archive"
ARCHIVE_VERSION = 1

ANALYSIS_OVERRIDE_ROOTS = {
    "dataset_id",
    "run_batch_id",
    "base_config",
    "feature_schema",
    "sweep",
    "base_overrides",
    "output",
}

SIMULATION_OVERRIDE_ROOTS = {
    "scale",
    "body",
    "flagella",
    "fluid",
    "motor",
    "potentials",
    "hook",
    "body_equiv_load",
    "run_tumble",
    "time",
    "output_sampling",
    "brownian",
    "render",
    "seed",
    "stiffness_scales",
}


def _merge_nested(dst: dict[str, Any], src: dict[str, Any]) -> dict[str, Any]:
    merged = dict(dst)
    for key, value in src.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_nested(merged[key], value)
        else:
            merged[key] = value
    return merged


def _dotted_to_nested(key: str, value: Any) -> dict[str, Any]:
    nested: dict[str, Any] = {}
    node = nested
    parts = key.split(".")
    for part in parts[:-1]:
        child = node.setdefault(part, {})
        if not isinstance(child, dict):
            raise ValueError(f"Override path conflict: {key}")
        node = child
    node[parts[-1]] = value
    return nested


def normalize_base_overrides(raw: dict[str, Any] | None) -> dict[str, Any]:
    """Normalize flat or nested analysis base_overrides to a nested dict."""

    normalized: dict[str, Any] = {}
    for key, value in (raw or {}).items():
        if "." in str(key):
            normalized = _merge_nested(normalized, _dotted_to_nested(str(key), value))
        elif isinstance(value, dict) and isinstance(normalized.get(key), dict):
            normalized[key] = _merge_nested(normalized[key], value)
        else:
            normalized[str(key)] = value
    return normalized


def _coerce_int_list(value: Any) -> list[int]:
    if isinstance(value, str):
        if value.strip() == "":
            return []
        return [int(part.strip()) for part in value.split(",")]
    if isinstance(value, (list, tuple)):
        return [int(part) for part in value]
    return [int(value)]


def _normalize_sweep(raw: dict[str, Any] | None) -> dict[str, Any]:
    sweep = dict(raw or {})
    if "n_flagella" in sweep:
        sweep["n_flagella"] = _coerce_int_list(sweep["n_flagella"])
    if "seeds" in sweep:
        sweep["seeds"] = _coerce_int_list(sweep["seeds"])
    return sweep


def apply_analysis_cli_overrides(
    analysis_config: dict[str, Any],
    cli_overrides: list[str] | None,
) -> dict[str, Any]:
    """Apply KEY=VALUE CLI overrides to an analysis config.

    In scripts.analysis, output.* targets the analysis output paths.
    Simulation shorthand such as time.* and motor.* is routed to base_overrides.
    """

    if not cli_overrides:
        effective = dict(analysis_config)
        effective["base_overrides"] = normalize_base_overrides(
            effective.get("base_overrides", {}) or {}
        )
        effective["sweep"] = _normalize_sweep(effective.get("sweep", {}) or {})
        return effective

    analysis_items: list[str] = []
    simulation_items: list[str] = []
    for raw in cli_overrides:
        if "=" not in raw:
            raise ValueError(f"Invalid override; expected KEY=VALUE: {raw}")
        key = raw.split("=", 1)[0].strip()
        root = key.split(".", 1)[0]
        if root in ANALYSIS_OVERRIDE_ROOTS:
            analysis_items.append(raw)
        elif root in SIMULATION_OVERRIDE_ROOTS:
            simulation_items.append(raw)
        else:
            raise ValueError(f"Unknown analysis override root: {root}")

    effective = dict(analysis_config)
    effective["base_overrides"] = normalize_base_overrides(
        effective.get("base_overrides", {}) or {}
    )

    analysis_nested = merge_overrides({}, analysis_items)
    explicit_base_overrides = normalize_base_overrides(
        analysis_nested.pop("base_overrides", {}) or {}
    )
    effective = _merge_nested(effective, analysis_nested)
    if explicit_base_overrides:
        effective["base_overrides"] = _merge_nested(
            normalize_base_overrides(effective.get("base_overrides", {}) or {}),
            explicit_base_overrides,
        )

    simulation_nested = merge_overrides({}, simulation_items)
    if simulation_nested:
        effective["base_overrides"] = _merge_nested(
            normalize_base_overrides(effective.get("base_overrides", {}) or {}),
            simulation_nested,
        )
    effective["sweep"] = _normalize_sweep(effective.get("sweep", {}) or {})

    return effective


def _stack_vectors(
    values: list[tuple[float, ...]],
    *,
    dtype: np.dtype[Any] | type[Any],
    width: int,
) -> np.ndarray:
    if not values:
        return np.zeros((0, width), dtype=dtype)
    return np.asarray(values, dtype=dtype)


def _stack_same_shape_arrays(values: list[np.ndarray]) -> np.ndarray:
    if not values:
        return np.zeros((0, 0, 3), dtype=float)
    first_shape = values[0].shape
    for value in values[1:]:
        if value.shape != first_shape:
            raise ValueError(
                "bead_positions_um shapes differ across states; cannot archive"
            )
    return np.stack(values, axis=0)


def _stack_int_sequences(values: list[tuple[int, ...]]) -> np.ndarray:
    if not values:
        return np.zeros((0, 0), dtype=np.int64)
    width = len(values[0])
    arr = np.empty((len(values), width), dtype=np.int64)
    for idx, seq in enumerate(values):
        if len(seq) != width:
            raise ValueError("Sequence lengths differ across states; cannot archive")
        arr[idx] = np.asarray(seq, dtype=np.int64)
    return arr


def save_state_archive(path: Path, states: list[SimulationState]) -> None:
    """Persist SimulationState rows for later replay/rendering."""

    path.parent.mkdir(parents=True, exist_ok=True)

    t = np.asarray([st.t for st in states], dtype=float)
    position_um = _stack_vectors(
        [st.position_um for st in states], dtype=float, width=3
    )
    quaternion = _stack_vectors([st.quaternion for st in states], dtype=float, width=4)
    velocity_um_s = _stack_vectors(
        [st.velocity_um_s for st in states], dtype=float, width=3
    )
    omega_rad_s = _stack_vectors(
        [st.omega_rad_s for st in states], dtype=float, width=3
    )
    bead_positions_um = _stack_same_shape_arrays(
        [np.asarray(st.bead_positions_um, dtype=float) for st in states]
    )
    flag_states = _stack_int_sequences([tuple(st.flag_states) for st in states])
    reverse_flagella = _stack_int_sequences(
        [tuple(st.reverse_flagella) for st in states]
    )

    np.savez_compressed(
        path,
        archive_format=np.asarray(ARCHIVE_FORMAT),
        archive_version=np.asarray(ARCHIVE_VERSION, dtype=np.int64),
        state_count=np.asarray(len(states), dtype=np.int64),
        t=t,
        position_um=position_um,
        quaternion=quaternion,
        velocity_um_s=velocity_um_s,
        omega_rad_s=omega_rad_s,
        bead_positions_um=bead_positions_um,
        flag_states=flag_states,
        reverse_flagella=reverse_flagella,
    )


def load_state_archive(path: Path) -> list[SimulationState]:
    """Load archived SimulationState rows for replay/rendering."""

    with np.load(path, allow_pickle=False) as data:
        archive_format = str(data["archive_format"].item())
        if archive_format != ARCHIVE_FORMAT:
            raise ValueError(f"Unsupported archive format: {archive_format}")
        archive_version = int(np.asarray(data["archive_version"]).item())
        if archive_version != ARCHIVE_VERSION:
            raise ValueError(f"Unsupported archive version: {archive_version}")

        t = np.asarray(data["t"], dtype=float)
        position_um = np.asarray(data["position_um"], dtype=float)
        quaternion = np.asarray(data["quaternion"], dtype=float)
        velocity_um_s = np.asarray(data["velocity_um_s"], dtype=float)
        omega_rad_s = np.asarray(data["omega_rad_s"], dtype=float)
        bead_positions_um = np.asarray(data["bead_positions_um"], dtype=float)
        flag_states = np.asarray(data["flag_states"], dtype=np.int64)
        reverse_flagella = np.asarray(data["reverse_flagella"], dtype=np.int64)

    state_count = int(t.shape[0])
    states: list[SimulationState] = []
    for idx in range(state_count):
        states.append(
            SimulationState(
                t=float(t[idx]),
                position_um=tuple(float(v) for v in position_um[idx]),
                quaternion=tuple(float(v) for v in quaternion[idx]),
                velocity_um_s=tuple(float(v) for v in velocity_um_s[idx]),
                omega_rad_s=tuple(float(v) for v in omega_rad_s[idx]),
                bead_positions_um=np.asarray(bead_positions_um[idx], dtype=float),
                flag_states=tuple(int(v) for v in flag_states[idx])
                if flag_states.ndim == 2
                else (),
                reverse_flagella=tuple(int(v) for v in reverse_flagella[idx])
                if reverse_flagella.ndim == 2
                else (),
            )
        )
    return states


def write_trajectory_csv(path: Path, states: list[SimulationState]) -> None:
    """Write a human-readable trajectory CSV for inspection."""

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "t",
        "x",
        "y",
        "z",
        "qx",
        "qy",
        "qz",
        "qw",
        "vx",
        "vy",
        "vz",
        "wx",
        "wy",
        "wz",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for st in states:
            writer.writerow(
                {
                    "t": st.t,
                    "x": st.position_um[0],
                    "y": st.position_um[1],
                    "z": st.position_um[2],
                    "qx": st.quaternion[0],
                    "qy": st.quaternion[1],
                    "qz": st.quaternion[2],
                    "qw": st.quaternion[3],
                    "vx": st.velocity_um_s[0],
                    "vy": st.velocity_um_s[1],
                    "vz": st.velocity_um_s[2],
                    "wx": st.omega_rad_s[0],
                    "wy": st.omega_rad_s[1],
                    "wz": st.omega_rad_s[2],
                }
            )


def archive_metadata(*, sample_id: str, config_path: str) -> dict[str, Any]:
    """Return a compact metadata block for raw output manifests."""

    return {
        "sample_id": sample_id,
        "config_path": config_path,
        "state_archive_format": ARCHIVE_FORMAT,
        "state_archive_version": ARCHIVE_VERSION,
    }
