"""Offline, same-state motor torque audit of completed Phase 2 archives.

This is a force diagnostic, not a simulation or an acceptance gate.  A
counterfactual body reaction is reported only after the archived-state
reconstruction agrees with the recorded diagnostic sample.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from sim_swim.analysis.phase2_replay import _build_cfg
from sim_swim.analysis.torque_weight_replay import reconstructed_segment_weights
from sim_swim.dynamics.forces import (
    _principal_axis_or_none,
    compute_root_torque_segment_couples_forces,
)
from sim_swim.sim.core import Simulator


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_samples(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No diagnostic samples in {path}")
    return rows


def _select_samples(
    rows: list[dict[str, str]], threshold: float
) -> dict[str, dict[str, str]]:
    def ratio(row: dict[str, str]) -> float:
        return float(row["diagnostic_motor_torque_balance_residual_ratio"])

    above = next((row for row in rows if ratio(row) > threshold), None)
    selected = {"sampled_maximum": max(rows, key=ratio), "final_observed": rows[-1]}
    if above is not None:
        selected["first_observed_above_threshold"] = above
    return selected


def _balance(
    positions_m: np.ndarray,
    forces: np.ndarray,
    body_indices: np.ndarray,
    flagella_indices: list[np.ndarray],
) -> dict[str, Any]:
    origin = positions_m.mean(axis=0)
    flag_indices = np.concatenate(flagella_indices)
    body_force = forces[body_indices].sum(axis=0)
    flag_force = forces[flag_indices].sum(axis=0)
    body_torque = np.cross(
        positions_m[body_indices] - origin, forces[body_indices]
    ).sum(axis=0)
    flag_torque = np.cross(
        positions_m[flag_indices] - origin, forces[flag_indices]
    ).sum(axis=0)
    torque_scale = np.linalg.norm(np.cross(positions_m - origin, forces), axis=1).sum()
    force_scale = np.linalg.norm(forces, axis=1).sum()
    return {
        "body_force_N": body_force.tolist(),
        "flag_force_N": flag_force.tolist(),
        "global_force_N": (body_force + flag_force).tolist(),
        "body_torque_Nm": body_torque.tolist(),
        "flag_torque_Nm": flag_torque.tolist(),
        "global_torque_Nm": (body_torque + flag_torque).tolist(),
        "global_force_residual_ratio": float(
            np.linalg.norm(body_force + flag_force) / max(force_scale, 1e-30)
        ),
        "global_torque_residual_ratio": float(
            np.linalg.norm(body_torque + flag_torque) / max(torque_scale, 1e-30)
        ),
        "torque_scale_Nm": float(torque_scale),
    }


def _per_flag(
    positions_m: np.ndarray,
    body_indices: np.ndarray,
    flagella_indices: list[np.ndarray],
    weights: list[np.ndarray],
    torque_Nm: float,
    full_vector_body_reaction: bool,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    origin = positions_m.mean(axis=0)
    total = np.zeros_like(positions_m)
    rows: list[dict[str, Any]] = []
    for flag_id, flag_indices in enumerate(flagella_indices):
        force, diagnostic = compute_root_torque_segment_couples_forces(
            positions_m,
            [flag_indices],
            body_indices,
            np.asarray([torque_Nm]),
            [weights[flag_id]],
            full_vector_body_reaction=full_vector_body_reaction,
        )
        if diagnostic.degenerate_axis_count:
            raise ValueError(f"Degenerate motor force in F{flag_id}")
        total += force
        axis = _principal_axis_or_none(positions_m[flag_indices])
        if axis is None:
            raise ValueError(f"Degenerate flagellar axis in F{flag_id}")
        if (
            float(
                np.dot(
                    axis, positions_m[flag_indices[-1]] - positions_m[flag_indices[0]]
                )
            )
            < 0
        ):
            axis = -axis
        flag_torque = np.cross(
            positions_m[flag_indices] - origin, force[flag_indices]
        ).sum(axis=0)
        body_torque = np.cross(
            positions_m[body_indices] - origin, force[body_indices]
        ).sum(axis=0)
        axial = float(np.dot(flag_torque, axis))
        transverse = flag_torque - axial * axis
        rows.append(
            {
                "flag_id": flag_id,
                "axis": axis.tolist(),
                "flag_torque_Nm": flag_torque.tolist(),
                "flag_axial_torque_Nm": axial,
                "flag_transverse_torque_Nm": transverse.tolist(),
                "flag_transverse_torque_norm_Nm": float(np.linalg.norm(transverse)),
                "body_reaction_torque_Nm": body_torque.tolist(),
                "body_reaction_axial_torque_Nm": float(np.dot(body_torque, axis)),
                "body_reaction_transverse_torque_Nm": (
                    body_torque - np.dot(body_torque, axis) * axis
                ).tolist(),
                "net_torque_Nm": (body_torque + flag_torque).tolist(),
            }
        )
    return total, rows


def _require_whole_force_agreement(
    positions_m: np.ndarray,
    body_indices: np.ndarray,
    flagella_indices: list[np.ndarray],
    weights: list[np.ndarray],
    torque_Nm: float,
    full_vector_body_reaction: bool,
    summed_per_flag_forces: np.ndarray,
) -> None:
    whole, _ = compute_root_torque_segment_couples_forces(
        positions_m,
        flagella_indices,
        body_indices,
        np.full(len(flagella_indices), torque_Nm),
        weights,
        full_vector_body_reaction=full_vector_body_reaction,
    )
    relative_error = float(
        np.linalg.norm(whole - summed_per_flag_forces)
        / max(float(np.linalg.norm(whole)), 1e-30)
    )
    if relative_error > 1e-10:
        raise ValueError(
            f"Per-flag motor forces do not reconstruct whole motor force: {relative_error}"
        )


def _observed_match(
    calculated: dict[str, Any], observed: dict[str, str], torque_Nm: float
) -> dict[str, Any]:
    errors: dict[str, float] = {}
    for side in ("body", "flag"):
        key = f"{side}_torque_Nm"
        stored = np.asarray(
            [
                float(observed[f"diagnostic_motor_net_torque_{side}_{axis}_Nm"])
                for axis in "xyz"
            ]
        )
        measured = np.asarray(calculated[key])
        errors[side] = float(
            np.linalg.norm(measured - stored)
            / max(abs(torque_Nm), float(np.linalg.norm(stored)), 1e-30)
        )
    ratio_error = abs(
        calculated["global_torque_residual_ratio"]
        - float(observed["diagnostic_motor_torque_balance_residual_ratio"])
    )
    force_ratio_error = abs(
        calculated["global_force_residual_ratio"]
        - float(observed["diagnostic_motor_force_balance_residual_ratio"])
    )
    matched = (
        max(errors.values()) <= 0.005
        and ratio_error <= 0.005
        and force_ratio_error <= 0.005
    )
    return {
        "matched": matched,
        "relative_torque_vector_errors": errors,
        "torque_ratio_absolute_error": ratio_error,
        "force_ratio_absolute_error": force_ratio_error,
        "recorded_torque_residual_ratio": float(
            observed["diagnostic_motor_torque_balance_residual_ratio"]
        ),
    }


def audit_completed_archives(
    evaluation_dir: Path, output_dir: Path, *, threshold: float = 0.02
) -> dict[str, Any]:
    replay_input = evaluation_dir / "replay_input"
    replay_manifest_path = replay_input / "run_manifest.json"
    replay_manifest = json.loads(replay_manifest_path.read_text(encoding="utf-8"))
    with (evaluation_dir / "summary.csv").open(encoding="utf-8", newline="") as handle:
        summary_rows = {row["condition_id"]: row for row in csv.DictReader(handle)}
    prepared = []
    for record in replay_manifest["conditions"]:
        condition_id = record["condition_id"]
        source_dir = Path(record["output_dir"])
        archive = source_dir / "state_archive.npz"
        summary = source_dir / "run_summary.json"
        samples_path = source_dir / "diagnostic_samples.csv"
        for path in (archive, summary, samples_path):
            if not path.is_file():
                raise FileNotFoundError(path)
        if _sha256(archive) != summary_rows[condition_id]["state_archive_sha256"]:
            raise ValueError(f"SHA-256 mismatch for {condition_id}: state_archive.npz")
        run_summary = json.loads(summary.read_text(encoding="utf-8"))
        if run_summary.get("execution", {}).get("status") != "completed":
            raise ValueError(
                f"{condition_id}: partial/incomplete source is not permitted"
            )
        samples = _read_samples(samples_path)
        selected = _select_samples(samples, threshold)
        with np.load(archive, allow_pickle=False) as data:
            times = np.asarray(data["t"], dtype=float)
            positions = np.asarray(data["bead_positions_um"], dtype=float) * 1e-6
        if len(times) == 0 or len(times) != len(positions):
            raise ValueError(f"{condition_id}: invalid completed archive")
        cfg = _build_cfg(
            base_cfg_path=Path(replay_manifest["base_config"]),
            condition_record=record,
            fps_out_3d=25,
        )
        if (
            cfg.motor.force_distribution != "root_torque_segment_couples"
            or cfg.motor.enable_switching
            or cfg.motor.torque_ramp_enabled
            or cfg.motor.body_reaction_full_vector
        ):
            raise ValueError(f"{condition_id}: unsupported motor reconstruction policy")
        simulator = Simulator(cfg)
        body_indices = np.flatnonzero(simulator.model.bead_is_body)
        flagella_indices = simulator.model.flagella_indices
        if positions.shape[1] != simulator.model.positions_m.shape[0]:
            raise ValueError(f"{condition_id}: geometry/archive bead count mismatch")
        chosen = {"initial": (0, None)}
        for name, sample in selected.items():
            sample_t = float(sample["diagnostic_t_s"])
            archive_index = int(np.argmin(np.abs(times - sample_t)))
            if abs(times[archive_index] - sample_t) > max(cfg.dt_s * 2, 1e-6):
                raise ValueError(
                    f"{condition_id}: no archive state near diagnostic {sample_t}"
                )
            chosen[name] = (archive_index, sample)
        prepared.append(
            (
                record,
                archive,
                samples_path,
                cfg,
                body_indices,
                flagella_indices,
                times,
                positions,
                chosen,
            )
        )

    # The local-twist field is deterministic under RUN/no-ramp/no-switching.
    # Reconstruct only selected archive times, once for identical motor configs.
    weight_cache: dict[tuple[Any, ...], dict[int, np.ndarray]] = {}
    for (
        _record,
        _archive,
        _samples,
        cfg,
        _body,
        flags,
        times,
        _positions,
        chosen,
    ) in prepared:
        key = (
            cfg.motor.torque_distribution_profile,
            len(flags[0]) - 1,
            cfg.dt_s,
            cfg.motor_torque_Nm,
        )
        weight_cache.setdefault(key, {})
        for archive_index, _ in chosen.values():
            weight_cache[key][int(round(times[archive_index] / cfg.dt_s))] = np.empty(0)
    for key, requested in weight_cache.items():
        profile, segment_count, dt_s, torque_Nm = key
        steps = sorted(requested)
        values = reconstructed_segment_weights(
            profile,
            segment_count,
            times_s=np.asarray(steps, dtype=float) * dt_s,
            dt_s=dt_s,
            torque_Nm=torque_Nm,
        )
        requested.update(zip(steps, values, strict=True))

    results: list[dict[str, Any]] = []
    for (
        record,
        archive,
        samples_path,
        cfg,
        body,
        flags,
        times,
        positions,
        chosen,
    ) in prepared:
        key = (
            cfg.motor.torque_distribution_profile,
            len(flags[0]) - 1,
            cfg.dt_s,
            cfg.motor_torque_Nm,
        )
        state_results: dict[str, Any] = {}
        for name, (archive_index, observed) in chosen.items():
            pos = positions[archive_index]
            weight = weight_cache[key][int(round(times[archive_index] / cfg.dt_s))]
            weights = [weight] * len(flags)
            current_forces, per_flag_current = _per_flag(
                pos, body, flags, weights, cfg.motor_torque_Nm, False
            )
            _require_whole_force_agreement(
                pos,
                body,
                flags,
                weights,
                cfg.motor_torque_Nm,
                False,
                current_forces,
            )
            calculated = _balance(pos, current_forces, body, flags)
            verification = (
                None
                if observed is None
                else _observed_match(
                    calculated, observed, cfg.motor_torque_Nm * len(flags)
                )
            )
            item: dict[str, Any] = {
                "archive_index": archive_index,
                "archive_t_s": float(times[archive_index]),
                "diagnostic_t_s": None
                if observed is None
                else float(observed["diagnostic_t_s"]),
                "recorded_match": verification,
                "reconstructed_current_axial_only": {
                    **calculated,
                    "per_flag": per_flag_current,
                },
                "full_vector_same_state": None,
                "interpretation": "initial state has no recorded motor diagnostic; counterfactual not validated"
                if observed is None
                else "recorded mismatch; counterfactual withheld",
            }
            if verification is not None and verification["matched"]:
                full_forces, per_flag_full = _per_flag(
                    pos, body, flags, weights, cfg.motor_torque_Nm, True
                )
                _require_whole_force_agreement(
                    pos,
                    body,
                    flags,
                    weights,
                    cfg.motor_torque_Nm,
                    True,
                    full_forces,
                )
                item["full_vector_same_state"] = {
                    **_balance(pos, full_forces, body, flags),
                    "per_flag": per_flag_full,
                }
                item["interpretation"] = (
                    "same archived geometry and torque weights; not a trajectory or stability prediction"
                )
            state_results[name] = item
        results.append(
            {
                "condition_id": record["condition_id"],
                "source_archive": str(archive.resolve()),
                "source_archive_sha256": _sha256(archive),
                "source_diagnostic_samples": str(samples_path.resolve()),
                "source_diagnostic_samples_sha256": _sha256(samples_path),
                "states": state_results,
            }
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "kind": "motor_torque_completed_archive_audit",
        "diagnostic_only": True,
        "threshold_unchanged": threshold,
        "force_model_unchanged": True,
        "replay_manifest": str(replay_manifest_path.resolve()),
        "replay_manifest_sha256": _sha256(replay_manifest_path),
        "condition_count": len(results),
        "conditions": results,
        "comparison_policy": "full_vector_body_reaction is compared only on the same saved state after matching recorded current-force diagnostics; no rerun or stability inference",
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    matched = sum(
        item["recorded_match"] is not None and item["recorded_match"]["matched"]
        for row in results
        for item in row["states"].values()
    )
    compared = sum(
        item["recorded_match"] is not None
        for row in results
        for item in row["states"].values()
    )
    (output_dir / "run.log").write_text(
        f"completed diagnostic-only motor torque audit\nconditions={len(results)}\nrecorded_matches={matched}/{compared}\n",
        encoding="utf-8",
    )
    return manifest


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    result = audit_completed_archives(args.evaluation_dir, args.output_dir)
    print(
        f"Audited {result['condition_count']} completed conditions: {args.output_dir}"
    )


if __name__ == "__main__":
    main()
