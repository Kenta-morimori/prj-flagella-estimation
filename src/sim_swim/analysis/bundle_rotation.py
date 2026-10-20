"""Reusable post-run bundle geometry and material-motion diagnostics.

Run common development_evaluation first. This report adds no acceptance gate
and never re-simulates archived states.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.dynamics.forces import _principal_axis_or_none
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.debug_summary import (
    _bundle_metrics,
    _flag_helix_bundle_radius_stats_um,
)
from sim_swim.sim.helix_axis import (
    angle_deg_between,
    estimate_body_axis,
    estimate_flag_helix_axis,
    helix_axis_alignment_metrics,
)
from sim_swim.sim.params import SimulationConfig


def section_spreads(points: np.ndarray, indices: list[np.ndarray]) -> dict[str, float]:
    """Cross-filament distances at matched normalized chain contour positions.

    These include root offsets and helical radius; they do not define touching
    or a topological bundle. No world-axis projection is used.
    """
    values = {}
    for fraction in (0.25, 0.5, 0.75, 1.0):
        samples = []
        for ids in indices:
            chain = points[ids[1:]]
            arc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(chain, axis=0), axis=1))]
            samples.append(
                [np.interp(fraction * arc[-1], arc, chain[:, k]) for k in range(3)]
            )
        samples = np.asarray(samples)
        radius = np.linalg.norm(samples - samples.mean(0), axis=1)
        values[f"contour_{fraction:.2f}_spread_um"] = float(
            np.sqrt(np.mean(radius**2)) * 1e6
        )
    return values


def motor_root_lateral_arms_um(
    points: np.ndarray, indices: list[np.ndarray]
) -> np.ndarray:
    """Saved-state first segment arms relative to the motor's fitted flag axis."""
    arms = np.full(len(indices), np.nan)
    for flag_id, ids in enumerate(indices):
        chain = points[ids]
        axis = _principal_axis_or_none(chain)
        if axis is None or len(chain) < 2:
            continue
        segment = chain[0] - chain[1]
        arms[flag_id] = np.linalg.norm(segment - np.dot(segment, axis) * axis) * 1e6
    return arms


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        # Per-flag fields differ across counts in a multi-condition campaign.
        fields = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def observation_windows(end_s: float) -> list[tuple[float, float]]:
    """Keep historical windows and cover the actual end of a longer archive."""
    end = round(end_s, 9)
    candidates = [
        (0.0, end),
        (0.0, 0.2),
        (0.2, 0.5),
        (0.5, 0.75),
        (0.75, 1.0),
        (max(0.0, end - 0.5), end),
        (max(0.0, end - 0.25), end),
    ]
    return list(dict.fromkeys((lo, hi) for lo, hi in candidates if 0 <= lo < hi <= end))


def completed_rotation_cycles(rows: list[dict], n_flags: int) -> tuple[list[dict], str]:
    """Compare full signed best-fit-spin turns, without asserting stationarity.

    Flag 0 is the reference clock. Reversals, missing phase or an invalid flag
    prevent a cycle report. No phase unwrapping, gap filling or partial cycles.
    """
    phase = np.array([r["flag0_phase_rad"] for r in rows], dtype=float)
    times = np.array([r["t_s"] for r in rows], dtype=float)
    if (
        not all(r[f"flag{i}_valid"] for r in rows for i in range(n_flags))
        or not np.isfinite(phase).all()
    ):
        return [], "unavailable: invalid rotation phase"
    direction = np.sign(phase[-1] - phase[0])
    progress = direction * (phase - phase[0])
    # A completed run can append its final state after the last regular archive
    # sample.  The two timestamps differ by roundoff while the phase is identical.
    # Exclude only that near-duplicate from the strictly monotone phase clock;
    # a real plateau or reversal over an archive interval remains unavailable.
    near_duplicate = (
        (np.diff(times) >= 0)
        & (np.diff(times) <= 1e-9)
        & (np.abs(np.diff(progress)) <= 1e-12)
    )
    keep = np.r_[True, ~near_duplicate]
    clock_times = times[keep]
    clock_progress = progress[keep]
    if direction == 0 or np.any(np.diff(clock_progress) <= 0):
        return [], "unavailable: non-monotone rotation"
    if np.max(np.diff(clock_progress)) >= np.pi / 2:
        return [], "unavailable: geometry sampling does not resolve a rotation cycle"
    turns = int(np.floor(clock_progress[-1] / (2 * np.pi)))
    if turns < 1:
        return [], "unavailable: less than one completed turn"
    boundaries = np.interp(
        np.arange(turns + 1) * 2 * np.pi, clock_progress, clock_times
    )
    records = []
    for i, (start, end) in enumerate(zip(boundaries[:-1], boundaries[1:]), 1):
        selected = [r for r in rows if start <= r["t_s"] <= end]
        record = {
            "condition_id": rows[0]["condition_id"],
            "cycle": i,
            "direction": int(direction),
            "start_s": float(start),
            "end_s": float(end),
            "period_s": float(end - start),
            "samples": len(selected),
        }
        for key in (
            "bundle_radius_mean_um",
            "helix_pair_angle_mean_deg",
            "bundle_rear_angle_deg",
            "contour_0.75_spread_um",
            "contour_1.00_spread_um",
            "speed_um_s",
            "forward_speed_um_s",
        ):
            if not all(key in row for row in rows):
                continue
            values = np.array([r[key] for r in rows])
            first, last = np.interp([start, end], times, values)
            record[key + "_start"] = float(first)
            record[key + "_end"] = float(last)
            record[key + "_boundary_change"] = float(last - first)
            record[key + "_mean"] = (
                float(np.mean([r[key] for r in selected])) if selected else float("nan")
            )
        if all("forward_speed_um_s" in row for row in rows):
            sample_times = np.r_[start, times[(times > start) & (times < end)], end]
            forward_speed = np.interp(
                sample_times,
                times,
                [r["forward_speed_um_s"] for r in rows],
            )
            record["forward_displacement_um"] = float(
                np.sum(
                    0.5
                    * (forward_speed[:-1] + forward_speed[1:])
                    * np.diff(sample_times)
                )
            )
        for flag in range(n_flags):
            values = np.array([r[f"flag{flag}_phase_rad"] for r in rows])
            first, last = np.interp([start, end], times, values)
            record[f"flag{flag}_net_turns"] = float((last - first) / (2 * np.pi))
        records.append(record)
    return records, "available: complete best-fit-spin turns; no stationarity gate"


def analyze(run_dir: Path, evaluation_dir: Path, output_dir: Path):
    # Compact summaries and common QC precede archives.
    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    qc = {
        r["condition_id"]: r
        for r in csv.DictReader((evaluation_dir / "summary.csv").open())
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    baseline = SimulationConfig.from_dict(load_yaml(Path(manifest["base_config"])))
    all_rows = []
    windows = []
    cycles = []
    provenance = []
    for record in manifest["conditions"]:
        cid = record["condition_id"]
        path = run_dir / "conditions" / cid
        summary = json.loads((path / "run_summary.json").read_text())
        if (
            qc[cid]["screen_status"] != "pass"
            or summary["execution"]["status"] != "completed"
        ):
            provenance.append(
                {"condition_id": cid, "status": "excluded: common QC/completion"}
            )
            continue
        cfg = baseline.with_overrides(record["config_overrides"])
        model = ModelBuilder(cfg).build()
        with np.load(path / "state_archive.npz", allow_pickle=False) as archive:
            times = archive["t"]
            positions = archive["bead_positions_um"] * 1e-6
        rotation = list(csv.DictReader((path / "rotation_diagnostics.csv").open()))
        rt = np.array([float(r["t_s"]) for r in rotation])
        selection = np.array([np.argmin(abs(rt - t)) for t in times])
        offset = float(np.max(abs(rt[selection] - times)))
        if offset > cfg.tau_s * cfg.dt_star + 1e-10:
            raise ValueError("rotation/archive timestamps do not match")
        phases = np.array(
            [
                [
                    float(rotation[j][f"flag{i}_phase_rad"])
                    for i in range(cfg.flagella.n_flagella)
                ]
                for j in selection
            ]
        )
        valid = np.array(
            [
                [
                    bool(int(rotation[j][f"flag{i}_valid"]))
                    for i in range(cfg.flagella.n_flagella)
                ]
                for j in selection
            ]
        )
        phases[~valid] = np.nan
        angular_speed = np.gradient(phases, times, axis=0)
        centers = positions[:, model.body_indices].mean(1)
        velocity = np.gradient(centers, times, axis=0)
        speed = np.linalg.norm(velocity, axis=1) * 1e6
        floor = max(float(speed.max()) * 1e-8, 1e-12)
        rows = []
        for j, (t, p) in enumerate(zip(times, positions)):
            body = estimate_body_axis(p, model.body_layer_indices, model.body_indices)
            helix = [
                estimate_flag_helix_axis(p, ids, i)
                for i, ids in enumerate(model.flagella_indices)
            ]
            alignment = helix_axis_alignment_metrics([h.axis for h in helix])
            bundle = _bundle_metrics(p, model, cfg, body.axis, body.rear_direction)
            radius = _flag_helix_bundle_radius_stats_um(p, model.flagella_indices)[0]
            root_arms = motor_root_lateral_arms_um(p, model.flagella_indices)
            finite_root_arms = root_arms[np.isfinite(root_arms)]
            angle = (
                angle_deg_between(velocity[j], body.axis)
                if speed[j] > floor
                else np.nan
            )
            row = {
                "condition_id": cid,
                "t_s": float(t),
                "bundle_radius_mean_um": radius,
                "bundle_rear_angle_deg": bundle["bundle_axis_vs_rear_angle_deg"],
                "helix_pair_angle_mean_deg": alignment.pair_angle_deg_mean,
                "helix_pair_angle_max_deg": alignment.pair_angle_deg_max,
                "speed_um_s": float(speed[j]),
                "forward_speed_um_s": float(np.dot(velocity[j], body.axis) * 1e6),
                "motor_root_min_lateral_arm_um": (
                    float(finite_root_arms.min()) if finite_root_arms.size else np.nan
                ),
                "travel_body_axis_angle_deg": min(angle, 180 - angle)
                if np.isfinite(angle)
                else np.nan,
                **section_spreads(p, model.flagella_indices),
            }
            for i in range(cfg.flagella.n_flagella):
                row[f"flag{i}_motor_root_lateral_arm_um"] = float(root_arms[i])
                row[f"flag{i}_phase_rad"] = float(phases[j, i])
                row[f"flag{i}_omega_rad_s"] = float(angular_speed[j, i])
                row[f"flag{i}_valid"] = bool(valid[j, i])
                row[f"flag{i}_fit_rms_max_m"] = float(
                    rotation[selection[j]][f"flag{i}_fit_rms_max_m"]
                )
                if i:
                    row[f"flag{i}_minus_flag0_phase_drift_rad"] = float(
                        phases[j, i] - phases[j, 0]
                    )
            rows.append(row)
        all_rows.extend(rows)
        for lo, hi in observation_windows(float(times[-1])):
            window = [r for r in rows if lo - 1e-9 <= r["t_s"] <= hi + 1e-9]
            item = {
                "condition_id": cid,
                "start_s": lo,
                "end_s": hi,
                "samples": len(window),
                "rotation_valid": all(
                    r[f"flag{i}_valid"]
                    for r in window
                    for i in range(cfg.flagella.n_flagella)
                ),
            }
            for key in [
                "bundle_radius_mean_um",
                "bundle_rear_angle_deg",
                "helix_pair_angle_mean_deg",
                "helix_pair_angle_max_deg",
                "travel_body_axis_angle_deg",
                "speed_um_s",
                "forward_speed_um_s",
            ]:
                data = np.array([r[key] for r in window])
                item[key + "_median"] = float(np.nanmedian(data))
                item[key + "_p95"] = float(np.nanpercentile(data, 95))
            sp = np.array([r["speed_um_s"] for r in window])
            item["speed_cv"] = float(np.std(sp) / max(np.mean(sp), 1e-30))
            item["forward_sample_fraction"] = float(
                np.mean([r["forward_speed_um_s"] > 0 for r in window])
            )
            root_arms = np.array([r["motor_root_min_lateral_arm_um"] for r in window])
            finite_root_arms = root_arms[np.isfinite(root_arms)]
            item["motor_root_min_lateral_arm_um_min"] = (
                float(finite_root_arms.min()) if finite_root_arms.size else np.nan
            )
            item["motor_root_min_lateral_arm_um_median"] = (
                float(np.median(finite_root_arms)) if finite_root_arms.size else np.nan
            )
            for i in range(cfg.flagella.n_flagella):
                phase = np.array([r[f"flag{i}_phase_rad"] for r in window])
                omega = np.array([r[f"flag{i}_omega_rad_s"] for r in window])
                item[f"flag{i}_net_turns"] = float((phase[-1] - phase[0]) / (2 * np.pi))
                item[f"flag{i}_omega_mean_rad_s"] = float(np.mean(omega))
                item[f"flag{i}_omega_std_rad_s"] = float(np.std(omega))
                if i:
                    relative = np.array(
                        [r[f"flag{i}_minus_flag0_phase_drift_rad"] for r in window]
                    )
                    item[f"flag{i}_phase_drift_range_rad"] = float(np.ptp(relative))
            windows.append(item)
        condition_cycles, cycle_status = completed_rotation_cycles(
            rows, cfg.flagella.n_flagella
        )
        cycles.extend(condition_cycles)
        provenance.append(
            {
                "condition_id": cid,
                "status": "analyzed",
                "config_overrides": record["config_overrides"],
                "state_count": len(times),
                "rotation_time_match_error_s": offset,
                "rotation_valid_fraction": float(valid.mean()),
                "cycle_report_status": cycle_status,
                "completed_reference_turns": len(condition_cycles),
                "whole_run_net_turns": (
                    (phases[-1] - phases[0]) / (2 * np.pi)
                ).tolist(),
            }
        )
    if not all_rows:
        raise ValueError("No completed QC-pass input")
    write_csv(output_dir / "bundle_rotation_timeseries.csv", all_rows)
    if windows:
        write_csv(output_dir / "bundle_rotation_windows.csv", windows)
    if cycles:
        write_csv(output_dir / "bundle_rotation_cycles.csv", cycles)
    else:
        (output_dir / "bundle_rotation_cycles.csv").unlink(missing_ok=True)
    _plot(all_rows, output_dir)
    (output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "analysis": "bundle_rotation",
                "analyzer_source_sha256": hashlib.sha256(
                    Path(__file__).read_bytes()
                ).hexdigest(),
                "run_manifest_sha256": hashlib.sha256(
                    (run_dir / "run_manifest.json").read_bytes()
                ).hexdigest(),
                "qc_summary_sha256": hashlib.sha256(
                    (evaluation_dir / "summary.csv").read_bytes()
                ).hexdigest(),
                "source_run": str(run_dir.resolve()),
                "qc_source": str(evaluation_dir.resolve()),
                "acceptance_gate": False,
                "phase_definition": "integrated internal-step best-fit axial spin relative to body; accumulated phase drift, not fitted absolute initial phase",
                "conditions": provenance,
            },
            indent=2,
        )
        + "\n"
    )


def _plot(rows, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ids = list(dict.fromkeys(r["condition_id"] for r in rows))
    figure, axes = plt.subplots(4, len(ids), figsize=(4 * len(ids), 12), squeeze=False)
    for col, cid in enumerate(ids):
        subset = [r for r in rows if r["condition_id"] == cid]
        t = [r["t_s"] for r in subset]
        for ax, key in zip(
            axes[:3, col],
            [
                "bundle_radius_mean_um",
                "helix_pair_angle_mean_deg",
                "bundle_rear_angle_deg",
            ],
        ):
            ax.plot(t, [r[key] for r in subset])
            ax.set_ylabel(key)
            ax.grid(alpha=0.3)
        for key in subset[0]:
            if key.startswith("flag") and key.endswith("_phase_rad"):
                axes[3, col].plot(t, [r[key] / (2 * np.pi) for r in subset], label=key)
        axes[0, col].set_title(cid)
        axes[3, col].set_ylabel("accumulated relative-body turns")
        axes[3, col].legend(fontsize=7)
        axes[3, col].set_xlabel("t [s]")
    figure.tight_layout()
    figure.savefig(output / "bundle_rotation.png", dpi=140)
    plt.close(figure)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--evaluation-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    analyze(args.run_dir, args.evaluation_dir, args.output_dir)


if __name__ == "__main__":
    main()
