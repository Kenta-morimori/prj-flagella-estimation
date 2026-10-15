"""Reusable post-run bundle geometry and material-motion diagnostics.

Run common development_evaluation first. This report adds no acceptance gate
and never re-simulates archived states.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from sim_swim.analysis.multi_run_campaign import load_yaml
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


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


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
                "travel_body_axis_angle_deg": min(angle, 180 - angle)
                if np.isfinite(angle)
                else np.nan,
                **section_spreads(p, model.flagella_indices),
            }
            for i in range(cfg.flagella.n_flagella):
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
        for lo, hi in [
            (0.0, float(times[-1])),
            (0.0, 0.2),
            (0.2, 0.5),
            (0.5, 0.75),
            (0.75, 1.0),
        ]:
            if hi > times[-1] + 1e-9:
                continue
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
            ]:
                data = np.array([r[key] for r in window])
                item[key + "_median"] = float(np.nanmedian(data))
                item[key + "_p95"] = float(np.nanpercentile(data, 95))
            sp = np.array([r["speed_um_s"] for r in window])
            item["speed_cv"] = float(np.std(sp) / max(np.mean(sp), 1e-30))
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
        provenance.append(
            {
                "condition_id": cid,
                "status": "analyzed",
                "config_overrides": record["config_overrides"],
                "state_count": len(times),
                "rotation_time_match_error_s": offset,
                "rotation_valid_fraction": float(valid.mean()),
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
    _plot(all_rows, output_dir)
    (output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "analysis": "bundle_rotation",
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
