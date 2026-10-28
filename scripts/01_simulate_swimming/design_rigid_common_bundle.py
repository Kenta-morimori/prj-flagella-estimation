"""Freeze the reproducible 2010 hex whole-chain rigid bundle coordinates."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import yaml

from sim_swim.model.builder import ModelBuilder
from sim_swim.model.initial_bundle import initial_bundle_force_metrics
from sim_swim.model.rigid_common_bundle import (
    MODE,
    PATTERNS,
    PERTURBATION,
    construct_chains,
)
from sim_swim.sim.params import SimulationConfig


def design(base_path: Path) -> dict:
    base = SimulationConfig.from_dict(yaml.safe_load(base_path.read_text()))
    records = []
    body = None
    for slots in PATTERNS:
        cfg = base.with_overrides(
            {
                "flagella": {
                    "n_flagella": len(slots),
                    "placement_mode": "explicit_center_layer_slots",
                    "attachment_slots": list(slots),
                },
                "seed": {"attach_seed": 0, "phase_seed": 0},
            }
        )
        model = ModelBuilder(cfg).build()
        this_body = model.positions_m[model.body_indices] / cfg.b_m
        if body is None:
            body = this_body
        else:
            np.testing.assert_allclose(this_body, body, rtol=0, atol=1e-12)
        originals = np.asarray(
            [model.positions_m[ids] / cfg.b_m for ids in model.flagella_indices]
        )
        attachments = this_body[model.flagella_attach_body_indices]
        variants = ("base", PERTURBATION) if slots == (0, 2, 4) else ("base",)
        for variant in variants:
            points, recipe = construct_chains(
                originals,
                attachments,
                slots,
                variant=variant,
                hook_length_over_b=cfg.hook.length_over_b,
            )
            records.append(
                {
                    "attachment_slots": list(slots),
                    "variant": variant,
                    "recipe": recipe,
                    "flagella_positions_over_b": points.tolist(),
                }
            )
    return {
        "version": 1,
        "kind": MODE,
        "units": "b",
        "source_profile": str(base_path),
        "source_implementation_commit": "82962a3cd1e8a467e968c7b8b84fb645021a58a4",
        "normal_bend_deg": base.potentials.bend.theta0_deg["normal"],
        "normal_torsion_deg": base.potentials.torsion.phi0_deg["normal"],
        "bead_diameter_over_b": 2 * base.bead_radius_m / base.b_m,
        "hook_length_over_b": base.hook.length_over_b,
        "body_positions_over_b": body.tolist(),
        "conditions": records,
    }


def audit(base_path: Path, state_path: Path) -> dict:
    base = SimulationConfig.from_dict(yaml.safe_load(base_path.read_text()))
    digest = hashlib.sha256(state_path.read_bytes()).hexdigest()
    records = []
    for slots in PATTERNS:
        for variant in ("base", PERTURBATION) if slots == (0, 2, 4) else ("base",):
            cfg = base.with_overrides(
                {
                    "flagella": {
                        "n_flagella": len(slots),
                        "placement_mode": "explicit_center_layer_slots",
                        "attachment_slots": list(slots),
                        "initial_bundle_mode": MODE,
                        "initial_bundle_variant": variant,
                        "initial_bundle_state_path": str(state_path),
                        "initial_bundle_state_sha256": digest,
                    },
                    "seed": {"attach_seed": 0, "phase_seed": 0},
                    "potentials": {
                        "spring_spring_repulsion": {"body_flagella_enabled": False}
                    },
                    "motor": {
                        "body_reaction_full_vector": True,
                        "body_reaction_support": "attach_one_ring",
                        "attach_frame_reaction": "energy_gradient",
                        "local_attach_frame_position_scale": 1.0,
                        "segment_torque_correction": "none",
                        "segment_torque_distribution": "fixed_weights",
                    },
                }
            )
            cfg = cfg.with_overrides(
                {"potentials": {"torsion": {"force_method": "analytic"}}}
            )
            model = ModelBuilder(cfg).build()
            records.append(
                {
                    "n_flagella": len(slots),
                    "attachment_slots": list(slots),
                    "variant": variant,
                    "geometry": model.initial_geometry_metadata,
                    "initial_force_and_energy": initial_bundle_force_metrics(
                        cfg, model
                    ),
                }
            )
    return {
        "source_implementation_commit": "82962a3cd1e8a467e968c7b8b84fb645021a58a4",
        "source_profile": str(base_path),
        "source_profile_sha256": hashlib.sha256(base_path.read_bytes()).hexdigest(),
        "state_sha256": digest,
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base", type=Path, default=Path("conf/sim_swim_2010_hex.yaml")
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path)
    args = parser.parse_args()
    content = json.dumps(design(args.base), indent=2, ensure_ascii=False) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(content, encoding="utf-8")
    if args.audit_dir:
        args.audit_dir.mkdir(parents=True, exist_ok=True)
        audit_path = args.audit_dir / "static_geometry_audit.json"
        audit_path.write_text(
            json.dumps(audit(args.base, args.output), indent=2) + "\n",
            encoding="utf-8",
        )
        (args.audit_dir / "run.log").write_text(
            "Generated and audited 13 base geometries plus one orientation perturbation.\n",
            encoding="utf-8",
        )
        (args.audit_dir / "manifest.json").write_text(
            json.dumps(
                {
                    "kind": "rigid_common_bundle_static_audit",
                    "created_at_jst": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
                    "state_sha256": hashlib.sha256(
                        args.output.read_bytes()
                    ).hexdigest(),
                    "audit_sha256": hashlib.sha256(audit_path.read_bytes()).hexdigest(),
                    "record_count": 14,
                    "status": "completed",
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
