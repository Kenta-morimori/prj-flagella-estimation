"""Issue #183 の reference torque 比較条件を展開する。"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.sim.params import TIME_SCALE_POLICIES, SimulationConfig

POLICIES = frozenset({"fixed-reference", "tracking-reference"})
TIME_BASES = frozenset({"same-real-time", "same-dimensionless-time"})


def _number_list(raw: Any, key: str) -> list[float]:
    values = [float(value) for value in list(raw or [])]
    if not values or any(value <= 0.0 for value in values):
        raise ValueError(f"{key} must be a non-empty list of positive values")
    return values


def build_plan(raw: dict[str, Any]) -> dict[str, Any]:
    """Return reproducible comparison conditions without running a simulation."""

    base_config = Path(str(raw["base_config"]))
    base = SimulationConfig.from_dict(load_yaml(base_config))
    base.validate_execution_supported()
    requested_scale_policy = raw.get("time_scale_policy")
    if requested_scale_policy is not None:
        requested_scale_policy = str(requested_scale_policy)
        if requested_scale_policy not in TIME_SCALE_POLICIES:
            raise ValueError(
                "time_scale_policy must be 'profile_default', 'reference_torque', or "
                "'legacy_fixed_tau_s_1'"
            )
    reference = float(raw["reference_torque_Nm"])
    if reference <= 0.0:
        raise ValueError("reference_torque_Nm must be positive")
    policies = [str(value) for value in raw.get("policies", [])]
    time_bases = [str(value) for value in raw.get("time_bases", [])]
    if not policies or set(policies) - POLICIES:
        raise ValueError(f"policies must use: {sorted(POLICIES)}")
    if not time_bases or set(time_bases) - TIME_BASES:
        raise ValueError(f"time_bases must use: {sorted(TIME_BASES)}")
    scales = _number_list(raw.get("motor_torque_scales"), "motor_torque_scales")
    output_base_dir = str(
        raw.get("output_base_dir") or "outputs/phase2_reference_torque"
    )
    durations = dict(raw.get("durations", {}) or {})
    real_time_s = float(durations.get("same_real_time_s", 0.0))
    dimensionless_tau = float(durations.get("same_dimensionless_tau", 0.0))
    if "same-real-time" in time_bases and real_time_s <= 0.0:
        raise ValueError("durations.same_real_time_s must be positive")
    if "same-dimensionless-time" in time_bases and dimensionless_tau <= 0.0:
        raise ValueError("durations.same_dimensionless_tau must be positive")

    conditions: list[dict[str, Any]] = []
    for policy in policies:
        for time_basis in time_bases:
            for scale in scales:
                motor_torque = reference * scale
                tracking = policy == "tracking-reference"
                overrides: dict[str, Any] = {
                    "motor": {
                        "torque_Nm": motor_torque,
                        "reference_torque_Nm": motor_torque if tracking else reference,
                        "torque_for_forces_override_Nm": 0.0 if tracking else reference,
                        "allow_reference_torque_mismatch": not tracking,
                    },
                    "time": {
                        **(
                            {"scale_policy": requested_scale_policy}
                            if requested_scale_policy is not None
                            else {}
                        ),
                        "duration": {
                            "value": real_time_s
                            if time_basis == "same-real-time"
                            else dimensionless_tau,
                            "unit": "s" if time_basis == "same-real-time" else "tau",
                        },
                    },
                }
                cfg = base.with_overrides(overrides)
                condition_id = f"{policy.replace('-', '_')}_{time_basis.replace('-', '_')}_motor_scale_{scale:g}"
                conditions.append(
                    {
                        "condition_id": condition_id,
                        "comparison_policy": policy,
                        "time_basis": time_basis,
                        "motor_torque_scale": scale,
                        "config_overrides": overrides,
                        "time": cfg.time_manifest(),
                        "equivalence": {
                            "motor_torque_Nm": cfg.motor_torque_Nm,
                            "reference_torque_Nm": cfg.reference_torque_Nm,
                            "torque_for_forces_Nm": cfg.torque_for_forces_Nm,
                            "physical_properties_fixed": not tracking,
                            "reference_linked_properties": tracking,
                            "same_real_time": time_basis == "same-real-time",
                            "same_dimensionless_time": time_basis
                            == "same-dimensionless-time",
                            "tau_tracks_reference_torque": cfg.time_scale_policy
                            == "reference_torque",
                        },
                        "execution_command": (
                            "uv run python scripts/01_simulate_swimming/01_simulate_swimming.py "
                            f"--config {base_config} "
                            f"motor.torque_Nm={motor_torque:.12g} "
                            f"motor.reference_torque_Nm={overrides['motor']['reference_torque_Nm']:.12g} "
                            f"motor.torque_for_forces_override_Nm={overrides['motor']['torque_for_forces_override_Nm']:.12g} "
                            f"motor.allow_reference_torque_mismatch={str(not tracking).lower()} "
                            f"time.scale_policy={cfg.time_scale_policy} "
                            f"time.duration.value={overrides['time']['duration']['value']:.12g} "
                            f"time.duration.unit={overrides['time']['duration']['unit']} "
                            f"output.base_dir={output_base_dir}/{condition_id}"
                        ),
                    }
                )
    return {
        "kind": "phase2_reference_torque_comparison_plan",
        "contract_version": 2,
        "base_config": str(base_config),
        "model_profile": base.model_profile_manifest(),
        "reference_torque_Nm": reference,
        "time_scale_policy": base.time_scale_policy
        if requested_scale_policy is None
        else requested_scale_policy,
        "output_base_dir": output_base_dir,
        "conditions": conditions,
        "interpretation_prohibitions": [
            "fixed-reference と tracking-reference を同じ物理系の感度差として混同しない。",
            "same-dimensionless-time の wall time / step 数を same-real-time の計算効率と比較しない。",
            "明示した legacy_fixed_tau_s_1 では reference torque を変えても tau_s は連動しない。",
            "この比較planは torque、dt_star、policy、dataset v2、supported status の正式採択を表さない。",
        ],
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    plan = build_plan(load_yaml(args.config))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(args.output)


if __name__ == "__main__":
    main()
