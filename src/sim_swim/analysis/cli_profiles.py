from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


BOOLEAN_KEYS = {
    "dry_run",
    "list_kind",
    "no_stop_on_shape_fail",
    "overwrite",
    "stop_on_shape_fail",
}

SWEEP_OVERRIDE_ALIASES = {
    "bundling_alignment": {
        "flagella.initial_helix_axis_from_rear_deg": "helix-axis-angles-deg",
        "flagella.n_flagella": "n-flagella",
        "motor.torque_Nm": "torques",
        "time.dt_star": "dt-star",
        "time.duration_s": "duration-s",
    },
    "hook_overstretch": {
        "flagella.n_flagella": "n-flagella",
        "motor.torque_Nm": "torque-nm",
        "seed.attach_seed": "attach-seed",
        "seed.phase_seed": "phase-seed",
        "time.dt_star": "dt-star",
        "time.duration_s": "duration-s",
    },
    "motor_scale": {
        "flagella.n_flagella": "n-flagella",
        "motor.torque_Nm": "torques",
        "time.dt_star": "dt-star",
        "time.duration_s": "duration",
    },
    "single_flagellum_torque": {
        "motor.torque_Nm": "torques",
        "time.dt_star": "dt-star",
        "time.duration_s": "duration",
    },
}


def load_profile(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError(f"profile must be a mapping: {path}")
    return data


def args_from_profile(profile: dict[str, Any]) -> list[str]:
    raw_args = profile.get("args", {})
    if raw_args is None:
        return []
    if not isinstance(raw_args, dict):
        raise ValueError("profile 'args' must be a mapping")

    out: list[str] = []
    for key, value in raw_args.items():
        option = "--" + str(key).replace("_", "-")
        if isinstance(value, bool):
            if value:
                out.append(option)
            continue
        if value is None:
            continue
        out.append(option)
        if isinstance(value, (list, tuple)):
            out.append(",".join(str(item) for item in value))
        else:
            out.append(str(value))
    return out


def split_config_key(argv: list[str]) -> tuple[Path | None, list[str]]:
    """Extract wrapper-level config=... while preserving legacy --config."""

    config_path: Path | None = None
    out: list[str] = []
    for raw in argv:
        if raw.startswith("config="):
            if config_path is not None:
                raise ValueError("config= specified multiple times")
            value = raw.split("=", 1)[1].strip()
            if not value:
                raise ValueError("config= requires a path")
            config_path = Path(value)
            continue
        out.append(raw)
    return config_path, out


def key_value_args_to_cli_args(
    raw_args: list[str],
    *,
    aliases: dict[str, str] | None = None,
) -> list[str]:
    """Convert KEY=VALUE wrapper overrides to argparse-style options."""

    aliases = aliases or {}
    out: list[str] = []
    for raw in raw_args:
        if raw.startswith("--") or "=" not in raw:
            out.append(raw)
            continue

        key, value = raw.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            raise ValueError(f"Invalid override; expected KEY=VALUE: {raw}")
        if key == "config":
            raise ValueError("config= must be specified before wrapper parsing")

        option_name = aliases.get(key)
        if option_name is None:
            if "." in key:
                out.append(raw)
                continue
            option_name = key.replace("_", "-")

        option = "--" + option_name
        if key in BOOLEAN_KEYS:
            normalized = value.lower()
            if normalized in {"1", "true", "yes", "y", "on"}:
                out.append(option)
                continue
            if normalized in {"0", "false", "no", "n", "off"}:
                raise ValueError(
                    f"{key}=false is not supported; omit the option instead."
                )
            raise ValueError(f"Invalid boolean value for {key}: {value}")

        out.extend([option, value])
    return out


def sweep_aliases(kind: str) -> dict[str, str]:
    return SWEEP_OVERRIDE_ALIASES.get(kind, {})
