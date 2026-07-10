from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]
PROFILE_DIR = REPO_ROOT / "conf" / "phase2_sweeps"

BOOLEAN_KEYS = {
    "describe_profile",
    "dry_run",
    "list_canonical_profiles",
    "list_kind",
    "list_profiles",
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
    "shape_stability_grid": {
        "flagella.initial_helix_axis_from_rear_deg": (
            "initial-helix-axis-from-rear-deg"
        ),
        "flagella.n_flagella": "n-flagella",
        "motor.torque_Nm": "torque-nm",
        "seed.attach_seed": "attach-seed",
        "seed.phase_seed": "phase-seed",
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


def _profile_path_str(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def _infer_profile_role(path: Path) -> str:
    return "heatmap" if "heatmap" in path.stem else "sweep"


def _load_profile_metadata(path: Path, profile: dict[str, Any]) -> dict[str, Any]:
    raw_metadata = profile.get("metadata") or {}
    if not isinstance(raw_metadata, dict):
        raise ValueError(f"profile metadata must be a mapping: {path}")

    role = str(raw_metadata.get("role") or _infer_profile_role(path)).strip()
    if role not in {"sweep", "heatmap"}:
        raise ValueError(f"Unsupported profile role {role!r}: {path}")

    canonical = bool(raw_metadata.get("canonical", True))
    alias_of = raw_metadata.get("alias_of")
    if alias_of is not None:
        alias_of = _profile_path_str(REPO_ROOT / str(alias_of))
        canonical = False

    metadata = {
        "role": role,
        "canonical": canonical,
        "alias_of": alias_of,
        "description": str(raw_metadata.get("description") or "").strip(),
        "recommended_heatmap_profile": raw_metadata.get("recommended_heatmap_profile"),
        "recommended_sweep_profile": raw_metadata.get("recommended_sweep_profile"),
    }
    for key in ("recommended_heatmap_profile", "recommended_sweep_profile"):
        value = metadata[key]
        if value is not None:
            metadata[key] = _profile_path_str(REPO_ROOT / str(value))
    return metadata


def load_profile_entry(path: Path) -> dict[str, Any]:
    profile = load_profile(path)
    entry = dict(profile)
    entry["path"] = _profile_path_str(path)
    entry["kind"] = str(profile.get("kind", "")).strip()
    entry["metadata"] = _load_profile_metadata(path, profile)
    return entry


def list_profile_entries(
    *, role: str | None = None, canonical_only: bool = False
) -> list[dict[str, Any]]:
    entries = [load_profile_entry(path) for path in sorted(PROFILE_DIR.glob("*.yaml"))]
    if role is not None:
        entries = [entry for entry in entries if entry["metadata"]["role"] == role]
    if canonical_only:
        entries = [entry for entry in entries if entry["metadata"]["canonical"]]
    return entries


def format_profile_listing(entries: list[dict[str, Any]]) -> list[str]:
    lines: list[str] = []
    for entry in entries:
        metadata = entry["metadata"]
        status = "canonical" if metadata["canonical"] else "alias"
        description = metadata["description"]
        line = f"{entry['path']}\t{entry['kind']}\t{status}"
        if description:
            line += f"\t{description}"
        lines.append(line)
    return lines


def format_profile_description(
    entry: dict[str, Any],
    catalog: list[dict[str, Any]] | None = None,
) -> list[str]:
    metadata = entry["metadata"]
    status = "canonical" if metadata["canonical"] else "alias"
    lines = [
        f"path: {entry['path']}",
        f"role: {metadata['role']}",
        f"kind: {entry['kind']}",
        f"status: {status}",
    ]
    if metadata["description"]:
        lines.append(f"description: {metadata['description']}")
    if metadata["alias_of"]:
        lines.append(f"alias_of: {metadata['alias_of']}")
    if metadata["recommended_heatmap_profile"]:
        lines.append(
            f"recommended_heatmap_profile: {metadata['recommended_heatmap_profile']}"
        )
    if metadata["recommended_sweep_profile"]:
        lines.append(
            f"recommended_sweep_profile: {metadata['recommended_sweep_profile']}"
        )
    if metadata["canonical"] and catalog is not None:
        aliases = [
            other["path"]
            for other in catalog
            if other["metadata"]["alias_of"] == entry["path"]
        ]
        if aliases:
            lines.append(f"aliases: {', '.join(sorted(aliases))}")
    return lines


def validate_profile_role(entry: dict[str, Any], expected_role: str) -> None:
    actual_role = entry["metadata"]["role"]
    if actual_role == expected_role:
        return
    script_name = "plot_heatmap.py" if actual_role == "heatmap" else "run_sweep.py"
    raise SystemExit(
        f"Profile {entry['path']!r} has role {actual_role!r}; use {script_name}."
    )


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
