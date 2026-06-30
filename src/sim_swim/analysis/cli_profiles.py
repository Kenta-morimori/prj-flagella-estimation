from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


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
