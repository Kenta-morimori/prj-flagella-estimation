"""Helpers for generic Phase 2 multi-run campaigns."""

from __future__ import annotations

import itertools
import re
from pathlib import Path
from typing import Any

import yaml

from sim_swim.analysis.flagella_count_behavior import normalize_base_overrides

CAMPAIGN_OVERRIDE_ROOTS = {
    "metadata",
    "base_config",
    "base_overrides",
    "sweep",
    "replay",
    "plot",
    "dataset",
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


def load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Campaign config must be a mapping: {path}")
    return data


def _merge_nested(dst: dict[str, Any], src: dict[str, Any]) -> dict[str, Any]:
    merged = dict(dst)
    for key, value in src.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_nested(merged[key], value)
        else:
            merged[key] = value
    return merged


def _coerce_cli_value(raw: str) -> Any:
    text = raw.strip()
    if text == "":
        return ""
    loaded = yaml.safe_load(text)
    if isinstance(loaded, str):
        lowered = loaded.lower()
        if lowered in {"true", "false"}:
            return lowered == "true"
        try:
            if any(ch in loaded for ch in (".", "e", "E")):
                return float(loaded)
            return int(loaded)
        except ValueError:
            return loaded
    return loaded


def _override_item_to_nested(raw: str) -> dict[str, Any]:
    if "=" not in raw:
        raise ValueError(f"Invalid override; expected KEY=VALUE: {raw}")
    key, value = raw.split("=", 1)
    node: dict[str, Any] = {}
    cursor = node
    parts = key.strip().split(".")
    if not all(parts):
        raise ValueError(f"Invalid override path: {raw}")
    for part in parts[:-1]:
        child: dict[str, Any] = {}
        cursor[part] = child
        cursor = child
    cursor[parts[-1]] = _coerce_cli_value(value)
    return node


def _normalize_axis(axis_name: str, raw_axis: dict[str, Any]) -> dict[str, Any]:
    key = str(raw_axis.get("key") or "").strip()
    if not key:
        raise ValueError(f"sweep.axes.{axis_name}.key is required")
    values = list(raw_axis.get("values") or [])
    if not values:
        raise ValueError(f"sweep.axes.{axis_name}.values must be non-empty")

    labels = list(raw_axis.get("labels") or [])
    if labels and len(labels) != len(values):
        raise ValueError(
            f"sweep.axes.{axis_name}.labels must match values length: "
            f"{len(labels)} != {len(values)}"
        )
    ids = list(raw_axis.get("ids") or [])
    if ids and len(ids) != len(values):
        raise ValueError(
            f"sweep.axes.{axis_name}.ids must match values length: "
            f"{len(ids)} != {len(values)}"
        )

    short_name = str(raw_axis.get("short_name") or axis_name).strip() or axis_name
    normalized_values = [_coerce_cli_value(str(value)) for value in values]
    normalized_labels = [
        str(label)
        for label in (
            labels or [format_axis_value(value) for value in normalized_values]
        )
    ]
    normalized_ids = [
        sanitize_token(str(item))
        for item in (
            ids
            or [
                f"{short_name}_{format_axis_value(value)}"
                for value in normalized_values
            ]
        )
    ]
    if len(set(normalized_ids)) != len(normalized_ids):
        raise ValueError(f"sweep.axes.{axis_name}.ids must be unique")

    return {
        "name": axis_name,
        "key": key,
        "short_name": short_name,
        "values": normalized_values,
        "labels": normalized_labels,
        "ids": normalized_ids,
    }


def normalize_campaign_config(raw_config: dict[str, Any]) -> dict[str, Any]:
    config = dict(raw_config)
    config["base_overrides"] = normalize_base_overrides(
        config.get("base_overrides", {}) or {}
    )
    sweep = dict(config.get("sweep", {}) or {})
    raw_axes = dict(sweep.get("axes", {}) or {})
    if not raw_axes:
        raise ValueError("sweep.axes is required")
    axes = [_normalize_axis(name, dict(axis)) for name, axis in raw_axes.items()]
    sweep["axes"] = axes
    config["sweep"] = sweep
    config["replay"] = dict(config.get("replay", {}) or {})
    config["plot"] = dict(config.get("plot", {}) or {})
    config["dataset"] = dict(config.get("dataset", {}) or {})
    config["output"] = dict(config.get("output", {}) or {})
    config["metadata"] = dict(config.get("metadata", {}) or {})
    return config


def apply_campaign_cli_overrides(
    campaign_config: dict[str, Any],
    cli_overrides: list[str] | None,
) -> dict[str, Any]:
    if not cli_overrides:
        return normalize_campaign_config(campaign_config)

    campaign_nested: dict[str, Any] = {}
    simulation_nested: dict[str, Any] = {}
    for raw in cli_overrides:
        if "=" not in raw:
            raise ValueError(f"Invalid override; expected KEY=VALUE: {raw}")
        key = raw.split("=", 1)[0].strip()
        root = key.split(".", 1)[0]
        if root in CAMPAIGN_OVERRIDE_ROOTS:
            campaign_nested = _merge_nested(
                campaign_nested, _override_item_to_nested(raw)
            )
        elif root in SIMULATION_OVERRIDE_ROOTS:
            simulation_nested = _merge_nested(
                simulation_nested, _override_item_to_nested(raw)
            )
        else:
            raise ValueError(f"Unknown campaign override root: {root}")

    effective = _merge_nested(dict(campaign_config), campaign_nested)
    effective["base_overrides"] = normalize_base_overrides(
        effective.get("base_overrides", {}) or {}
    )
    if simulation_nested:
        effective["base_overrides"] = _merge_nested(
            effective["base_overrides"], simulation_nested
        )
    return normalize_campaign_config(effective)


def dotted_override(key: str, value: Any) -> dict[str, Any]:
    node: dict[str, Any] = {}
    cursor = node
    parts = key.split(".")
    for part in parts[:-1]:
        child: dict[str, Any] = {}
        cursor[part] = child
        cursor = child
    cursor[parts[-1]] = value
    return node


def format_axis_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if value is None:
        return "null"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def sanitize_token(raw: str) -> str:
    token = raw.strip().replace(".", "p")
    token = re.sub(r"[^A-Za-z0-9_-]+", "_", token)
    token = re.sub(r"_+", "_", token).strip("_")
    return token or "value"


def axis_value_lookup(axis: dict[str, Any]) -> dict[str, tuple[Any, str, str, int]]:
    return {
        axis["ids"][index]: (
            axis["values"][index],
            axis["labels"][index],
            axis["ids"][index],
            index,
        )
        for index in range(len(axis["values"]))
    }


def build_campaign_conditions(config: dict[str, Any]) -> list[dict[str, Any]]:
    axes = list(config["sweep"]["axes"])
    layout = dict(config.get("replay", {}) or {})
    row_axis_name = str(layout.get("row_axis") or "").strip() or None
    col_axis_name = str(layout.get("col_axis") or "").strip() or None
    axis_by_name = {axis["name"]: axis for axis in axes}
    if row_axis_name and row_axis_name not in axis_by_name:
        raise ValueError(f"Unknown render.row_axis: {row_axis_name}")
    if col_axis_name and col_axis_name not in axis_by_name:
        raise ValueError(f"Unknown render.col_axis: {col_axis_name}")

    value_maps = {axis["name"]: axis_value_lookup(axis) for axis in axes}
    choices = [list(value_maps[axis["name"]].values()) for axis in axes]

    conditions: list[dict[str, Any]] = []
    for condition_index, combo in enumerate(itertools.product(*choices)):
        nested_override: dict[str, Any] = {}
        axis_values: dict[str, Any] = {}
        axis_labels: dict[str, str] = {}
        axis_ids: dict[str, str] = {}
        axis_order: dict[str, int] = {}
        label_parts: list[str] = []
        id_parts: list[str] = []
        for axis, (value, label, value_id, order_index) in zip(axes, combo):
            nested_override = _merge_nested(
                nested_override, dotted_override(axis["key"], value)
            )
            axis_values[axis["name"]] = value
            axis_labels[axis["name"]] = label
            axis_ids[axis["name"]] = value_id
            axis_order[axis["name"]] = order_index
            label_parts.append(f"{axis['short_name']}={label}")
            id_parts.append(value_id)
        condition = {
            "condition_index": condition_index,
            "condition_id": "__".join(id_parts),
            "condition_label": ", ".join(label_parts),
            "axis_values": axis_values,
            "axis_labels": axis_labels,
            "axis_ids": axis_ids,
            "axis_order": axis_order,
            "config_overrides": _merge_nested(
                config.get("base_overrides", {}) or {}, nested_override
            ),
        }
        if row_axis_name:
            condition["grid_row_axis"] = row_axis_name
            condition["grid_row_index"] = axis_order[row_axis_name]
            condition["grid_row_label"] = axis_labels[row_axis_name]
        if col_axis_name:
            condition["grid_col_axis"] = col_axis_name
            condition["grid_col_index"] = axis_order[col_axis_name]
            condition["grid_col_label"] = axis_labels[col_axis_name]
        conditions.append(condition)
    return conditions


def campaign_axes_metadata(config: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "name": axis["name"],
            "key": axis["key"],
            "short_name": axis["short_name"],
            "labels": list(axis["labels"]),
            "ids": list(axis["ids"]),
            "values": list(axis["values"]),
        }
        for axis in config["sweep"]["axes"]
    ]


def summary_axis_fields(condition: dict[str, Any]) -> dict[str, Any]:
    row = {
        "condition_index": condition["condition_index"],
        "condition_label": condition["condition_label"],
    }
    for axis_name, value in condition["axis_values"].items():
        row[f"axis_{axis_name}_value"] = value
        row[f"axis_{axis_name}_label"] = condition["axis_labels"][axis_name]
        row[f"axis_{axis_name}_id"] = condition["axis_ids"][axis_name]
        row[f"axis_{axis_name}_index"] = condition["axis_order"][axis_name]
    for field in (
        "grid_row_axis",
        "grid_row_index",
        "grid_row_label",
        "grid_col_axis",
        "grid_col_index",
        "grid_col_label",
    ):
        if field in condition:
            row[field] = condition[field]
    return row
