from __future__ import annotations

import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SCHEMA_PATH = ROOT / "schemas/phase3_clip_metadata.schema.json"
EXAMPLE_PATH = ROOT / "examples/phase3/clip_metadata_minimal.json"


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _check_required(schema: dict[str, Any], instance: dict[str, Any]) -> None:
    missing = set(schema.get("required", [])) - set(instance)
    assert not missing


def _conditional_for_label_source(
    schema: dict[str, Any], label_source: str
) -> dict[str, Any]:
    for conditional in schema["allOf"]:
        label_source_schema = (
            conditional.get("if", {})
            .get("properties", {})
            .get("labels", {})
            .get("properties", {})
            .get("label_source", {})
        )
        if label_source_schema.get("const") == label_source:
            return conditional
        if label_source in label_source_schema.get("enum", []):
            return conditional
    raise AssertionError(f"missing conditional for label_source={label_source}")


def test_phase3_clip_metadata_schema_keeps_common_required_contract() -> None:
    schema = _load_json(SCHEMA_PATH)

    assert schema["properties"]["schema_version"]["const"] == "phase3_clip_metadata/v0"
    assert set(schema["required"]) == {
        "schema_version",
        "dataset_id",
        "source_video",
        "processing_mode",
        "provenance",
        "track",
        "clip",
        "normalization",
        "frames",
        "labels",
        "qc",
    }
    assert set(schema["properties"]["processing_mode"]["enum"]) == {
        "detection_tracking",
        "gt_passthrough",
    }
    assert "group_key" in schema["properties"]["track"]["required"]
    assert "clip_id" not in schema["properties"]["track"]["required"]
    assert schema["properties"]["frames"]["minItems"] == 1


def test_phase3_minimal_fixture_matches_schema_required_fields() -> None:
    schema = _load_json(SCHEMA_PATH)
    example = _load_json(EXAMPLE_PATH)

    _check_required(schema, example)
    for object_key in (
        "source_video",
        "provenance",
        "track",
        "clip",
        "normalization",
        "labels",
        "qc",
    ):
        _check_required(
            schema["properties"][object_key],
            example[object_key],
        )
    _check_required(
        schema["properties"]["frames"]["items"],
        example["frames"][0],
    )

    assert example["processing_mode"] == "gt_passthrough"
    assert example["source_video"]["source_kind"] == "phase2_pseudo"
    assert example["labels"] == {
        "n_flagella": 1,
        "label_source": "phase2_gt",
    }
    assert example["track"]["group_key"] != example["clip"]["clip_id"]
    assert example["track"]["group_key"].endswith(str(example["provenance"]["run_id"]))
    assert example["clip"]["frame_count"] >= len(example["frames"])


def test_phase3_real_video_label_state_is_representable() -> None:
    schema = _load_json(SCHEMA_PATH)
    label_schema = schema["properties"]["labels"]["properties"]
    source_kind_enum = schema["properties"]["source_video"]["properties"][
        "source_kind"
    ]["enum"]

    assert "real_microscopy" in source_kind_enum
    assert "unavailable" in label_schema["label_source"]["enum"]
    assert "null" in label_schema["n_flagella"]["type"]


def test_phase3_real_video_audit_fields_are_representable() -> None:
    schema = _load_json(SCHEMA_PATH)
    source_video_properties = schema["properties"]["source_video"]["properties"]
    frame_properties = schema["properties"]["frames"]["items"]["properties"]
    qc_properties = schema["properties"]["qc"]["properties"]

    for field_name in (
        "frame_count",
        "codec_fourcc",
        "file_size_bytes",
    ):
        assert field_name in source_video_properties

    for field_name in (
        "body_axis_angle_rad",
        "body_length_px",
        "body_width_px",
        "detection_confidence",
    ):
        assert field_name in frame_properties

    assert frame_properties["detection_confidence"]["minimum"] == 0
    assert frame_properties["detection_confidence"]["maximum"] == 1

    for field_name in (
        "detection_confidence_min",
        "tracking_gap_count",
        "notes",
    ):
        assert field_name in qc_properties


def test_phase3_gt_passthrough_requires_non_null_run_id() -> None:
    schema = _load_json(SCHEMA_PATH)
    gt_passthrough_condition = next(
        conditional
        for conditional in schema["allOf"]
        if conditional["if"]["properties"].get("processing_mode", {}).get("const")
        == "gt_passthrough"
    )
    pseudo_source_condition = next(
        conditional
        for conditional in schema["allOf"]
        if "source_video" in conditional["if"].get("properties", {})
    )

    for conditional in (gt_passthrough_condition, pseudo_source_condition):
        provenance_then = conditional["then"]["properties"]["provenance"]
        assert "run_id" in provenance_then["required"]
        assert provenance_then["properties"]["run_id"] == {
            "type": "string",
            "minLength": 1,
        }


def test_phase3_label_source_and_n_flagella_are_conditionally_consistent() -> None:
    schema = _load_json(SCHEMA_PATH)

    unavailable_condition = _conditional_for_label_source(schema, "unavailable")
    unavailable_n_flagella = unavailable_condition["then"]["properties"]["labels"][
        "properties"
    ]["n_flagella"]
    assert unavailable_n_flagella == {"type": "null"}

    phase2_gt_condition = _conditional_for_label_source(schema, "phase2_gt")
    labeled_n_flagella = phase2_gt_condition["then"]["properties"]["labels"][
        "properties"
    ]["n_flagella"]
    assert labeled_n_flagella == {
        "type": "integer",
        "minimum": 0,
    }
    assert _conditional_for_label_source(schema, "manual") is phase2_gt_condition
