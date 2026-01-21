from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class DataConfig:
    video_path: Path
    fps: float


@dataclass(frozen=True)
class OutputConfig:
    base_dir: Path


@dataclass(frozen=True)
class TrackingConfig:
    max_link_distance: float


@dataclass(frozen=True)
class PreprocessConfig:
    method: str  # "tophat" | "bg_subtract" | "none"
    kernel_size: int


@dataclass(frozen=True)
class ThresholdConfig:
    method: str  # "otsu" | "adaptive"
    invert: bool
    block_size: int


@dataclass(frozen=True)
class FilterConfig:
    min_area_px: float
    max_area_px: float | None
    max_area_frac: float | None
    reject_border_touch: bool


@dataclass(frozen=True)
class DetectionConfig:
    preprocess: PreprocessConfig
    threshold: ThresholdConfig
    filter: FilterConfig


@dataclass(frozen=True)
class ButtEstimationConfig:
    smooth_window: int
    freeze_speed_thresh: float
    features: list[str]


@dataclass(frozen=True)
class SaveConfig:
    contour: bool


@dataclass(frozen=True)
class TrackingButtConfig:
    detection: DetectionConfig
    tracking: TrackingConfig
    butt_estimation: ButtEstimationConfig
    save: SaveConfig


@dataclass(frozen=True)
class Config:
    data: DataConfig
    output: OutputConfig
    tracking_butt: TrackingButtConfig


def _get(dict_obj: dict[str, Any], key: str, default: Any) -> Any:
    value = dict_obj.get(key, default)
    return value if value is not None else default


def load_config(path: Path) -> Config:
    raw: dict[str, Any] = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}

    data_raw = raw.get("data", {}) or {}
    video_path = data_raw.get("video_path") or data_raw.get("data_dir") or "data/sample1.mp4"
    data_cfg = DataConfig(
        video_path=Path(video_path),
        fps=float(_get(data_raw, "fps", 0.0)),
    )

    output_raw = raw.get("output", {}) or {}
    output_cfg = OutputConfig(base_dir=Path(_get(output_raw, "base_dir", "outputs")))

    tracking_raw = raw.get("tracking_butt", {}) or {}

    detection_raw = tracking_raw.get("detection", {}) or {}
    preprocess_raw = detection_raw.get("preprocess", {}) or {}
    preprocess_cfg = PreprocessConfig(
        method=str(_get(preprocess_raw, "method", "tophat")),
        kernel_size=int(_get(preprocess_raw, "kernel_size", 31)),
    )
    threshold_raw = detection_raw.get("threshold", {}) or {}
    threshold_cfg = ThresholdConfig(
        method=str(_get(threshold_raw, "method", "otsu")),
        invert=bool(_get(threshold_raw, "invert", False)),
        block_size=int(_get(threshold_raw, "block_size", 35)),
    )
    filter_raw = detection_raw.get("filter", {}) or {}
    max_area_px = _get(filter_raw, "max_area_px", None)
    filter_cfg = FilterConfig(
        min_area_px=float(_get(filter_raw, "min_area_px", 0.0)),
        max_area_px=float(max_area_px) if max_area_px not in (None, "") else None,
        max_area_frac=float(_get(filter_raw, "max_area_frac", 0.02)),
        reject_border_touch=bool(_get(filter_raw, "reject_border_touch", True)),
    )
    detection_cfg = DetectionConfig(
        preprocess=preprocess_cfg, threshold=threshold_cfg, filter=filter_cfg
    )

    tracking_cfg = TrackingConfig(
        max_link_distance=float(
            _get(tracking_raw.get("tracking", {}) or {}, "max_link_distance", 30.0)
        )
    )

    butt_raw = tracking_raw.get("butt_estimation", {}) or {}
    butt_cfg = ButtEstimationConfig(
        smooth_window=int(_get(butt_raw, "smooth_window", 5)),
        freeze_speed_thresh=float(_get(butt_raw, "freeze_speed_thresh", 0.5)),
        features=list(_get(butt_raw, "features", [])),
    )

    save_raw = tracking_raw.get("save", {}) or {}
    save_cfg = SaveConfig(contour=bool(_get(save_raw, "contour", False)))

    tracking_butt_cfg = TrackingButtConfig(
        detection=detection_cfg,
        tracking=tracking_cfg,
        butt_estimation=butt_cfg,
        save=save_cfg,
    )

    return Config(data=data_cfg, output=output_cfg, tracking_butt=tracking_butt_cfg)


def with_save_contour(config: Config, enabled: bool) -> Config:
    updated_save = replace(config.tracking_butt.save, contour=enabled)
    updated_tb = replace(config.tracking_butt, save=updated_save)
    return replace(config, tracking_butt=updated_tb)
