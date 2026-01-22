from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class DataConfig:
    video_path: Path
    fps: float
    bac_short_axis_length_um: float
    px_per_um: float


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
    max_minor_factor: float | None
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
class OverlayConfig:
    draw_history: bool
    history_length: int
    hide_history_after: int


@dataclass(frozen=True)
class TrackingButtConfig:
    detection: DetectionConfig
    tracking: TrackingConfig
    butt_estimation: ButtEstimationConfig
    save: SaveConfig
    overlay: OverlayConfig


@dataclass(frozen=True)
class Config:
    data: DataConfig
    output: OutputConfig
    tracking_butt: TrackingButtConfig


def _get(dict_obj: dict[str, Any], key: str, default: Any) -> Any:
    """辞書から値を取得し、None の場合はデフォルトを返す。

    Args:
        dict_obj: 参照する辞書。
        key: 取得するキー。
        default: デフォルト値。

    Returns:
        Any: 取得した値またはデフォルト。
    """
    value = dict_obj.get(key, default)
    return value if value is not None else default


def load_config(path: Path) -> Config:
    """YAML設定を読み込み、型付きConfigを構築する。

    Args:
        path: 設定ファイルパス。

    Returns:
        Config: 型付き設定オブジェクト。
    """
    raw: dict[str, Any] = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}

    data_raw = raw.get("data", {}) or {}
    video_path = data_raw.get("video_path") or "data/sample1.mp4"
    px_per_um_raw = _get(data_raw, "px_per_um", None)
    px2um_raw = _get(data_raw, "px2um", None)
    um_per_px_raw = _get(data_raw, "um_per_px", None)
    px_per_um: float
    if px_per_um_raw not in (None, ""):
        px_per_um = float(px_per_um_raw)
    elif px2um_raw not in (None, ""):
        px_per_um = 1.0 / float(px2um_raw)
    elif um_per_px_raw not in (None, ""):
        px_per_um = 1.0 / float(um_per_px_raw)
    else:
        px_per_um = 1.0
    data_cfg = DataConfig(
        video_path=Path(video_path),
        fps=float(_get(data_raw, "fps", 0.0)),
        bac_short_axis_length_um=float(_get(data_raw, "bac_short_axis_length_um", 1.0)),
        px_per_um=px_per_um,
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
    max_minor_factor = _get(filter_raw, "max_minor_factor", 3.0)
    filter_cfg = FilterConfig(
        min_area_px=float(_get(filter_raw, "min_area_px", 5.0)),
        max_area_px=float(max_area_px) if max_area_px not in (None, "") else None,
        max_area_frac=float(_get(filter_raw, "max_area_frac", 0.05)),
        max_minor_factor=float(max_minor_factor)
        if max_minor_factor not in (None, "")
        else None,
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
    overlay_raw = tracking_raw.get("overlay", {}) or {}
    overlay_cfg = OverlayConfig(
        draw_history=bool(_get(overlay_raw, "draw_history", True)),
        history_length=int(_get(overlay_raw, "history_length", 30)),
        hide_history_after=int(_get(overlay_raw, "hide_history_after", 60)),
    )

    tracking_butt_cfg = TrackingButtConfig(
        detection=detection_cfg,
        tracking=tracking_cfg,
        butt_estimation=butt_cfg,
        save=save_cfg,
        overlay=overlay_cfg,
    )

    return Config(data=data_cfg, output=output_cfg, tracking_butt=tracking_butt_cfg)


def with_save_contour(config: Config, enabled: bool) -> Config:
    """save.contour を上書きした新しい Config を返す。

    Args:
        config: 元の設定。
        enabled: 輪郭保存を有効にするか。

    Returns:
        Config: 保存設定を変更したコピー。
    """
    updated_save = replace(config.tracking_butt.save, contour=enabled)
    updated_tb = replace(config.tracking_butt, save=updated_save)
    return replace(config, tracking_butt=updated_tb)


def apply_overrides(config: Config, overrides: dict[str, Any]) -> Config:
    """key=value の上書きを Config に適用する。

    Args:
        config: 元の設定。
        overrides: 上書き辞書（例: {"data": {"video_path": ...}}）。

    Returns:
        Config: 上書き適用後の設定。
    """
    cfg = config
    data_over = overrides.get("data", {})
    if data_over:
        kwargs: dict[str, Any] = {}
        if "video_path" in data_over:
            kwargs["video_path"] = Path(data_over["video_path"])
        if "fps" in data_over:
            kwargs["fps"] = float(data_over["fps"])
        if "bac_short_axis_length_um" in data_over:
            kwargs["bac_short_axis_length_um"] = float(
                data_over["bac_short_axis_length_um"]
            )
        px_per_um_override: float | None = None
        if "px_per_um" in data_over and data_over["px_per_um"] not in (None, ""):
            px_per_um_override = float(data_over["px_per_um"])
        elif "px2um" in data_over and data_over["px2um"] not in (None, ""):
            px_per_um_override = 1.0 / float(data_over["px2um"])
        elif "um_per_px" in data_over and data_over["um_per_px"] not in (None, ""):
            px_per_um_override = 1.0 / float(data_over["um_per_px"])
        if px_per_um_override is not None:
            kwargs["px_per_um"] = px_per_um_override
        if kwargs:
            cfg = replace(cfg, data=replace(cfg.data, **kwargs))
    return cfg
