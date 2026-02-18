"""シミュレーション用のパラメータ定義（Phase2用）。"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict


def _get(raw: Dict[str, Any], key: str, default: Any) -> Any:
    """辞書から値を取得し、None/空文字の場合はデフォルトを返す。"""
    val = raw.get(key, default)
    return default if val in (None, "") else val


@dataclass(frozen=True)
class BodyParams:
    """菌体（spherocylinder）の幾何パラメータ。"""

    length_total_um: float = 3.0
    diameter_um: float = 0.8

    @property
    def length_cyl_um(self) -> float:
        """円柱部の長さ（全長から両端の半球2個の直径を引いた値）。"""
        return max(0.0, self.length_total_um - self.diameter_um)


@dataclass(frozen=True)
class FlagellumParams:
    """べん毛の幾何・駆動パラメータ。"""

    n_flagella: int = 5  # デフォルト本数
    length_um: float = 12.0
    pitch_um: float = 2.3
    radius_um: float = 0.2
    filament_diameter_um: float = 0.02
    placement_mode: str = "uniform"  # uniform / biased など
    motor_freq_hz: float = 100.0
    helix_step_um: float = 0.2  # 描画用の離散間隔


@dataclass(frozen=True)
class EnvironmentParams:
    """媒質条件と外力スイッチ。"""

    viscosity_mpas: float = 1.0
    temperature_k: float = 298.0
    include_brownian: bool = True
    include_gravity: bool = False


@dataclass(frozen=True)
class TimeParams:
    """時間刻みと出力レート設定。"""

    fps_out: float = 50.0
    dt_sim: float = 0.0005
    duration_s: float = 5.0

    @property
    def dt_out(self) -> float:
        return 1.0 / self.fps_out if self.fps_out > 0 else 0.02


@dataclass(frozen=True)
class RenderParams:
    """2D投影・デバッグ描画設定。"""

    image_size_px: int = 256
    pixel_size_um: float = 0.203
    render_flagella: bool = False
    flagella_linewidth_px: float = 3.0


@dataclass(frozen=True)
class SeedParams:
    """乱数 seed をまとめて管理する。"""

    global_seed: int = 0


@dataclass(frozen=True)
class SimulationConfig:
    """シミュレーション全体の設定。"""

    body: BodyParams
    flagella: FlagellumParams
    env: EnvironmentParams
    time: TimeParams
    render: RenderParams
    seed: SeedParams

    @staticmethod
    def from_dict(raw: Dict[str, Any]) -> "SimulationConfig":
        """辞書（YAML読込後）からデフォルトを補完して構築する。"""

        body_raw = raw.get("body", {}) or {}
        body = BodyParams(
            length_total_um=float(_get(body_raw, "length_total_um", 3.0)),
            diameter_um=float(_get(body_raw, "diameter_um", 0.8)),
        )

        flag_raw = raw.get("flagella", {}) or {}
        flagella = FlagellumParams(
            n_flagella=int(_get(flag_raw, "n_flagella", FlagellumParams().n_flagella)),
            length_um=float(_get(flag_raw, "length_um", 12.0)),
            pitch_um=float(_get(flag_raw, "pitch_um", 2.3)),
            radius_um=float(_get(flag_raw, "radius_um", 0.2)),
            filament_diameter_um=float(_get(flag_raw, "filament_diameter_um", 0.02)),
            placement_mode=str(_get(flag_raw, "placement_mode", "uniform")),
            motor_freq_hz=float(_get(flag_raw, "motor_freq_hz", 100.0)),
        )

        env_raw = raw.get("environment", {}) or {}
        env = EnvironmentParams(
            viscosity_mpas=float(_get(env_raw, "viscosity_mpas", 1.0)),
            temperature_k=float(_get(env_raw, "temperature_k", 298.0)),
            include_brownian=bool(_get(env_raw, "include_brownian", True)),
            include_gravity=bool(_get(env_raw, "include_gravity", False)),
        )

        time_raw = raw.get("time", {}) or {}
        fps_out = float(_get(time_raw, "fps_out", 50.0))
        time_params = TimeParams(
            fps_out=fps_out,
            dt_sim=float(_get(time_raw, "dt_sim", 0.0005)),
            duration_s=float(_get(time_raw, "duration_s", 5.0)),
        )

        render_raw = raw.get("render", {}) or {}
        render = RenderParams(
            image_size_px=int(_get(render_raw, "image_size_px", 256)),
            pixel_size_um=float(_get(render_raw, "pixel_size_um", 0.203)),
            render_flagella=bool(_get(render_raw, "render_flagella", False)),
            flagella_linewidth_px=float(_get(render_raw, "flagella_linewidth_px", 3.0)),
        )

        seed_raw = raw.get("seed", {}) or {}
        seed = SeedParams(global_seed=int(_get(seed_raw, "global_seed", 0)))

        return SimulationConfig(
            body=body,
            flagella=flagella,
            env=env,
            time=time_params,
            render=render,
            seed=seed,
        )

    def with_overrides(self, overrides: Dict[str, Any]) -> "SimulationConfig":
        """ネスト辞書の上書きから新しい Config を生成する。"""

        cfg_dict: Dict[str, Any] = {
            "body": vars(self.body),
            "flagella": vars(self.flagella),
            "environment": vars(self.env),
            "time": vars(self.time),
            "render": vars(self.render),
            "seed": vars(self.seed),
        }

        # ネスト辞書を更新
        for top_key, sub in overrides.items():
            if not isinstance(sub, dict):
                cfg_dict[top_key] = sub
                continue
            current = cfg_dict.get(top_key, {}) or {}
            current.update(sub)
            cfg_dict[top_key] = current

        return SimulationConfig.from_dict(cfg_dict)


def merge_overrides(
    overrides: Dict[str, Any], items: list[str] | None
) -> Dict[str, Any]:
    """`a.b.c=val` 形式をネスト辞書へ畳み込む。"""

    if not items:
        return overrides
    merged = {**overrides}
    for raw in items:
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        parts = key.split(".")
        node = merged
        for p in parts[:-1]:
            node = node.setdefault(p, {})  # type: ignore[assignment]
        node[parts[-1]] = value
    return merged
