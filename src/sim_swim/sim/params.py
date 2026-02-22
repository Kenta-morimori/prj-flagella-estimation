"""Phase2 シミュレーション設定の定義。"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


def _get(raw: dict[str, Any], key: str, default: Any) -> Any:
    """辞書から値を取得し、None/空文字はデフォルトへ置換する。"""
    val = raw.get(key, default)
    return default if val in (None, "") else val


def _state_angles(
    raw: dict[str, Any] | None, default: dict[str, float]
) -> dict[str, float]:
    src = raw or {}
    return {
        "normal": float(_get(src, "normal", default["normal"])),
        "semicoiled": float(_get(src, "semicoiled", default["semicoiled"])),
        "curly1": float(_get(src, "curly1", default["curly1"])),
    }


@dataclass(frozen=True)
class DiscretizationParams:
    """離散化設定。"""

    ds_um: float = 0.8


@dataclass(frozen=True)
class BodyParams:
    """菌体形状設定。"""

    length_total_um: float = 3.0
    diameter_um: float = 0.8
    bond_L_um: float | None = None


@dataclass(frozen=True)
class FlagellumParams:
    """べん毛設定。"""

    n_flagella: int = 5
    length_um: float = 8.0
    bond_L_um: float | None = None
    pitch_um: float = 2.3
    radius_um: float = 0.2
    filament_diameter_um: float = 0.02
    placement_mode: str = "uniform"
    helix_step_um: float = 0.2


@dataclass(frozen=True)
class ScaleParams:
    """ビーズスケール。"""

    bead_radius_a_um: float = 0.15


@dataclass(frozen=True)
class FluidParams:
    """流体設定。"""

    viscosity_Pa_s: float = 1.0e-3


@dataclass(frozen=True)
class MotorParams:
    """モータートルク設定。"""

    torque_Nm: float = 4.0e-18
    reverse_n_flagella: int = 1


@dataclass(frozen=True)
class SpringPotentialParams:
    """Fraenkel型springポテンシャル設定。"""

    H: float = 1.0e-4
    s_um: float = 0.15


@dataclass(frozen=True)
class BendPotentialParams:
    """曲げポテンシャル設定。"""

    kb: float = 2.0e-19
    theta0_deg: dict[str, float] | None = None


@dataclass(frozen=True)
class TorsionPotentialParams:
    """ねじりポテンシャル設定。"""

    kt: float = 2.0e-19
    phi0_deg: dict[str, float] | None = None


@dataclass(frozen=True)
class SpringSpringRepulsionParams:
    """spring-spring 反発設定。"""

    A_ss: float = 2.0e-19
    a_ss_um: float = 0.1
    cutoff_um: float = 0.6


@dataclass(frozen=True)
class PotentialsParams:
    """ポテンシャル設定。"""

    spring: SpringPotentialParams
    bend: BendPotentialParams
    torsion: TorsionPotentialParams
    spring_spring_repulsion: SpringSpringRepulsionParams


@dataclass(frozen=True)
class HookParams:
    """フック設定。"""

    enabled: bool = True
    kb: float = 8.0e-20
    threshold_deg: float = 90.0


@dataclass(frozen=True)
class RunTumbleParams:
    """run-and-tumble の時系列設定。"""

    run_tau: float = 0.20
    tumble_tau: float = 0.08
    semicoiled_tau: float = 0.03
    curly1_tau: float = 0.03


@dataclass(frozen=True)
class BrownianParams:
    """Brownianノイズ設定。"""

    enabled: bool = False
    temperature_K: float = 298.0
    method: str = "cholesky"
    jitter: float = 1.0e-20


@dataclass(frozen=True)
class TimeParams:
    """時間刻み設定。"""

    dt: float = 0.002
    fps_out: float = 25.0
    duration_s: float = 0.5

    @property
    def dt_out(self) -> float:
        return 1.0 / self.fps_out if self.fps_out > 0 else self.dt


@dataclass(frozen=True)
class RenderParams:
    """レンダ設定。"""

    image_size_px: int = 256
    pixel_size_um: float = 0.203
    render_flagella: bool = False
    flagella_linewidth_px: float = 3.0


@dataclass(frozen=True)
class SeedParams:
    """乱数seed設定。"""

    global_seed: int = 0


@dataclass(frozen=True)
class OutputParams:
    """出力設定。"""

    base_dir: str = "outputs"


@dataclass(frozen=True)
class SimulationConfig:
    """Phase2シミュレーション設定。"""

    discretization: DiscretizationParams
    body: BodyParams
    flagella: FlagellumParams
    scale: ScaleParams
    fluid: FluidParams
    motor: MotorParams
    potentials: PotentialsParams
    hook: HookParams
    run_tumble: RunTumbleParams
    brownian: BrownianParams
    time: TimeParams
    render: RenderParams
    seed: SeedParams
    output: OutputParams

    @staticmethod
    def from_dict(raw: dict[str, Any]) -> "SimulationConfig":
        """辞書（YAML読込後）から設定を構築する。"""

        discretization_raw = raw.get("discretization", {}) or {}
        discretization = DiscretizationParams(
            ds_um=float(_get(discretization_raw, "ds_um", 0.8))
        )

        body_raw = raw.get("body", {}) or {}
        body_bond = _get(body_raw, "bond_L_um", None)
        body = BodyParams(
            length_total_um=float(_get(body_raw, "length_total_um", 3.0)),
            diameter_um=float(_get(body_raw, "diameter_um", 0.8)),
            bond_L_um=None if body_bond is None else float(body_bond),
        )

        flag_raw = raw.get("flagella", {}) or {}
        flag_bond = _get(flag_raw, "bond_L_um", None)
        flagella = FlagellumParams(
            n_flagella=int(_get(flag_raw, "n_flagella", 5)),
            length_um=float(_get(flag_raw, "length_um", 8.0)),
            bond_L_um=None if flag_bond is None else float(flag_bond),
            pitch_um=float(_get(flag_raw, "pitch_um", 2.3)),
            radius_um=float(_get(flag_raw, "radius_um", 0.2)),
            filament_diameter_um=float(_get(flag_raw, "filament_diameter_um", 0.02)),
            placement_mode=str(_get(flag_raw, "placement_mode", "uniform")),
            helix_step_um=float(_get(flag_raw, "helix_step_um", 0.2)),
        )

        scale_raw = raw.get("scale", {}) or {}
        scale = ScaleParams(
            bead_radius_a_um=float(_get(scale_raw, "bead_radius_a_um", 0.15))
        )

        fluid_raw = raw.get("fluid", {}) or {}
        fluid = FluidParams(
            viscosity_Pa_s=float(_get(fluid_raw, "viscosity_Pa_s", 1e-3))
        )

        motor_raw = raw.get("motor", {}) or {}
        motor = MotorParams(
            torque_Nm=float(_get(motor_raw, "torque_Nm", 4e-18)),
            reverse_n_flagella=int(_get(motor_raw, "reverse_n_flagella", 1)),
        )

        spring_raw = (raw.get("potentials", {}) or {}).get("spring", {}) or {}
        spring = SpringPotentialParams(
            H=float(_get(spring_raw, "H", 1e-4)),
            s_um=float(_get(spring_raw, "s_um", 0.15)),
        )

        bend_raw = (raw.get("potentials", {}) or {}).get("bend", {}) or {}
        bend = BendPotentialParams(
            kb=float(_get(bend_raw, "kb", 2e-19)),
            theta0_deg=_state_angles(
                bend_raw.get("theta0_deg"),
                {"normal": 25.0, "semicoiled": 55.0, "curly1": 75.0},
            ),
        )

        torsion_raw = (raw.get("potentials", {}) or {}).get("torsion", {}) or {}
        torsion = TorsionPotentialParams(
            kt=float(_get(torsion_raw, "kt", 2e-19)),
            phi0_deg=_state_angles(
                torsion_raw.get("phi0_deg"),
                {"normal": 15.0, "semicoiled": 95.0, "curly1": 145.0},
            ),
        )

        ssr_raw = (raw.get("potentials", {}) or {}).get(
            "spring_spring_repulsion", {}
        ) or {}
        spring_spring_repulsion = SpringSpringRepulsionParams(
            A_ss=float(_get(ssr_raw, "A_ss", 2e-19)),
            a_ss_um=float(_get(ssr_raw, "a_ss_um", 0.1)),
            cutoff_um=float(_get(ssr_raw, "cutoff_um", 0.6)),
        )

        potentials = PotentialsParams(
            spring=spring,
            bend=bend,
            torsion=torsion,
            spring_spring_repulsion=spring_spring_repulsion,
        )

        hook_raw = raw.get("hook", {}) or {}
        hook = HookParams(
            enabled=bool(_get(hook_raw, "enabled", True)),
            kb=float(_get(hook_raw, "kb", 8e-20)),
            threshold_deg=float(_get(hook_raw, "threshold_deg", 90.0)),
        )

        rt_raw = raw.get("run_tumble", {}) or {}
        run_tumble = RunTumbleParams(
            run_tau=float(_get(rt_raw, "run_tau", 0.20)),
            tumble_tau=float(_get(rt_raw, "tumble_tau", 0.08)),
            semicoiled_tau=float(_get(rt_raw, "semicoiled_tau", 0.03)),
            curly1_tau=float(_get(rt_raw, "curly1_tau", 0.03)),
        )

        brown_raw = raw.get("brownian", {}) or {}
        brownian = BrownianParams(
            enabled=bool(_get(brown_raw, "enabled", False)),
            temperature_K=float(_get(brown_raw, "temperature_K", 298.0)),
            method=str(_get(brown_raw, "method", "cholesky")),
            jitter=float(_get(brown_raw, "jitter", 1e-20)),
        )

        time_raw = raw.get("time", {}) or {}
        time = TimeParams(
            dt=float(_get(time_raw, "dt", 0.002)),
            fps_out=float(_get(time_raw, "fps_out", 25.0)),
            duration_s=float(_get(time_raw, "duration_s", 0.5)),
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

        output_raw = raw.get("output", {}) or {}
        output = OutputParams(base_dir=str(_get(output_raw, "base_dir", "outputs")))

        return SimulationConfig(
            discretization=discretization,
            body=body,
            flagella=flagella,
            scale=scale,
            fluid=fluid,
            motor=motor,
            potentials=potentials,
            hook=hook,
            run_tumble=run_tumble,
            brownian=brownian,
            time=time,
            render=render,
            seed=seed,
            output=output,
        )

    def with_overrides(self, overrides: dict[str, Any]) -> "SimulationConfig":
        """ネスト辞書の上書きから新しい設定を生成する。"""

        merged = asdict(self)

        def _merge(dst: dict[str, Any], src: dict[str, Any]) -> None:
            for key, value in src.items():
                if isinstance(value, dict) and isinstance(dst.get(key), dict):
                    _merge(dst[key], value)
                else:
                    dst[key] = value

        _merge(merged, overrides)
        return SimulationConfig.from_dict(merged)


def _coerce_value(text: str) -> Any:
    """CLI override 文字列を bool/int/float/None へ変換する。"""

    val = text.strip()
    low = val.lower()
    if low in {"null", "none"}:
        return None
    if low in {"true", "false"}:
        return low == "true"

    try:
        if any(ch in val for ch in (".", "e", "E")):
            return float(val)
        return int(val)
    except ValueError:
        return val


def merge_overrides(
    overrides: dict[str, Any], items: list[str] | None
) -> dict[str, Any]:
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
        node[parts[-1]] = _coerce_value(value)
    return merged
