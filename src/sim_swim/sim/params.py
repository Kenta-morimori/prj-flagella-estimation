"""Phase2 シミュレーション設定の定義。"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import math
from typing import Any

K_B = 1.380649e-23
DT_STAR_TARGET = 1.0e-3


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


def _isclose(
    a: float, b: float, rel_tol: float = 1.0e-12, abs_tol: float = 1.0e-12
) -> bool:
    return math.isclose(float(a), float(b), rel_tol=rel_tol, abs_tol=abs_tol)


@dataclass(frozen=True)
class ScaleParams:
    """スケール設定。"""

    b_um: float = 1.0
    bead_radius_a_over_b: float = 0.1


@dataclass(frozen=True)
class BodyPrismParams:
    """菌体のn角柱離散化設定。"""

    n_prism: int = 3
    dz_over_b: float = 0.5
    radius_over_b: float = 0.5
    axis: str = "x"


@dataclass(frozen=True)
class BodyParams:
    """菌体設定。"""

    prism: BodyPrismParams
    length_total_um: float = 2.0


@dataclass(frozen=True)
class FlagellaDiscretizationParams:
    """べん毛離散化設定。"""

    ds_over_b: float = 0.58


@dataclass(frozen=True)
class FlagellaHelixInitParams:
    """べん毛の初期ヘリックス形状設定。"""

    radius_over_b: float = 0.25
    pitch_over_b: float = 2.5


@dataclass(frozen=True)
class FlagellumParams:
    """べん毛設定。"""

    n_flagella: int = 3
    placement_mode: str = "uniform"
    discretization: FlagellaDiscretizationParams = field(
        default_factory=FlagellaDiscretizationParams
    )
    bond_L_over_b: float = 0.58
    length_over_b: float = 5.8
    helix_init: FlagellaHelixInitParams = field(default_factory=FlagellaHelixInitParams)


@dataclass(frozen=True)
class FluidParams:
    """流体設定。"""

    viscosity_Pa_s: float = 1.0e-3


@dataclass(frozen=True)
class MotorParams:
    """モータ設定。"""

    torque_Nm: float = 4.0e-18
    reverse_n_flagella: int = 1
    enable_switching: bool = False


@dataclass(frozen=True)
class SpringPotentialParams:
    """springポテンシャル設定。"""

    H_over_T_over_b: float = 10.0
    s: float = 0.1


@dataclass(frozen=True)
class BendPotentialParams:
    """bendingポテンシャル設定。"""

    kb_over_T: float = 20.0
    theta0_deg: dict[str, float] | None = None


@dataclass(frozen=True)
class TorsionPotentialParams:
    """torsionポテンシャル設定。"""

    kt_over_T: float = 10.0
    phi0_deg: dict[str, float] | None = None


@dataclass(frozen=True)
class SpringSpringRepulsionParams:
    """spring-spring反発ポテンシャル設定。"""

    A_ss_over_T: float = 1.0
    a_ss_over_b: float = 0.2
    cutoff_over_b: float = 0.2


@dataclass(frozen=True)
class PotentialsParams:
    """ポテンシャル設定。"""

    spring: SpringPotentialParams
    bend: BendPotentialParams
    torsion: TorsionPotentialParams
    spring_spring_repulsion: SpringSpringRepulsionParams


@dataclass(frozen=True)
class HookParams:
    """hook設定。"""

    enabled: bool = True
    threshold_deg: float = 90.0
    kb_over_T: float = 20.0


@dataclass(frozen=True)
class RunTumbleParams:
    """run-and-tumble設定。"""

    run_tau: float = 1200.0
    tumble_tau: float = 800.0
    semicoiled_tau: float = 400.0
    curly1_tau: float = 400.0


@dataclass(frozen=True)
class BrownianParams:
    """Brownian項設定。"""

    enabled: bool = False
    temperature_K: float = 298.0
    method: str = "cholesky"
    jitter: float = 1.0e-20


@dataclass(frozen=True)
class ProjectionParams:
    """拘束投影設定。"""

    enable_flagella_template_projection: bool = True
    enable_flagella_chain_length_projection_when_template_off: bool = True


@dataclass(frozen=True)
class StiffnessScaleParams:
    """剛性スケール設定。"""

    body: float = 200.0
    flag_bend: float = 300.0
    flag_torsion: float = 300.0


@dataclass(frozen=True)
class LocalHelixParams:
    """基部近傍の局所ヘリックス拘束設定。"""

    enabled: bool = False
    n_local: int = 4
    k_radius_over_torque: float = 0.5
    k_phase_over_torque: float = 0.5


@dataclass(frozen=True)
class CollapseDiagnosticsParams:
    """collapse診断設定。"""

    enabled: bool = True
    write_every_step: bool = True
    max_flagella_points: int = 4
    collapse_distance_um: float = 0.15
    collapse_consecutive_steps: int = 3


@dataclass(frozen=True)
class DiagnosticsParams:
    """診断設定。"""

    collapse: CollapseDiagnosticsParams = field(
        default_factory=CollapseDiagnosticsParams
    )


@dataclass(frozen=True)
class TimeParams:
    """時間設定。"""

    duration_s: float = 0.1
    dt_s: float = 1.0e-3


@dataclass(frozen=True)
class OutputSamplingParams:
    """出力サンプリング設定。"""

    out_all_steps_3d: bool = True
    fps_out_2d: float = 25.0


@dataclass(frozen=True)
class RenderParams:
    """描画設定。"""

    image_size_px: int = 256
    pixel_size_um: float = 0.203
    flagella_linewidth_px: float = 2.0
    render_flagella: bool = True
    render_flagella_2d: bool = False

    save_frames_3d: bool = True
    follow_camera_3d: bool = True
    view_range_um: float = 8.0
    timestamp_3d: bool = True
    timestamp_fmt: str = "t = {t:.3f} s"
    label_flagella: bool = True

    follow_camera_2d: bool = False
    center_body_in_2d: bool = True
    save_frames_2d: bool = True


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

    scale: ScaleParams
    body: BodyParams
    flagella: FlagellumParams
    fluid: FluidParams
    motor: MotorParams
    potentials: PotentialsParams
    hook: HookParams
    run_tumble: RunTumbleParams
    time: TimeParams
    output_sampling: OutputSamplingParams
    brownian: BrownianParams
    projection: ProjectionParams
    stiffness_scales: StiffnessScaleParams
    local_helix: LocalHelixParams
    diagnostics: DiagnosticsParams
    render: RenderParams
    seed: SeedParams
    output: OutputParams

    @property
    def b_m(self) -> float:
        return max(self.scale.b_um, 1e-9) * 1e-6

    @property
    def viscosity_Pa_s(self) -> float:
        return max(self.fluid.viscosity_Pa_s, 1e-12)

    @property
    def input_torque_Nm(self) -> float:
        return float(self.motor.torque_Nm)

    @property
    def use_eta_b3_torque(self) -> bool:
        return math.isclose(self.input_torque_Nm, -1.0, rel_tol=0.0, abs_tol=1e-12)

    @property
    def is_motor_off_torque(self) -> bool:
        return math.isclose(self.input_torque_Nm, 0.0, rel_tol=0.0, abs_tol=1e-30)

    @property
    def torque_eta_b3_Nm(self) -> float:
        return self.viscosity_Pa_s * (self.b_m**3)

    @property
    def torque_scale_Nm(self) -> float:
        """τ計算に使うスケールトルク絶対値。"""
        if self.use_eta_b3_torque or self.is_motor_off_torque:
            return self.torque_eta_b3_Nm
        return abs(self.input_torque_Nm)

    @property
    def torque_for_forces_Nm(self) -> float:
        """力学パラメータ（剛性/反発）のスケーリングに使うトルク絶対値。"""
        if self.use_eta_b3_torque or self.is_motor_off_torque:
            return self.torque_eta_b3_Nm
        return abs(self.motor_torque_Nm)

    @property
    def motor_torque_Nm(self) -> float:
        """モータ力計算で使う符号付きトルク。"""
        if self.use_eta_b3_torque:
            return self.torque_eta_b3_Nm
        if self.is_motor_off_torque:
            return 0.0
        return self.input_torque_Nm

    @property
    def torque_Nm(self) -> float:
        """後方互換: モータ力計算で使う符号付きトルク。"""
        return self.motor_torque_Nm

    @property
    def bead_radius_m(self) -> float:
        return self.scale.bead_radius_a_over_b * self.b_m

    @property
    def tau_s(self) -> float:
        """内部計算で使う時間スケールτ[s]。Phase2では常に1固定。"""
        return 1.0

    @property
    def output_dt_s(self) -> float:
        """設定上の出力間隔[s]。"""
        return max(self.time.dt_s, 0.0)

    @property
    def dt_s(self) -> float:
        """内部計算ステップ幅Δt[s]。"""
        return DT_STAR_TARGET

    @property
    def dt_star(self) -> float:
        return self.dt_s / max(self.tau_s, 1e-30)

    @property
    def duration_star(self) -> float:
        return max(self.time.duration_s, 0.0) / max(self.tau_s, 1e-30)

    @property
    def total_steps(self) -> int:
        return max(1, int(math.ceil(self.duration_star / max(self.dt_star, 1e-30))))

    def validate_time_scaling(
        self,
        target_dt_star: float = DT_STAR_TARGET,
        rel_tol: float = 1e-12,
    ) -> None:
        """内部計算で `dt_star=1e-3` が成立しているか検証する。"""

        if not math.isclose(
            self.dt_star, target_dt_star, rel_tol=rel_tol, abs_tol=1e-18
        ):
            raise ValueError(
                "内部時間刻みが論文条件に一致しない。"
                f" dt_s={self.dt_s:.12e}, tau_s={self.tau_s:.12e},"
                f" dt_star={self.dt_star:.12e}, target_dt_star={target_dt_star:.12e}; "
                "dt_star(=Δt/τ) must be 1e-3"
            )

    def compute_body_n_layers(
        self,
        rel_tol: float = 1e-12,
        abs_tol: float = 1e-12,
    ) -> int:
        """`length_total_um` と `dz_over_b` から n_layers を算出する。"""

        dz_over_b = float(self.body.prism.dz_over_b)
        if dz_over_b <= 0.0:
            raise ValueError("body.prism.dz_over_b は正である必要がある。")

        l_over_b = float(self.body.length_total_um) / max(self.scale.b_um, 1e-12)
        if l_over_b <= 0.0:
            raise ValueError("body.length_total_um は正である必要がある。")

        ratio = l_over_b / dz_over_b
        n_intervals = int(round(ratio))
        if n_intervals < 1 or not math.isclose(
            ratio, float(n_intervals), rel_tol=rel_tol, abs_tol=abs_tol
        ):
            raise ValueError(
                "body.length_total_um / scale.b_um は body.prism.dz_over_b の整数倍である必要がある。"
                f" L_over_b={l_over_b:.12e}, dz_over_b={dz_over_b:.12e}, ratio={ratio:.12e}"
            )

        return n_intervals + 1

    def paper_reference_checks(self) -> list[tuple[str, bool, str, str]]:
        """論文条件との一致チェックを返す。"""

        checks: list[tuple[str, bool, str, str]] = []

        def add(name: str, actual: Any, expected: Any) -> None:
            if isinstance(actual, (int, float)) and isinstance(expected, (int, float)):
                ok = _isclose(float(actual), float(expected))
            elif isinstance(actual, dict) and isinstance(expected, dict):
                ok = set(actual.keys()) == set(expected.keys()) and all(
                    _isclose(float(actual[k]), float(expected[k])) for k in expected
                )
            else:
                ok = actual == expected
            checks.append((name, ok, repr(actual), repr(expected)))

        add("flagella.n_flagella", self.flagella.n_flagella, 3)
        add("scale.bead_radius_a_over_b", self.scale.bead_radius_a_over_b, 0.1)
        add(
            "potentials.spring.H_over_T_over_b",
            self.potentials.spring.H_over_T_over_b,
            10.0,
        )
        add("potentials.spring.s", self.potentials.spring.s, 0.1)
        add("potentials.bend.kb_over_T", self.potentials.bend.kb_over_T, 20.0)
        add("potentials.torsion.kt_over_T", self.potentials.torsion.kt_over_T, 10.0)
        add(
            "potentials.bend.theta0_deg",
            self.potentials.bend.theta0_deg
            or {"normal": 142.0, "semicoiled": 90.0, "curly1": 105.0},
            {"normal": 142.0, "semicoiled": 90.0, "curly1": 105.0},
        )
        add(
            "potentials.torsion.phi0_deg",
            self.potentials.torsion.phi0_deg
            or {"normal": -60.0, "semicoiled": 65.0, "curly1": 120.0},
            {"normal": -60.0, "semicoiled": 65.0, "curly1": 120.0},
        )
        add("run_tumble.run_tau", self.run_tumble.run_tau, 1200.0)
        add("run_tumble.tumble_tau", self.run_tumble.tumble_tau, 800.0)
        add("run_tumble.semicoiled_tau", self.run_tumble.semicoiled_tau, 400.0)
        add("run_tumble.curly1_tau", self.run_tumble.curly1_tau, 400.0)

        return checks

    @staticmethod
    def from_dict(raw: dict[str, Any]) -> "SimulationConfig":
        """辞書（YAML読込後）から設定を構築する。"""

        brown_raw = raw.get("brownian", {}) or {}
        brownian = BrownianParams(
            enabled=bool(_get(brown_raw, "enabled", False)),
            temperature_K=float(_get(brown_raw, "temperature_K", 298.0)),
            method=str(_get(brown_raw, "method", "cholesky")),
            jitter=float(_get(brown_raw, "jitter", 1e-20)),
        )

        projection_raw = raw.get("projection", {}) or {}
        projection = ProjectionParams(
            enable_flagella_template_projection=bool(
                _get(projection_raw, "enable_flagella_template_projection", True)
            ),
            enable_flagella_chain_length_projection_when_template_off=bool(
                _get(
                    projection_raw,
                    "enable_flagella_chain_length_projection_when_template_off",
                    True,
                )
            ),
        )

        stiffness_raw = raw.get("stiffness_scales", {}) or {}
        stiffness_scales = StiffnessScaleParams(
            body=float(_get(stiffness_raw, "body", 200.0)),
            flag_bend=float(_get(stiffness_raw, "flag_bend", 300.0)),
            flag_torsion=float(_get(stiffness_raw, "flag_torsion", 300.0)),
        )

        local_helix_raw = raw.get("local_helix", {}) or {}
        local_helix = LocalHelixParams(
            enabled=bool(_get(local_helix_raw, "enabled", False)),
            n_local=int(_get(local_helix_raw, "n_local", 4)),
            k_radius_over_torque=float(
                _get(local_helix_raw, "k_radius_over_torque", 0.5)
            ),
            k_phase_over_torque=float(
                _get(local_helix_raw, "k_phase_over_torque", 0.5)
            ),
        )

        diagnostics_raw = raw.get("diagnostics", {}) or {}
        collapse_raw = diagnostics_raw.get("collapse", {}) or {}
        diagnostics = DiagnosticsParams(
            collapse=CollapseDiagnosticsParams(
                enabled=bool(_get(collapse_raw, "enabled", True)),
                write_every_step=bool(_get(collapse_raw, "write_every_step", True)),
                max_flagella_points=int(_get(collapse_raw, "max_flagella_points", 4)),
                collapse_distance_um=float(
                    _get(collapse_raw, "collapse_distance_um", 0.15)
                ),
                collapse_consecutive_steps=int(
                    _get(collapse_raw, "collapse_consecutive_steps", 3)
                ),
            )
        )

        scale_raw = raw.get("scale", {}) or {}
        b_um = float(_get(scale_raw, "b_um", 1.0))
        bead_radius_a_over_b = scale_raw.get("bead_radius_a_over_b")
        if bead_radius_a_over_b is None:
            bead_radius_um = scale_raw.get("bead_radius_a_um")
            if bead_radius_um is not None:
                bead_radius_a_over_b = float(bead_radius_um) / max(b_um, 1e-12)
            else:
                bead_radius_a_over_b = 0.1
        scale = ScaleParams(
            b_um=b_um,
            bead_radius_a_over_b=float(bead_radius_a_over_b),
        )

        body_raw = raw.get("body", {}) or {}
        prism_raw = body_raw.get("prism", {}) or {}
        diameter_um = body_raw.get("diameter_um", None)
        radius_over_b = prism_raw.get("radius_over_b", None)
        if radius_over_b is None and diameter_um is not None:
            radius_over_b = (float(diameter_um) * 0.5) / max(scale.b_um, 1e-12)

        prism = BodyPrismParams(
            n_prism=int(_get(prism_raw, "n_prism", 3)),
            dz_over_b=float(_get(prism_raw, "dz_over_b", 0.5)),
            radius_over_b=float(radius_over_b if radius_over_b is not None else 0.5),
            axis=str(_get(prism_raw, "axis", "x")),
        )

        body = BodyParams(
            prism=prism,
            length_total_um=float(_get(body_raw, "length_total_um", 2.0 * scale.b_um)),
        )

        flag_raw = raw.get("flagella", {}) or {}
        flag_dis_raw = flag_raw.get("discretization", {}) or {}
        ds_over_b = flag_dis_raw.get("ds_over_b")
        if ds_over_b is None:
            ds_um = flag_dis_raw.get(
                "ds_um", raw.get("discretization", {}).get("ds_um")
            )
            if ds_um is not None:
                ds_over_b = float(ds_um) / max(scale.b_um, 1e-12)
            else:
                ds_over_b = 0.58

        bond_L_over_b = flag_raw.get("bond_L_over_b")
        if bond_L_over_b is None:
            bond_um = flag_raw.get("bond_L_um")
            if bond_um is not None:
                bond_L_over_b = float(bond_um) / max(scale.b_um, 1e-12)
            else:
                bond_L_over_b = float(ds_over_b)

        length_over_b = flag_raw.get("length_over_b")
        if length_over_b is None:
            length_um = flag_raw.get("length_um")
            if length_um is not None:
                length_over_b = float(length_um) / max(scale.b_um, 1e-12)
            else:
                length_over_b = 5.8

        helix_raw = flag_raw.get("helix_init", {}) or {}
        radius_over_b_h = helix_raw.get("radius_over_b")
        if radius_over_b_h is None:
            old_r_um = flag_raw.get("radius_um")
            radius_over_b_h = (
                float(old_r_um) / max(scale.b_um, 1e-12)
                if old_r_um is not None
                else 0.25
            )
        pitch_over_b_h = helix_raw.get("pitch_over_b")
        if pitch_over_b_h is None:
            old_pitch_um = flag_raw.get("pitch_um")
            pitch_over_b_h = (
                float(old_pitch_um) / max(scale.b_um, 1e-12)
                if old_pitch_um is not None
                else 2.5
            )

        flagella = FlagellumParams(
            n_flagella=int(_get(flag_raw, "n_flagella", 3)),
            placement_mode=str(_get(flag_raw, "placement_mode", "uniform")),
            discretization=FlagellaDiscretizationParams(ds_over_b=float(ds_over_b)),
            bond_L_over_b=float(bond_L_over_b),
            length_over_b=float(length_over_b),
            helix_init=FlagellaHelixInitParams(
                radius_over_b=float(radius_over_b_h),
                pitch_over_b=float(pitch_over_b_h),
            ),
        )

        fluid_raw = raw.get("fluid", {}) or {}
        fluid = FluidParams(
            viscosity_Pa_s=float(_get(fluid_raw, "viscosity_Pa_s", 1e-3))
        )

        motor_raw = raw.get("motor", {}) or {}
        motor = MotorParams(
            torque_Nm=float(_get(motor_raw, "torque_Nm", 4e-18)),
            reverse_n_flagella=int(_get(motor_raw, "reverse_n_flagella", 1)),
            enable_switching=bool(_get(motor_raw, "enable_switching", False)),
        )

        thermal = K_B * max(brownian.temperature_K, 1e-9)
        b_m = max(scale.b_um, 1e-9) * 1e-6

        spring_raw = (raw.get("potentials", {}) or {}).get("spring", {}) or {}
        h_over = spring_raw.get("H_over_T_over_b")
        if h_over is None and "H" in spring_raw:
            h_over = float(spring_raw["H"]) * b_m / max(thermal, 1e-30)
        s_val = spring_raw.get("s")
        if s_val is None and "s_um" in spring_raw:
            s_val = float(spring_raw["s_um"]) / max(scale.b_um, 1e-12)
        spring = SpringPotentialParams(
            H_over_T_over_b=float(h_over if h_over is not None else 10.0),
            s=float(s_val if s_val is not None else 0.1),
        )

        bend_raw = (raw.get("potentials", {}) or {}).get("bend", {}) or {}
        kb_over = bend_raw.get("kb_over_T")
        if kb_over is None and "kb" in bend_raw:
            kb_over = float(bend_raw["kb"]) / max(thermal, 1e-30)
        bend = BendPotentialParams(
            kb_over_T=float(kb_over if kb_over is not None else 20.0),
            theta0_deg=_state_angles(
                bend_raw.get("theta0_deg"),
                {"normal": 142.0, "semicoiled": 90.0, "curly1": 105.0},
            ),
        )

        torsion_raw = (raw.get("potentials", {}) or {}).get("torsion", {}) or {}
        kt_over = torsion_raw.get("kt_over_T")
        if kt_over is None and "kt" in torsion_raw:
            kt_over = float(torsion_raw["kt"]) / max(thermal, 1e-30)
        torsion = TorsionPotentialParams(
            kt_over_T=float(kt_over if kt_over is not None else 10.0),
            phi0_deg=_state_angles(
                torsion_raw.get("phi0_deg"),
                {"normal": -60.0, "semicoiled": 65.0, "curly1": 120.0},
            ),
        )

        rep_raw = (raw.get("potentials", {}) or {}).get(
            "spring_spring_repulsion", {}
        ) or {}
        a_over = rep_raw.get("A_ss_over_T")
        if a_over is None and "A_ss" in rep_raw:
            a_over = float(rep_raw["A_ss"]) / max(thermal, 1e-30)
        a_len = rep_raw.get("a_ss_over_b")
        if a_len is None and "a_ss_um" in rep_raw:
            a_len = float(rep_raw["a_ss_um"]) / max(scale.b_um, 1e-12)
        cutoff = rep_raw.get("cutoff_over_b")
        if cutoff is None and "cutoff_um" in rep_raw:
            cutoff = float(rep_raw["cutoff_um"]) / max(scale.b_um, 1e-12)

        repulsion = SpringSpringRepulsionParams(
            A_ss_over_T=float(a_over if a_over is not None else 1.0),
            a_ss_over_b=float(a_len if a_len is not None else 0.2),
            cutoff_over_b=float(cutoff if cutoff is not None else 0.2),
        )

        potentials = PotentialsParams(
            spring=spring,
            bend=bend,
            torsion=torsion,
            spring_spring_repulsion=repulsion,
        )

        hook_raw = raw.get("hook", {}) or {}
        kb_hook_over = hook_raw.get("kb_over_T")
        if kb_hook_over is None and "kb" in hook_raw:
            kb_hook_over = float(hook_raw["kb"]) / max(thermal, 1e-30)
        hook = HookParams(
            enabled=bool(_get(hook_raw, "enabled", True)),
            threshold_deg=float(_get(hook_raw, "threshold_deg", 90.0)),
            kb_over_T=float(kb_hook_over if kb_hook_over is not None else 20.0),
        )

        rt_raw = raw.get("run_tumble", {}) or {}
        run_tumble = RunTumbleParams(
            run_tau=float(_get(rt_raw, "run_tau", 1200.0)),
            tumble_tau=float(_get(rt_raw, "tumble_tau", 800.0)),
            semicoiled_tau=float(_get(rt_raw, "semicoiled_tau", 400.0)),
            curly1_tau=float(_get(rt_raw, "curly1_tau", 400.0)),
        )

        time_raw = raw.get("time", {}) or {}
        if "dt_over_tau" in time_raw:
            raise ValueError(
                "time.dt_over_tau は入力キーとして廃止。time.dt_s を使用してください。"
            )

        dt_s = time_raw.get("dt_s")
        if dt_s is None:
            raise ValueError(
                "time.dt_s は必須です。出力間隔として設定してください（内部計算のΔtは1e-3固定）。"
            )

        time = TimeParams(
            duration_s=float(_get(time_raw, "duration_s", 0.1)),
            dt_s=float(dt_s),
        )

        out_sample_raw = raw.get("output_sampling", {}) or {}
        old_fps = (raw.get("time", {}) or {}).get("fps_out")
        output_sampling = OutputSamplingParams(
            out_all_steps_3d=bool(_get(out_sample_raw, "out_all_steps_3d", True)),
            fps_out_2d=float(_get(out_sample_raw, "fps_out_2d", old_fps or 25.0)),
        )

        render_raw = raw.get("render", {}) or {}
        render = RenderParams(
            image_size_px=int(_get(render_raw, "image_size_px", 256)),
            pixel_size_um=float(_get(render_raw, "pixel_size_um", 0.203)),
            flagella_linewidth_px=float(_get(render_raw, "flagella_linewidth_px", 2.0)),
            render_flagella=bool(_get(render_raw, "render_flagella", True)),
            render_flagella_2d=bool(_get(render_raw, "render_flagella_2d", False)),
            save_frames_3d=bool(_get(render_raw, "save_frames_3d", True)),
            follow_camera_3d=bool(_get(render_raw, "follow_camera_3d", True)),
            view_range_um=float(_get(render_raw, "view_range_um", 8.0)),
            timestamp_3d=bool(_get(render_raw, "timestamp_3d", True)),
            timestamp_fmt=str(_get(render_raw, "timestamp_fmt", "t = {t:.3f} s")),
            label_flagella=bool(_get(render_raw, "label_flagella", True)),
            follow_camera_2d=bool(_get(render_raw, "follow_camera_2d", False)),
            center_body_in_2d=bool(_get(render_raw, "center_body_in_2d", True)),
            save_frames_2d=bool(_get(render_raw, "save_frames_2d", True)),
        )

        seed_raw = raw.get("seed", {}) or {}
        seed = SeedParams(global_seed=int(_get(seed_raw, "global_seed", 0)))

        output_raw = raw.get("output", {}) or {}
        output = OutputParams(base_dir=str(_get(output_raw, "base_dir", "outputs")))

        return SimulationConfig(
            scale=scale,
            body=body,
            flagella=flagella,
            fluid=fluid,
            motor=motor,
            potentials=potentials,
            hook=hook,
            run_tumble=run_tumble,
            time=time,
            output_sampling=output_sampling,
            brownian=brownian,
            projection=projection,
            stiffness_scales=stiffness_scales,
            local_helix=local_helix,
            diagnostics=diagnostics,
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
