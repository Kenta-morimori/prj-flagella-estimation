"""Issue #163 spring formulation comparison artifacts and default selection."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sim_swim.dynamics.forces import compute_spring_forces

FORMULATIONS = ("legacy", "fene_fraenkel")
FORMULATION_COLUMN = "axis_formulation_value"


def force_extension_rows(
    *,
    rest_lengths_over_b: Iterable[float] = (0.58, 1.0),
    h_over_t_over_b: float = 10.0,
    s: float = 0.1,
    samples: int = 199,
    torque_nm: float = 2.5e-20,
    b_um: float = 1.0,
) -> list[dict[str, float | str]]:
    """Return nondimensional force-extension curves for both formulations.

    Values are produced through the production ``compute_spring_forces`` API.
    Positive values pull the left bead right under extension; negative values
    push it left under compression. The sample range stays just inside both
    singular boundaries.
    """

    if not math.isfinite(h_over_t_over_b) or h_over_t_over_b <= 0.0:
        raise ValueError("h_over_t_over_b must be finite and positive")
    if not math.isfinite(s) or not 0.0 < s < 1.0:
        raise ValueError("s must satisfy 0 < s < 1")
    if samples < 3:
        raise ValueError("samples must be at least 3")
    if not math.isfinite(torque_nm) or torque_nm <= 0.0:
        raise ValueError("torque_nm must be finite and positive")
    if not math.isfinite(b_um) or b_um <= 0.0:
        raise ValueError("b_um must be finite and positive")

    rest_length_values = [float(value) for value in rest_lengths_over_b]
    if not rest_length_values:
        raise ValueError("rest_lengths_over_b must be non-empty")

    rows: list[dict[str, float | str]] = []
    b_m = b_um * 1.0e-6
    force_scale = torque_nm / b_m
    h_const = h_over_t_over_b * force_scale
    spring_pairs = np.array([[0, 1]], dtype=int)
    for rest in rest_length_values:
        if not math.isfinite(rest) or rest <= 0.0:
            raise ValueError("rest lengths must be finite and positive")
        relative_extensions = np.linspace(-0.99 * s, 0.99 * s, samples)
        for formulation in FORMULATIONS:
            for q in relative_extensions:
                delta_over_b = float(q) * rest
                distance_m = (rest + delta_over_b) * b_m
                forces = compute_spring_forces(
                    np.array([[0.0, 0.0, 0.0], [distance_m, 0.0, 0.0]]),
                    spring_pairs,
                    np.array([rest * b_m]),
                    h_const=h_const,
                    s=s,
                    b_m=b_m,
                    formulation=formulation,
                )
                force = float(forces[0, 0])
                rows.append(
                    {
                        "formulation": formulation,
                        "rest_length_over_b": rest,
                        "relative_extension_q": float(q),
                        "extension_over_b": delta_over_b,
                        "distance_over_b": rest + delta_over_b,
                        "force_N": force,
                        "force_over_t_over_b": force / force_scale,
                    }
                )
    return rows


def write_force_extension_artifacts(
    output_dir: Path,
    *,
    rest_lengths_over_b: Iterable[float] = (0.58, 1.0),
    h_over_t_over_b: float = 10.0,
    s: float = 0.1,
    samples: int = 199,
) -> tuple[Path, Path]:
    rows = force_extension_rows(
        rest_lengths_over_b=rest_lengths_over_b,
        h_over_t_over_b=h_over_t_over_b,
        s=s,
        samples=samples,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "force_extension.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    rest_values = list(dict.fromkeys(float(row["rest_length_over_b"]) for row in rows))
    figure, axes_grid = plt.subplots(
        1,
        len(rest_values),
        figsize=(5 * len(rest_values), 4),
        sharey=False,
        squeeze=False,
    )
    axes = axes_grid[0]
    for axis, rest in zip(axes, rest_values, strict=True):
        for formulation in FORMULATIONS:
            selected = [
                row
                for row in rows
                if row["formulation"] == formulation
                and row["rest_length_over_b"] == rest
            ]
            axis.plot(
                [float(row["relative_extension_q"]) for row in selected],
                [float(row["force_N"]) for row in selected],
                label=formulation,
            )
        axis.axhline(0.0, color="black", linewidth=0.7)
        axis.axvline(0.0, color="black", linewidth=0.7)
        axis.set_title(f"L/b = {rest:g}")
        axis.set_xlabel("relative extension q = (d-L)/L")
        axis.set_ylabel("production spring force on left bead [N]")
        axis.grid(alpha=0.25)
    axes[0].legend()
    figure.tight_layout()
    png_path = output_dir / "force_extension.png"
    figure.savefig(png_path, dpi=160)
    plt.close(figure)
    return csv_path, png_path


def _read_summary(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"summary.csv is empty: {path}")
    return rows


def _parse_bool(value: str | None) -> bool | None:
    text = str(value or "").strip().lower()
    if text in {"1", "true", "yes", "pass"}:
        return True
    if text in {"0", "false", "no", "fail"}:
        return False
    return None


def _formulation(row: dict[str, str]) -> str:
    for key in (FORMULATION_COLUMN, "axis_formulation_label", "formulation"):
        value = str(row.get(key, "")).strip()
        if value:
            return value
    return ""


def _evaluate_row(row: dict[str, str] | None) -> dict[str, Any]:
    reasons: list[str] = []
    if row is None:
        return {"pass": False, "reasons": ["summary row is missing"]}

    nonbody = _parse_bool(row.get("final_shape_pass_nonbody"))
    body = _parse_bool(row.get("body_shape_pass"))
    if nonbody is not True:
        reasons.append("final_shape_pass_nonbody is not true")
    if body is not True:
        reasons.append("body_shape_pass is not true")
    if str(row.get("first_fail_t_s", "")).strip():
        reasons.append("first_fail_t_s is present")

    final_t_s: float | None = None
    duration_s: float | None = None
    for field in ("duration_s", "final_t_s"):
        value = str(row.get(field, "")).strip()
        if not value:
            reasons.append(f"{field} is missing")
            continue
        try:
            parsed = float(value)
        except ValueError:
            reasons.append(f"{field} is not numeric")
            continue
        if not math.isfinite(parsed):
            reasons.append(f"{field} is non-finite")
            continue
        if field == "duration_s":
            duration_s = parsed
        else:
            final_t_s = parsed
    if (
        final_t_s is not None
        and duration_s is not None
        and final_t_s < duration_s - max(1.0e-12, duration_s * 1.0e-9)
    ):
        reasons.append("final_t_s did not reach duration_s")

    finite_metric_fields = {
        "max_flag_bond_rel_err": ("max_flag_bond_rel_err",),
        "max_hook_len_rel_err": (
            "max_hook_len_rel_err",
            "hook_len_rel_err_max",
        ),
    }
    for field, aliases in finite_metric_fields.items():
        value = next(
            (
                str(row.get(alias, "")).strip()
                for alias in aliases
                if str(row.get(alias, "")).strip()
            ),
            "",
        )
        if value:
            try:
                if not math.isfinite(float(value)):
                    reasons.append(f"{field} is non-finite")
            except ValueError:
                reasons.append(f"{field} is not numeric")
    return {"pass": not reasons, "reasons": reasons}


def decide_default(
    motor_off_summary: Path,
    motor_on_summary: Path,
) -> dict[str, Any]:
    """Apply the accepted #163 default-selection rule to two campaign summaries."""

    evaluations: dict[str, dict[str, Any]] = {name: {} for name in FORMULATIONS}
    inputs = {
        "motor_off": motor_off_summary,
        "motor_on": motor_on_summary,
    }
    for run_name, path in inputs.items():
        rows = _read_summary(path)
        by_formulation: dict[str, dict[str, str]] = {}
        for row in rows:
            formulation = _formulation(row)
            if formulation in by_formulation:
                raise ValueError(f"duplicate {formulation} row in {path}")
            by_formulation[formulation] = row
        for formulation in FORMULATIONS:
            evaluations[formulation][run_name] = _evaluate_row(
                by_formulation.get(formulation)
            )

    for formulation in FORMULATIONS:
        evaluations[formulation]["overall_pass"] = all(
            evaluations[formulation][run_name]["pass"] for run_name in inputs
        )
    fene_pass = evaluations["fene_fraenkel"]["overall_pass"]
    selected = "fene_fraenkel" if fene_pass else "legacy"
    reason = (
        "fene_fraenkel passed strict body and non-body gates in both runs"
        if fene_pass
        else "fene_fraenkel did not pass both runs; retain legacy compatibility"
    )
    return {
        "decision_rule": "prefer fene_fraenkel only when it passes both runs",
        "selected_default": selected,
        "reason": reason,
        "inputs": {name: str(path) for name, path in inputs.items()},
        "evaluations": evaluations,
    }


def write_decision_artifacts(
    output_dir: Path,
    *,
    motor_off_summary: Path,
    motor_on_summary: Path,
) -> tuple[Path, Path, dict[str, Any]]:
    decision = decide_default(motor_off_summary, motor_on_summary)
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "default_decision.json"
    json_path.write_text(
        json.dumps(decision, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    lines = [
        "# Spring formulation default decision",
        "",
        f"- Selected default: `{decision['selected_default']}`",
        f"- Rule: {decision['decision_rule']}",
        f"- Reason: {decision['reason']}",
        "",
        "| formulation | motor-off | motor-on | overall |",
        "|---|---:|---:|---:|",
    ]
    for formulation in FORMULATIONS:
        evaluation = decision["evaluations"][formulation]
        lines.append(
            f"| `{formulation}` | {evaluation['motor_off']['pass']} | "
            f"{evaluation['motor_on']['pass']} | {evaluation['overall_pass']} |"
        )
        for run_name in ("motor_off", "motor_on"):
            for reason in evaluation[run_name]["reasons"]:
                lines.append(f"\n- `{formulation}` / `{run_name}`: {reason}")
    markdown_path = output_dir / "default_decision.md"
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, markdown_path, decision
