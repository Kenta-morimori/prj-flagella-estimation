"""Bounded, every-step aggregation for compact Phase 2 outputs."""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping
from zoneinfo import ZoneInfo

from sim_swim.sim.body_shape_gate import (
    BODY_BEND_ERR_MAX_DEG_LIMIT,
    BODY_CENTERLINE_DEV_UM_MAX_LIMIT,
    BODY_SPRING_STRETCH_RATIO_MAX_LIMIT,
    BODY_TRIANGLE_AREA_RATIO_MIN_LIMIT,
)


class OnlineRunSummary:
    """Aggregate diagnostic rows without retaining their time series."""

    def __init__(self, *, expected_steps: int) -> None:
        self.expected_steps = expected_steps
        self.count = 0
        self.last_t_s: float | None = None
        self.extrema: dict[str, dict[str, float | None]] = {}
        self.gates = {
            key: self._gate() for key in ("finite", "shape_nonbody", "shape_body")
        }
        self._initial_body_area: float | None = None

    @staticmethod
    def _gate() -> dict[str, Any]:
        return {
            "status": "available",
            "final_pass": None,
            "any_fail": False,
            "first_observed_fail_t_s": None,
            "last_observed_fail_t_s": None,
            "observed_fail_sample_count": 0,
            "first_failure_category": None,
            "first_failure_target": None,
        }

    def record(self, row: Mapping[str, Any]) -> None:
        self.count += 1
        t_s = float(row.get("t_s", float("nan")))
        self.last_t_s = t_s if math.isfinite(t_s) else self.last_t_s
        for key, value in row.items():
            if isinstance(value, bool):
                continue
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            item = self.extrema.setdefault(
                key, {"min": None, "max": None, "final": None, "finite": True}
            )
            item["final"] = number if math.isfinite(number) else None
            item["finite"] = bool(item["finite"]) and math.isfinite(number)
            if math.isfinite(number):
                item["min"] = (
                    number if item["min"] is None else min(float(item["min"]), number)
                )
                item["max"] = (
                    number if item["max"] is None else max(float(item["max"]), number)
                )
        self._record_gate(
            "finite", bool(row.get("finite_pass", False)), t_s, "finite", "all"
        )

    def record_body(self, row: Mapping[str, Any]) -> None:
        """Apply the existing body gate to every step without emitting a CSV."""

        t_s = float(row.get("t_s", float("nan")))
        area = float(row.get("body_triangle_area_min", float("nan")))
        if self._initial_body_area is None and math.isfinite(area) and area > 0.0:
            self._initial_body_area = area
        stretch = float(row.get("body_spring_max_stretch_ratio", float("nan")))
        bend = float(row.get("body_bend_max_error_deg", float("nan")))
        center = float(row.get("body_centerline_max_deviation_um", float("nan")))
        ratio = (
            area / self._initial_body_area if self._initial_body_area else float("nan")
        )
        values = (stretch, bend, center, ratio)
        if not all(math.isfinite(value) for value in values):
            passed, category = False, "body_nonfinite"
        elif stretch > BODY_SPRING_STRETCH_RATIO_MAX_LIMIT:
            passed, category = False, "body_spring"
        elif bend > BODY_BEND_ERR_MAX_DEG_LIMIT:
            passed, category = False, "body_bend"
        elif center > BODY_CENTERLINE_DEV_UM_MAX_LIMIT:
            passed, category = False, "body_centerline"
        elif ratio < BODY_TRIANGLE_AREA_RATIO_MIN_LIMIT:
            passed, category = False, "body_area"
        else:
            passed, category = True, "none"
        for key, value in row.items():
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            item = self.extrema.setdefault(
                key, {"min": None, "max": None, "final": None, "finite": True}
            )
            item["final"] = number if math.isfinite(number) else None
            item["finite"] = bool(item["finite"]) and math.isfinite(number)
            if math.isfinite(number):
                item["min"] = (
                    number if item["min"] is None else min(float(item["min"]), number)
                )
                item["max"] = (
                    number if item["max"] is None else max(float(item["max"]), number)
                )
        self._record_gate("shape_body", passed, t_s, category, "body")

    def _record_gate(
        self, name: str, passed: bool, t_s: float, category: str, target: str
    ) -> None:
        gate = self.gates[name]
        gate["final_pass"] = passed
        if not passed:
            gate["any_fail"] = True
            gate["observed_fail_sample_count"] += 1
            gate["last_observed_fail_t_s"] = t_s
            if gate["first_observed_fail_t_s"] is None:
                gate["first_observed_fail_t_s"] = t_s
                gate["first_failure_category"] = category
                gate["first_failure_target"] = target

    def document(
        self,
        *,
        run_dir: Path,
        completed: bool,
        reason: str,
        time_manifest: Mapping[str, Any],
        policy: str,
    ) -> dict[str, Any]:
        return {
            "schema_version": "1.0",
            "kind": "phase2_run_summary",
            "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
            "input": {
                "run_dir": str(run_dir),
                "step_summary_csv": None,
                "body_constraint_diagnostics_csv": None,
            },
            "execution": {
                "status": "completed" if completed else "partial",
                "row_count": self.count,
                "step_indices_contiguous_from_zero": True,
                "observed_first_t_s": 0.0 if self.count else None,
                "observed_final_t_s": self.last_t_s,
                "expected_total_steps": self.expected_steps,
                "expected_final_step_summary_t_s": time_manifest.get(
                    "final_step_summary_t_s"
                ),
                "reason": reason,
            },
            "sampling": {
                "step_summary_row_count": 0,
                "time_spacing_s": {"min_s": None, "median_s": None, "max_s": None},
                "episode_definition": "every internal step aggregated online",
                "persistent_observed_min_consecutive_fail_samples": 1,
                "episode_storage_limit_per_gate": 1,
                "output_policy": policy,
            },
            "gates": self.gates,
            "all_step_metrics": self.extrema,
            "extrema": {
                key: {"value": val["max"], "t_s": None}
                for key, val in self.extrema.items()
            },
            "source_file_size_bytes": {
                "step_summary_csv": None,
                "body_constraint_diagnostics_csv": None,
            },
        }

    def write(
        self,
        *,
        run_dir: Path,
        completed: bool,
        reason: str,
        time_manifest: Mapping[str, Any],
        policy: str,
    ) -> Path:
        path = run_dir / "run_summary.json"
        path.write_text(
            json.dumps(
                self.document(
                    run_dir=run_dir,
                    completed=completed,
                    reason=reason,
                    time_manifest=time_manifest,
                    policy=policy,
                ),
                ensure_ascii=False,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        return path
