"""Bounded, every-step aggregation for compact Phase 2 outputs."""

from __future__ import annotations

import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping
from zoneinfo import ZoneInfo


class OnlineRunSummary:
    """Aggregate diagnostic rows without retaining their time series."""

    def __init__(self, *, expected_steps: int) -> None:
        self.expected_steps = expected_steps
        self.count = 0
        self.last_t_s: float | None = None
        self.extrema: dict[str, dict[str, float | None]] = {}
        self.gates = {key: self._gate() for key in ("finite", "shape_nonbody")}

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
        self._record_gate(
            "shape_nonbody",
            bool(row.get("shape_pass_nonbody", False)),
            t_s,
            str(row.get("first_fail_category_nonbody") or "shape_nonbody"),
            str(row.get("first_fail_flag_id_nonbody") or "nonbody"),
        )

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
        body_unavailable = {
            "status": "unavailable",
            "reason": "not recorded by compact output policy",
        }
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
            "gates": {**self.gates, "shape_body": body_unavailable},
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
