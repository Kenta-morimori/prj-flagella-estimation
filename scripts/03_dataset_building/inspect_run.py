#!/usr/bin/env python3
"""Inspect a bounded time interval in an existing run's step summary."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.run_summary import resolve_run_dir


def _episode_window(run_dir: Path, gate: str, episode: int) -> tuple[float, float]:
    summary = json.loads((run_dir / "run_summary.json").read_text(encoding="utf-8"))
    episodes = summary.get("gates", {}).get(gate, {}).get("episodes", [])
    if not 1 <= episode <= len(episodes):
        raise ValueError(
            f"episode must be between 1 and {len(episodes)} for gate={gate}"
        )
    selected = episodes[episode - 1]
    if selected.get("start_t_s") is None or selected.get("end_t_s") is None:
        raise ValueError("selected episode has no finite time bounds")
    return float(selected["start_t_s"]), float(selected["end_t_s"])


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--columns", required=True, help="Comma-separated columns; max 12."
    )
    parser.add_argument("--gate")
    parser.add_argument("--episode", type=int)
    parser.add_argument("--start-t-s", type=float)
    parser.add_argument("--end-t-s", type=float)
    parser.add_argument("--max-rows", type=int, default=100)
    args = parser.parse_args(argv)
    if (args.gate is None) != (args.episode is None):
        parser.error("--gate and --episode must be used together")
    run_dir = resolve_run_dir(args.run_dir)
    if args.gate:
        start, end = _episode_window(run_dir, args.gate, args.episode)
    elif args.start_t_s is not None and args.end_t_s is not None:
        start, end = args.start_t_s, args.end_t_s
    else:
        parser.error("provide --gate/--episode or both --start-t-s and --end-t-s")
    columns = [value.strip() for value in args.columns.split(",") if value.strip()]
    if end < start or not 1 <= args.max_rows <= 1000 or not 1 <= len(columns) <= 12:
        parser.error("invalid time window, row limit, or columns")
    with (run_dir / "step_summary.csv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        missing = [
            column for column in columns if column not in set(reader.fieldnames or [])
        ]
        if missing:
            parser.error(f"unknown columns: {', '.join(missing)}")
        rows = []
        matching = 0
        for row in reader:
            if start <= float(row["t_s"]) <= end:
                matching += 1
                if len(rows) < args.max_rows:
                    rows.append({column: row[column] for column in columns})
    print(
        json.dumps(
            {
                "run_dir": str(run_dir),
                "time_window_s": [start, end],
                "columns": columns,
                "matching_row_count": matching,
                "returned_row_count": len(rows),
                "truncated": matching > len(rows),
                "rows": rows,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
