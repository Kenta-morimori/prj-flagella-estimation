"""Bounded inspection CLI for a selected window of ``step_summary.csv``."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from sim_swim.analysis.run_summary import resolve_run_dir


def _episode_window(run_dir: Path, gate: str, episode: int) -> tuple[float, float]:
    summary_path = run_dir / "run_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError("run_summary.json is required for --gate/--episode")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    episodes = summary.get("gates", {}).get(gate, {}).get("episodes", [])
    if episode < 1 or episode > len(episodes):
        raise ValueError(
            f"episode must be between 1 and {len(episodes)} for gate={gate}"
        )
    selected = episodes[episode - 1]
    start, end = selected.get("start_t_s"), selected.get("end_t_s")
    if start is None or end is None:
        raise ValueError("selected episode has no finite time bounds")
    return float(start), float(end)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
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
    if args.gate is not None:
        start, end = _episode_window(
            resolve_run_dir(args.input_dir), args.gate, args.episode
        )
    elif args.start_t_s is not None and args.end_t_s is not None:
        start, end = args.start_t_s, args.end_t_s
    else:
        parser.error("provide --gate/--episode or both --start-t-s and --end-t-s")
    if end < start or not 1 <= args.max_rows <= 1000:
        parser.error("time window must be ordered and --max-rows must be 1..1000")
    columns = [value.strip() for value in args.columns.split(",") if value.strip()]
    if not columns or len(columns) > 12:
        parser.error("--columns must contain 1..12 columns")
    run_dir = resolve_run_dir(args.input_dir)
    selected: list[dict[str, str]] = []
    match_count = 0
    with (run_dir / "step_summary.csv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle)
        available = set(reader.fieldnames or [])
        missing = [column for column in columns if column not in available]
        if missing:
            parser.error(f"unknown columns: {', '.join(missing)}")
        for row in reader:
            try:
                t_s = float(row["t_s"])
            except (KeyError, ValueError):
                continue
            if start <= t_s <= end:
                match_count += 1
                if len(selected) < args.max_rows:
                    selected.append({column: row[column] for column in columns})
    print(
        json.dumps(
            {
                "run_dir": str(run_dir),
                "time_window_s": [start, end],
                "columns": columns,
                "matching_row_count": match_count,
                "returned_row_count": len(selected),
                "truncated": match_count > len(selected),
                "rows": selected,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
