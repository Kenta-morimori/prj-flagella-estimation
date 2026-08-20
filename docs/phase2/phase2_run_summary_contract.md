# Phase 2 `run_summary.json` contract

## Purpose

`run_summary.json` is the compact, formal first-read artifact for one Phase 2 simulation or campaign condition. It aggregates existing diagnostics and does not change the physical model, gate definitions, thresholds, or `step_summary.csv`.

Issue #186 の `output.policy: compact` では、全内部 step の診断・strict QC をオンライン集約する一方、全step CSV は書かない。`all_step_metrics` は各数値列の min/max/final/finite、gate は first failure の時刻・category・target を持つ。compact archive は物理時間一様で、標準 `output.archive_interval_s: 0.001`（0.5 s なら約501 state）である。これは 25--1000 fps replay、約100 Hz 回転の1周期約10点、`dt_star` 間の同じ物理時間解像度を両立するためである。archive より高い fps は補間せず拒否する。未定義の将来の全step指標は compact archive だけから完全再構成できないため、必要時は短時間 debug policy を使う。

## Location and reading order

- Single simulation: `<run>/sim/run_summary.json`
- Sweep / multi-run: `<campaign>/<condition>/run_summary.json`

Read `manifest.json` and `run_summary.json` first. Inspect raw diagnostics only through the bounded CLI below; do not load a complete `step_summary.csv` into Codex context during routine analysis.

```bash
uv run python scripts/03_dataset_building/inspect_run.py \
  --run-dir outputs/.../condition_001 \
  --gate shape_nonbody --episode 1 \
  --columns t_s,shape_pass_nonbody,first_fail_category_nonbody,hook_len_rel_err_max
```

The CLI requires a time window or a `gate`/`episode`, requires 1--12 columns, and returns at most 1,000 rows (100 by default).

## Schema version 1.0

Machine-readable top-level validation contract: `schemas/phase2_run_summary.schema.json`.

- `execution`: `completed`, `partial`, or `unknown`. `unknown` means an older output lacks enough manifest time metadata; it does not mean pass or fail.
- `sampling`: observed temporal spacing and the episode policy. Durations are physical seconds but are bounded by observed samples; they do not claim an unobserved threshold-crossing time.
- `gates.finite` and `gates.shape_nonbody`: existing per-step gate values, without reinterpretation.
- `gates.shape_body`: `available` only when body diagnostics exist. Missing diagnostics are `unavailable`, never treated as a body failure or pass.
- `episodes`: maximal consecutive observed fail samples. `start_t_s` / `end_t_s` are the first / last observed fail times; `next_observed_pass_t_s` records recovery when observed. `observed_duration_s` is therefore a sampled span, not an unobserved threshold-crossing duration. `persistent_observed` means at least three consecutive observations, not a new physical acceptance threshold.
- `extrema`: selected existing diagnostic maxima and their observed times.

Each gate stores at most 32 episodes. When that limit is exceeded, it retains early, late, and long episodes and records the omitted count. No per-step time series is copied into JSON. The generated file is capped at 64 KiB; generation fails rather than silently emitting an unbounded artifact.

## Generation and compatibility

New simulations generate the file automatically after their diagnostics are written. Existing outputs can be summarized without re-simulation:

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind run-summary \
  --input-dir outputs/.../condition_001
```

The command does not modify source diagnostics. Replacing an existing summary requires `--overwrite`.
