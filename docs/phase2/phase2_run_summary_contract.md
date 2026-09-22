# Phase 2 `run_summary.json` contract

## Purpose

`run_summary.json` is the compact, formal first-read artifact for one Phase 2 simulation or campaign condition. It aggregates existing diagnostics and does not change the physical model, gate definitions, thresholds, or `step_summary.csv`.

Issue #186 の `output.policy: compact` では、全内部 step の診断・strict QC をオンライン集約する一方、全step CSV は書かない。`all_step_metrics` は各数値列の min/max/final/finite、gate は first failure の時刻・category・target を持つ。compact archive は物理時間一様で、標準 `output.archive_interval_s: 0.001`（0.5 s なら約501 state）である。これは 25--1000 fps replay、約100 Hz 回転の1周期約10点、`dt_star` 間の同じ物理時間解像度を両立するためである。archive より高い fps は補間せず拒否する。未定義の将来の全step指標は compact archive だけから完全再構成できないため、必要時は短時間 debug policy を使う。

## compact heartbeat と partial evidence

`generic_multi_run` が compact condition を実行すると、既定で
`output.checkpoint_interval_steps: 2500` ごと（および終端時）に次を同じ
condition directoryへ原子的に更新する。これは resume 用stateではなく、停止・例外時にも
定量・定性診断に残す評価証跡である。

- `progress.json`: status、completed/total internal steps、`t_star`、`t_s`、wall time、steps/s、online QC extrema、first failure。
- `diagnostic_samples.csv`: checkpoint境界のraw diagnostic/body diagnostic行（各列は`diagnostic_` / `body_`接頭辞で衝突なく保存）。平均・標準偏差・分位点・移動窓・任意の`tau`範囲は、ここから後処理で計算する。onlineの平均値を正本にしてはならない。
- `state_archive.partial.npz` と `trajectory.partial.csv`: checkpoint時点までのreplay入力。最終の`state_archive.npz` / `trajectory.csv`とは区別する。

SIGTERM/SIGINTは次のinternal step境界で安全停止し、`run_summary.json`の
`execution.status=partial`、`performance.json`、conditionのpartial evidence、campaignの
`campaign_completion.json`を残してexit code 130で終了する。この状態では`summary.csv`と
`run_manifest.json`を作らず、parallel aggregateも生成しない。正常完了時は既存final artifact
contractを維持し、`progress.json.status=completed`へ確定する。debug policyは全step CSV互換を
維持し、heartbeat artifactを追加しない。

partial evidenceをreplayするには、先に
`python -m sim_swim.analysis.partial_generic_multi_run --include-partial-checkpoint`で
analysis-only inputを作り、さらにreplay側で`--allow-partial`を明示する。映像には`PARTIAL`を
記録し、dataset採択、profile昇格、canonical判定の入力に使用しない。

partial evidenceは通常の3D+2D replayだけを対象とする。hydrodynamics flow archiveは正常完了時だけ
`hydro_archive.npz`として保存し、checkpointでは保存しない。そのため`--flow-overlay` / hydrodynamics
replayはpartial manifestを明示的に拒否する。partial flow replayが必要になった場合は、保存量・QC・
provenanceを別Issueで設計する。

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
