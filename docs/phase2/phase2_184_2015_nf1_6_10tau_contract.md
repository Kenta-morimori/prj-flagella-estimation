# Phase 2 Issue #184: 2015 project nf1–6 10τ stability screen

cs10の起動・同期・解析手順は`docs/phase2/phase2_184_2015_nf1_6_cs10_runbook.md`を正本とする。

2015 projectの`2.5e-20 N m` per-flagellum torqueにおける、長時間の形状・遊泳診断である。
これはdiagnostic-only campaignである。dataset採択、canonical torque選定、supported profile昇格、
Phase 3 handoffは行わない。#61のpitch QCの採否見直しも本campaignでは行わない。

| item | value |
| --- | --- |
| profile | `conf/sim_swim_2015.yaml` project |
| torque / scale policy | `2.5e-20 N m`; motor = reference = force torque; `reference_torque` |
| conditions | `n_flagella=1,2,3,4,5,6`; `seeded_surface`; attach/phase seed `0` |
| integration | `dt_star=1e-5`, `duration_tau=10` |
| motion | RUN fixed、switchingなし、Brownian OFF |
| execution | cs10、`cs10_qualified` 3 worker、6 isolated shards |

開始前に、再解析済み#61の3 torque・1τ decision JSONとsummary CSVが揃っていなければならない。
現行の3/3 strict FAIL（pitch / motor torque residual）は既知の診断前提としてmanifestへ保存する。
jobは同一campaign root・3 condition・全strict FAILをauditし、不在・不整合・別campaignならsimulationを起動しない。

`seeded_surface`はseed 0で決定論的に中心層を優先する。nf1–3は中心層の先頭slotから、nf4–6は中心層の6 slotを順に使う。各conditionについてstrict QC、最初のfailure criterion / 時刻 / step、wall time、steps/s、
body/flagella motion、3d+2d replayを保存する。見積りは約5.4日/condition、3 workerの二波で約11日である。

2026-09-08の`seeded_center_layer` jobはnf4/nf5が`n_prism=6`を割り切らず開始前geometry構成で失敗した。
nf1–3と完走するnf6は保持するが、clean seeded-surface campaignのaggregate、定量比較、profile判断には使わない。
以後のdry-runは全conditionを`ModelBuilder`でgeometry-only検証し、attachment topologyをjob manifestに記録する。

## User-run procedure

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/issue184_2015_nf1_6_10tau/job.yaml --dry-run
.venv-cs10/bin/python scripts/cs10/queue.py enqueue \
  --job-yaml conf/phase2_parallel/issue184_2015_nf1_6_10tau/job.yaml
```

1 reservationが6 shardを3 workerで二波に分け、全conditionとaggregateが終端後にActions通知を一度だけ送る。
停止・失敗時は完了artifactをコピー・混在・再利用せず、別のdated output rootでclean jobとして扱う。
完了後にmanifest、summary、compact diagnostics、state archive、trajectory、replayをlocal同期してSHA-256とQCを照合する。

集約済みclean campaignだけを次で定量評価する。partial / failed shardや旧`seeded_center_layer` outputは入力にできない。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py \
  --analysis-kind issue184-2015-nf1-6 \
  --run-root <campaign-root> \
  --output-dir <campaign-root>/analysis/issue184_2015_nf1_6
```
