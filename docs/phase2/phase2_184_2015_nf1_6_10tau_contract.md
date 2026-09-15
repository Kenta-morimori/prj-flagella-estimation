# Phase 2 Issue #184: 2015 project nf1–6 10τ stability screen

cs10の起動・同期・解析手順は`docs/phase2/phase2_184_2015_nf1_6_cs10_runbook.md`を正本とする。

2015 projectの`2.5e-20 N m` per-flagellum torqueにおける、長時間の形状・遊泳診断である。
これはdiagnostic-only campaignである。dataset採択、canonical torque選定、supported profile昇格、
Phase 3 handoffは行わない。#61のpitch QCの採否見直しも本campaignでは行わない。

| item | value |
| --- | --- |
| profile | `conf/sim_swim_2015.yaml` project |
| torque / scale policy | `2.5e-20 N m`; motor = reference = force torque; `reference_torque` |
| conditions | 今回の補完は`n_flagella=4,5`のみ、`seeded_surface`、attach/phase seed `0`。nf1–3/nf6は既存`seeded_center_layer` artifactを保持する。 |
| integration | `dt_star=1e-5`, `duration_tau=10` |
| motion | RUN fixed、switchingなし、Brownian OFF |
| execution | 今回限定でcs10 direct 2 worker・nf4/nf5 isolated shards・Actions通知なし |

開始前に、再解析済み#61の3 torque・1τ decision JSONとsummary CSVが揃っていなければならない。
現行の3/3 strict FAIL（pitch / motor torque residual）は既知の診断前提としてmanifestへ保存する。
jobは同一campaign root・3 condition・全strict FAILをauditし、不在・不整合・別campaignならsimulationを起動しない。

`seeded_surface`はseed 0で決定論的に中心層を優先する。今回のnf4/nf5は中心層の先頭slotから構成する。各conditionについてstrict QC、最初のfailure criterion / 時刻 / step、wall time、steps/s、
body/flagella motion、3d+2d replayを保存する。見積りは約5.4日/condition、2 worker並列で約6日である。

2026-09-08の`seeded_center_layer` jobはnf4/nf5が`n_prism=6`を割り切らず開始前geometry構成で失敗した。
nf1–3と完走するnf6は保持する。今回だけは、nf4/nf5の`seeded_surface`補完結果と合わせたtopology混在の暫定横断表へ明示provenance付きで載せられるが、clean seeded-surface campaignのaggregate、dataset/profile判断、canonical選定には使わない。
以後のdry-runは全conditionを`ModelBuilder`でgeometry-only検証し、attachment topologyをjob manifestに記録する。

## 今回限定の direct procedure

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  --config conf/phase2_parallel/issue184_2015_nf4_nf5_10tau/job.yaml --dry-run
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  --config conf/phase2_parallel/issue184_2015_nf4_nf5_10tau/job.yaml \
  --max-workers 2 --output-root <new-nas-root>
```

このdirect jobはActions通知を送らない。停止済みnf1/nf2 direct rootは保持・隔離し、評価に使わない。nf4/nf5完走後、旧nf1–3/nf6と新nf4/nf5を暫定横断表へ記録する。
次回のclean nf1–6 screenは同一topology・single reservation・単発通知を必須とする。

集約済みclean campaignだけを次で正式評価する。partial / failed shardや旧`seeded_center_layer` outputは入力にできない。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py \
  --analysis-kind issue184-2015-nf1-6 \
  --run-root <campaign-root> \
  --output-dir <campaign-root>/analysis/issue184_2015_nf1_6
```

今回の暫定横断表は別CLIで明示したchild rootだけを入力にする。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py \
  --analysis-kind issue184-2015-nf1-6-provisional \
  --source nf01=<old-nf01-run-root> --source nf02=<old-nf02-run-root> \
  --source nf03=<old-nf03-run-root> --source nf04=<new-nf04-run-root> \
  --source nf05=<new-nf05-run-root> --source nf06=<old-nf06-run-root> \
  --output-dir <local-output>/analysis/issue184_provisional
```
