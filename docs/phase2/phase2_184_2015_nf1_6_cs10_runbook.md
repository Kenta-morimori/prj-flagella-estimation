# Issue #184: 2015 project nf1–6 10τ cs10 runbook

このcampaignは`execution:cs10`のUser-run diagnostic-only実行である。6条件は
`seeded_surface`・attach/phase seed `0`に統一し、3 workerの二波で実行する。
全shardとaggregateが終端になったときだけ、queue reservation単位でActions通知を一度送る。

## 現行failed job

2026-09-08の`seeded_center_layer` reservationはnf4/nf5が`n_prism=6`の割り切れ制約で
開始前に失敗した。nf1–3と完走するnf6は保持するが、clean campaignへコピー・再利用・比較しない。
このreservationは停止しない。

## 実行前

cs10のclean worktreeで、以下のdry-runにより6 condition、worker数、分離output、
tracking-reference、seed、およびattachment body bead / layer / slotを確認する。
geometry preflightが失敗したらqueueへ投入しない。

```bash
export PATH=/home/people/Ktakemori/.local/bin:/usr/local/bin:/usr/bin
cd ~/src/prj-flagella-estimation
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/issue184_2015_nf1_6_10tau/job.yaml --dry-run
```

## ユーザー実行

```bash
.venv-cs10/bin/python scripts/cs10/queue.py enqueue \
  --branch codex/issue-61-2015-10tau-stability \
  --config conf/phase2_parallel/issue184_2015_nf1_6_10tau/job.yaml
```

queue dispatcherは別途稼働しているものを使用する。停止・再開が必要なら、reservation IDを
確認してからユーザー承認済みのqueue手順を使う。完了済みchild artifactを別reservationへ
コピーしない。

## 完了後

`job_manifest.json`が`succeeded`、`aggregation.status=completed`、6 conditionを確認してから
campaign manifest、summary、compact diagnostics、state archive、trajectoryをローカル同期する。
件数・SHA-256・QCを照合後、以下を実行する。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py \
  --analysis-kind issue184-2015-nf1-6 \
  --run-root <campaign-root> \
  --output-dir <campaign-root>/analysis/issue184_2015_nf1_6

uv run python scripts/03_dataset_building/replay_dataset.py \
  --run-dir <campaign-root> --output-dir <campaign-root>/analysis/replay \
  --view 3d+2d --mode both --camera-3d fixed --camera-2d fixed \
  --view-range-mode campaign-envelope --fps-out-3d 25 --fps-out-2d 25
```

`issue184_decision.json`がPASSでも、dataset採択、profile昇格、canonical selection、Phase 3 handoffは行わない。
