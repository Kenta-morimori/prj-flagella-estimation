# Issue #184: 2015 project nf1–6 10τ cs10 runbook

このcampaignは`execution:cs10`のdiagnostic-only実行である。今回に限り、未実行だったnf4/nf5を
`seeded_surface`・attach/phase seed `0`で2 worker同時実行する。Actions通知は送らない。

## 現行failed job

2026-09-08の`seeded_center_layer` reservationはnf4/nf5が`n_prism=6`の割り切れ制約で
開始前に失敗した。nf1–3と完走するnf6は保持する。今回停止したnf1/nf2 direct rootは隔離し、
暫定表にも使わない。nf1–3/nf6と新nf4/nf5だけは、topology混在を明記した暫定横断表に載せられる。
このreservationは停止しない。

## 実行前

cs10のworktreeで、以下のdry-runによりnf4/nf5、worker数、分離output、
tracking-reference、seed、およびattachment body bead / layer / slotを確認する。
geometry preflightが失敗したらqueueへ投入しない。

```bash
export PATH=/home/people/Ktakemori/.local/bin:/usr/local/bin:/usr/bin
cd ~/src/prj-flagella-estimation
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  --config conf/phase2_parallel/issue184_2015_nf4_nf5_10tau/job.yaml --dry-run
```

## ユーザー実行

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  --config conf/phase2_parallel/issue184_2015_nf4_nf5_10tau/job.yaml \
  --max-workers 2 --output-root <new-nas-root>
```

今回のdirect runはqueue reservationを作らず、Actions通知も送らない。旧nf6は継続する。
完了済みchild artifactを別rootへコピーしない。

## 完了後

`job_manifest.json`が`succeeded`、`aggregation.status=completed`、nf4/nf5を確認してから
campaign manifest、summary、compact diagnostics、state archive、trajectoryをローカル同期する。
件数・SHA-256・QCを照合後、以下を実行する。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py \
  --analysis-kind issue184-2015-nf1-6-provisional \
  --source nf01=<old-nf01-run-root> --source nf02=<old-nf02-run-root> \
  --source nf03=<old-nf03-run-root> --source nf04=<new-nf04-run-root> \
  --source nf05=<new-nf05-run-root> --source nf06=<old-nf06-run-root> \
  --output-dir <local-output>/analysis/issue184_provisional

uv run python scripts/03_dataset_building/replay_dataset.py \
  --run-dir <campaign-root> --output-dir <campaign-root>/analysis/replay \
  --view 3d+2d --mode both --camera-3d fixed --camera-2d fixed \
  --view-range-mode campaign-envelope --fps-out-3d 25 --fps-out-2d 25
```

暫定`issue184_provisional_decision.json`は常に`status=provisional`であり、dataset採択、profile昇格、canonical selection、Phase 3 handoffは行わない。
