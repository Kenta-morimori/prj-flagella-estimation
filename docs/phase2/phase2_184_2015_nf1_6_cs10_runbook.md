# Issue #184: 2015 project 計算費用probe cs10 runbook

本runbookはPR #242のmerge前に、2015 projectの計算費用を実測するためのもの。物理的な採択や10τの安定性評価ではない。Issue #61の既存3条件はstrict FAILであり、pitch QCの見直しは保留する。

## 固定契約

- `conf/phase2_parallel/issue184_2015_runtime_probe_0p01s/job.yaml`を1予約として使い、`seeded_surface`・seed 0・`2.5e-20 N m`・tracking-reference・`dt_star=1e-5`・RUN固定・Brownian OFFのnf1–6を独立shardで実行する。
- 各条件は0.01実秒=0.25τ=25,000 internal steps。configでは浮動小数点の端数による25,001 stepを避けるため`0.25 tau`として指定する。`cs10_qualified`、最大3 worker、compact checkpoint 2,500 steps、output分離。全6 shardとaggregateの終端確定後にActions通知を予約単位で最大1回送る。
- 既存の2010 hex予約6には触れない。新probeの後にもdispatcherを`--once`で終え、予約6を自動起動しない。

## 旧runの扱い

reservation 5の旧`seeded_center_layer` jobはnf1–3のみ完走。nf4/nf5はgeometry構成に失敗し、nf6は未完走である。直接実行中のnf4/nf5 `seeded_surface`も未完走である。queueをpauseし、対象PID・process group・child状態を再確認してからreservation 5をcancel、direct nf4/nf5を停止する。旧reservation 5の`cancelled`通知1件は許容する。

完走済みnf1–3、停止済みdirect nf1/nf2、および親jobのmanifestsは保持する。削除するのは対象・未完走を記録したnf4/nf5/nf6のchild出力のみで、削除対象と回復可能性を実施後に報告する。旧artifactを新probeへコピー・混在させない。

## 起動前と通知

PR #242の固定commitで隔離worktreeを用意する。checkoutを変更する前に、旧run停止と対象出力の処理を完了する。dry-runとgeometry/output preflightで6 condition、attachment body bead/layer/slot、3 worker、分離pathを確認する。`gh auth status`とbranch上の`cs10-queue-notify.yml`を読み取り確認し、認証情報は記録・同期しない。

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  --config conf/phase2_parallel/issue184_2015_runtime_probe_0p01s/job.yaml --dry-run
```

queueにprobeを予約6より高priorityでenqueueする。dispatcherは環境変数`CS10_RUNTIME_PYTHON=/home/people/Ktakemori/src/prj-flagella-estimation/.venv-cs10/bin/python`と`CS10_QUEUE_NOTIFICATION_REF=codex/issue-61-2015-10tau-stability`を設定し、`scripts/cs10/queue.py run --once`で起動する。実行用Pythonは共有venvを指すが、コードとconfigは固定commitのclean worktreeから読む。通知refの既定は`main`であり、PR branch上workflowを使う今回だけoverrideする。失敗・取消でも当該予約の終端時に1通知を試み、再試行は自動では行わない。

## 終了後

`job_manifest.json`の`succeeded`、`aggregation.status=completed`、campaign completion、6つのcondition symlink/summary/performanceを確認する。必要なmanifest、summary、performance、compact diagnosticsとreplay用archiveをローカルへ同期し、件数・SHA-256・QCを照合する。cs10 operational logやcredentialsは同期しない。

```bash
uv run python scripts/01_simulate_swimming/estimate_runtime.py \
  --job-root <local-job-root> --target-duration-s 0.5 \
  --conditions nf01,nf02,nf03,nf04,nf05,nf06 \
  --historical-performance nf01=<local-old-nf01-performance.json> \
  --historical-performance nf02=<local-old-nf02-performance.json> \
  --historical-performance nf03=<local-old-nf03-performance.json> \
  --output-dir <local-job-root>/analysis/runtime_projection
```

`runtime_projection.csv/json`に各条件の実測wall time・steps/s、0.5秒=12.5τ=1,250,000 stepsへの50倍外挿、3-worker makespan、nf1–3の既存10τ実測との比、不確かさ、provenanceを記録する。速度測定の成功はstrict QC PASSではない。6条件が揃わなければ採否を判断しない。結果表をユーザーに提示し、2015を採用候補から外すか確認してからPR #242へ判断を記録する。
