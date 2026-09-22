# Issue #245: 2010 hex 2 s long-duration runbook

## Scope

`2010_hex_project` pending candidateを、固定`T=2.5e-20 N m / flagellum`、`dt_star=1e-3`、RUN固定、Brownian/switching OFF、archive interval `0.001 s`で2.0 s（reference torqueで50τ、50,000 steps）評価する。canonical model、dataset、ML policy、遊泳特徴量解析は本Issueの対象外である。

対象は18条件である。`n=1..5`は`attach_seed=0,1,2`・`phase_seed=0`、`n=6`は`attach_seed=0`・`phase_seed=0,1,2`を使う。n=6は6 slot全占有の`full_ring_rotation_equivalent`であり、対向する2べん毛はhex固有の制約なので一般的な付着配置の多様性の根拠には使わない。

Heavy/runtime execution targetは`cs10_user_run`（label: `execution:cs10`）である。Codexはcs10への接続、tmux起動、job開始・停止を行わない。3 workersのscreen実測から本番の暫定見積りは約11時間であり、qualification後の実測wall timeで更新する。

## Qualification

予約・起動の前に、ユーザーはreview済みPR commit、clean worktree、18 shard dry-run、NAS容量、
queue/GitHub認証、既存reservationとの衝突なしを確認する。以下は予約を作成しないpre-reservation
checklistである。

```bash
cd ~/src/prj-flagella-estimation
git fetch origin
git show --quiet --format='%H %s' origin/codex/issue-245-2010-hex-long-duration
git status --porcelain
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  --config conf/phase2_parallel/issue245_2010_hex_long_duration/qualification_job.yaml \
  --dry-run
df -h /net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs
gh auth status --hostname github.com
.venv-cs10/bin/python scripts/cs10/queue.py status
```

上記が成立し、**ユーザーがqualification reservationを明示許可した後だけ**、同一18 shardを`0.001 s`に短縮して予約する。

```bash
cd ~/src/prj-flagella-estimation
git fetch origin
.venv-cs10/bin/python scripts/cs10/queue.py enqueue \
  --branch origin/codex/issue-245-2010-hex-long-duration \
  --config conf/phase2_parallel/issue245_2010_hex_long_duration/qualification_job.yaml \
  --priority 0
```

```bash
.venv-cs10/bin/python scripts/cs10/queue.py status
```

本番へ進める条件は18/18の`job=succeeded`、`failed_configs=[]`、aggregate=`completed`、campaignの`run_summary_count=18`、全condition strict PASS、archive/manifest/summaryのSHA-256検証成功である。`progress.json`、`diagnostic_samples.csv`、`state_archive.partial.npz`はstrict failure時のdiagnostics/replay専用であり、completed archiveや採択根拠として使わない。strict failureが1件でもあれば、本番、特徴量評価、採択へ進まずdiagnostics/replayレビューで停止する。

## 2 s campaign, sync, and aggregation

上の条件、同期済みrequired artifactのhash検証、ユーザーの明示許可後だけ、本番reservationを作成できる。

```bash
.venv-cs10/bin/python scripts/cs10/queue.py enqueue \
  --branch origin/codex/issue-245-2010-hex-long-duration \
  --config conf/phase2_parallel/issue245_2010_hex_long_duration/job.yaml \
  --priority 0
```

完了したcampaignをMacへ同期する際は、再解析に必要な各conditionの`state_archive.npz`、`run_summary.json`、`performance.json`、root manifestとsummaryを保持し、SHA-256を検証する。

```bash
.venv/bin/python scripts/cs10/sync_reference_from_cs10.py \
  --host cs10 --remote-dir "$CAMPAIGN" \
  --local-dir outputs/YYYY-MM-DD/HHMMSS/issue245_2010_hex_2s

.venv/bin/python scripts/03_dataset_building/analyze_dataset.py \
  --analysis-kind model-development-evaluation \
  --config conf/phase2_multi_run/2010_hex_project_long_duration_2s_issue245.yaml \
  --run-dir outputs/YYYY-MM-DD/HHMMSS/issue245_2010_hex_2s \
  --output-dir outputs/YYYY-MM-DD/HHMMSS/model_development_evaluation \
  --render-replay
```

集約器は`manifest.json`、`run.log`、統合summary、`window_qc.csv`、first failure、wall time、steps/s、artifact SHA-256、およびn別・seed別固定camera 3D/2D replayを作る。hook angleは値とfirst failureを保存するdiagnostic-onlyであり、finite/body/hook length/flag/motorだけがstrict PASS/FAILである。
