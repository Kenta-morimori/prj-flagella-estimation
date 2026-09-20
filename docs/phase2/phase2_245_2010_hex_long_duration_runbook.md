# Issue #245: 2010 hex 2 s long-duration runbook

## Scope

`2010_hex_project` pending candidateを、固定`T=2.5e-20 N m / flagellum`、`dt_star=1e-3`、RUN固定、Brownian/switching OFF、archive interval `0.001 s`で2.0 s（reference torqueで50τ、50,000 steps）評価する。canonical model、dataset、ML policy、遊泳特徴量解析は本Issueの対象外である。

対象は18条件である。`n=1..5`は`attach_seed=0,1,2`・`phase_seed=0`、`n=6`は`attach_seed=0`・`phase_seed=0,1,2`を使う。n=6は6 slot全占有の`full_ring_rotation_equivalent`であり、対向する2べん毛はhex固有の制約なので一般的な付着配置の多様性の根拠には使わない。

Heavy/runtime execution targetは`cs10_user_run`（label: `execution:cs10`）である。Codexはcs10への接続、tmux起動、job開始・停止を行わない。3 workersのscreen実測から本番の暫定見積りは約11時間であり、qualification後の実測wall timeで更新する。

## Qualification

ユーザーはPRのcommitを取得後、同一18 shardを`0.001 s`に短縮して実行する。

```bash
cd ~/src/prj-flagella-estimation
git pull --ff-only origin codex/issue-245-2010-hex-long-duration
git rev-parse --short HEAD

.venv-cs10/bin/python scripts/cs10/parallel_tmux.py start \
  --config conf/phase2_parallel/issue245_2010_hex_long_duration/qualification_job.yaml \
  --session issue245-qualification --label issue245_2010_hex_qualification
```

```bash
.venv-cs10/bin/python scripts/cs10/parallel_tmux.py status \
  --control-dir outputs/YYYY-MM-DD/HHMMSS/cs10_parallel/issue245_2010_hex_qualification
```

本番へ進める条件は18/18の`job=succeeded`、`failed_configs=[]`、aggregate=`completed`、campaignの`run_summary_count=18`、全condition strict PASS、archive/manifest/summaryのSHA-256検証成功である。strict failureが1件でもあれば、本番、特徴量評価、採択へ進まずdiagnostics/replayレビューで停止する。

## 2 s campaign, sync, and aggregation

上の条件とユーザーの明示許可後だけ、本番を開始する。

```bash
.venv-cs10/bin/python scripts/cs10/parallel_tmux.py start \
  --config conf/phase2_parallel/issue245_2010_hex_long_duration/job.yaml \
  --session issue245-2s --label issue245_2010_hex_2s
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
