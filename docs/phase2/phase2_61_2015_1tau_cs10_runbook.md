# Issue #61: 2015 project 1τ cs10 実行runbook

`execution:cs10`のUser-run campaignである。対象は独立した3 torque条件`1e-21`、`2.5e-20`、`1e-19 N m`で、Mac serial見積りは20時間超である。本runbookの正本は`cs10_qualified` parallel jobであり、3 workerとaggregateの全完了後にqueue reservation単位でActions通知を一度だけ行う。conditionごと・queue空状態では通知しない。

## 実行前

cs10上でStage A（および存在する場合はTask D）のsource manifestを確認し、次のJSONを作成する。`manifest_path`は実在するmanifestへの絶対path、`source_run_root`はそのrun rootにする。

```json
[
  {"label": "2015 Stage A", "source_run_root": "/path/to/stage-a", "manifest_path": "/path/to/stage-a/manifest.json"}
]
```

Task Dが未実行ならentryを追加しない。既存manifestがcs10に存在しない場合はevidence JSONを指定せずに実行し、campaign manifestの空配列とdecisionの`none_available_at_run_start`を保持する。存在しないoutput、推測したSHA-256、Stage A/Task D outputのコピーは許可しない。

まず3条件を確認する（simulationは起動しない）。

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/issue61_2015_1tau/job.yaml --dry-run
```

## ユーザー実行

2026-09-05のserial選択はlauncher制約による運用不備であり、分析根拠にはしない。この判断は将来のserial例外を作らない。新規launchはqueueの1 reservationとして行い、各workerを次のjob YAMLから起動する。

```bash
.venv-cs10/bin/python scripts/cs10/queue.py enqueue \
  --job-yaml conf/phase2_parallel/issue61_2015_1tau/job.yaml
```

失敗・中断時はshard artifactを保持する。完了conditionをコピーして混在させず、新しいdated output rootでcleanなparallel jobを実行する。

## 完了済み shard の再集約

2026-09-06のparallel runは3 shardが完了済みであり、child内部の単一condition IDがすべて`project`だったため、初回aggregateだけが重複として失敗した。child outputは変更・コピーせず、aggregate側でtask IDと実torqueからcanonical IDを付ける。以下はsimulationを起動しない。

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/issue61_2015_1tau/job.yaml \
  --output-root /net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs/2026-09-06/220250/parallel/issue61-2015-1tau__3654804140d9 \
  --aggregate-existing
```

成功時のみ`job_manifest.json`は`succeeded`、`aggregation.status`は`completed`となる。`campaign/conditions/`には`project_torque_1em21`、`project_torque_2p5em20`、`project_torque_1em19`の3 symlinkを作る。統合manifestとsummaryにはcanonical IDと、不変のchild ID・child pathをprovenanceとして記録する。

## 実行後

cs10からローカルへ、campaign manifest、summary、performance、各conditionのcompact summary、state archive、trajectory、必要なQC/evidenceを同期する。件数・SHA-256・QC recordを照合し、operational logやcredentialは同期しない。

simulationを再起動せず、次を実行する。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind issue61-2015-1tau \
  --run-root <campaign-root> \
  --output-dir <campaign-root>/analysis/issue61_2015_1tau
```

最初に`issue61_decision.json`、次に`issue61_summary.csv`を確認する。`status=fail`ならそのcriterionをIssue #61へ記録し、#184へのhandoffやprofile昇格を行わない。

その後、campaignの3 conditionを`3d+2d`でreplayする。これはstate archiveからのrenderであり、simulationを再起動しない。

```bash
.venv-cs10/bin/python scripts/03_dataset_building/replay_dataset.py \
  --run-dir <campaign-root> --output-dir <campaign-root>/analysis/replay \
  --view 3d+2d --camera-3d fixed --camera-2d fixed \
  --view-range-mode campaign-envelope --fps-out-3d 25 --fps-out-2d 25
```
