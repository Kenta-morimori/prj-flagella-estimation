# Issue #61: 2015 project 1τ cs10 実行runbook

`execution:cs10`のUser-run campaignである。Codexは予約、tmux、queue、開始・停止を行わない。独立3条件で、Mac serial見積りは20時間超である。現在のserial runはparallel job完了時の単発Actions通知を実装する前に開始したため、通知対象ではない。

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
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_issue61_1tau_tracking_stability.yaml dry_run=true
```

## ユーザー実行

このrunはparallel launcherが同一sweep profileへの異なるoverrideを表現できないことを理由にserialで開始した。3 torqueは物理的に独立であり、この理由は今後のserial例外として認めない。以後は、先にcondition shardとparallel aggregateを実装・dry-runし、`cs10_qualified` parallel jobで開始する。1 parallel jobは3 condition workerの終了とaggregate完了後に最終stateをActionsへ一度だけ通知し、conditionごと・queue空状態では通知しない。現在のrunは停止・再起動しない。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_issue61_1tau_tracking_stability.yaml
```

失敗・中断時はcampaign root、`run.log`、`manifest.json`、`run_manifest.json`、conditionごとの`run_summary.json`を保持する。再開は完了conditionをコピーせず、新しいdated output rootで未完了conditionを再実行する。

## 実行後

cs10からローカルへ、campaign manifest、summary、performance、各conditionのcompact summary、state archive、trajectory、必要なQC/evidenceを同期する。件数・SHA-256・QC recordを照合し、operational logやcredentialは同期しない。

simulationを再起動せず、次を実行する。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind issue61-2015-1tau \
  --run-root <campaign-root> \
  --output-dir <campaign-root>/analysis/issue61_2015_1tau
```

最初に`issue61_decision.json`、次に`issue61_summary.csv`を確認する。`status=fail`ならそのcriterionをIssue #61へ記録し、#184へのhandoffやprofile昇格を行わない。
