# Issue #203 uniform torque profile runbook

## Scope

このcampaignは `root_torque_segment_couples + uniform` と既存 `diffusive` の同一seed比較である。2010 project default、dataset v2、`dt_star=1e-3`の収束性は変更・確定しない。後者は #200 の責務である。

## cs10 execution

恒常的なruntime、Git更新、tmux、helperの操作は
[`docs/codex/cs10_runbook.md`](../codex/cs10_runbook.md)を正本とする。cs10では
`.venv-cs10/bin/python`を使い、`uv run`、直列の`run_multi_run.py`直接起動、手書きの
marker scriptを使わない。

このcampaignは #207 の independent-condition launcher を使う。27条件を最大8 worker
（8 + 8 + 8 + 3）、各 child process の BLAS thread=1で実行する。途中停止した
`outputs/2026-08-23/213359/phase2_issue203_uniform` は canonical comparison に使わない。

まず本番と同じ27 shardで `0.001 s` qualification を実行する。

```bash
.venv-cs10/bin/python scripts/cs10/parallel_tmux.py start \
  --config conf/phase2_parallel/issue203_uniform_torque_profile/qualification_job.yaml \
  --session issue203-qualification --label issue203_uniform_qualification
```

出力された`control_dir`で状態を確認する。PASS条件は job `succeeded`、`failed_configs=[]`、
aggregate `completed`、campaign `run_summary_count=27`である。

```bash
.venv-cs10/bin/python scripts/cs10/parallel_tmux.py status \
  --control-dir outputs/YYYY-MM-DD/HHMMSS/cs10_parallel/issue203_uniform_qualification
```

qualification後に本番を起動する。本番の所要時間は同じworker policyの実測から再見積もる。

```bash
.venv-cs10/bin/python scripts/cs10/parallel_tmux.py start \
  --config conf/phase2_parallel/issue203_uniform_torque_profile/job.yaml \
  --session issue203 --label issue203_uniform
```

`status`の marker、job、campaign がすべて成功であることを確認してから転送・paired解析へ進む。

## Transfer and Mac analysis

最初は `manifest.json`、`run_manifest.json`、`summary.csv`、`run.log`、`campaign_completion.json`、`user_exit_marker.json`、各conditionの`run_summary.json`、失敗時だけ`stderr.log`と`failure_record.json`を転送する。`step_summary.csv`と`state_archive.npz`は最小転送から除外する。

全27本のMac replayはUser指定のため、第2段階で全conditionの`state_archive.npz`も転送する。paired解析はarchiveのbead座標からaxis指標を再計算するため、`flag_helix_axis_diagnostics.csv`の転送は不要である。

```bash
mkdir -p outputs/YYYY-MM-DD/HHMMSS
# `parallel_tmux.py status` の launch.output_root をここへ設定する。
REMOTE_JOB=~/src/prj-flagella-estimation/outputs/YYYY-MM-DD/HHMMSS/parallel/issue203_uniform__UUID
REMOTE="$REMOTE_JOB/campaign"
LOCAL=outputs/YYYY-MM-DD/HHMMSS/phase2_issue203_uniform
mkdir -p "$LOCAL"
scp Ktakemori@cs10:"$REMOTE_JOB"/job_manifest.json "$LOCAL/"
scp Ktakemori@cs10:"$REMOTE"/{manifest.json,run_manifest.json,summary.csv,run.log,campaign_completion.json} "$LOCAL/"
for n in 01 02 03; do for a in 000 001 002; do for p in 000 001 002; do
  id="as${a}__ps${p}__nf${n}"; mkdir -p "$LOCAL/$id"
  scp Ktakemori@cs10:"$REMOTE/conditions/$id/run_summary.json" "$LOCAL/$id/"
done; done; done
uv run python scripts/03_dataset_building/analyze_issue203_torque_profiles.py \
  --config conf/phase2_analysis/issue203_uniform_paired_comparison.yaml \
  --uniform-run-dir "$LOCAL" --output-dir "$LOCAL/analysis/paired_comparison" --overwrite
```

全27本のMac replayを行う場合だけ、最小転送の後にarchiveを追加する。

```bash
for n in 01 02 03; do for a in 000 001 002; do for p in 000 001 002; do
  id="as${a}__ps${p}__nf${n}"
  scp Ktakemori@cs10:"$REMOTE/conditions/$id/state_archive.npz" "$LOCAL/$id/"
done; done; done
for n in 01 02 03; do for a in 000 001 002; do for p in 000 001 002; do
  uv run python scripts/03_dataset_building/render_issue203_composite_replay.py \
    --run-dir "$LOCAL" \
    --condition-id "as${a}__ps${p}__nf${n}" \
    --output-dir "$LOCAL/analysis/composite"
done; done; done
```

失敗時だけ、該当 child の `stderr.log` と condition の `failure_record.json` を追加転送する。
`job_manifest.json` の失敗 record から child path を確定してから実行する。

```bash
scp Ktakemori@cs10:"$REMOTE_JOB/children/NNN_CONDITION/stderr.log" "$LOCAL/failures/"
scp Ktakemori@cs10:"$REMOTE/conditions/CONDITION/failure_record.json" "$LOCAL/failures/"
```

比較結果の `paired_conditions.csv`、`paired_aggregate.csv`、`paired_delta.png` をレビューし、missing / unpaired / non-passをPASSとして扱わない。profile採否はUser承認後に #200 へ渡す。
