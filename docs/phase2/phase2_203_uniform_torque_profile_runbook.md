# Issue #203 uniform torque profile runbook

## Scope

このcampaignは `root_torque_segment_couples + uniform` と既存 `diffusive` の同一seed比較である。2010 project default、dataset v2、`dt_star=1e-3`の収束性は変更・確定しない。後者は #200 の責務である。

## cs10 execution

cs10へCodexは接続しない。Userは次でtmuxを開始する。

```bash
ssh -t Ktakemori@cs10 'cd ~/src/prj-flagella-estimation && tmux new -As issue203'
```

cs10ではMac用`.venv`や`uv run`を使わない。Issue #208のwheel-only runtime
`.venv-cs10`を準備し、以降のPython commandには必ず`.venv-cs10/bin/python`を使う。

```bash
export CS10_PYTHON="$(readlink -f .venv/bin/python)"
export CS10_VENV="$PWD/.venv-cs10"
unset C_INCLUDE_PATH LIBRARY_PATH
bash scripts/cs10/setup_environment.sh
.venv-cs10/bin/python -c 'import matplotlib, numpy, yaml; print("cs10 runtime OK")'
```

このcampaignは #207 の independent-condition launcher を使う。27条件を最大8 worker
（8 + 8 + 8 + 3）で実行し、各子 process の BLAS thread は1に固定する。直列の
`run_multi_run.py` を直接起動してはならない。途中停止した
`outputs/2026-08-23/213359/phase2_issue203_uniform` は canonical comparison に使わない。

まず本番と同じ27 shardで `0.001 s` の parallel qualification を実施する。

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/issue203_uniform_torque_profile/qualification_job.yaml \
  dry_run=true

RUN_ID=$(TZ=Asia/Tokyo date +%Y-%m-%d/%H%M%S)
MARKER_DIR="outputs/${RUN_ID}/phase2_issue203_uniform_parallel_qualification"
mkdir -p "$MARKER_DIR"
JOB_ROOT=$(.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/issue203_uniform_torque_profile/qualification_job.yaml)
status=$?
printf '{"exit_code": %s, "job_root": "%s"}\n' "$status" "$JOB_ROOT" \
  > "$MARKER_DIR/user_exit_marker.json"
printf '%s\n' "$JOB_ROOT" > "$MARKER_DIR/job_root.txt"
exit "$status"
```

qualification が PASS（`job_manifest.json.status=succeeded`、`failed_configs=[]`、
`campaign/campaign_completion.json.status=completed`、27 summaries）なら、同様に本番を起動する。
本番の実時間は条件依存である。直列の約44時間という旧見積りは適用しない。8 worker の理想下限は
約5.5時間で、I/O・condition差を含めた保守的な予約枠は **8--12時間** とする。

```bash
RUN_ID=$(TZ=Asia/Tokyo date +%Y-%m-%d/%H%M%S)
MARKER_DIR="outputs/${RUN_ID}/phase2_issue203_uniform_parallel"
mkdir -p "$MARKER_DIR"
JOB_ROOT=$(.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/issue203_uniform_torque_profile/job.yaml)
status=$?
printf '{"exit_code": %s, "job_root": "%s"}\n' "$status" "$JOB_ROOT" \
  > "$MARKER_DIR/user_exit_marker.json"
printf '%s\n' "$JOB_ROOT" > "$MARKER_DIR/job_root.txt"
exit "$status"
```

完了確認（Macから）は stdout ではなく marker、parallel job、canonical campaign を確認する。

```bash
ssh Ktakemori@cs10 'cd ~/src/prj-flagella-estimation && \
  cat outputs/YYYY-MM-DD/HHMMSS/phase2_issue203_uniform_parallel/user_exit_marker.json && \
  JOB_ROOT=$(cat outputs/YYYY-MM-DD/HHMMSS/phase2_issue203_uniform_parallel/job_root.txt) && \
  jq '{status,failed_configs,aggregation}' "$JOB_ROOT/job_manifest.json" && \
  cat "$JOB_ROOT/campaign/campaign_completion.json" && \
  test -f "$JOB_ROOT/campaign/summary.csv" && \
  find -L "$JOB_ROOT/campaign/conditions" -name run_summary.json | wc -l'
```

## Transfer and Mac analysis

最初は `manifest.json`、`run_manifest.json`、`summary.csv`、`run.log`、`campaign_completion.json`、`user_exit_marker.json`、各conditionの`run_summary.json`、失敗時だけ`stderr.log`と`failure_record.json`を転送する。`step_summary.csv`と`state_archive.npz`は最小転送から除外する。

全27本のMac replayはUser指定のため、第2段階で全conditionの`state_archive.npz`も転送する。paired解析はarchiveのbead座標からaxis指標を再計算するため、`flag_helix_axis_diagnostics.csv`の転送は不要である。

```bash
mkdir -p outputs/YYYY-MM-DD/HHMMSS
REMOTE_JOB=~/src/prj-flagella-estimation/outputs/YYYY-MM-DD/HHMMSS/parallel/issue203_uniform_torque_profile__UUID
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
