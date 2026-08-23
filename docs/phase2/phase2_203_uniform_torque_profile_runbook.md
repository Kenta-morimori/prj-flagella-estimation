# Issue #203 uniform torque profile runbook

## Scope

このcampaignは `root_torque_segment_couples + uniform` と既存 `diffusive` の同一seed比較である。2010 project default、dataset v2、`dt_star=1e-3`の収束性は変更・確定しない。後者は #200 の責務である。

## cs10 execution

cs10へCodexは接続しない。Userは次でtmuxを開始する。

```bash
ssh -t Ktakemori@cs10 'cd ~/src/prj-flagella-estimation && tmux new -As issue203'
```

tmux内でJST run IDを決め、実行する。cs10 qualificationの短縮runは約8.6 steps/sだったため、50,001 steps × 27条件の逐次実行は概算44時間である。これは実測ではなく保守的な開始目安である。

```bash
RUN_ID=$(TZ=Asia/Tokyo date +%Y-%m-%d/%H%M%S)
RUN_DIR="outputs/${RUN_ID}/phase2_issue203_uniform"
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_uniform_torque_profile_2s_issue203.yaml \
  output.base_dir="$RUN_DIR" output.timestamp_subdir=false
status=$?
python - <<PY
import json, pathlib
p = pathlib.Path("$RUN_DIR") / "user_exit_marker.json"
p.write_text(json.dumps({"exit_code": $status}) + "\\n")
PY
exit "$status"
```

完了確認（Macから）は、stdoutではなくmarker、campaign completion、summary、27件のrun summaryで行う。

```bash
ssh Ktakemori@cs10 'cd ~/src/prj-flagella-estimation && \
  cat outputs/YYYY-MM-DD/HHMMSS/phase2_issue203_uniform/user_exit_marker.json && \
  cat outputs/YYYY-MM-DD/HHMMSS/phase2_issue203_uniform/campaign_completion.json && \
  test -f outputs/YYYY-MM-DD/HHMMSS/phase2_issue203_uniform/summary.csv && \
  find outputs/YYYY-MM-DD/HHMMSS/phase2_issue203_uniform -name run_summary.json | wc -l'
```

## Transfer and Mac analysis

最初は `manifest.json`、`run_manifest.json`、`summary.csv`、`run.log`、`campaign_completion.json`、`user_exit_marker.json`、各conditionの`run_summary.json`、失敗時だけ`stderr.log`と`failure_record.json`を転送する。`step_summary.csv`と`state_archive.npz`は最小転送から除外する。

全27本のMac replayはUser指定のため、第2段階で全conditionの`state_archive.npz`も転送する。paired解析はarchiveのbead座標からaxis指標を再計算するため、`flag_helix_axis_diagnostics.csv`の転送は不要である。

```bash
mkdir -p outputs/YYYY-MM-DD/HHMMSS
REMOTE=~/src/prj-flagella-estimation/outputs/YYYY-MM-DD/HHMMSS/phase2_issue203_uniform
LOCAL=outputs/YYYY-MM-DD/HHMMSS/phase2_issue203_uniform
mkdir -p "$LOCAL"
scp Ktakemori@cs10:"$REMOTE"/{manifest.json,run_manifest.json,summary.csv,run.log,campaign_completion.json,user_exit_marker.json} "$LOCAL/"
for n in 01 02 03; do for a in 000 001 002; do for p in 000 001 002; do
  id="as${a}__ps${p}__nf${n}"; mkdir -p "$LOCAL/$id"
  scp Ktakemori@cs10:"$REMOTE/$id/run_summary.json" "$LOCAL/$id/"
done; done; done
uv run python scripts/03_dataset_building/analyze_issue203_torque_profiles.py \
  --config conf/phase2_analysis/issue203_uniform_paired_comparison.yaml \
  --uniform-run-dir "$LOCAL" --output-dir "$LOCAL/analysis/paired_comparison" --overwrite
```

全27本のMac replayを行う場合だけ、最小転送の後にarchiveを追加する。

```bash
for n in 01 02 03; do for a in 000 001 002; do for p in 000 001 002; do
  id="as${a}__ps${p}__nf${n}"
  scp Ktakemori@cs10:"$REMOTE/$id/state_archive.npz" "$LOCAL/$id/"
done; done; done
for n in 01 02 03; do for a in 000 001 002; do for p in 000 001 002; do
  uv run python scripts/03_dataset_building/render_issue203_composite_replay.py \
    --run-dir "$LOCAL" \
    --condition-id "as${a}__ps${p}__nf${n}" \
    --output-dir "$LOCAL/analysis/composite"
done; done; done
```

比較結果の `paired_conditions.csv`、`paired_aggregate.csv`、`paired_delta.png` をレビューし、missing / unpaired / non-passをPASSとして扱わない。profile採否はUser承認後に #200 へ渡す。
