# Issue #203 uniform torque profile runbook

## Scope

このcampaignは `root_torque_segment_couples + uniform` と既存 `diffusive` の同一seed比較である。2010 project default、dataset v2、`dt_star=1e-3`の収束性は変更・確定しない。後者は #200 の責務である。

## Canonical artifact locations and video compatibility

検証済み #203 uniform artifact は、日付別の一時出力ではなく、比較対象と並べて次へ置く。

```bash
UNIFORM=outputs/phase2_multi_run/flagella_count_behavior_v1_r2/reference/2010_project_uniform_tau_linked_2s_nf1_3_as3_ps3_2026-08-23
DIFFUSIVE=outputs/phase2_multi_run/flagella_count_behavior_v1_r2/reference/2010_project_tau_linked_2s_nf1_4_as3_ps3_2026-08-18
```

canonical composite collection は各 reference の `analysis/composite` である。動画・manifest は
`nfNN_composite_grid.*` とし、各gridは attach seed × phase seed の9条件をsubplot表示する。
H.264 High profile (`yuv420p`, `faststart`) のみを残す。logical `condition_id` は従来の
`asNNN__psNNN__nfNN` のままとする。renderer は全経路でFFmpeg `libx264` を必須とし、
OpenCV/`mp4v` fallback は行わない。`composite_h264` は canonical output 名として使わない。

```bash
uv run python scripts/03_dataset_building/transcode_mp4_h264.py \
  --input-dir "$UNIFORM/analysis/composite" \
  --output-dir "$UNIFORM/analysis/.composite_h264_staging" --overwrite
```

過去artifactをMacで一括移行する場合も同じ既存CLIを使う。`--replace` は各MP4を staging で
H.264/duration 検証してから置換し、参照可能な全 manifest のcodec情報も更新する。PNGは既に
PowerPointへ貼り付け可能な形式なので変換しない。

```bash
uv run python scripts/03_dataset_building/transcode_mp4_h264.py \
  --input-dir outputs --output-dir /tmp/phase2_h264_staging \
  --replace --overwrite
```

## cs10 execution

恒常的なruntime、Git更新、tmux、helperの操作は
[`docs/codex/cs10_runbook.md`](../codex/cs10_runbook.md)を正本とする。cs10では
`.venv-cs10/bin/python`を使い、`uv run`、直列の`run_multi_run.py`直接起動、手書きの
marker scriptを使わない。

このcampaignは #207 の independent-condition launcher を使う。27条件を最大8 worker
（8 + 8 + 8 + 3）、各 child process の BLAS thread=1で実行する。途中停止した
`outputs/2026-08-23/213359/phase2_issue203_uniform` は canonical comparison に使わない。

初回またはFFmpeg未導入時は、render前にこの既存setupを実行する。これは固定版の
ユーザー領域 static FFmpeg (`~/.local/opt`) をSHA-256検証して導入し、`libx264`を検査する。

```bash
bash scripts/cs10/setup_environment.sh
~/.local/bin/ffmpeg -hide_banner -encoders | grep libx264
```

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

## cs10 analysis and replay

2 s archiveはcs10に残し、motion feature解析・grid replay・composite replayもcs10で行う。
aggregate campaignの`conditions/<condition_id>`はchild artifactへの相対symlinkであり、
`run_manifest.json`と同じcondition IDを正本として解決する。archiveをMacへ一括転送しない。

```bash
CAMPAIGN=outputs/YYYY-MM-DD/HHMMSS/parallel/issue203_uniform_torque_profile__UUID/campaign

# Reference #204 と同じ 3D/2D時系列・0.25/0.5/1.0 s window解析。
.venv-cs10/bin/python scripts/03_dataset_building/analyze_motion_features.py \
  --config conf/phase2_analysis/issue203_uniform_motion_feature_study.yaml \
  run_dir="$CAMPAIGN" output_dir="$CAMPAIGN/analysis/motion_features" overwrite=true

# 27条件の3D/2D replay grid。長時間renderはtmux内で実行する。
.venv-cs10/bin/python scripts/03_dataset_building/replay_dataset.py \
  --run-dir "$CAMPAIGN" --output-dir "$CAMPAIGN/analysis/replay" \
  --view 3d+2d --max-panels-per-grid 9 --overwrite

# rendererがH.264を必須化しているため、生成後に形式だけを確認する。
~/.local/bin/ffprobe -v error -select_streams v:0 \
  -show_entries stream=codec_name,pix_fmt -of default=nokey=1:noprint_wrappers=1 \
  "$CAMPAIGN/analysis/replay/grid_swim3d.mp4"

# n=1/2/3ごとに9条件を3×3 subplotで表示する、3本の3D +
# local-segment nominal torque-weight composite replay。
mkdir -p "$CAMPAIGN/analysis/composite"
for n in 1 2 3; do
  .venv-cs10/bin/python scripts/03_dataset_building/render_issue203_composite_replay.py \
    --run-dir "$CAMPAIGN" --n-flagella "$n" \
    --output-dir "$CAMPAIGN/analysis/composite"
done
```

通常のreplay gridでも、`--show-torque-weight-panels`を指定すれば各3D panelの横に
べん毛別のlocal-segment nominal torque weightを表示できる。uniformは時間不変、
diffusiveは保存時刻までのlocal-twist update再構成であり、realized force / torqueの
表示ではない。

```bash
.venv-cs10/bin/python scripts/03_dataset_building/replay_dataset.py \
  --run-dir "$CAMPAIGN" --output-dir "$CAMPAIGN/analysis/replay_with_weights" \
  --view 3d --show-torque-weight-panels --overwrite
```

`motion_features/manifest.json`、`replay/manifest.json`、各
`composite/nf??_composite_grid_manifest.json`と、MP4数（grid replayはページ数、compositeは3）を
確認する。strict non-PASSを通常の
比較plotへ混入させず、condition CSVとmanifestには残す。

diffusive reference も比較の対照表示として n=1--3 の27条件を同じ composite CLI で生成する。
`diffusive` panel は実装と同じ local-twist update から各保存時刻の segment weight を再構成する
（uniform の time-invariant panel とは区別される）。これは archive を読む render であり、cs10
では tmux 内で、Mac では archive を移送済みならローカルで実行する。

```bash
mkdir -p "$DIFFUSIVE/analysis/composite"
for n in 1 2 3; do
  uv run python scripts/03_dataset_building/render_issue203_composite_replay.py \
    --run-dir "$DIFFUSIVE" --n-flagella "$n" \
    --output-dir "$DIFFUSIVE/analysis/composite"
done
```

## n=3 profile × dt contact diagnostic

uniform n=3の2 failure seed (`as001__ps000__nf03`, `as002__ps002__nf03`)と
PASS control (`as000__ps000__nf03`)を診断する。既存uniform / diffusiveの`dt_star=1e-3`
6 artifactは再利用するため、cs10で新規実行するのは`dt_star=3e-4,1e-4`の12 shardだけである。
各runは2.0 s、RUN fixed、Brownian OFFで、compact outputでも全internal stepのfinite / nonbody /
body QCを集約する。これはprofile・dt・datasetの採択や排除力parameter変更を行うcampaignではない。

uniformの小dt 6条件を通常3D gridとして確認する場合も、既存archiveから再描画するだけでよい。
simulationは再実行しない。

```bash
.venv-cs10/bin/python scripts/03_dataset_building/replay_dataset.py \
  --run-dir "$CAMPAIGN" --output-dir "$CAMPAIGN/analysis/replay_uniform_dt" \
  --view 3d --max-panels-per-grid 6 --overwrite \
  --condition-id nf03__as000__ps000__uniform__dt3e-4 \
  --condition-id nf03__as000__ps000__uniform__dt1e-4 \
  --condition-id nf03__as001__ps000__uniform__dt3e-4 \
  --condition-id nf03__as001__ps000__uniform__dt1e-4 \
  --condition-id nf03__as002__ps002__uniform__dt3e-4 \
  --condition-id nf03__as002__ps002__uniform__dt1e-4
```

PRを最新化してからtmux管理下で開始する。

```bash
cd ~/src/prj-flagella-estimation
git pull --ff-only origin codex/issue-203-uniform-torque-profile
git rev-parse --short HEAD

.venv-cs10/bin/python scripts/cs10/parallel_tmux.py start \
  --config conf/phase2_parallel/issue203_torque_profile_dt_contact/job.yaml \
  --session issue203-dt-contact --label issue203_dt_contact
```

出力された`control_dir`を用い、job / aggregate / 12 run summary がすべて完了することを確認する。

```bash
.venv-cs10/bin/python scripts/cs10/parallel_tmux.py status \
  --control-dir outputs/YYYY-MM-DD/HHMMSS/cs10_parallel/issue203_dt_contact
```

目安は8 workerで2.5--4時間である。完了後はcs10上で解析する。`$CAMPAIGN`には上記statusで示される
parallel jobの`campaign_root`を指定する。

```bash
CAMPAIGN=outputs/YYYY-MM-DD/HHMMSS/parallel/issue203_torque_profile_dt_contact__UUID/campaign
.venv-cs10/bin/python scripts/03_dataset_building/analyze_torque_profile_dt_contact.py \
  extract --run-dir "$CAMPAIGN" --source new_diagnostic \
  --seed-case nf03__as000__ps000 --seed-case nf03__as001__ps000 \
  --seed-case nf03__as002__ps002 \
  --output-dir "$CAMPAIGN/analysis/contact_stability_fragment" --overwrite
```

この処理はarchiveをcs10上で読む。生成するfragmentはCSVと`fragment_manifest.json`だけで、
source artifactのSHA-256を記録する。`state_archive.npz`はcs10に残す。

## Transfer and paired analysis

contact診断でcs10からMacへ転送するのはfragment CSVとmanifestだけである。archive、
`step_summary.csv`、動画は転送しない。

```bash
mkdir -p outputs/YYYY-MM-DD/HHMMSS
# `parallel_tmux.py status` の launch.output_root をここへ設定する。
REMOTE_JOB=~/src/prj-flagella-estimation/outputs/YYYY-MM-DD/HHMMSS/parallel/issue203_uniform_torque_profile__UUID
REMOTE="$REMOTE_JOB/campaign"
LOCAL=outputs/YYYY-MM-DD/HHMMSS/phase2_issue203_uniform
mkdir -p "$LOCAL"
mkdir -p "$LOCAL/contact_fragments/cs10"
scp Ktakemori@cs10:"$REMOTE/analysis/contact_stability_fragment"/{contact_stability_fragment.csv,fragment_manifest.json} "$LOCAL/contact_fragments/cs10/"
```

Macでは既存reference archiveから各3条件のfragmentを作る。`$UNIFORM` と
`$DIFFUSIVE` は既存reference campaign rootを指定する。

```bash
uv run python scripts/03_dataset_building/analyze_torque_profile_dt_contact.py \
  extract --run-dir "$UNIFORM" --source reused_reference --profile uniform --dt-star 1e-3 \
  --seed-case nf03__as000__ps000 --seed-case nf03__as001__ps000 --seed-case nf03__as002__ps002 \
  --output-dir "$LOCAL/contact_fragments/uniform_reference" --overwrite
uv run python scripts/03_dataset_building/analyze_torque_profile_dt_contact.py \
  extract --run-dir "$DIFFUSIVE" --source reused_reference --profile diffusive --dt-star 1e-3 \
  --seed-case nf03__as000__ps000 --seed-case nf03__as001__ps000 --seed-case nf03__as002__ps002 \
  --output-dir "$LOCAL/contact_fragments/diffusive_reference" --overwrite
uv run python scripts/03_dataset_building/analyze_torque_profile_dt_contact.py \
  combine --config conf/phase2_analysis/torque_profile_dt_contact_diagnostic.yaml \
  --fragment-dir "$LOCAL/contact_fragments/uniform_reference" \
  --fragment-dir "$LOCAL/contact_fragments/diffusive_reference" \
  --fragment-dir "$LOCAL/contact_fragments/cs10" --overwrite
```

`combine` は3 seed × uniform/diffusive × `1e-3/3e-4/1e-4`の18条件を厳格に要求する。
欠損・重複・profile/dt/seed/source不整合、またはCSV/condition provenanceのSHA-256不一致は失敗する。
成功時のcanonical出力は`$UNIFORM/analysis/contact_stability_dt/`のCSV、図、manifestである。

cs10で生成した解析・動画をMacへ回収する場合は、次を追加転送する。

```bash
scp -r Ktakemori@cs10:"$REMOTE/analysis/motion_features" "$LOCAL/analysis/"
scp -r Ktakemori@cs10:"$REMOTE/analysis/replay" "$LOCAL/analysis/"
scp -r Ktakemori@cs10:"$REMOTE/analysis/composite" "$LOCAL/analysis/"
```

失敗時だけ、該当 child の `stderr.log` と condition の `failure_record.json` を追加転送する。
`job_manifest.json` の失敗 record から child path を確定してから実行する。

```bash
scp Ktakemori@cs10:"$REMOTE_JOB/children/NNN_CONDITION/stderr.log" "$LOCAL/failures/"
scp Ktakemori@cs10:"$REMOTE/conditions/CONDITION/failure_record.json" "$LOCAL/failures/"
```

比較結果の `paired_conditions.csv`、`paired_aggregate.csv`、`paired_delta.png` をレビューし、missing / unpaired / non-passをPASSとして扱わない。profile採否はUser承認後に #200 へ渡す。
