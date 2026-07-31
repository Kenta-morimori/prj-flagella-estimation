# scripts

このディレクトリには、リポジトリで直接実行する CLI を置いています。

## 01_simulate_swimming

Phase 2 の 3D 遊泳シミュレーション、条件 sweep、sweep 結果の heatmap 生成を扱います。

### 単発シミュレーション

```bash
uv run python -m scripts.01_simulate_swimming
```

設定は `conf/sim_swim.yaml` を読みます。一時的な変更や設定ファイル指定は `KEY=VALUE` 形式で指定します。

```bash
uv run python -m scripts.01_simulate_swimming \
  config=conf/sim_swim.yaml \
  time.duration_s=0.05 \
  motor.torque_Nm=2.0e-20 \
  time.dt_star=1.0e-4
```

`--config`、`--duration-s`、`--fps-out` は既存互換用に残しています。新規コマンド例では `config=...`、`time.duration_s=...`、`output_sampling.fps_out_2d=...` を使います。

主な出力は `outputs/YYYY-MM-DD/HHMMSS/` 配下に作成され、`manifest.json` と `run.log` に実行条件が記録されます。

### Sweep

開発用・診断用の task-specific sweep は `run_sweep.py` を使います。対象 profile は `conf/phase2_sweeps/` に置きます。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/shape_stability_grid.yaml \
  dry_run=true \
  sample_limit=3
```

実行時は `[1/3] shape_stability_grid ...` のように進捗が標準出力へ表示されます。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/shape_stability_grid.yaml \
  mode=first-second-grid \
  time.duration_s=0.001 \
  motor.torque_Nm=0 \
  first_second_spring_scales=1 \
  output_dir=/private/tmp/phase2_smoke \
  overwrite=true
```

profile の既定値は `KEY=VALUE` で上書きできます。sweep の標準 summary は出力先の `summary.csv` です。`--config` などの legacy option 形式も互換用に残しています。

### Multi-Run

ユーザが複数条件を一度で実行し、その結果から replay / plot へ進む入口は `run_multi_run.py` です。profile は `conf/phase2_multi_run/` に置き、1つの config を run / plot / replay で共用します。campaign 単位で `run.log`、`manifest.json`、`run_manifest.json`、`summary.csv` を残します。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml \
  dry_run=true
```

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml \
  sweep.axes.torque.values=[1.5e-20,2.0e-20,2.5e-20] \
  time.duration_s=1.0 \
  overwrite=true
```

`output.timestamp_subdir=false` の multi-run profile では、`output.base_dir` がそのまま run root になります。再実行で同じ run root を置き換える場合だけ `overwrite=true` を明示します。

主な profile:

| profile | 用途 |
| --- | --- |
| `conf/phase2_sweeps/motor_scale.yaml` | motor-local scale の sweep |
| `conf/phase2_sweeps/single_flagellum_torque.yaml` | single flagellum torque 条件評価 |
| `conf/phase2_sweeps/bundling_alignment.yaml` | 複数べん毛の helix axis alignment 評価 |
| `conf/phase2_sweeps/shape_stability_grid.yaml` | hook / proximal flagellum を含む shape stability grid |
| `conf/phase2_sweeps/torque_distribution_grid.yaml` | #97 用 torque distribution 2x2 比較 |
| `conf/phase2_sweeps/hook_overstretch.yaml` | 旧名互換 profile |

新規の user-facing 実行では `shape_stability_grid.yaml` を正本として使います。
`hook_overstretch.yaml` は historical alias であり，既存メモや過去コマンドの互換用です。

利用可能な sweep profile を CLI から確認する場合:

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  list_canonical_profiles=true
```

個別 profile の `role` / `canonical` / 推奨 heatmap を確認する場合:

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/shape_stability_grid.yaml \
  describe_profile=true
```

### Heatmap

sweep summary から heatmap を作る場合は `plot_heatmap.py` を使います。

```bash
uv run python scripts/01_simulate_swimming/plot_heatmap.py \
  config=conf/phase2_sweeps/shape_stability_heatmap.yaml \
  summary_csv=/private/tmp/phase2_smoke/summary.csv \
  mode=first-second-grid
```

heatmap profile は出力先を固定しません。`output_dir` を省略すると、`summary_csv` と同じ directory の `plots/` へ出力します。明示した場合はその出力先を使います。
`shape_stability_heatmap.yaml` は `mode=position-only-grid` などの実行時 override で対象 grid を切り替えられます。

generic multi-run の summary plot も `plot_heatmap.py` から行います。同じ config をそのまま使います。`plot.default_y_axis` が未設定の profile では heatmap ではなく 1 軸 line plot を出します。

```bash
uv run python scripts/01_simulate_swimming/plot_heatmap.py \
  config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml
```

`output.timestamp_subdir=false` の multi-run profile では、`summary_csv` / `run_dir` を省略すると `output.base_dir/summary.csv` を読み、`output.base_dir/plots/` へ出力します。

主な profile:

| profile | 用途 |
| --- | --- |
| `conf/phase2_sweeps/motor_scale_heatmap.yaml` | motor scale sweep の collapse heatmap |
| `conf/phase2_sweeps/dt_star_torque_heatmap.yaml` | `dt_star` x torque heatmap |
| `conf/phase2_sweeps/local_scale_mode_heatmap.yaml` | local scale mode x torque heatmap |
| `conf/phase2_sweeps/shape_stability_heatmap.yaml` | shape stability grid heatmap |
| `conf/phase2_multi_run/latest_model_torque_shape_stability.yaml` | generic multi-run summary plot / replay metadata |
| `conf/phase2_sweeps/hook_overstretch_heatmap.yaml` | 旧名互換 heatmap profile |

heatmap も `shape_stability_heatmap.yaml` を正本として使います。
`hook_overstretch_heatmap.yaml` は historical alias です。

利用可能な heatmap profile を CLI から確認する場合:

```bash
uv run python scripts/01_simulate_swimming/plot_heatmap.py \
  list_canonical_profiles=true
```

### Replay Render

既存 sweep 出力の `summary.csv`、`run_manifest.json`、各 condition directory の `state_archive.npz` から、再シミュレーションなしで比較 plot / 3D replay を生成する場合は `render_shape_stability_grid_replay.py` を使います。generic multi-run 出力でも同じ CLI を使います。

```bash
uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py \
  config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml \
  overwrite=true
```

`--mode plot-only` は metrics CSV / PNG のみ、`--mode render-only` は 3D grid movie のみ、`--mode both` は両方を生成します。
`output.timestamp_subdir=false` の multi-run profile では、`run_dir` / `input_dir` を省略すると `output.base_dir` を読み、`output.base_dir/replay/` へ出力します。legacy 互換として `summary_csv=...` や `--input-dir ... --output-dir ...` も引き続き使えます。

## 02_phase2_analysis

Phase 2.8 の RUN 固定べん毛数差分析用 dataset 作成 CLI です。raw 出力は `scripts/01_simulate_swimming/run_multi_run.py` で作り、同じ `conf/phase2_multi_run/*.yaml` を dataset builder / heatmap / replay でも使います。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml \
  dry_run=true sample_limit=3
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/02_phase2_analysis/build_dataset.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
```

標準設定:

| config | 用途 |
| --- | --- |
| `conf/phase2_multi_run/flagella_count_behavior_v0.yaml` | Issue #71 の診断用 dataset v0。run / heatmap / replay / dataset 作成で共通に使う |

Issue #71 の診断用 dataset v0 は次で実行します。36 sample の長時間 run なので、まず `dry_run=true` や `sample_limit=1` で確認してください。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/01_simulate_swimming/plot_heatmap.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/02_phase2_analysis/build_dataset.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/02_phase2_analysis/plot_distributions.py \
  --dataset-id v0
```

raw condition は `step_summary.csv`、`trajectory.csv`、`state_archive.npz` を保存します。dataset directory から 3D / 2D render を作る場合は次を使います。

```bash
uv run python scripts/02_phase2_analysis/render_sample.py \
  --dataset-dir outputs/phase2_analysis/flagella_count_behavior/datasets/v0
```

dataset の分布 plot は次で生成します。

```bash
uv run python scripts/02_phase2_analysis/plot_distributions.py \
  --dataset-id v0
```

dataset v1 の `n_flagella=1,2,3` について，2D投影後の運動特徴量と grouped baseline を確認する場合は次を使います。

```bash
uv run python scripts/02_phase2_analysis/analyze_2d_separability.py \
  --dataset-dir outputs/phase2_analysis/flagella_count_behavior/datasets/v1 \
  --output-dir outputs/phase2_analysis/flagella_count_behavior/datasets/v1/analysis/projection_2d \
  --overwrite
```

## 03_phase3

### v1 r1 3秒runから0.5秒clip datasetを作成

Issue #159 の Phase 3 common clip dataset は、既存の3秒 state archive から生成します。長時間の Phase 2 simulation は再実行しません。canonical `.npy` clip は、べん毛を描画せず、菌体のみを `body_capsule_rigid_v1` のtracking cropで描画します。菌体変形はML入力へ反映せず、QC確認対象として扱います。Phase 3 configでは、出力fpsは `output_sampling.fps_out`、画像条件は `render.image_size_px` / `render.pixel_size_um` を正本とします。

probe は1 classあたり1 runだけ処理します。probe 出力は実行ごとの timestamp path に置き、最終 dataset v1 には使いません。

```bash
uv run python scripts/03_phase3/build_clip_dataset.py \
  config=conf/phase3/gt_passthrough_v1_r1_duration_3s_clips.yaml \
  filters.max_per_class=1
```

final build は27 runすべてを処理し、Phase 3 common clip dataset v1 の canonical path として `outputs/phase3_common_clip/datasets/v1/` に出力します。期待値は5 clips/run、合計135 clips、独立group 27です。

```bash
uv run python scripts/03_phase3/build_clip_dataset.py \
  config=conf/phase3/gt_passthrough_v1_r1_duration_3s_clips.yaml \
  output_dir=outputs/phase3_common_clip/datasets/v1
```

`outputs/YYYY-MM-DD/HHMMSS/phase3_v1_r1_clip_dataset/` はstaging / probe用途の実行ログ付き出力です。最終採択済みの Phase 4 入力 dataset は `outputs/phase3_common_clip/datasets/v1/` を参照します。`data/` は外部入力や手元データ置き場として残し、今回の生成済み共通clip dataset正本には使いません。

QC確認対象は `dataset_summary.csv`、`qc_summary.csv`、`split_summary.csv`、`clip_metadata.jsonl` です。`n_flagella=3` の first-fail run では、first-failを含むwindowとそれ以降のwindowが `qc_label=diagnostic` / `training_candidate=false` になります。ここでの `diagnostic` は、Phase 2 shape QC の `first_fail_t_s` を含む、またはそれ以後のためtrainingから除外した診断用clipという意味であり、全frameで形状崩壊が目視できることを意味しません。early clipは削除せず、Phase 4側の `freeze.warmup_s=0.0/0.5/1.0` で選択します。
v1 r1 では source summary の `use_for_ml_candidate=false` run も除外せず、clip artifact を作成したうえで window QC により diagnostic-only として扱います。

MP4 grid replay は source state archive から同じclip windowを再構成し、defaultでは各panelに3D source replayと2D body-only renderを横並びで表示します。3D側はPhase 2共通3D rendererを使い、べん毛を描画します。2D側はcanonical ML入力と同じ菌体のみのcapsule renderです。titleには `n_flagella`、`run_id`、`clip_index`、time band、QC label が入り、3D status textには first-fail時刻と before/after first-fail status が入ります。

`--max-clips` を省略すると、filterに一致する全clipを処理します。1本のMP4に描画する上限はdefault 12 clipで、対象が12 clipを超える場合は `3d_2d_grid_001.mp4`、`3d_2d_grid_002.mp4` のように分割します。この上限は `--clips-per-video` で変更できます。`--max-clips` は全体の処理件数を明示的に制限したいprobe用途に使います。

```bash
uv run python scripts/03_phase3/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1
uv run python scripts/03_phase3/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1 \
  --qc-label diagnostic
uv run python scripts/03_phase3/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1 \
  --training-candidate true \
  --n-flagella 3 \
  --max-clips 12
```

MP4 panelは必要に応じて選択できます。defaultは `--panel-layout 3d+2d` です。

```bash
uv run python scripts/03_phase3/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1 \
  --panel-layout 3d
uv run python scripts/03_phase3/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1 \
  --panel-layout 2d
```

従来の `.npy` tensor contact sheetを確認したい場合は `--mode contact-sheet` を付けます。

```bash
uv run python scripts/03_phase3/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1 \
  --mode contact-sheet \
  --max-clips 12
```

定性確認点は、early/middle/late clipでtracking cropが破綻していないこと、2D側にflagellaが描画されていないこと、菌体がcapsule状で不自然に伸縮しないこと、3Dと2Dの姿勢・時刻対応が破綻していないこと、pre-first-failとdiagnosticの境界表示がmetadataと一致すること、class/run/time bandが期待どおり見分けられることです。

Phase 4 freeze audit は同じ Phase 3 dataset を直接読みます。`warmup_s` は比較条件として切り替えます。

```bash
uv run python scripts/04_phase4/audit_dataset_freeze.py \
  config=conf/phase4/dataset_freeze_v1_r1.yaml \
  dataset_dir=outputs/phase3_common_clip/datasets/v1
uv run python scripts/04_phase4/audit_dataset_freeze.py \
  config=conf/phase4/dataset_freeze_v1_r1.yaml \
  dataset_dir=outputs/phase3_common_clip/datasets/v1 \
  freeze.warmup_s=0.5
uv run python scripts/04_phase4/audit_dataset_freeze.py \
  config=conf/phase4/dataset_freeze_v1_r1.yaml \
  dataset_dir=outputs/phase3_common_clip/datasets/v1 \
  freeze.warmup_s=1.0
```

## 04_phase4

### Duration / Seed Study

Issue #129 のサンプル時間長・run数検討では、3秒raw run、dataset revision、duration / seed解析を順番に実行します。条件だけを確認する場合は `probe`、full factorial 27 runを実行する場合は `full` を指定します。

```bash
./scripts/04_phase4/run_duration_seed_study.sh probe
caffeinate -i ./scripts/04_phase4/run_duration_seed_study.sh full
```

raw run、dataset revision、duration / seed解析は、実行ごとのJST timestampを持つ同じcampaign rootの `dataset/`、`analysis/` 配下へ保存されます。

```text
outputs/YYYY-MM-DD/HHMMSS/phase4_duration_seed_study/
├── dataset/v1_r1_duration_3s/
├── analysis/phase4_duration_seed_study/
└── replay/
```

同一特徴量のtime plotは `0.25 / 0.5 / 1.0 s` で共通y軸を使います。seed heatmapはdurationを行、`n_flagella` を列とする1枚にまとめ、durationごとに独立したcolor scaleを使います。`x` と `!` はそれぞれQC fail windowとQC fail runを示します。

3秒raw runがすでに存在する場合は、`RUN_DIR`にcampaign rootを指定し、シミュレーションを再実行せずdataset revisionと解析結果だけを生成できます。以下は既存のcanonical campaignを再解析する例です。

```bash
RUN_DIR=outputs/phase2_multi_run/flagella_count_duration_3s_r1
uv run python scripts/02_phase2_analysis/build_dataset.py \
  config=conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml \
  run_dir="${RUN_DIR}" \
  dataset.output_dir="${RUN_DIR}/dataset/v1_r1_duration_3s"
uv run python scripts/04_phase4/evaluate_duration_seed_study.py \
  config=conf/phase4/duration_seed_study_v1_r1.yaml \
  dataset_dir="${RUN_DIR}/dataset/v1_r1_duration_3s" \
  output_dir="${RUN_DIR}/analysis/phase4_duration_seed_study"
```

3秒raw runを再シミュレーションせず3D比較renderする場合:

```bash
uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py \
  config=conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml \
  run_dir="${RUN_DIR}" \
  mode=render-only \
  max_panels_per_grid=9 \
  overwrite=true
```

27条件は最大9条件ずつ複数pageへ出力されます。各パネルにはconditionとともに `PASS`、または `FAIL <category>@<first_fail_t_s>` が表示されます。
