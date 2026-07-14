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
  --dataset-id fc_v0_nf1_2_3_6_as3_ps3_dur1p0
```

raw condition は `step_summary.csv`、`trajectory.csv`、`state_archive.npz` を保存します。dataset directory から 3D / 2D render を作る場合は次を使います。

```bash
uv run python scripts/02_phase2_analysis/render_sample.py \
  --dataset-dir outputs/phase2_analysis/flagella_count_behavior/datasets/fc_v0_nf1_2_3_6_as3_ps3_dur1p0
```

dataset の分布 plot は次で生成します。

```bash
uv run python scripts/02_phase2_analysis/plot_distributions.py \
  --dataset-id fc_v0_nf1_2_3_6_as3_ps3_dur1p0
```
