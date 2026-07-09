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

複数条件の sweep は `run_sweep.py` を使います。条件セットは `conf/phase2_sweeps/` の YAML profile で選びます。

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

主な profile:

| profile | 用途 |
| --- | --- |
| `conf/phase2_sweeps/motor_scale.yaml` | motor-local scale の sweep |
| `conf/phase2_sweeps/single_flagellum_torque.yaml` | single flagellum torque 条件評価 |
| `conf/phase2_sweeps/bundling_alignment.yaml` | 複数べん毛の helix axis alignment 評価 |
| `conf/phase2_sweeps/shape_stability_grid.yaml` | hook / proximal flagellum を含む shape stability grid |
| `conf/phase2_sweeps/torque_distribution_grid.yaml` | #97 用 torque distribution 2x2 比較 |

### Heatmap

sweep summary から heatmap を作る場合は `plot_heatmap.py` を使います。

```bash
uv run python scripts/01_simulate_swimming/plot_heatmap.py \
  config=conf/phase2_sweeps/shape_stability_heatmap.yaml \
  summary_csv=/private/tmp/phase2_smoke/summary.csv \
  mode=first-second-grid
```

heatmap profile は出力先を固定しません。`output_dir` を省略すると、`summary_csv` と同じ directory の `plots/` へ出力します。明示した場合はその出力先を使います。

主な profile:

| profile | 用途 |
| --- | --- |
| `conf/phase2_sweeps/motor_scale_heatmap.yaml` | motor scale sweep の collapse heatmap |
| `conf/phase2_sweeps/dt_star_torque_heatmap.yaml` | `dt_star` x torque heatmap |
| `conf/phase2_sweeps/local_scale_mode_heatmap.yaml` | local scale mode x torque heatmap |
| `conf/phase2_sweeps/shape_stability_heatmap.yaml` | shape stability grid heatmap |

### Replay Render

既存 sweep 出力の `summary.csv`、`run_manifest.json`、各 condition directory の `state_archive.npz` から、再シミュレーションなしで比較 plot / 3D replay を生成する場合は `render_shape_stability_grid_replay.py` を使います。

```bash
uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py \
  --input-dir outputs/phase2_103/stage_c_lateral_position_only_dur0p6 \
  --mode both \
  --output-dir /private/tmp/phase2_103_lateral_replay \
  --overwrite
```

`--mode plot-only` は metrics CSV / PNG のみ、`--mode render-only` は 3D grid movie のみ、`--mode both` は両方を生成します。

## 02_phase2_analysis

Phase 2.8 の RUN 固定べん毛数差分析用 batch / dataset 作成 CLI です。

```bash
uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py --dry-run --sample-limit 3
uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py
uv run python scripts/02_phase2_analysis/build_flagella_count_behavior_dataset.py
```

標準設定:

| config | 用途 |
| --- | --- |
| `conf/phase2_analysis/flagella_count_behavior_dataset.yaml` | 標準 dataset |
| `conf/phase2_analysis/flagella_count_behavior_dataset_center_prefix.yaml` | center-priority 付着点 seed 比較 |

raw sample は `step_summary.csv`、`trajectory.csv`、`state_archive.npz` を保存します。後から 3D / 2D render を作る場合は次を使います。

```bash
uv run python scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py \
  --dataset-dir outputs/phase2_analysis/flagella_count_behavior/datasets/fc_nf1_2_3_6_as3_ps3_dur0p5
```

dataset の分布 plot は次で生成します。

```bash
uv run python scripts/02_phase2_analysis/plot_flagella_count_behavior_distributions.py \
  --dataset-id fc_nf1_2_3_6_as3_ps3_dur0p5
```
