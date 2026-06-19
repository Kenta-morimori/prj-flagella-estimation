# scripts ディレクトリ

## 01_simulate_swimming/
Phase2 実行系を集約したディレクトリです。

- `01_simulate_swimming/01_simulate_swimming.py`:
  3D遊泳シミュレーションを行い、2D投影と可視化成果物を出力するCLI
- `01_simulate_swimming/run_motor_scale_sweep.py`:
  torque×local scale の sweep を実行するCLI
- `01_simulate_swimming/plot_motor_scale_collapse_heatmap.py`:
  sweep 結果CSVから heatmap を生成するCLI

実行例:
```
uv run python -m scripts.01_simulate_swimming
uv run python -m scripts.01_simulate_swimming.run_motor_scale_sweep --help
uv run python -m scripts.01_simulate_swimming.plot_motor_scale_collapse_heatmap --help
```

Phase 2.6 の single flagellum 螺旋維持 gate では、`run_motor_scale_sweep` の `--dt-star` と `--stub-mode full_flagella` を使って内部刻みと full flagellum 条件を明示します。

### `scripts.01_simulate_swimming` のCLI指定方法

- `duration_s` やトルク、`dt_star` など設定キーは `KEY=VALUE` 形式の override で指定します。

例:
```bash
# 実行時間を 0.05s に変更
uv run python -m scripts.01_simulate_swimming time.duration_s=0.05

# トルクを 3.0e-21 N*m に変更
uv run python -m scripts.01_simulate_swimming motor.torque_Nm=3.0e-21

# 実行時間・トルク・内部刻みを同時指定
uv run python -m scripts.01_simulate_swimming \
  time.duration_s=0.05 \
  motor.torque_Nm=3.0e-21 \
  time.dt_star=1.0e-4
```

`scripts.plot_motor_scale_collapse_heatmap` は Phase2 実行系の責務境界を明確にするため削除済みです。heatmap 生成は `scripts.01_simulate_swimming.plot_motor_scale_collapse_heatmap` を使用してください。

## 02_phase2_analysis/
Phase 2.8 の RUN 固定べん毛数差分析用 batch / dataset 作成CLIです。

- `02_phase2_analysis/run_flagella_count_behavior_sweep.py`:
  `n_flagella = 1, 2, 3, 6` × `seed = 0` の 4 samples をまとめて実行し、sample ごとの raw output と `run_manifest.json` を出力する CLI
- `02_phase2_analysis/build_flagella_count_behavior_dataset.py`:
  `run_manifest.json` から `summary.csv`、`qc_summary.csv`、`timeseries/<sample_id>.csv`、`dataset_manifest.json` を作成する CLI
- `02_phase2_analysis/render_flagella_count_behavior_sample.py`:
  raw sample archive から 3D/2D render を後出し生成する CLI

標準 4 samples の実行例:

```bash
uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py
uv run python scripts/02_phase2_analysis/build_flagella_count_behavior_dataset.py
```

raw sample からの再描画:

```bash
uv run python scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py \
  --sample-dir outputs/phase2_analysis/flagella_count_behavior/runs/fc_nf1_2_3_6_seed1_dur0p5/samples/nf01_seed000 \
  --output-dir outputs/phase2_analysis/flagella_count_behavior/replays/nf01_seed000
```

補助オプション:

```bash
# 実行対象を一部に絞って確認する
uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py --dry-run --sample-limit 2

# 既存 sample / dataset を上書きして再生成する
uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py --overwrite
uv run python scripts/02_phase2_analysis/build_flagella_count_behavior_dataset.py --overwrite

# sweep 実行時に analysis 設定と simulation 設定を override する
uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py \
  dataset_id=fc_nf1_2_3_6_seed1_dur1p0 \
  run_batch_id=fc_nf1_2_3_6_seed1_dur1p0 \
  output.run_batch_dir=outputs/phase2_analysis/flagella_count_behavior/runs/fc_nf1_2_3_6_seed1_dur1p0 \
  time.duration_s=1.0

# 上記 run から dataset を作る
uv run python scripts/02_phase2_analysis/build_flagella_count_behavior_dataset.py \
  dataset_id=fc_nf1_2_3_6_seed1_dur1p0 \
  output.run_batch_dir=outputs/phase2_analysis/flagella_count_behavior/runs/fc_nf1_2_3_6_seed1_dur1p0 \
  output.dataset_dir=outputs/phase2_analysis/flagella_count_behavior/datasets/fc_nf1_2_3_6_seed1_dur1p0
```

`02_phase2_analysis` の override は `KEY=VALUE` 形式です。`dataset_id`、`run_batch_id`、`output.run_batch_dir`、`output.dataset_dir` は Phase2 analysis 側の設定として扱います。`time.duration_s`、`time.dt_star`、`motor.torque_Nm`、`render.*` などは simulation 設定の省略形として扱い、各 sample の `base_overrides` に反映します。simulation 側の `output.base_dir` を変えたい場合は `base_overrides.output.base_dir=...` と明示してください。

`runs/<run_batch_id>/samples/<sample_id>/raw/` には、`step_summary.csv` に加えて `trajectory.csv` と `state_archive.npz` を残します。`state_archive.npz` は後から 3D / 2D render を再生成するための状態保存です。

デフォルト設定は `conf/phase2_analysis/flagella_count_behavior_dataset.yaml` を参照してください。標準出力先は以下です。

- `outputs/phase2_analysis/flagella_count_behavior/runs/fc_nf1_2_3_6_seed1_dur0p5/`
- `outputs/phase2_analysis/flagella_count_behavior/datasets/fc_nf1_2_3_6_seed1_dur0p5/`

実行時に使われた analysis 設定は `analysis_config_used.yaml` と各 manifest の `effective_analysis_config` に記録されます。
