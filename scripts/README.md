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
- `02_phase2_analysis/plot_flagella_count_behavior_distributions.py`:
  dataset の `summary.csv` から本数別特徴量分布、QC集計、NaN集計、簡易スクリーニングを出力する CLI

標準 4 samples の実行例:

```bash
uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py
uv run python scripts/02_phase2_analysis/build_flagella_count_behavior_dataset.py
```

標準条件では override 指定は不要です。`KEY=VALUE` override は、YAMLを複製せずに `duration_s` などを一時的に変えたい場合の補助機能です。

軽量 sweep を試す場合は、診断と archive state を間引く fast config を使います。

```bash
uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py \
  --config conf/phase2_analysis/flagella_count_behavior_dataset_fast.yaml
```

`runner.step_summary_stride` は `step_summary.csv` と `flag_helix_axis_diagnostics.csv` の記録間隔、`runner.state_stride` は `trajectory.csv` と `state_archive.npz` の保存間隔です。標準 config はどちらも実質 `1` で、従来通り全 step を保存します。
`--stop-on-shape-fail` は shape fail を全記録stepで確認する前提のため、`runner.step_summary_stride=1` のときだけ使用できます。

raw sample からの再描画:

```bash
uv run python scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py \
  --sample-dir outputs/phase2_analysis/flagella_count_behavior/runs/fc_nf1_2_3_6_seed1_dur0p5/samples/nf01_seed000 \
  --output-dir outputs/phase2_analysis/flagella_count_behavior/replays/nf01_seed000
```

replay render の3D出力はデフォルトで `output_sampling.out_all_steps_3d=false` として扱い、`--fps-out-3d` で間引きfpsを指定します。全archive stateを3D描画したい場合のみ `--out-all-steps-3d` を指定してください。2D側は `--fps-out-2d` で指定できます。
3D/2Dのmp4はVSCodeなどで再生しやすいH.264 (`avc1` / `H264`) を優先し、環境にencoderがない場合は従来の `mp4v` にfallbackします。実際に使われたcodecは `manifest.json` の `render_video` に記録されます。

dataset からの特徴量分布可視化:

```bash
uv run python scripts/02_phase2_analysis/plot_flagella_count_behavior_distributions.py \
  --dataset-id fc_nf1_2_3_6_seed1_dur0p5
```

出力は dataset directory 配下の `plots/distributions/`、`plots/qc/`、`analysis/` に保存されます。
既存の分析出力を再生成する場合のみ `--overwrite` を指定してください。

補助オプション:

```bash
# 実行対象を一部に絞って確認する
uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py --dry-run --sample-limit 2

# 既存 sample / dataset を上書きして再生成する
uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py --overwrite
uv run python scripts/02_phase2_analysis/build_flagella_count_behavior_dataset.py --overwrite

# 条件を変える場合は dataset_id / run_batch_id / output path も変える
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
`runner.step_summary_stride`、`runner.state_stride`、`runner.flush_interval_steps`、`runner.sample_order` は Phase2 analysis runner 側の設定です。`runner.sample_order=interleave_n_flagella` を指定すると、seed 条件ごとに `n_flagella` を混ぜて実行します。

既存の `run_batch_dir` / `sample_id` に `step_summary.csv` がある場合、runner は保存済み sample config と、raw内容に影響する runner stride (`runner.step_summary_stride` / `runner.state_stride`) が今回の設定と一致するときだけ既存rawを再利用します。条件が異なる場合は、古いrawと新しいmanifest metadataの混在を避けるため停止します。同じ出力先で条件を変えて再生成する場合は `--overwrite` を指定してください。

`runs/<run_batch_id>/samples/<sample_id>/raw/` には、`step_summary.csv` に加えて `trajectory.csv` と `state_archive.npz` を残します。`state_archive.npz` は後から 3D / 2D render を再生成するための状態保存です。

デフォルト設定は `conf/phase2_analysis/flagella_count_behavior_dataset.yaml` を参照してください。標準出力先は以下です。

- `outputs/phase2_analysis/flagella_count_behavior/runs/fc_nf1_2_3_6_seed1_dur0p5/`
- `outputs/phase2_analysis/flagella_count_behavior/datasets/fc_nf1_2_3_6_seed1_dur0p5/`

実行時に使われた analysis 設定は `analysis_config_used.yaml` と各 manifest の `effective_analysis_config` に記録されます。
