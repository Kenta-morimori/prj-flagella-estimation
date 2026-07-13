# Phase 2 multi-run profiles

この directory には，Phase 2 の user-facing な複数条件実行 profile YAML を置く。
開発・診断用の task-specific sweep は `conf/phase2_sweeps/` を使い，既存 config に対する複数条件実行と結果の比較は `conf/phase2_multi_run/` を使う。

## Canonical profiles

現時点の標準 profile は次である。

| config | 用途 |
| --- | --- |
| `conf/phase2_multi_run/latest_model_torque_shape_stability.yaml` | 最新モデルの torque 複数条件 shape stability 比較 |
| `conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml` | Issue #71 の RUN 固定べん毛本数差 diagnostic dataset v0 |

どちらも run / plot / replay の設定を 1 枚にまとめる。
`flagella_count_behavior_diagnostic.yaml` はさらに `dataset:` section を持ち，dataset 作成にも同じ config を使う。

## 標準実行コマンド

multi-run，plot，replay は 3 コマンドに分けて実行する。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml overwrite=true

uv run python scripts/01_simulate_swimming/plot_heatmap.py config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml

uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml overwrite=true
```

`run_multi_run.py` は simulation の複数条件実行だけを担当する。plot と replay は，生成済みの `summary.csv` と各 condition の出力を別コマンドで読む。

Issue #71 の diagnostic dataset は，同じ profile から dataset 作成と分布分析まで続ける。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py config=conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml overwrite=true

uv run python scripts/01_simulate_swimming/plot_heatmap.py config=conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml

uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py config=conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml overwrite=true

uv run python scripts/02_phase2_analysis/build_dataset.py config=conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml overwrite=true

uv run python scripts/02_phase2_analysis/plot_distributions.py --dataset-id fc_nf1_2_3_6_as3_ps3_torque2p0_dur0p5 --overwrite
```

36 sample の本実行は長時間 run として扱う。確認だけなら `dry_run=true sample_limit=5` を使う。

## 出力先

`output.base_dir` は campaign 全体の run root である。
`latest_model_torque_shape_stability.yaml` では次を固定 run root にする。

```yaml
output:
  base_dir: outputs/phase2_multi_run/latest_model_torque_shape_stability
  timestamp_subdir: false
```

`output.timestamp_subdir: false` の profile では，`plot_heatmap.py config=...` と `render_shape_stability_grid_replay.py config=...` が `output.base_dir` を参照する。
そのため，通常は `summary_csv=...` や `run_dir=...` を毎回指定しなくてよい。

代表的な出力は次の構成になる。

```text
outputs/phase2_multi_run/latest_model_torque_shape_stability/
  summary.csv
  run_manifest.json
  <condition_id>/
    step_summary.csv
    trajectory.csv
    state_archive.npz
    run.log
    manifest.json
  plots/
    plot_data.csv
    <metric>_vs_<axis>.png
  replay/
    ...
```

`plots/` と `replay/` は，それぞれ plot / replay コマンドを実行したときに作られる。`run_multi_run.py` だけでは作られない。

`flagella_count_behavior_diagnostic.yaml` では run root と dataset 出力先を分ける。

```yaml
output:
  base_dir: outputs/phase2_multi_run/flagella_count_behavior_diagnostic
  timestamp_subdir: false

dataset:
  dataset_id: fc_nf1_2_3_6_as3_ps3_torque2p0_dur0p5
  output_dir: outputs/phase2_analysis/flagella_count_behavior/datasets/fc_nf1_2_3_6_as3_ps3_torque2p0_dur0p5
```

## plot 設定

`plot` section は，`plot_heatmap.py` がどの metric をどの軸で可視化するかを決める。

```yaml
plot:
  default_x_axis: torque
  default_y_axis: null
  metrics:
    - first_fail_t_s
    - hook_len_rel_err_max
    - max_flag_bond_rel_err
    - axis_center_to_body_roll_ratio_mean
```

`plot.metrics` に列挙した metric が PNG 出力対象になる。
各 metric は `summary.csv` の列名として存在している必要がある。
正規化済みの入力行は `plots/plot_data.csv` にも保存される。

`plot.default_y_axis: null` の場合，条件軸が 1 つなので heatmap ではなく line plot を出す。

```text
plots/first_fail_t_s_vs_torque.png
plots/hook_len_rel_err_max_vs_torque.png
```

2 軸 heatmap を出したい場合は，`sweep.axes` に 2 つ目の軸を追加し，`plot.default_y_axis` にその軸名を指定する。

```yaml
sweep:
  axes:
    torque:
      key: motor.torque_Nm
      short_name: torque
      values: [1.5e-20, 2.0e-20, 2.5e-20]
    viscosity:
      key: fluid.viscosity_Pa_s
      short_name: viscosity
      values: [0.001, 0.002, 0.003]

plot:
  default_x_axis: torque
  default_y_axis: viscosity
  metrics:
    - first_fail_t_s
```

この場合は次のような heatmap PNG を出す。

```text
plots/first_fail_t_s_heatmap.png
```

条件軸が 3 つ以上ある profile では，表示しない軸を `plot.filter_axes` で固定する。
`flagella_count_behavior_diagnostic.yaml` は `n_flagella x attach_seed` の heatmap を描くため，`phase_seed` を固定している。

```yaml
plot:
  default_x_axis: n_flagella
  default_y_axis: attach_seed
  filter_axes:
    phase_seed: "0"
```

`filter_axes` を指定しないまま 3 軸以上の profile を plot しようとすると，同じ heatmap cell へ複数条件が重なるためエラーにする。

## dataset 設定

`dataset` section は，`scripts/02_phase2_analysis/build_dataset.py` が読む。
`run_multi_run.py` の `run_manifest.json` にある `conditions` を dataset sample に変換し，`summary.csv`，`qc_summary.csv`，`timeseries/<sample_id>.csv`，`dataset_manifest.json` を生成する。

```yaml
dataset:
  dataset_id: fc_nf1_2_3_6_as3_ps3_torque2p0_dur0p5
  feature_schema: conf/phase2_analysis/flagella_count_behavior_features.yaml
  output_dir: outputs/phase2_analysis/flagella_count_behavior/datasets/fc_nf1_2_3_6_as3_ps3_torque2p0_dur0p5
  sample_id_template: "nf{n_flagella:02d}_as{attach_seed:03d}_ps{phase_seed:03d}"
  timeseries_sampling: all_steps
```

`dataset.dataset_id` や `dataset.output_dir` は CLI override できる。

```bash
uv run python scripts/02_phase2_analysis/build_dataset.py \
  config=conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml \
  dataset.dataset_id=my_diagnostic \
  dataset.output_dir=/private/tmp/my_diagnostic_dataset \
  overwrite=true
```

## replay 設定

`replay` section は，`render_shape_stability_grid_replay.py` の出力を決める。

```yaml
replay:
  mode: both
  fps_out_3d: 10.0
  output_subdir: replay
```

`output_subdir` は `output.base_dir` からの相対 directory である。
上の例では replay 出力先は次になる。

```text
outputs/phase2_multi_run/latest_model_torque_shape_stability/replay/
```

replay は `summary.csv` と各 condition の `state_archive.npz` / `trajectory.csv` を読む。
replay 出力種別の細かい on/off 指定や script 名の整理は，後続タスクで扱う。

## CLI override

profile の値は `KEY=VALUE` 形式で上書きできる。

例: 実行時間だけ 0.2 s に短縮する。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml \
  base_overrides.time.duration_s=0.2 \
  overwrite=true
```

例: torque 条件を差し替える。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml \
  sweep.axes.torque.values=1.0e-20,1.5e-20,2.0e-20 \
  overwrite=true
```

CLI override で条件軸を差し替えた場合でも，`summary.csv` には `axis_<name>_index` と `axis_<name>_label` が保存される。
plot はこの summary 側の軸 label を優先して読む。
