# Phase 2 multi-run profiles

この directory には，Phase 2 の user-facing な複数条件実行 profile YAML を置く。
開発・診断用の task-specific sweep は `conf/phase2_sweeps/` を使い，既存 config に対する複数条件実行と結果の比較は `conf/phase2_multi_run/` を使う。

## Canonical profile

現時点の標準 profile は次である。

```text
conf/phase2_multi_run/latest_model_torque_shape_stability.yaml
```

この profile は，最新モデルで `motor.torque_Nm` を複数条件に振り，shape stability の定量 plot と replay を作るための設定を 1 枚にまとめる。

## 標準実行コマンド

multi-run，plot，replay は 3 コマンドに分けて実行する。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml overwrite=true

uv run python scripts/01_simulate_swimming/plot_heatmap.py config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml

uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml overwrite=true
```

`run_multi_run.py` は simulation の複数条件実行だけを担当する。plot と replay は，生成済みの `summary.csv` と各 condition の出力を別コマンドで読む。

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
