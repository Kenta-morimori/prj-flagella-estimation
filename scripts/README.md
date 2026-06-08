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

## 02_detect_bac.py
Phase3: 動画から菌体検出と個体クリップ生成を行うCLIの雛形です。

## 03_train_evaluate.py
Phase4: べん毛本数推定モデルの学習・評価を行うCLIの雛形です。

## 10_overlay_flagella.py
Phase5+: 推定結果を元動画へ重畳表示するCLIの雛形です。
