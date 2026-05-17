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

互換性のため、以下の旧エントリも利用可能です（内部で新モジュールへ委譲）:
- `scripts.run_motor_scale_sweep`
- `scripts.plot_motor_scale_collapse_heatmap`

## 02_detect_bac.py
Phase3: 動画から菌体検出と個体クリップ生成を行うCLIの雛形です。

## 03_train_evaluate.py
Phase4: べん毛本数推定モデルの学習・評価を行うCLIの雛形です。

## 10_overlay_flagella.py
Phase5+: 推定結果を元動画へ重畳表示するCLIの雛形です。
