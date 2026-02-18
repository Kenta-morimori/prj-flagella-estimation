# scripts ディレクトリ

## 01_simulate_swimming.py
Phase2: 3D遊泳シミュレーションを行い、2D投影と可視化成果物を出力するCLIです。
パラメータは`conf/sim_swim.yaml`をデフォルトで参照しています。別のconfigを参照したい場合は、以下の実行コマンドに`--config [configのパス]`を追加ください。

実行例:
```
uv run python -m scripts.01_simulate_swimming
```

パラメータを修正した実行方法
```
uv run python -m scripts.01_simulate_swimming flagella.n_flagella=6
```

## 02_detect_bac.py
Phase3: 動画から菌体検出と個体クリップ生成を行うCLIの雛形です。

## 03_train_evaluate.py
Phase4: べん毛本数推定モデルの学習・評価を行うCLIの雛形です。

## 10_overlay_flagella.py
Phase5+: 推定結果を元動画へ重畳表示するCLIの雛形です。
