# プロジェクト概要
顕微鏡動画からバクテリアを検知・追跡し、姿勢推定とお尻点推定を行うパイプラインです。Phase1では重心トラッキングとデバッグ用オーバーレイ生成が中心です。

# セットアップ
- Python 3.11 推奨。uv を利用して依存関係を解決します。
- 初回: `uv sync` で依存をインストール。
- フック有効化（必要に応じて）: `git config core.hooksPath .githooks`

# 実行例
```
uv run python -m scripts.01_tracking_butt --config conf/config.yaml data.video_path=データパス
```
- `data.video_path` は CLI から `key=value` で上書き可能。
- 反転は初期フレームで自動判定し全フレーム固定。`threshold.invert` を config で明示的に固定してもよい。
- `tracking_butt.save.contour` は `--save-contour` で有効化。

# Phase2 (simulation → projection) 実行例
```
uv run python -m scripts.10_simulate_swimming --config conf/sim_swim.yaml --duration-s 5.0 --fps-out 50 --render-flagella
```
- 上記の他、`flagella.n_flagella=8` のような `key=value` overrides も併用可能。
- 出力（Phase2）: `sim/trajectory.csv`, `render/swim3d.mp4`, `render/swim3d_final.png`, `render2d/projection.mp4`, `render2d/frames/*.png`, `run.log`, `manifest.json`
- **tracking/* は出力しません（Phase2のデータ生成では不要）**

# 出力
- `outputs/YYYY-MM-DD/HHMMSS/` 配下に `run.log`, `manifest.json`, `tracking/track.csv`, `butt.json`, `qc.json`, `overlay.mp4` などを保存。
- 初期フレームの検知スナップショット `tracking/initial_detection.png` を出力。

# テスト
- すべて: `uv run pytest`
- スタイル: `uv run ruff check`
- 合成トラッキング検証は `tests/test_synthetic_tracking.py` を参照。
