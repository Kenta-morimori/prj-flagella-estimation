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

# 出力
- `outputs/YYYY-MM-DD/HHMMSS/` 配下に `run.log`, `manifest.json`, `tracking/track.csv`, `butt.json`, `qc.json`, `overlay.mp4` などを保存。
- 初期フレームの検知スナップショット `tracking/initial_detection.png` を出力。

# テスト
- すべて: `uv run pytest`
- スタイル: `uv run ruff check`
- 合成トラッキング検証は `tests/test_synthetic_tracking.py` を参照。

