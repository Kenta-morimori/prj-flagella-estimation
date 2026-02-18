# プロジェクト概要
顕微鏡動画からバクテリアのべん毛本数を推定するための開発リポジトリです。

## セットアップ
- Python 3.11 推奨
- 依存インストール: `uv sync`

## Phase2 実行
`uv run python -m scripts.01_simulate_swimming --config conf/sim_swim.yaml`

実行後、`outputs/YYYY-MM-DD/HHMMSS/` 配下に `run.log` と `manifest.json` を含む成果物が生成されます。
