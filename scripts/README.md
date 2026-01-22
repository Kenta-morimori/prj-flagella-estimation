# scripts ディレクトリについて

## 01_tracking_butt.py
バクテリアの検知・追跡・姿勢推定・お尻推定を実行するエントリーポイントです。

### 実行例
```
uv run python -m scripts.01_tracking_butt --config conf/config.yaml data.video_path=data/your_video.avi
```

### 主なオプション
- `--config`: 設定ファイルパス（省略時は `conf/config.yaml`）。
- `--save-contour`: 指定すると輪郭 `.npy` を `tracking/contours/` に保存。
- `key=value` 形式で config を上書き可能（例: `data.fps=30`, `data.um_per_px=0.2`）。

### 前提
- 実行時は `PYTHONPATH=src` が通っているか、`uv run` でモジュールを解決してください。
- Gitワークツリーがクリーンでないと実行エラーになります（manifestにcommit情報を残すため）。

