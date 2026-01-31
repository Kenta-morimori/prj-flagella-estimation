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

---

## 10_simulate_swimming.py (Phase2: simulation → projection)
剛体菌体＋螺旋べん毛のシミュレーションを行い、3D/2Dの可視化と軌跡を出力します。tracking/* は出力しません。

### 実行例
```
uv run python -m scripts.10_simulate_swimming --config conf/sim_swim.yaml --duration-s 5.0 --fps-out 50 --render-flagella
```

### 主なオプション
- `--config`: 設定ファイルパス（デフォルト `conf/sim_swim.yaml`）。
- `--duration-s`: time.duration_s の上書き（秒）。
- `--fps-out`: time.fps_out の上書き（フレームレート）。
- `--render-flagella`: べん毛のデバッグ描画を有効化。
- `key=value` 形式で config を上書き可能（例: `flagella.n_flagella=8`）。

### 出力
- `sim/trajectory.csv`
- `render/swim3d.mp4`, `render/swim3d_final.png`（Matplotlib 3Dグリッド表示）
- `render2d/projection.mp4`, `render2d/frames/*.png`
- `run.log`, `manifest.json`

### 前提
- 実行時は `PYTHONPATH=src` が通っているか、`uv run` でモジュールを解決してください。
- Gitワークツリーがクリーンでないと実行エラーになります（manifestにcommit情報を残すため）。
