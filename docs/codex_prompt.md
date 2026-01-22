````markdown
# Phase1 Codex Prompt（tracking → butt point）

## 目的
顕微鏡動画から、菌体の検知・追跡（track_id付与）・姿勢推定（楕円近似）・お尻点（butt point）推定を行い、デバッグ用の重畳動画と解析用の成果物（CSV/JSON）を出力する。  
Phase1では、べん毛の長さ・本数推定は行わない（擬似べん毛の重畳は行ってよい）。

---

## 改善用追加指示（重心トラッキング強化）
- 遊泳する複数バクテリアを安定トラッキングし、重心ベースでID/フレームをcsv出力・重畳描画する改善を行う。お尻推定は実装してもよいが、評価対象は重心トラッキングのみ（将来評価する可能性を明記）。
- 初期フレームで菌体は背景と反転しないと決め打ち（撮影中に全体フォーカスのズレは想定しない）。極性が逆転するケースはトラッキング不可扱いとし、全フレーム同じ極性で処理する（初期フレームで一度だけ極性判定）。
- `um_per_px` と `bac_short_axis_length_um` から得る期待短径を使い、巨大成分や太すぎる成分を除外。境界接触や面積上限も組み合わせて画面サイズ級の誤検知を避ける。
- 複数オブジェクトは重ならない前提だが、近接時にもIDが飛ばないよう距離ゲーティングで対応付けする。
- 合成データによる自動テストを追加する: 512x512 のフレームで黒寄りグレーの楕円（短軸 = `bac_short_axis_length_um` を `um_per_px` でピクセル換算）がジグザグ移動するAVI。200フレーム、速度5–10px/frame 相当、縦に領域を分けて約10体が重ならず遊泳し、長短軸比は最大10倍程度までばらつかせる。`data/test_*.avi` を生成し、重心トラッキング継続とID保持を確認する。
- 暗い菌体を失わないため、前処理はまず `bg_subtract` または `none` を優先して検討し（tophatは暗い菌体を消しやすい）、カーネルは背景ムラより小さめ（例: 11〜17）にする。微粒子を抑えるため `min_area_px` を8〜10に上げ、誤検知を避けたい場合は `max_area_frac` を0.05程度まで調整する。
- すべての関数に日本語docstringを付与し、概要と Args/Returns を記載する。
- コミット前に ruff / pytest を自動実行し、CIを通す。

---

## 実行条件（必須）
### 1) コミットベース実行の強制
- 実行時、Gitが以下の状態でなければエラーで終了すること。
  - Git管理下である
  - HEADが存在する（初回コミット済み）
  - `git status --porcelain` が空（未コミット変更なし）
- 実行できる場合、commit id（short/long）をログとmanifestに必ず記録すること。

### 2) 出力ディレクトリ規則
- 実行結果は必ず `outputs/YYYY-MM-DD/HHMMSS/` に保存すること（JST）。
- logファイルも同ディレクトリに保存すること（例：`run.log`）。

### 3) 既存の実行コンテキストを利用する
- `flagella_estimation.core.run_context.init_run()` を必ず使用し、出力先・ログ・manifest作成を統一すること。
- 既存構造・既存テストを壊さないこと（CIが通ること）。

---

## 入力
- 動画ファイル（例：mp4/avi）
- 設定ファイル：`conf/config.yaml`（YAML）
- 実行スクリプト：`scripts/01_tracking_butt.py`

---

## 出力（変更禁止）
`outputs/YYYY-MM-DD/HHMMSS/` 配下に最低限以下を生成すること。

### 必須
- `run.log`
- `manifest.json`
- `tracking/track.csv`
- `tracking/butt.json`
- `tracking/qc.json`
- `tracking/overlay.mp4`

### オプション（輪郭保存が有効な場合のみ）
- `tracking/contours/` 以下に輪郭データを保存する（詳細は後述）

---

## 追加の要求（config管理）
### 尻推定で使用する特徴量の列挙
- 尻推定に利用する特徴量は **configに列挙**すること。
- 列挙された特徴量のみを計算し、`track.csv` に列として保存すること。
- 将来のべん毛推定（Phase2/3）でも同名特徴量を流用する可能性が高い旨をコメントとして残すこと（実装は不要）。

---

## 輪郭保存（デフォルトfalse＋CLIで上書き）
- `tracking_butt.save.contour` は **デフォルトfalse** とする。
- `scripts/01_tracking_butt.py` に `--save-contour` を追加し、指定時のみ実行中に `tracking_butt.save.contour=true` に上書きして保存を有効化すること。
- 保存形式は `.npy` を推奨する（OpenCV contour相当の (N,2) 座標配列）。

ファイル命名規則（推奨・固定）：
- `tracking/contours/frame_{frame:06d}_track_{track_id:04d}.npy`

---

## 座標系と表現（固定）
- 画像座標：xは右向き、yは下向き（OpenCV標準）。
- `cx, cy` はピクセル座標（float可）。
- `theta` はラジアンで統一する。`theta=0` は +x 方向。
- `major, minor` はピクセル長。
- `vx, vy` は px/frame（差分でよい）。

---

## 実装範囲（Phase1でやること）
1. 各フレームから菌体候補を検知する（古典的手法のみ）
2. フレーム間で最近傍対応付けによりトラッキングし `track_id` を付与する
3. 検知結果から姿勢（theta、major/minor）を推定する
4. お尻点を推定し、時間窓で安定化する
5. QC指標を計算して保存する
6. overlay動画を生成する（菌体輪郭/楕円・track_id・お尻点・擬似べん毛）

---

## 実装制約（必須）
- 外部モデルのダウンロード、学習済みモデルの導入は禁止（OpenCV等の古典処理のみ）。
- 既存のディレクトリ構造・既存の出力名は変更しないこと。
- `ruff` / `pytest` を通すこと（既存のCIを壊さない）。
- importは `src/flagella_estimation/` 配下を前提とすること。
- 関数docstringは日本語で概要と Args/Returns を記載すること。

---

## 検知（detection）の最小仕様
### 推奨処理（最小・壊れにくさ優先）
- grayscale → 軽い平滑化 → 二値化（Otsuを基本。初期フレームで一度だけ極性を判定して固定。極性が逆転するものはトラッキング不可） → 連結成分 → 面積フィルタ → 楕円フィット（可能なら）

### 検知で保持する情報（最低限）
- `frame`
- `cx, cy`
- `theta, major, minor`
- `area`
- `bbox`（x,y,w,h）（overlay・デバッグ用）

※ 楕円フィット不能な場合は、その候補を捨てるか、`theta/major/minor` を欠損扱いにしてよい。ただし後段が破綻しないよう `is_valid` 等で扱うこと。

---

## トラッキング（tracking）の最小仕様
- 最近傍 + 距離ゲーティング（`max_link_distance`）で対応付けすること。
- 欠損は許容してよい（Phase1では簡易でよい）。
- 速度 `vx, vy` は `cx,cy` の差分から計算してよい。

---

## お尻点推定（butt estimation）の最小仕様（固定）
### 1) お尻候補点の計算
- 長軸方向の単位ベクトル `u=(cos(theta), sin(theta))` を用い、
  - `p1 = c + (major/2) * u`
  - `p2 = c - (major/2) * u`
  を候補端点とする。

### 2) どちらが尻かの決定（速度整合）
- 速度 `v=(vx,vy)` を用い、原則として **進行方向と反対側の端点を尻**とする。
- 具体例：`dot(v, p1-c)` と `dot(v, p2-c)` を比較し、より小さい方（より逆向き）を尻とする。
- `speed < freeze_speed_thresh` の場合、`frozen=true` とし、尻点は前フレームの値を保持する。

### 3) 時間窓での安定化
- `smooth_window` フレームの窓で、尻の選択（p1側/p2側）の符号を多数決し、反転チラつきを抑えること。

### 4) べん毛方向（重畳用）
- `flagella_dir = normalize(butt_point - center)` を基本とし、尻点から `flagella_dir` 方向へ描画すること。

---

## 特徴量（configで列挙し、列として保存）
### 方針
- configに列挙された特徴量のみ計算し、`track.csv` に保存すること。
- Phase1で最低限実装してよい候補（例）：
  - `speed`
  - `vel_axis_dot`（`dot(v,u)`）
  - `heading_change`（速度方向の角度差分。実装できる範囲でよい）
- 未実装の特徴量名が列挙された場合は、エラーにせず警告ログを残し、その列は出力しない（またはNaNで埋める）。どちらかに統一すること。

---

## 成果物フォーマット（固定）

### `tracking/track.csv`（必須）
- 行：`frame, track_id` 単位
- 最低カラム（必須）：
  - `frame, track_id, cx, cy, theta, major, minor, vx, vy, is_valid`
- 追加カラム（config列挙の特徴量）：
  - 例：`speed, vel_axis_dot, ...`

### `tracking/butt.json`（必須）
- 構造：`track_id -> frame -> {x, y, conf, frozen, flagella_dir_x, flagella_dir_y}`
- `conf` は簡易でよい（例：`abs(dot(v,u))` の正規化、またはspeedに基づくスコア）。

### `tracking/qc.json`（必須）
- 構造：`track_id -> {frames, missing_rate, speed_spike_rate, ...}`
- Phase1では最低限 `frames` と `missing_rate` があればよい。

### `tracking/overlay.mp4`（必須）
- 描画要素（最低限）：
  - track_id
  - 中心点
  - お尻点
  - 楕円（または長軸線）
  - 擬似べん毛（尻点から `flagella_dir` に線分/曲線）
- 擬似べん毛は固定の仮パラメータでよい（Nや長さはPhase1では推定しない）。
- 菌体の短辺長(`bac_short_axis_length_um`)を右下にスケールバーとして追加する．想定している長さであるかを確認するため．

---

## config案（必要キー）
`conf/config.yaml` に最低限以下を用意し、コード側はこのキーを参照すること。

例：
```yaml
data:
  video_path: "data/sample.mp4"
  fps: 0

output:
  base_dir: "outputs"

tracking_butt:
  tracking:
    max_link_distance: 30.0
  butt_estimation:
    smooth_window: 5
    freeze_speed_thresh: 0.5
    features:
      - speed
      - vel_axis_dot
  save:
    contour: false
````

---

## CLI仕様（必須）

`scripts/01_tracking_butt.py` に以下を用意すること。

* `--config conf/config.yaml`（デフォルトありでよい）
* `--save-contour`（指定時に `tracking_butt.save.contour=true` に上書き）

---

## DoD（完了条件）

* 未コミット変更がある状態では実行が停止し、コミットを促すエラーになる。
* コミット済み状態で実行すると `outputs/YYYY-MM-DD/HHMMSS/` が生成され、必須成果物が揃う。
* `run.log` と `manifest.json` に commit id が残る。
* overlayで track_id とお尻点が極端にチラつかない（完全一致は不要）。
* CI（ruff/pytest）が通る。

---

## Codexへの指示（出力形式）

* 変更・追加したファイル一覧（パス）を先に列挙すること。
* その後、変更した各ファイルの最終版全文を提示すること。
* 既存ファイルの削除や大規模な構造変更は禁止（必要がある場合は理由と代替案を提示すること）。
