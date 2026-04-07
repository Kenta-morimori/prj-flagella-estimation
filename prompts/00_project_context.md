# Project Context（全体像・目的）

## 概要・目的
サルモネラなどのバクテリアは、螺旋状の繊維構造であるべん毛を複数本有し、べん毛が束になりながら回転することで推進力を得て遊泳する。べん毛本数は運動様式と関係する。
本プロジェクトの目的は、サルモネラの遊泳を計測した顕微鏡動画像データから、個体ごとのべん毛本数を推定することである。
予測には機械学習を用いる。学習データは、3D物理シミュレーションを2D投影した擬似顕微鏡像データを生成して用いる。

## フェーズ構成
- Phase1: セットアップ
- Phase2: 物理シミュレーション（3D→2D投影データ生成）
- Phase3: 菌体検出（動画→個体クリップ生成）
- Phase4: 画像認識モデル開発（学習・評価・推定）
- Phase5+: 推定結果の可視化

## ディレクトリと責務（最重要）
- `scripts/` は実行用（CLI入口）であり、設定読み込み・入出力・実行順の制御のみを担う。アルゴリズム/処理本体は `src/` に置く。
- `src/` は再利用可能な実装本体（ライブラリ）。
- `conf/` は設定ファイル（YAML）。
- `outputs/` は実行成果物。各実行ごとに `outputs/<run_id>/` を作成し、`run.log` と `manifest.json` を必ず残す。

---

## MVP制約（必須）
MVPでは以下のみ許容し、それ以外は ValueError とする。
- `body.n_prism = 3`
- `0 <= flagella.n_flagella <= 9`

---

## ログと出力（必須）
- `step_summary_full.csv` は生成しない（廃止）。
- `step_summary.csv` のみに統一し、`outputs/<run_id>/sim/step_summary.csv` に出力する。
- `step_summary.csv` は **1 step = 1 row** の時系列ログであり、解析に必要な情報はすべてこのCSVへ集約する。

## 実行不可時の運用（必須）
- テスト・シミュレーション・コマンドが環境要因で実行できない場合でも、作業差分は必ず commit/push する。
- commit メッセージまたはPR本文に、未実行項目と理由（例: 環境制約、依存不足、権限不足）を明記する。
- 実行結果の代替として、実施済みの静的確認（差分確認、型/linters の結果など）を残す。

---

## テスト方針（必須）
手動確認（動画目視、CSVを人が読む等）に依存しない。
重要な物理制約・形状制約・長時間安定性は **pytest（tests/配下）**で自動検証し、CIで常時実行する。

### 長時間（multi-step）テストの必須化
「初期は安定だが途中で崩れる」失敗モードを取り逃さないため、1 step テストは不十分。
以下の検証は **N=1000〜2000 step** を必須とする。

### テストは原則 “ライブラリ層を直接呼ぶ”
subprocess で scripts を叩く方式は遅く壊れやすい。
原則として `src/` の ModelBuilder / DynamicsEngine / simulate loop を直接呼び、Python内でメトリクスを計算して assert する。
（初期段階の回帰防止として subprocess + step_summary.csv を併用するのは可。将来的にライブラリ直接呼び出しへ置換。）

---

## べん毛・フック実装要件（論文準拠MVP）
目的：トルク増大時に「螺旋が潰れて団子になる」を防ぎ、論文の “nearly rigid” な螺旋を維持する。

### 1) べん毛の幾何（論文固定値）
- Helix pitch（軸方向ピッチ）: **2.5 b**
- Helix diameter: **0.5 b**（半径 R=0.25 b）
- 隣接ビーズ距離（3D距離＝平衡ばね長）: **0.58 b**

### 2) べん毛の角度条件（論文固定値）
- 平衡 bending 角（normal）: **142°**
- 平衡 torsion 角（normal）: **−60°**
（semicoiled / curly1 への遷移を扱う場合は Table 1 の値へ切替）
※MVPではまず normal の固定値で「螺旋を潰さない」ことを優先する。

### 3) 初期形状の生成方針（必須）
初期形状は「初期角＝平衡角」が成り立つように生成する。
つまり、べん毛の各要素で
- bending angle ≈ 142°
- torsion angle ≈ −60°
を満たす点列を生成し、その上で
- 隣接3D距離 ≈ 0.58b
- pitch ≈ 2.5b
- radius ≈ 0.25b
も満たすことをチェックする。

### 4) 形状チェック（tests/で必須）
初期生成直後に、べん毛について以下を自動検証する：
- `bond_len_flag_intra` が 0.58b 近傍（mean±2%、min/max±5%）
- helix pitch が 2.5b 近傍（推定値が許容内）
- helix radius が 0.25b 近傍（推定値が許容内）
- bending/torsion の平均角が論文固定値近傍（例：±数度）

### 5) フック要件（MVP）
- フック長（attach–flag0）: **0.25 b**
- フックは3点（attach, flag0, flag1）で扱い、角度条件付きで曲げ拘束（交差防止）を入れる。
- フック長は long-step（N>=2000）で維持されることを tests/ で担保する。

---

## べん毛基部接線の「後方整列」要件（MVP・必須）
目的：初期状態で flagellum が菌体の後方へ伸びるようにし、初期条件の不自然さ（前方へ伸びる等）を除去する。

### 定義（固定）
- body長軸 `body_axis`：
  - bodyビーズ座標のPCA第1主成分（単位ベクトル）を用いる。
- 後方方向 `rear_dir`：
  - `rear_dir = -body_axis`
- べん毛基部接線 `tangent0`：
  - `tangent0 = normalize(flag1_pos - flag0_pos)`
- 整列条件：
  - `angle(tangent0, rear_dir) <= 10°`（初期状態で満たす）

### 実装方針（最小破壊）
- べん毛のローカル座標系の軸 `axis_u` を `rear_dir` に一致させる。
- 初期点列生成（ヘリックス点列）は既存の処理を使いつつ、「生成座標を `axis_u/axis_v/axis_w` で回転」させる。

### 自動検証（tests/で必須）
- 初期形状生成直後（step0）の全べん毛について、`angle(tangent0, rear_dir) <= 10°` をassertする。

---

## Motor switching（★追加・決定）
目的：run状態（polymorph state固定）を全時間で維持できるようにする。

- 追加パラメータ：`motor.enable_switching`（bool）
  - デフォルト：false（run固定）
  - false のとき：run/tumble更新（多形切替）を行わない
- `step_summary.csv` には `flag_state_changed`（bool）を追加し、
  `motor.enable_switching=false` のとき常に False であることを multi-step pytest で担保する。