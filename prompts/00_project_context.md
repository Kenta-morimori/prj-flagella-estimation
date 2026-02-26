# Project Context（全体像・目的）

## 概要・目的
サルモネラなどのバクテリアは、螺旋状の繊維構造であるべん毛を複数本有し、べん毛が束になりながら回転することで推進力を得て遊泳する。べん毛本数は運動様式と関係する。

本プロジェクトの目的は、サルモネラの遊泳を計測した顕微鏡動画像データから、個体ごとのべん毛本数を推定することである。
予測には機械学習を用いる。学習データは、3D物理シミュレーションを2D投影した擬似顕微鏡像データを生成して用いる。
Phase3でYOLO等により菌体を検出・クリップ化し、Phase4でCNN/ViT等により学習・推定する。

## フェーズ構成
- Phase1: セットアップ（環境構築・実行基盤整備）※前提作業として扱う
- Phase2: 物理シミュレーション（3D→2D投影データ生成）
- Phase3: 菌体検出（動画→個体クリップ生成）
- Phase4: 画像認識モデル開発（学習・評価・推定）
- Phase5+: 推定結果の可視化（元動画への重畳など）

## ディレクトリと責務（最重要）
- `scripts/` は実行用（CLI入口）であり、設定読み込み・入出力・実行順の制御のみを担う。アルゴリズム/処理本体は `src/` に置く。
- `src/` は再利用可能な実装本体（ライブラリ）を置く。
- `conf/` は設定ファイル（YAML）を置く。
- `outputs/` は実行成果物を保存する。各実行ごとに `outputs/YYYY-MM-DD/HHMMSS/` を作成し、`run.log` と `manifest.json` を必ず残す。

## scripts（入口）
- `scripts/01_simulate_swimming.py`（Phase2）: 3Dシミュレーション結果を2Dへ投影し、擬似顕微鏡像（動画/画像列）とメタ情報を生成する。
- `scripts/02_detect_bac.py`（Phase3）: 動画像から菌体を検出し、個体ごとのクリップを生成する。
- `scripts/03_train_evaluate.py`（Phase4）: べん毛本数推定モデルを学習・評価する。
- `scripts/10_overlay_flagella.py`（Phase5+）: 元動画へ推定結果を重畳し可視化する。

## src（実装本体）
- `src/sim_swim/`: Phase2の実装
- `src/flagella_estimation/`: Phase3-4の実装

---

## MVP制約（必須）
MVPでは以下のみ許容し、それ以外は ValueError とする。
- `body.n_prism = 3`
- `0 <= flagella.n_flagella <= 3`

---

## ログと出力（必須）
- `step_summary_full.csv` は生成しない（廃止）。
- `step_summary.csv` のみに統一し、`outputs/<run_id>/sim/step_summary.csv` に出力する（`sim_debug/` は使わない）。
- `step_summary.csv` は **1 step = 1 row** の時系列ログであり、解析に必要な情報はすべてこのCSVへ集約する。

---

## テスト方針（必須）
手動確認（動画目視、CSVを人が読む、など）に依存しない。
重要な物理制約・形状制約・長時間安定性は **pytest（tests/配下）** で自動検証し、CIで常時実行する。

### 1) 長時間（multi-step）テストの必須化
「初期は安定だが途中で崩れる」失敗モードを取り逃さないため、1 step テストは不十分。
以下の検証は **N=1000〜2000 step** を必須とする。

### 2) テストは原則 “ライブラリ層を直接呼ぶ”
subprocess で scripts を叩く方式は遅く壊れやすい。
原則として `src/` の ModelBuilder / DynamicsEngine / simulate loop を直接呼び、Python内でメトリクスを計算して assert する。
（初期段階の回帰防止として subprocess + step_summary.csv を併用するのは可。将来的にライブラリ直接呼び出しへ置換。）

---

## べん毛・フック実装要件（論文準拠MVP）
（目的：フック長とべん毛ボンド長を長ステップで維持し、安定に回転できるモデルを作る）

### 1) べん毛（Flagellum）幾何（MVP）
- Helix pitch（軸方向ピッチ）: 2.5 b
- Helix diameter: 0.5 b（半径 R=0.25 b）
- 隣接ビーズの平衡ばね長（3D距離）: 0.58 b

注意：
- 0.58b は「隣接ビーズの3D距離（ばね自然長）」であり、「軸方向刻み Δx=0.58b」とは限らない。
- よって、ヘリックス点列は「隣接3D距離=0.58b」を満たすように生成する（Δxを解いて決める）。

### 2) べん毛点列生成（必須仕様）
R=0.25b, p=2.5b, L=0.58b として、Δxを
d(Δx) = sqrt(Δx^2 + (2R sin(pi Δx / p))^2) = L
で解き、以下で点列を生成する：
x_n = x0 + nΔx
y_n = R cos(2π x_n/p + φ0)
z_n = ± R sin(2π x_n/p + φ0)

### 3) べん毛に課す拘束（MVP）
- spring（距離）：隣接 (i,i+1)
- bending（角度）：(i-1,i,i+1)
- torsion（ねじれ）：(i-1,i,i+1,i+2)
- 排除体積/反発：近づきすぎ・交差を防ぐ（必要なら seg–seg 距離で評価）

### 4) フック（Hook）要件（MVP）
- フック長（attach–flag0）: 0.25 b
- フックは3点（attach, flag0, flag1）で扱い、角度条件付きで曲げ拘束（交差防止）を入れる。
- フック長および隣接べん毛ボンド長（flag intra）は **multi-step で維持されること**を tests/ で担保する。