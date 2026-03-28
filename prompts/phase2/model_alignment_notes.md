# model_alignment_notes

## 2026-03-26 (after 2-15 Phase1-3)

### 論文モデルと同じ点
- 力学の主骨格は spring, bending, torsion, hook, repulsion, motor force, hydrodynamics で構成。
- template projection を OFF にでき、局所ポテンシャル主体の比較実験が可能。
- repulsion は `potentials.spring_spring_repulsion.*` として設定・比較可能。
- stiffness scale は config から指定可能で、default は従来値を維持。

### 論文モデルと異なる点
- post-step projection が存在する。
: body rigid projection
: hook distance projection
: template OFF 時の chain-length projection（任意）
- stiffness scale は論文の直接記述ではなく、実装上のスケーリング係数として導入。

### 差分の理由
- template OFF での発散/凝集を比較するために、数値安定性と比較可能性を確保する必要があるため。
- ただし template 全体再投影は拘束が強いため、OFF 系では chain-length のような弱い拘束で挙動を観察する方針。

### 差分の位置づけ
- body rigid projection: 数値安定化用（暫定仕様）
- hook distance projection: 数値安定化用（暫定仕様）
- template projection: 比較実験用スイッチ（恒久的に保持予定）
- chain-length projection when template off: 比較実験用（暫定仕様）
- configurable stiffness scales: 比較実験用（恒久仕様）

### 設定キー補足
- stiffness の実キーは `stiffness_scales.body`, `stiffness_scales.flag_bend`, `stiffness_scales.flag_torsion`。
- 一部メモで記載の `dynamics.*` とは異なるため、CLI override は `stiffness_scales.*` を使う。

### 現時点の課題
- template OFF + chain-length ON でも即時凝集が残るケースがある。
- repulsion 強化や stiffness 緩和のみでは、現時点の観測では十分な改善に至っていない。
- 基部近傍の shape 保持を弱く補助する local-helix constraint の評価が次フェーズ課題。

## 2026-03-26 (after 2-15 Phase4 diagnostics)

### 論文モデルと同じ点
- 力学更新そのもの（spring, bending, torsion, hook, repulsion, motor, hydrodynamics）は変更していない。
- diagnostics は step 後の状態を読み取るだけで、力計算や積分更新に作用しない。

### 論文モデルと異なる点
- `collapse_diagnostics.csv` と `collapse_summary.csv` を追加し、collapse onset を定量観測する出力経路を導入した。
- collapse 判定を `min_interflagella_distance_um < 0.15` の 3 step 連続で判定する実装ルールを導入した。

### 差分の理由
- OFF 系の即時凝集/発散で「何が先に壊れるか」を比較可能にするため。
- 物理拘束を追加する前に、崩壊の開始指標を時系列で可視化するため。

### 差分の位置づけ
- collapse diagnostics: 比較実験用（数値観測用。力学への副作用なし）

### 現時点の課題
- diagnostics で onset の兆候は見えるが、凝集を解消する機構そのものは未導入。
- 既存 local-helix を含む拘束群との相互作用の整理が必要。

## 2026-03-28 (after 2-15 Phase5)

### Phase1-4 の要約
- template OFF 系の拘束経路を整理し、repulsion/stiffness を CLI 比較可能にした。
- collapse diagnostics を導入し、崩壊開始は距離指標の急減が先行することを確認できる状態にした。

### diagnostics から見えた collapse onset の支配量
- `min_interflagella_distance_um` の急減。
- `min_body_distance_um` の急減。

### Phase5 で導入した basal_link 仕様
- `basal_link` 設定を追加し、`enabled=true` のとき `body bead - first flag bead` の向きを post-step projection で補正。
- 補正は body 側 bead を固定し、first flag bead のみ更新。
- リンク長は hook 接続の spring 平衡長を維持。
- 補正順は `hook 距離拘束 -> basal_link 方向固定 -> chain-length projection`。

### Phase5 で導入した steric_exclusion 仕様
- `steric_exclusion` 設定を追加し、WCA 型 bead-bead 反発を力項として導入。
- 反発は total force の末尾に加算し、post-step projection では適用しない。
- 対象 pair は事前計算し、step ごとに再構築しない。
- 対象は `flagella-body` と `異なる flagella 間`、除外は `same-flagellum` `body-body` `hook-neighbor`。

### 論文モデルとの一致点
- hook は `body bead - first flag bead - second flag bead` の 3 点系として扱う。
- 非貫通性を強化する意図で反発機構を導入する点は論文の物理的意図と整合。

### 論文モデルとの相違点
- `body bead - first flag bead` を固定長かつ body 長軸に垂直方向へ補正する点。
- segment-segment repulsion に加えて bead-bead steric exclusion を追加する点。

### 差分の位置づけ
- basal_link 方向固定: 比較実験用 / 数値安定化用（恒久仕様ではない）。
- bead-bead steric exclusion: 比較実験用 / 数値安定化用（恒久仕様ではない）。
