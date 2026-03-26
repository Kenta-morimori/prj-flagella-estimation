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
