# Phase2 / model_alignment_notes

## 目的
この文書は、**物理モデル実装**について、参照論文モデルと現状実装の整合・非整合を明確に管理するためのノートである。

具体的には、各項目について以下を整理する。

- 論文モデルではどう定義されているか
- 現状実装ではどう実装されているか
- 両者の差分は何か
- その差分をどのタスクで修正するか
- 現時点で修正済みか未修正か

この文書では、まず**論文モデルの整理**を行い、その後に**現状実装との差分**および**修正履歴**を追記する。

---

## 対象論文
- Watari, N., Larson, R. G.
- *The Hydrodynamics of a Run-and-Tumble Bacterium Propelled by Polymorphic Helical Flagella*
- Biophysical Journal 98, 12–17 (2010)
- DOI: `10.1016/j.bpj.2009.09.044`

---

## この文書の読み方
各項目は以下の順で記載する。

1. **論文モデル**
2. **現状実装（main 時点）**
3. **差分評価**
4. **対応タスク**
5. **状態**
6. **更新履歴**

状態は以下のいずれかを使う。

- `一致`
- `概ね一致`
- `部分一致`
- `不一致`
- `未実装`
- `要確認`

---

# 1. 論文モデルの整理

## 1.1 全体像
論文モデルは、**cell body + multiple flagella** を bead-spring model として離散化し、低 Reynolds 数下で、各 bead に働く力と Rotne-Prager-Yamakawa (RPY) hydrodynamic interaction によって時間発展させるモデルである。

主要な構成要素は以下である。

- cell body
- multiple flagella
- flexible hook
- motor torque
- spring potential
- bending potential
- torsion potential
- spring-spring repulsion
- hydrodynamic interaction (RPY)

Brownian motion は、論文の基本シミュレーションでは無視される。

---

## 1.2 物理モデルの主要要件（論文）
以下は、論文モデルの実装上重要な要件である。

### A. Body geometry
- body は bead で離散化される
- 図1では 15 beads で表現される
- 形状は長軸 `2b`、短軸 `b` 相当の E. coli body を近似する
- 論文では「shape of a body ... are kept nearly rigid with the employed potentials」とされる
- したがって body は**厳密剛体ではなく、ポテンシャルによりほぼ剛体として保たれる**

### B. Flagellum geometry
- 各 flagellum も bead で離散化される
- 図1では 1 本あたり 15 beads
- normal state の flagellum は左巻き helix
- 代表パラメータ
  - helix pitch: `2.5 b`
  - helix diameter: `0.5 b`
  - helix radius: `0.25 b`
  - adjacent bond equilibrium length: `0.58 b`

### C. Hook
- hook は **1 body bead + 先頭2つの flagellum beads** の 3 点系
- hook bending angle が 90° を超えるときは bending potential をかけない
- 90° 以下でのみ反発的に効かせ、flagellum が body に食い込むのを防ぐ
- hook は flexible であり、flagella の向きは hydrodynamics によって決まる

### D. Motor torque
- motor torque は hook 周辺 3 beads 近傍へ分配される
- flagellum と body への torque balance を保つように与えられる
- motor 自体の詳細機械構造は直接は解かない

### E. Shape maintenance
- 隣接 bead 間距離は spring potential
- bending は 3 bead angle potential
- torsion は 4 bead torsion potential
- spring-spring repulsion で spring crossing を防ぐ
- これらにより body / flagella の shape は**nearly rigid**に保たれる
- ただし、**明示的な剛体再投影やテンプレート再投影を前提としていない**

### F. Hydrodynamics
- low Reynolds number
- bead 間 hydrodynamics は RPY tensor
- 更新式は mobility × total force による overdamped update
- Brownian motion は基本モデルでは無視

### G. Run / tumble / polymorph
- run では全 flagella が同方向回転し bundle を形成
- tumble では少なくとも 1 本の motor を reverse
- reverse された flagellum は normal → semicoiled → curly1 の polymorphic transformation を起こす

### H. Time scaling
- 論文は `b`, `η`, `T` を基本単位にとり
- 時間スケールは `τ = η b^3 / T`
- time step は `Δt = 10^-3`（無次元時間）

---

# 2. 現状実装との差分管理

---

## 2.1 Body geometry

### 論文モデル
- body は bead-spring により離散化される
- shape は potentials により nearly rigid に維持される
- 厳密剛体であるとは書かれていない

### 現状実装（main 時点）
- `builder.py` では body を三角柱状に離散化している
- `body_layer_indices`, `body_ring_edges`, `body_vertical_edges` を持つ
- body は `engine.py` の `_project_body_rigid()` で毎 step rigid projection される

### 差分評価
- 論文は **nearly rigid by potentials**
- 現状実装は **explicit rigid projection**
- よって、body の剛性の入れ方は論文より強い

### 対応タスク
- 必要なら後続タスクで整理
- 直近では差分記録のみ

### 状態
- `部分一致`

### 更新履歴
- 初版作成時点で記録

---

## 2.2 Flagellum geometry

### 論文モデル
- normal helix の代表値
  - pitch `2.5 b`
  - diameter `0.5 b`
  - radius `0.25 b`
  - adjacent bond length `0.58 b`
- 1 本あたり 15 beads

### 現状実装（main 時点）
- `params.py` の既定値は
  - `helix_init.radius_over_b = 0.25`
  - `helix_init.pitch_over_b = 2.5`
  - `bond_L_over_b = 0.58`
- builder 側で helix を離散生成している
- bead 数は flagella length と bond length から導出している

### 差分評価
- 基本幾何は概ね論文に整合
- bead 数の決め方は論文の固定 15 beads と完全には一致しない場合がある

### 対応タスク
- 必要なら離散化仕様の明確化タスクで整理

### 状態
- `概ね一致`

### 更新履歴
- 初版作成時点で記録

---

## 2.3 Attach position / basal link

### 論文モデル
- hook は body 表面の attach 点から flagellum に接続される
- Fig.1, Fig.2 から、基部は body 局所幾何に対して定義される
- body bead - first flag bead の基部リンクは、body 表面から外向きに置かれる構造として解釈される
- hook の自由度は body bead - first flag bead - second flag bead の 3 点系で表される

### 現状実装（main 時点）
- `builder.py` では `body bead - first flag bead` を body rear direction ベースで初期化している
- 現状の基部リンクは、attach body bead が属する断面の
  - `layer centroid -> attach vertex`
  方向では定義されていない
- `engine.py` でも basal link の方向拘束は、局所断面中心→頂点方向ではなく、実質的に body 軸垂直化 / 基部軸ベースで扱われている

### 差分評価
- **論文モデルに対して不一致**
- 現状は、基部リンク方向が body 局所断面幾何に整合していない

### 対応タスク
- `2-16`: basal link を局所断面中心→attach頂点方向の固定長リンクへ修正

### 状態
- `不一致`

### 更新履歴
- 初版作成時点で main の差分を記録
- `2-16` 実装後に更新予定

---

## 2.4 Basal link length

### 論文モデル
- Fig.1 では hook 基部長として `0.25 b` が示される
- 基部リンクは hook の一部として固定的な長さスケールを持つ

### 現状実装（main 時点）
- `builder.py` では hook length `0.25 b` を使っている
- spring rest length も初期配置から導出される

### 差分評価
- 長さスケール自体は概ね整合
- ただし、方向定義がずれているため、長さだけ一致していても論文整合とは言えない

### 対応タスク
- `2-16`: 方向と合わせて固定長リンクとして整理

### 状態
- `部分一致`

### 更新履歴
- 初版作成時点で記録

---

## 2.5 Hook の自由度

### 論文モデル
- hook は 3 点系
- flexible hook
- bending potential は 90° 以下でのみ作用
- flagella の向きは hydrodynamics により決まる

### 現状実装（main 時点）
- `hook_triplets` を持つ
- `threshold_deg = 90.0`
- `compute_hook_forces()` により hook angle 制約をかける
- 一方で、post-step の template projection が flag 全体を戻すため、hook の自由度に対して間接的な強拘束が入っている

### 差分評価
- hook 3 点系・90° しきい値は整合
- しかし template projection により、hook が担うべき自由度を一部上書きしている可能性がある

### 対応タスク
- `2-16`: basal link 補正を first bead のみに限定し、second bead 以降への直接拘束を避ける
- 必要なら後続で hook の自由度評価を追記

### 状態
- `部分一致`

### 更新履歴
- 初版作成時点で記録

---

## 2.6 Helix initial direction

### 論文モデル
- run 状態では flagella は後方へ伸び、bundle を形成する
- 初期条件の詳細手順は簡潔だが、Fig.3 の Normal(run) は body 後方への束化形状を示す

### 現状実装（main 時点）
- `builder.py` では helix 初期接線が rear direction を向くように構成している
- `tests/test_model_builder.py` にも base tangent rear 向きの検証がある

### 差分評価
- 初期後方向きという点では概ね整合
- ただし first bead 自体の置き方が局所半径方向でないため、基部条件込みでは不完全

### 対応タスク
- `2-16`: first bead を局所半径方向に置いた上で、second bead 以降の helix 軸を後方へ向ける

### 状態
- `部分一致`

### 更新履歴
- 初版作成時点で記録

---

## 2.7 Flagellar shape maintenance

### 論文モデル
- flagellum shape は spring / bending / torsion / spring-spring repulsion により nearly rigid に維持される
- 明示的なテンプレート再投影は論文本文の主要実装には現れない

### 現状実装（main 時点）
- spring / bending / torsion / repulsion を持つ
- さらに
  - `chain-length projection`
  - `template projection`
  を post-step 補正として持つ
- `_project_flagella_template()` は flag 全体を基部 frame から再投影する

### 差分評価
- **論文モデルより強い人工拘束**
- 特に template projection は、hydrodynamics や elastic response ではなく、flag shape を直接戻している

### 対応タスク
- `2-16`: basal link 補正と template projection の責務分離を開始
- 後続タスクで template projection の要否・縮小・撤去を整理

### 状態
- `不一致`

### 更新履歴
- 初版作成時点で記録

---

## 2.8 Chain-length projection

### 論文モデル
- 隣接 bead 距離は spring potential により保たれる
- 明示的な距離投影は論文の基本モデルには書かれていない

### 現状実装（main 時点）
- flag chain に対して明示的な chain-length projection を行っている

### 差分評価
- 数値安定化としては有効
- ただし論文モデルの「potential-based nearly rigid」より強い

### 対応タスク
- 差分として継続記録
- 必要なら後続で検討

### 状態
- `不一致`

### 更新履歴
- 初版作成時点で記録

---

## 2.9 Hydrodynamics

### 論文モデル
- bead 間 hydrodynamics は RPY tensor
- overdamped update
- Brownian motion は無視

### 現状実装（main 時点）
- `compute_rpy_mobility()`
- mobility × force による update
- Brownian は default off

### 差分評価
- 基本方針は整合

### 対応タスク
- 直近なし

### 状態
- `一致`

### 更新履歴
- 初版作成時点で記録

---

## 2.10 Run / tumble / polymorph switching

### 論文モデル
- tumble 時に少なくとも 1 本の motor reverse
- semicoiled → curly1 の polymorphic transformation
- run and tumble trajectory を再現する

### 現状実装（main 時点）
- switching 機構自体はある
- ただし main の標準運用では `motor.enable_switching = false`
- 現在は run 固定をデフォルトとしている

### 差分評価
- 機構としては部分的に整合
- main の標準挙動は論文 full model とは一致しない

### 対応タスク
- run 固定は Phase2 の比較のための運用
- full run/tumble 再現タスク時に再整理

### 状態
- `部分一致`

### 更新履歴
- 初版作成時点で記録

---

## 2.11 Time scaling

### 論文モデル
- `τ = η b^3 / T`
- `Δt = 10^-3`（無次元時間）

### 現状実装（main 時点）
- `dt_star = 1e-3` を固定
- `tau_s = 1.0` として内部時間を簡略化
- 物理量から τ を都度再構成する形ではない

### 差分評価
- 数値条件 `dt* = 1e-3` は一致
- 物理スケーリングとしては簡略化あり

### 対応タスク
- 差分として記録
- 必要なら後続で physical τ の扱いを再検討

### 状態
- `部分一致`

### 更新履歴
- 初版作成時点で記録

---

## 2.12 Wall effects

### 論文モデル
- 壁近傍の wall-corrected RPY
- bead-wall short-range repulsion
- 壁近傍での clockwise swimming

### 現状実装（main 時点）
- free-space RPY が基本
- wall effect は main では未実装

### 差分評価
- 論文 full scope に対して未実装

### 対応タスク
- 将来の別タスク

### 状態
- `未実装`

### 更新履歴
- 初版作成時点で記録

---

# 3. 直近の重要修正対象

## 3.1 最優先
- basal link direction を論文整合にする
- `body bead - first flag bead` を局所断面中心→attach頂点方向の固定長リンクにする
- second bead 以降を basal link 補正の直接対象にしない
- template projection と basal link 補正の責務を分離する

## 3.2 次点
- template projection の必要性と影響の評価
- chain-length projection の役割整理
- body rigid projection の論文整合性整理

---

# 4. タスク別更新欄

## Task 2-16
### 目的
- basal link 条件を論文整合へ修正する

### 実装前の状態
- attach 局所半径方向に未整合
- body 軸垂直化 / 基部軸ベースで扱われている
- first bead のみを対象にした局所半径方向固定ではない

### 実装後の記録欄
- [ ] `body bead - first flag bead` が固定長になった
- [ ] `body bead - first flag bead` が attach layer centroid → attach vertex 方向になった
- [ ] 初期配置で上記を満たす
- [ ] step 後も上記を満たす
- [ ] second bead 以降は basal link 補正の直接対象ではない
- [ ] tests 追加済み
- [ ] 本文の該当節を更新済み

### メモ
- 実装後にここへ結果を追記する

---

# 5. 要確認事項

## 5.1 Body 断面幾何の厳密一致
- 論文 Fig.1 の body 幾何と、現状三角柱離散化の厳密対応は要確認

## 5.2 Motor torque distribution の厳密一致
- 現状実装が Fig.2 の torque distribution をどこまで忠実に再現しているかは、必要なら追加確認する

---

# 6. 今後の更新ルール
- 新しい実装差分が見つかったら、まず対応する節に追記する
- 新タスクで修正した場合は
  - `差分評価`
  - `対応タスク`
  - `状態`
  - `更新履歴`
  を更新する
- 可能な限り、「論文モデル」「現状実装」「差分」を混ぜずに書く

---

# 7. 初版作成時の要約
初版時点で、論文モデルとの主要差分は以下である。

1. basal link direction が局所断面中心→頂点方向になっていない
2. template projection により flag 全体を直接再投影している
3. chain-length projection により距離拘束を明示投影している
4. body rigid projection により body が厳密剛体化されている
5. wall effect は未実装
6. run/tumble は機構を持つが、標準運用は run 固定である

このうち、直近で最優先なのは **basal link direction の修正**である。