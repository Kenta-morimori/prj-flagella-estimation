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

## Current status (2026-04-07)
- body static は安定
- minimal basal stub static（`stub_mode=minimal_basal_stub`, `motor.torque_Nm=0`）は安定
- full static（`stub_mode=full_flagella`, `motor.torque_Nm=0`）は finite 完走だが幾何崩壊
- 現時点の最優先残タスクは full flagella 静的安定
- `n_flagella=1`, motor ON 安定化は次段（次PR）の主対象
- 解釈として、body 単独は主因候補から後退し、full chain / torsion / chain-maintenance 側が未解決で残っている

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

### 現状実装（2-16 実装後）
- `builder.py` で各 attach body bead について、所属 layer の
  - `layer centroid -> attach vertex`
  方向を算出し、`first flag bead = attach + 0.25b * radial_unit` として初期化する
- `engine.py` に `_project_basal_link_direction()` を追加し、各 step の body layer 幾何から
  - `current layer centroid -> current attach vertex`
  を再計算して `body bead - first flag bead` を固定長・固定方向へ戻す
- basal link 補正が直接書き換える対象は `first flag bead` のみ

### 差分評価
- basal link の方向定義は、論文解釈に沿って局所断面幾何ベースへ改善した
- 断面中心→attach 頂点方向を step ごとに再評価するため、時間発展中も初期定義と整合する

### 対応タスク
- `2-16`: basal link を局所断面中心→attach頂点方向の固定長リンクへ修正

### 状態
- `概ね一致`

### 更新履歴
- 初版作成時点で main の差分を記録
- 2-16 で局所断面中心→attach頂点方向の固定長 basal link へ修正

---

## 2.4 Basal link length

### 論文モデル
- Fig.1 では hook 基部長として `0.25 b` が示される
- 基部リンクは hook の一部として固定的な長さスケールを持つ

### 現状実装（2-16 実装後）
- `builder.py` では hook length `0.25 b` を使い、局所半径方向に first bead を配置する
- `spring_rest_lengths_m` に基づく basal 長を `engine.py` 側で保持し、step ごとに固定長へ復元する

### 差分評価
- 長さスケールは初期・step 後とも `0.25 b` の固定長として維持される
- 方向定義も同時に局所幾何へ整合したため、基部リンクとしての解釈が改善した

### 対応タスク
- `2-16`: 方向と合わせて固定長リンクとして整理

### 状態
- `概ね一致`

### 更新履歴
- 初版作成時点で記録
- 2-16 で固定長維持を step 後補正に明示化

---

## 2.5 Hook の自由度

### 論文モデル
- hook は 3 点系
- flexible hook
- bending potential は 90° 以下でのみ作用
- flagella の向きは hydrodynamics により決まる

### 現状実装（2-16 実装後）
- `hook_triplets` を持つ
- `threshold_deg = 90.0`
- `compute_hook_forces()` により hook angle 制約をかける
- post-step では
  - `hook distance projection`
  - `basal_link direction fix (first bead のみ)`
  - `chain-length projection`
  - `template projection`
  を順に適用している
- `second bead` 以降は basal link 補正の直接対象ではない（ただし他拘束の影響は受ける）

### 差分評価
- hook 3 点系・90° しきい値は整合
- basal link 補正の直接対象を first bead のみに限定できた
- 一方で template projection が残るため、hook 自由度に対する間接拘束は依然として存在する

### 対応タスク
- `2-16`: basal link 補正を first bead のみに限定し、second bead 以降への直接拘束を避ける
- 必要なら後続で hook の自由度評価を追記

### 状態
- `部分一致`

### 更新履歴
- 初版作成時点で記録
- 2-16 で basal link 補正の直接ターゲットを first bead のみに限定

---

## 2.6 Helix initial direction

### 論文モデル
- run 状態では flagella は後方へ伸び、bundle を形成する
- 初期条件の詳細手順は簡潔だが、Fig.3 の Normal(run) は body 後方への束化形状を示す

### 現状実装（2-16 実装後）
- `builder.py` で first bead を局所半径方向へ置きつつ、`first-second` 接線が rear direction を向くように構成している
- `tests/test_model_builder.py` の base tangent rear 検証は維持されている

### 差分評価
- first bead の局所半径方向配置と helix 後方向き初期条件を両立できた
- 基部条件込みで run 初期形状の解釈が改善した

### 対応タスク
- `2-16`: first bead を局所半径方向に置いた上で、second bead 以降の helix 軸を後方へ向ける

### 状態
- `概ね一致`

### 更新履歴
- 初版作成時点で記録
- 2-16 で first bead 局所半径方向配置と rear 向き接線条件を同時に満たす実装へ更新

---

## 2.7 Flagellar shape maintenance

### 論文モデル
- flagellum shape は spring / bending / torsion / spring-spring repulsion により nearly rigid に維持される
- 明示的なテンプレート再投影は論文本文の主要実装には現れない

### 現状実装（2-16 実装後）
- spring / bending / torsion / repulsion を持つ
- さらに
  - `chain-length projection`
  - `template projection`
  を post-step 補正として持つ
- `_project_basal_link_direction()`（first bead のみ）を追加し、
  basal link 補正と shape maintenance の責務を関数として分離した
- `_project_flagella_template()` は引き続き flag 全体を基部 frame から再投影する

### 差分評価
- **論文モデルより強い人工拘束**
- 特に template projection は、hydrodynamics や elastic response ではなく、flag shape を直接戻している
- ただし 2-16 で basal link 補正責務は first bead の固定長・固定方向に分離された

### 対応タスク
- `2-16`: basal link 補正と template projection の責務分離を開始
- 後続タスクで template projection の要否・縮小・撤去を整理

### 状態
- `不一致`

### 更新履歴
- 初版作成時点で記録
- 2-16 で basal link 補正と template/chain-length 補正を処理順として分離

---

## 2.8 Chain-length projection

### 論文モデル
- 隣接 bead 距離は spring potential により保たれる
- 明示的な距離投影は論文の基本モデルには書かれていない

### 現状実装（2-16 実装後）
- flag chain に対して明示的な chain-length projection を行っている
- post-step では basal link 方向固定の後段に chain-length projection を適用し、
  `second bead` 以降の長さ拘束を担わせている

### 差分評価
- 数値安定化としては有効
- ただし論文モデルの「potential-based nearly rigid」より強い
- 2-16 後は、basal link 補正（first bead）との責務分離が明確になった

### 対応タスク
- 差分として継続記録
- 必要なら後続で検討

### 状態
- `不一致`

### 更新履歴
- 初版作成時点で記録
- 2-16 で chain-length projection の役割を「basal 補正後の鎖長維持」に明確化

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
- [x] `body bead - first flag bead` が固定長になった
- [x] `body bead - first flag bead` が attach layer centroid → attach vertex 方向になった
- [x] 初期配置で上記を満たす
- [x] step 後も上記を満たす
- [x] second bead 以降は basal link 補正の直接対象ではない
- [x] tests 追加済み
- [x] 本文の該当節を更新済み

### メモ
- 実装:
  - `builder.py` で first bead 初期位置を局所半径方向へ修正

## Task 2-29
### 目的
- body rigid projection を常時適用から config 切替へ変更し、projection OFF でも body-only 形状安定性を検証可能にする

### 実装後の記録欄
- [x] `projection.enable_body_rigid_projection` を追加（default=true, 既存互換）
- [x] `enable_body_rigid_projection=false` で `engine.step()` の body rigid projection をスキップ
- [x] `n_flagella=0` の body-only テストを追加（短時間静止 / 長時間静止 / 一様外力並進 / 回転couple）
- [x] `body_constraint_diagnostics.csv` を追加し、force成分を含めて projection OFF の崩れ始めを追跡可能にした
- [x] `body_constraint_local_diagnostics.csv` を追加し、spring pair / bending triplet / layer face の局所異常を追跡可能にした
- [x] 三角柱側面長方形の side-brace を各面両対角にし、layer 間 6 本を topology test で検証可能にした
- [x] body-only の attach_proxy（center layer 3頂点への局所集中荷重）テストを追加した
- [x] 3本べん毛条件でも `body_constraint_local_diagnostics.csv` を出力し、body崩壊原因候補をCSVから説明可能にした

### Body geometry / Body rigidity treatment
- 論文は body を potentials により nearly rigid に保つ想定であり、明示的な rigid projection は主要構成ではない
- 現状実装は rigid projection を残すが、config で OFF 検証を可能化した
- default を ON のままにして既存 run 互換性は維持した

### 現在の論文との差分
- body rigid projection は論文モデル外の補助拘束として残っている
- ただし今回、OFF 検証系を常設したため、body 側を potential 主体で評価できる状態になった
- body 側の不足していた構造拘束は、側面長方形の片対角ではなく両対角が必要だった点にあり、各 layer 間 6 本へ修正した
- attach_proxy は body-only 検証のための proxy 負荷であり、論文の flagellum 力学そのものを置換する実装ではない
- エンジン側では `body_stiffness_scale` により body 関連ポテンシャルの剛性スケールを調整しており、論文の単一スケール想定からの差分として扱う
- 計算安定性とコストの理由から、実装では body–body セグメント間の短距離反発（排除）相互作用をスキップしているが、論文モデルでは body セグメント間の排除効果が存在する前提であるため、この点も論文との差分として明示的に管理する

### 第1ビーズへの含意
- 本タスクの実装対象は body のみであり、第1ビーズの新規拘束は導入していない
- 一方で、body 側が projection OFF で安定することを確認できれば、今後の basal link / 第1ビーズ安定化を切り分けて議論する前提が整う

### 今後の課題
- 次PRで flagellum projection 側の影響を切り分ける
- body 側の安定が前提として確立した後に、flagellum 側のショットノイズ仮説を検証する
- body-only diagnostics の force成分と局所CSVで得た「どの拘束が先に悪化するか」の知見を、次PRの flagellum projection 縮退設計に接続する
- 3本べん毛条件で body が崩れる場合は、`step_summary.csv` の `F_motor_mean_body` / `F_hook_mean_body` と
  `body_constraint_diagnostics.csv`, `body_constraint_local_diagnostics.csv` を併読して原因候補を説明する
  - `engine.py` に `_project_basal_link_direction()` を追加
  - post-step を `hook distance -> basal direction -> chain-length -> template` に整理
- 追加 tests:
  - `tests/test_model_builder.py::test_basal_link_initial_direction_matches_local_layer_radial`
  - `tests/test_simulation.py::test_basal_link_initial_geometry_matches_local_radial`
  - `tests/test_simulation.py::test_basal_link_geometry_after_one_step`
  - `tests/test_simulation.py::test_basal_link_geometry_after_multi_step`
  - `tests/test_simulation.py::test_basal_projection_directly_moves_first_bead_only`

## Task 2-33
### 目的
- hook / flagella projection を分解して個別 ON/OFF 可能にし、projection なし挙動を比較可能にする

### 実装後の記録欄
- [x] `projection.enable_hook_length_projection` を追加（default=true）
- [x] `projection.enable_basal_link_direction_projection` を追加（default=true）
- [x] `projection.enable_flagella_chain_length_projection` を追加（default=true）
- [x] `projection.enable_flagella_template_projection` を追加（default=true）
- [x] `engine.step()` の post-step projection を個別分岐へ変更
- [x] projection OFF 時に同等の hidden correction を追加していない
- [x] `step_summary.csv` に projection 状態列を追加
- [x] `step_summary.csv` に hook / flagella の over_b・rel_err 指標を追加
- [x] hook angle 観測列（mean/min/max/err）を追加
- [x] 1 flagellum, torque=0, projection all OFF で短時間完走・有限値を確認するテストを追加
- [x] motor ON, projection all OFF で比較用 diagnostics が出るテストを追加
- [x] `flagella.init_mode`（`legacy_radius_pitch` / `paper_table1`）を追加
- [x] `flagella.n_beads_per_flagellum` を追加し、初期離散化の source-of-truth を明示化
- [x] `paper_table1` で初期形状を `bond_L_over_b` + `theta0.normal` + `phi0.normal` + bead数 から生成
- [x] `legacy_radius_pitch` を後方互換モードとして維持
- [x] `sim/initial_geometry_summary.json` を追加し、初期幾何を保存
- [x] manifest に `initial_geometry_summary_json` を追記
- [x] 初期形状の target 整合性（bend/torsion）テストを追加
- [x] torsion finite-difference 幅を `potentials.torsion.fd_eps_over_b` として設定可能化
- [x] `step_summary.csv` に `torsion_fd_eps_m`, `torsion_fd_eps_over_b` を追加
- [x] paper mode + n=1 + torque=0 + all projection OFF で fd_eps sweep の比較テストを追加

### Hook / flagella projection の扱い（論文整合）
- hook は角度固定を常時強制するのではなく、まず `attach-first` 距離を spring 主体で扱う方針を維持
- basal link direction は論文より強い拘束になり得るため、独立に ON/OFF 可能にした
- flagella chain/template projection も独立に ON/OFF 可能にし、寄与を分離できるようにした
- 既定値はすべて true とし、既存 run 互換を維持した

### 初期べん毛生成の扱い（論文整合）
- 旧実装では `helix_init.radius/pitch` と `theta0/phi0` が独立で、初期 mismatch を生み得た
- 新実装では `paper_table1` モードを追加し、`bond_L/theta0.normal/phi0.normal/n_beads` を source-of-truth に統一した
- `radius/pitch` は paper モードでは派生量として扱い、初期 summary に記録する運用へ変更した
- `legacy_radius_pitch` は既存互換・比較用として残した

### 現在の論文との差分
- body rigid projection, basal direction projection, chain-length projection, template projection は、いずれも論文の potential 主体モデルより強い拘束になり得る
- ただし本タスクで各拘束を独立に無効化できるようになり、差分を観測ベースで切り分け可能になった
- 初期化については paper モードで source-of-truth を一意化したが、legacy モードは比較互換のため残存している

### Torsion force 検証（Phase A）
- `compute_torsion_forces()` の有限差分幅を固定値から設定値へ変更し、`fd_eps_over_b` を sweep 可能にした
- `fd_eps_over_b = 0.1 / 0.001 / 0.0001` の比較をテストで実施し、step 0 の `F_torsion_mean_flag` が縮小方向になることを確認した
- sweep 条件は `step_summary.csv` の `torsion_fd_eps_over_b` に残るため、後追い解析可能になった

### stiffness scale の扱い結論（Phase A 検証後）
- 本タスクでは Phase A を先行し、まず `fd_eps` 側の離散化誤差寄与を検証した
- 結果として torsion 支配の一因が `fd_eps` 側にあることを確認できたため、`flag_torsion_stiffness_scale` / `flag_bend_stiffness_scale` の default は本PRでは変更しない
- 上記 scale は現時点で repo 独自ノブとして扱い、paper mode の既定値見直しは `fd_eps` 調整後の追加比較タスクで最終判断する

### Torsion force 実装の有限差分幅（Phase A 完了結果）
- **`fd_eps_over_b = 0.1` では step 0 から幾何崩壊を観測**
- **`fd_eps_over_b = 1e-3` および `1e-4` では, `motor.torque_Nm=0` / `n_flagella=1` / all projection OFF 条件で静止平衡をほぼ維持**
- 現在の current candidate は **`fd_eps_over_b = 1e-3`** （1e-4 との改善は顕著でなかったため）
- この改善により，projection の有無そのものや paper mode における初期 shape mismatch が主因ではなく，torsion force 計算側の離散化誤差が支配的であったことが強く示唆される
- 現在設定値として，この `fd_eps_over_b = 1e-3` を適用し，次の Phase B (`motor ON, n_flagella=1`) での検証を実施予定

### Phase B 検証結果（motor ON）
- 条件 `paper_table1` + `n_flagella=1` + `all projection OFF` + `fd_eps_over_b=1e-3` で，
  `motor.torque_Nm = 4e-18 / 1e-18 / 5e-20` を比較した
- いずれのトルクでも基部付近から不安定化し，motor ON 問題は `fd_eps` 調整だけでは解決しなかった
- したがって，次の主因候補は `compute_motor_forces()` の torque-to-force split，または motor 印加時の basal 幾何拘束不足である

### Phase C 実装（motor diagnostics 追加）
- `step_summary.csv` に motor/basal 切り分け用の列を追加した
  - `motor_ra_len_um`, `motor_rb_len_um`
  - `motor_Ta_norm`, `motor_Tb_norm`
  - `motor_Fa_norm`, `motor_Fb_norm`
  - `motor_axis_vs_rear_direction_angle_deg`
  - `motor_attach_force_norm`, `motor_first_force_norm`, `motor_second_force_norm`
- 上記は `compute_motor_forces()` と `StepDiagnostics` を拡張して毎ステップ出力される
- 既存の `motor_degenerate_axis_count`, `motor_split_rank_deficient_count`, `motor_bond_length_clipped_count` と併せて，
  motor force split と basal 幾何のどちらが支配要因かをログから追跡できるようになった

### Phase D 比較結果（torque=5e-20, n_flagella=1）
- 比較条件: `paper_table1`, `fd_eps_over_b=1e-3`
- 4ケース比較:
  - all projection OFF
  - `hook_length` のみ ON
  - `hook_length + basal_link_direction` ON
  - all projection ON
- `step_summary.csv` last step の代表値:
  - all OFF: `hook_len_mean_over_b=37.98`, `flag_bond_len_mean_over_b=92.81`
  - hook only: `hook_len_mean_over_b=0.25`, `flag_bond_len_mean_over_b=82.81`
  - hook+basal: `hook_len_mean_over_b=0.25`, `flag_bond_len_mean_over_b=94.69`
  - all ON: `hook_len_mean_over_b=0.25`, `flag_bond_len_mean_over_b=0.58`
- 解釈:
  - hook長の崩れは `hook_length` 投影で抑えられる
  - しかし basal まで戻しても flagella chain 長の崩れは抑えられない
  - 以前の仮説「`hook_length + basal_link_direction` があれば十分」は，Phase D の比較結果では支持されなかった
  - 少なくとも `torque=5e-20`, `n_flagella=1`, `fd_eps_over_b=1e-3` 条件では，flagella chain / template 側の拘束がないと幾何破綻を防げなかった
  - この結果は，少なくとも `basal` 最小拘束不足だけでは説明が不十分で，
    `compute_motor_forces()` の torque-to-force split 実装（Phase E）を優先検証すべきことを示す

### Phase E の優先検証観点（具体化）
- `Ta + Tb = Ttot` の split が妥当か
- `f = (T × r) / ||r||^2` による force couple が短い basal link で過大になっていないか
- attach / first / second bead への分配が局所伸長を直接作っていないか

### Phase E 実施結果（2026-03-31, torque=5e-20, n=1, fd_eps=1e-3）
- 実行ケース:
  - `all_off`: `outputs/inspect_motor_diag/phaseE/all_off/2026-03-31/205653/sim/step_summary.csv`
  - `hook_basal`: `outputs/inspect_motor_diag/phaseE/hook_basal/2026-03-31/205711/sim/step_summary.csv`
  - `all_on`: `outputs/inspect_motor_diag/phaseE/all_on/2026-03-31/205730/sim/step_summary.csv`
- 観点1 (`Ta + Tb = Ttot`):
  - `motor_split_residual_norm max` は `5.98e-29` / `1.14e-26` / `1.09e-28`（all_off / hook_basal / all_on）
  - split residual は機械誤差レベルで，torque split 自体の整合性は維持される
- 観点2（短リンク過大化）:
  - `motor_attach_force_norm` と `motor_first_force_norm` は hook_basal / all_on でピーク `~2.0e-13`
  - `motor_second_force_norm` は all_on で `~1e-22`，hook_basal でも `~1e-15` 程度
  - force-aware split で second bead 側の過大分配は緩和できたが，基部2点の大荷重は残存
- 観点3（局所伸長）:
  - hook_basal の `flag_bond_rel_err_max last=3.62e+02`（崩壊継続）
  - all_on の `flag_bond_rel_err_max last=5.48e-16`（ほぼゼロ）
  - chain/template 拘束なしでは，split改良後も幾何維持不能
- 結論:
  - Phase E で検証した3点のうち，split 整合性は問題なし，second bead 分配は改善
  - それでも `hook_basal` は破綻するため，主因は split 式の単純不整合ではなく，motor ON 時の形状保持拘束（chain/template）欠如側が支配的

### Phase F 実装（2026-04-01, body-only motor equivalent load）
- `SimulationConfig` に `body_equiv_load` 設定を追加した
  - `enabled`
  - `mode` (`none` / `pure_couple` / `pure_couple_distributed` / `attach_proxy_local` / `distributed_rear_load`)
  - `target_torque_Nm`
  - `target_force_N`
  - `attach_region_id`
- `DynamicsEngine` に body-only 等価荷重生成を実装し，`total_forces` へ加算するようにした
- `step_summary.csv` に以下を追加した
  - `F_body_equiv_load_mean`
  - `F_body_equiv_load_max`
  - `body_equiv_load_mode`
  - `body_equiv_load_target_torque_Nm`
  - `body_equiv_load_target_force_N`
  - `body_equiv_attach_region_id`
- テスト追加/更新:
  - `tests/test_params.py`: `body_equiv_load` の default/override 検証
  - `tests/test_simulation.py`: `step_summary.csv` 新列の存在検証
  - `tests/test_body_stability.py`: body-only で `pure_couple` ON/OFF の診断値差分検証

### Phase F body-only proxy の位置づけと結果整理（2026-04-01 追記）
- `n_flagella=0` の body-only 条件は，論文 full model の再現ではなく，body を主因候補から切り分けるための proxy test として位置づける
- `attach_proxy_local`，`pure_couple`，`pure_couple_distributed` はそのための補助負荷モードであり，本番モデルを置換する目的ではない
- body-only short-time 比較では，baseline と `attach_proxy_local` は body diagnostics 上ほぼ無変形だった
- 旧 `pure_couple` では局所ねじれが発生し，body spring / bend diagnostics が相対的に悪化した
- 作用点・分配方法を見直した `pure_couple_distributed` では，short-time body-only 条件で baseline / attach-proxy と同程度まで改善し，nearly rigid motion に近い応答を示した
- したがって，旧 pure-couple で観測された悪化は body 一般の不安定性よりも，load distribution 起因である可能性が高い

### Phase F の論文整合性メモ（2026-04-01 追記）
- `pure_couple_distributed` は論文 full model そのものではない
- ただし body-only 切り分け proxy としては，旧 `pure_couple` より妥当な負荷近似として扱う
- body に対して望ましい挙動は，局所変形の増幅ではなく，形状をほぼ保った並進 / 回転応答である
- low Reynolds number / overdamped update では，pure couple は合力 0・合トルクのみとなるため，並進が小さく回転中心の応答になること自体は不自然ではない

### Phase F のパラメータ差分メモ（2026-04-01 追記）
- `body_stiffness_scale` は論文の一次パラメータではなく repo 独自の補助ノブである
- 本タスクでは `body_stiffness_scale` の設定化・調整は行っていない
- `dt_star` は Phase2 前提どおり内部 `1e-3` 固定であり，本タスクの調整対象にはしていない
- `fd_eps_over_b=1e-3` は論文固定値ではなく，torsion finite-difference 実装の current candidate として継続運用する

### 注意
- この段階の `body_equiv_load` は切り分け実験用であり，論文モデルの本番拘束ではない
- projection は引き続き診断専用であり，最終解として採用しない方針を維持する

### 今後の課題
- 2x2（hook系 ON/OFF × flagella系 ON/OFF）に加え、body projection ON/OFF を含む比較を継続する
- body-only 側の proxy は `pure_couple_distributed` を採用し，次段は `n_flagella=1`・`motor.torque_Nm=5e-20`・all projection OFF 条件へ戻す
- basal〜first / second bead の local stretch / bend / torsion / repulsion diagnostics を強化し，projection なし安定化に必要な potential / integrator 側不足を特定する
- `flag_bend_stiffness_scale` / `flag_torsion_stiffness_scale` はまだ先に触らない
- `n_flagella=3` 条件では hard gate より diagnostics による原因説明を優先する
- hook angle は現時点で観測列のみ運用し、論文根拠に基づく厳密しきい値は後続で再定義する
- Stage 1/2/3/4 の2x2比較は、`initial_geometry_summary.json` の初期誤差と併読して解釈する

### Phase G 実装（2026-04-01, local basal breakdown diagnostics & minimal_basal_stub mode）
- **minimal_basal_stub モード追加**:
  - `FlagellumParams` に `stub_mode` パラメータ追加（`none` / `minimal_basal_stub` / `full_flagella`）
  - builder で `stub_mode == "minimal_basal_stub"` の場合、n_beads_per_flagellum を 3 に制限
  - attach / first / second の最小基部構成を実装定義化
  - initial_geometry_summary.json に stub_mode を記録
  
- **局所 diagnostics（Phase G-2, G-3） - 既に実装済み**:
  - bond長・相対誤差：`local_attach_first_len_over_b`, `local_attach_first_rel_err`, `local_first_second_len_over_b`, `local_first_second_rel_err`, `local_second_third_len_over_b`, `local_second_third_rel_err`
  - bend/torsion角度・誤差：`local_basal_bend_angle_deg`, `local_basal_bend_err_deg`, `local_first_torsion_angle_deg`, `local_first_torsion_err_deg`
  - 力分解：`local_F_spring_attach_first`, `local_F_spring_first_second`, `local_F_spring_second_third`, `local_F_bend_basal`, `local_F_torsion_first`, `local_F_motor_attach`, `local_F_motor_first`, `local_F_motor_second`, `local_F_repulsion_basal_region`
  
- **Stage 6 3段階比較の構造確立**:
  - body-only（n_flagella=0 + body_equiv_load with pure_couple_distributed）
  - minimal basal stub（n_flagella=1 + stub_mode=minimal_basal_stub）
  - full flagella（n_flagella=1 + stub_mode=full_flagella）
  - すべて all projection OFF で実行可能に

- **テスト検証**:
  - 既存テスト 58 パス（ruff format/check も完全パス）
  - minimal_basal_stub 導入による既存テストへの影響なし

- **Stage 6 実行結果（2026-04-02）**:
  - body-only baseline: `pos_all_finite=True`、`any_nan=False`、`any_inf=False`、`body_equiv_load_mode=pure_couple_distributed`
  - minimal_basal_stub: `stub_mode=minimal_basal_stub` が config に反映され、`flag_intra_count=2`、`hook_count=1`、`local_torsion_first_err_deg=nan`（3 bead stub では torsion 未定義）、`flag_bond_rel_err_max=198.317...`、`local_attach_first_rel_err=351.446...`
  - full_flagella: `stub_mode=full_flagella`、`flag_intra_count=10`、`hook_count=1`、`local_first_torsion_err_deg=56.437...`、`flag_bond_rel_err_max=269.496...`、`local_attach_first_rel_err=150.935...`
  - 3 条件とも 0.01 s の短時間窓では finite を維持したため、first-fail 自体はまだ観測されていない
  - ただし minimal_basal_stub は full_flagella と異なる局所診断を示し、3-bead stub と full chain を切り分ける比較軸として機能することを確認した
  - 現時点では basal/hook 近傍の局所誤差と full chain / torsion 側の追加異常を二層で保持し、first-fail 判定は `local_attach_first_rel_err` を一次指標として整理する

- **Stage 6 追試（0.05 s, all projection OFF, 2026-04-07）**:
  - body-only / minimal_basal / full_flagella の3条件とも `pos_all_finite=True` を維持し、NaN/Inf は発生しなかった
  - minimal_basal では `local_attach_first_rel_err` が step 2 で 100 を超過し、最大 354.917（step 8）
  - full_flagella では `flag_bond_rel_err_max` と `local_first_torsion_err_deg` が step 1 で立ち上がり（それぞれ最大 278.045, 83.144）、chain/torsion 側の追加異常が早期に観測された
  - first-fail の一次判定軸は引き続き `local_attach_first_rel_err`、full 側の補助判定は `flag_bond_rel_err_max` と `local_first_torsion_err_deg` を優先する

### Phase G 実装の位置づけ
- projection 比較から系構成比較へシフト
- minimal_basal_stub か full flagella で先に壊れるかで basal/hook vs helix chain 側を切り分け
- body-only proxy は`pure_couple_distributed` で確定
- 局所 diagnostics により，最初に破綻するモード（stretch/bend/torsion/repulsion/motor）を観測可能化
- 最終的な first-fail 判定は `local_attach_first_rel_err` を軸にし，full flagella では `local_first_torsion_err_deg` と `flag_bond_rel_err_max` を補助して読む

### Phase H 方針転換（2026-04-07, projection 全削除 + Phase1/2 ゲート再定義）
- `projection` 実装を engine/params/conf/step_summary から削除した
  - 削除対象: body rigid / hook length / basal link direction / chain length / template projection
  - `step_summary.csv` から `projection_*_enabled` 列を削除
- 安定化の手段を potential + motor/load handling + integrator に限定した
- `body_equiv_load` は切り分け用途として維持した
- 今回PRのゲートは Phase1/2 までに限定した
  - Phase1: `n_flagella=0`, `motor.torque_Nm=0` で finite 完走
  - Phase2: `n_flagella=1`, `stub_mode=minimal_basal_stub`, `motor.torque_Nm=0` で finite 完走
- Codex 提案の暫定閾値（Phase2 静的）
  - `max(local_attach_first_rel_err) < 0.05`
  - `max(local_first_second_rel_err) < 0.05`
  - 根拠: 0.01 s 試行で実測が 1e-4 オーダーだったため、5% は十分に保守的な静的 gate

---

# 5. 要確認事項

## 5.1 Body 断面幾何の厳密一致
- 論文 Fig.1 の body 幾何と、現状三角柱離散化の厳密対応は要確認

## 5.2 Motor torque distribution の厳密一致
- 現状実装が Fig.2 の torque distribution をどこまで忠実に再現しているかは、必要なら追加確認する

---

# 6. Diagnostics 管理ルール

## 6.1 Diagnostics 列変更時の docs 同時更新ルール

**ルール**: diagnostics CSV / JSON の列を追加・削除・改名した場合は、同一 PR で必ず `docs/phase2/sim_diagnostics.md` を更新すること。

**対象**:
- `src/sim_swim/sim/debug_summary.py` の出力列定義
- `src/sim_swim/sim/core.py` の diagnostics 集計

**更新内容には以下を含めること**:
- 新規列の名前・意味・主要用途フェーズ
- 削除列の理由
- 改名列の旧名・新名・理由

**チェックリスト**:
- [ ] Diagnostics 列の追加/削除/改名がある
- [ ] `docs/phase2/sim_diagnostics.md` が更新されている
- [ ] 新規列説明に「主にどの phase で重要か」が記載されている

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

---

# 8. 2026-04-08 追記（Issue #37 / PR #42 継続）

## 8.1 変更内容（motor split / force-couple）
- `src/sim_swim/dynamics/forces.py` の `compute_motor_forces()` を更新。
- 既存の torque couple 分配（attach-first / first-second）は維持しつつ、
  torque を変えない軸方向 preload を追加。
- preload は link 方向の内力（合力ゼロ・追加トルクほぼゼロ）として与える。
  - attach-first preload scale: `-0.01`
  - first-second preload scale: `+0.01`

## 8.2 目的
- finite 完走だけでなく、Phase B/C の局所誤差を下げる。
- 特に Phase B (`minimal_basal_stub + motor on`) の
  `local_attach_first_rel_err`, `local_first_second_rel_err`, `flag_bond_rel_err_max`
  を改善する。

## 8.3 before/after（canonical torque = 4e-21, 0.1 s）
- Phase B before (baseline2):
  - `local_attach_first_rel_err`: `8746.056%`
  - `local_first_second_rel_err`: `4126.332%`
  - `flag_bond_rel_err_max`: `4126.332%`
- Phase B after (split final):
  - `local_attach_first_rel_err`: `2674.163%`
  - `local_first_second_rel_err`: `1159.254%`
  - `flag_bond_rel_err_max`: `1172.083%`

- Phase C before (baseline2):
  - `local_attach_first_rel_err`: `15740.976%`
  - `local_first_second_rel_err`: `5300.143%`
  - `flag_bond_rel_err_max`: `15612.983%`
- Phase C after (split final):
  - `local_attach_first_rel_err`: `16024.432%`（悪化）
  - `local_first_second_rel_err`: `5214.044%`（改善）
  - `flag_bond_rel_err_max`: `15611.907%`（微改善）

## 8.4 評価
- Phase B は明確な改善を確認。
- Phase C は一部改善（first-second / bond）だが、attach-first が未改善。
- よって Phase C は引き続き残課題。

## 8.5 現時点での未反映差分（2026-04-08 時点の明示）
- **single-flagella 前提の差分**:
  - 現在の canonical 検証は `n_flagella=1` を前提にしており、
    論文の multi-flagella + bundle 形成条件とは未整合。
- **motor split 実装差分**:
  - `compute_motor_forces()` には repo 独自の preload 調整が含まれ、
    論文由来の式をそのまま実装したものではない。
- **torque 設定の位置づけ差分**:
  - canonical torque `4e-21` は比較・安定化の debug 条件であり、
    論文整合の最終 torque 値を確定した状態ではない。
- **run/tumble 運用差分（現フェーズ）**:
  - run/tumble switching は引き続き OFF 運用であり、論文 full model との差分として継続管理する。
- **wall effect 差分（現フェーズ）**:
  - wall effect は今回も未導入で、near-wall 挙動の再現は別タスクで管理する。
- **Phase C 未達の主因候補**:
  - `full_flagella + motor on` の Phase C は未達であり、現時点で論文から最も遠い主因候補は
    `compute_motor_forces()` 近傍（attach-first 側の力分配）にある。

## 8.6 attach-first preload 再設計の追試（2026-04-08, PR #42 継続）

### 8.6.1 変更内容
- `src/sim_swim/dynamics/forces.py` の `compute_motor_forces()` で、
  attach-first 側 preload 係数を `-0.01` から `+0.01` に変更した。
- first-second 側 preload 係数は `+0.01` のまま固定し、
  torque (`4e-21`)・local constraint は変更していない。

### 8.6.2 目的
- Phase C (`full_flagella + motor on`) の `0.1 s` hard gate で、
  主指標 `local_attach_first_rel_err` を優先改善する。

### 8.6.3 主要比較（Phase C, 0.1 s）
- baseline2:
  - `local_attach_first_rel_err`: `15777.011%`
  - `local_first_second_rel_err`: `5187.382%`
  - `flag_bond_rel_err_max`: `15614.998%`
  - `local_first_torsion_err_deg`: `20.580`
- after_split_final（前回候補）:
  - `local_attach_first_rel_err`: `15999.423%`
  - `local_first_second_rel_err`: `5080.634%`
  - `flag_bond_rel_err_max`: `15601.708%`
  - `local_first_torsion_err_deg`: `17.810`
- after_attach_plus001（今回候補）:
  - `local_attach_first_rel_err`: `15712.137%`（改善）
  - `local_first_second_rel_err`: `5256.111%`（悪化）
  - `flag_bond_rel_err_max`: `15607.222%`（微悪化）
  - `local_first_torsion_err_deg`: `15.192`（改善）

### 8.6.4 評価
- Phase C の主指標 `local_attach_first_rel_err` は、
  前回候補 (`15999.423%`) から改善し、baseline2 (`15777.011%`) よりも低下した。
- ただし `local_first_second_rel_err` と `flag_bond_rel_err_max` は前回候補比で悪化し、
  Phase C の総合 DoD は未達のまま。

## 8.7 body axis 依存の重み関数による再調整（2026-04-08, PR #42 継続）

### 8.7.1 変更内容
- `compute_motor_forces()` の attach-first preload を、定数ではなく
  `motor_axis_vs_rear_direction_angle_deg` 相当の局所状態に依存する重み関数へ変更。
- `DynamicsEngine` から body 軸単位ベクトルを渡し、
  attach-first 側 preload は body 軸に対する motor 軸角で切り替える。
- first-second 側 preload は `+0.01` のまま維持し、
  torque (`4e-21`)・local constraint は固定したままにした。

### 8.7.2 目的
- Phase C の `0.1 s` hard gate で attach-first を維持しつつ、
  previous after 比で `local_first_second_rel_err` と `flag_bond_rel_err_max` の悪化を抑える。

### 8.7.3 主要比較（Phase C, 0.1 s）
- baseline2:
  - `local_attach_first_rel_err`: `15777.011%`
  - `local_first_second_rel_err`: `5187.382%`
  - `flag_bond_rel_err_max`: `15614.998%`
  - `local_first_torsion_err_deg`: `20.580`
- after_split_final（前回候補）:
  - `local_attach_first_rel_err`: `15999.423%`
  - `local_first_second_rel_err`: `5080.634%`
  - `flag_bond_rel_err_max`: `15601.708%`
  - `local_first_torsion_err_deg`: `17.810`
- after_axis_weight_final（今回候補）:
  - `local_attach_first_rel_err`: `15890.223%`（改善）
  - `local_first_second_rel_err`: `5041.721%`（改善）
  - `flag_bond_rel_err_max`: `15595.971%`（改善）
  - `local_first_torsion_err_deg`: `18.521`（改善）

### 8.7.4 評価
- `motor_axis_vs_rear_direction_angle_deg` を使う重み関数により、
  previous after 比で 0.1 s の `local_attach_first_rel_err` / `local_first_second_rel_err` / `flag_bond_rel_err_max` がすべて改善した。
- したがって Phase C は CSV 実値ベースでは「条件付き達成候補」に到達した。
- ただし Phase B 側の副作用は別途残るため、比較表は 0.1 s hard gate の改善事実と、短時間全体の副作用を分けて読む。

## 8.8 midpoint = 65 試行の記録（2026-04-08, PR #42 継続）

### 8.8.1 変更内容
- body-axis-angle 重み関数の midpoint を `60.0` から `65.0` に一度だけ変更して試行した。
- 他の条件（high 側係数、first-second preload、torque、local constraint）は固定した。

### 8.8.2 0.1 s 比較（mid60 -> mid65）
- Phase B:
  - `local_attach_first_rel_err`: `2649.919%` -> `1815.522%`
  - `local_first_second_rel_err`: `1097.286%` -> `1140.906%`
  - `flag_bond_rel_err_max`: `1224.663%` -> `1160.827%`
- Phase C:
  - `local_attach_first_rel_err`: `15950.680%` -> `16431.162%`（悪化）
  - `local_first_second_rel_err`: `5163.968%` -> `5122.239%`（改善）
  - `flag_bond_rel_err_max`: `15606.292%` -> `15633.602%`（悪化）

### 8.8.3 評価
- midpoint=65 は Phase B には効くが、Phase C の hard gate を壊すため不採用。
- 現行の midpoint=60 を維持し、Phase C 0.1 s の 3 指標改善を優先する。
- 副作用として Phase B 指標が大きく悪化するため、
  本変更は「Phase C attach-first 改善の方向性確認」段階として扱う。

---

## 8.9 Phase B absolute hard gate 課題整理（2026-04-09, PR #42 Phase B 対策）

### 8.9.1 現状：Phase B (0.1 s) は numerical finite だが shape-preserving finite 未達

**Phase B hard gate チェックリスト（0.1 s）**:

```
✓ pos_all_finite                           = True              (threshold: True)
✓ any_nan                                  = False             (threshold: False)
✓ any_inf                                  = False             (threshold: False)
✗ local_attach_first_rel_err               = 2649.92%          (threshold: 0.5%)
✗ local_first_second_rel_err               = 1097.29%          (threshold: 0.5%)
✗ local_second_third_rel_err               = 1224.66%          (threshold: 0.5%)
✗ flag_bond_rel_err_max                    = 1224.66%          (threshold: 0.5%)
✗ local_basal_bend_err_deg                 = 83.23°            (threshold: 45.0°)
✓ flag_state_changed                       = False             (threshold: False)
✓ motor_degenerate_axis_count              = 0                 (threshold: 0)
✓ motor_split_rank_deficient_count         = 0                 (threshold: 0)
✓ motor_bond_length_clipped_count          = 0                 (threshold: 0)
```

### 8.9.2 局所伸長の実態

- `attach-first` 実測値: **6.875b** （target: 0.58b） → **11.85倍拡張**
- `first-second` 実測値: **6.944b** （target: 0.58b） → **11.97倍拡張**
- `second-third` 実測値: **7.683b** （target: 0.58b） → **13.25倍拡張**

このレベルの過度拡張は、minimal_basal_stub （3 beads only） の構成で motor torque （4e-21 Nm） が加わったときに、
**spring hard limit （s = 0.1b）に達して、spring force がそれ以上の反発力を提供できず**、
単なる preload / weighting の微調整では対応不可の領域。

### 8.9.3 試行と結果

**Trial 1-5（2026-04-09）**: preload / weighting 調整による改善試行
- Trial 1: preload 増加 （-0.02, +0.02） → 悪化（61.7%）
- Trial 2: Both low/high compression → 悪化（約100倍）
- Trial 3: preload 削減 （-0.005, +0.008） → 若干改善（42.5%→26.5%相当）だが未達
- Trial 4: Torque split bias （w_a_aggressive = 1.5） → 悪化（54.3%）
- Trial 5: Direct torque fraction （attach = 0.25） → 悪化（約100倍）

### 8.9.4 物理的解釈

問題の根本：minimial_basal_stub の 3 bead system に対して、motor torque が「大きすぎる」わけではなく、
**spring potential の硬度設定（H_over_T_over_b = 10.0, s = 0.1b）が「too soft」**である可能性

- Spring force: `dU/dd = H * x / (1 - (x/s)^2)^2` で、x が s に近づくと指数増加
- x = 0.1b でもう hard limit に達し、それ以上の伸長を抑止できない
- Motor force の結果として 6.875b 伸長が生じても、spring はもう一定力を超えられない

### 8.9.5 今後の検討方向（このPR での実装予定）

**Option A**: preload / weighting のさらなる調整 （現在進行中）

**Option B**: Spring potential の見直し （conf パラメータ調整）
- `H_over_T_over_b`: 10.0 → より大きい値へ？
- `s`: 0.1b → より大きい値へ？
- ただし、ユーザー指示「local constraint の変更は最小限」に抵触する可能性

**Option C**: Motor torque の見直し
- ユーザー指示「torque をさらに小さくしない」により、現在値 4e-21 Nm は維持方針
- ただし「さらに」の含意の確認が必要

**Option D**: Physics-level redesign
- Motor force split logic の根本的見直し
- OR attach-first preload 以外の、別の stabilization mechanism の導入
- これは「禁止事項」の影響を受ける（7.5b, wall effect, run/tumble 等）のため検討困難

### 8.9.6 まとめ

- Phase B 0.1s は currently **numerical finite でありながら shape-preserving finite 未達**
- 伸長倍率（11-13倍）は preload/weighting のみでは対応困難
- 次ステップ: Option A の追加試行、Option B の necessity evaluation、および report 作成

---

## 8.10 Motor Force Split Logic 根本レビュー（2026-04-09, PR #42 Phase B 本質原因調査）

### 8.10.1 現行 Split ロジック（Weighted Minimum Norm）

**定式化**:
```
Constraint: Ta + Tb = Ttot, Ta·ra = 0, Tb·rb = 0
Objective:  minimize ||Ta||²/w_a + ||Tb||²/w_b
Weight:     w_a = 1/|ra|², w_b = 1/|rb|²
```

**Rationale**: 短いボンド（r_a, r_b）ほど同じトルク下で過大力が生じるため、short-link bias を入れる。

**Problem in minimal_basal_stub**:
- Both |ra| ≈ |rb| ≈ 0.25–0.3 μm （同等の短さ）
- Weight becoming nearly equal despite formula
- Result: torque split nearly uniform → Ta ≈ Tb ≈ Ttot/2
- Consequence: both bonds undergo equivalent overstretching (11–13x)

### 8.10.2 改善案検討と実装結果

**Option 1: Weight Asymmetry Amplification** (Trials 6–7)
- w_a = 5.0/|ra|² (5x penalty) → attach-first 1992% (baseline 1880% より悪化)
- w_a = 10.0/|ra|² (10x penalty) → attach-first 2777% (さらに悪化)
- Conclusion: Weight adjustment だけでは limit （むしろ悪化）

**Option 2: Explicit Torque Fraction** (Trials α=0.25, α=0.1)
- Direct split: Ta_mag = α|Ttot|, Tb_mag = (1-α)|Ttot|
- Maintain perpendicularity: Ta·ra = 0, Tb·rb = 0 by construction
- Result: attach-first 19100% (14x worse than baseline) → FAILED
- Interpretation: explicit fraction 割当は direction を distort する

**Option 3: Revert to Standard Weighted Minimum Norm**
- Back to w_a = 1/|ra|², w_b = 1/|rb|²
- Combined with spring hardness increase (H=20, s=0.2)
- Final result: attach-first 1880% ← BEST among all trials

### 8.10.3 最終状況（Spring hardness increase + Original weighted minimum norm）

**Phase B Hard Gate Check (0.1s)**:

| Metric | Value | Target | Status |
|--------|-------|--------|--------|
| local_attach_first_rel_err | 1880.30% | 0.5% | ✗ (3760x over) |
| local_first_second_rel_err | 1008.11% | 0.5% | ✗ (2016x over) |
| local_second_third_rel_err | 1248.94% | 0.5% | ✗ (2498x over) |
| flag_bond_rel_err_max | 1248.94% | 0.5% | ✗ (2498x over) |
| local_basal_bend_err_deg | 47.76° | 45.0° | ✗ (1.06x over) |
| motor_degenerate_axis_count | 0 | 0 | ✓ |
| motor_split_rank_deficient_count | 0 | 0 | ✓ |
| motor_bond_length_clipped_count | 0 | 0 | ✓ |

**Phase C Hard Gate Check (0.1s)** (reference):

| Metric | Value | Target | Status |
|--------|-------|--------|--------|
| local_attach_first_rel_err | 16164.01% | 1.0% | ✗ (16164x over) |
| local_first_second_rel_err | 5096.94% | 1.0% | ✗ (5097x over) |
| local_second_third_rel_err | 1635.02% | 1.0% | ✗ (1635x over) |
| flag_bond_rel_err_max | 15548.35% | 1.0% | ✗ (15548x over) |
| local_basal_bend_err_deg | 58.96° | 45.0° | ✗ (1.31x over) |

### 8.10.4 根本原因の解釈

**Necessary but Insufficient Condition**:
- Current Weighted Minimum Norm → theoretically sound principle
- Spring hardness increase → necessary for any improvement attempt
- BUT: still 3760x away from Phase B target (absolute improvement marginal)

**Fundamental Limitation**:
minimal_basal_stub (3 beads, ~0.6 μm total) under motor torque (4e-21 Nm) produces:
- Bond extension 1880% = 18.8x stretch
- Spring hard limit (s = 0.2b ~ 0.1 μm in physics) reached → no further force resistance
- Bond length already 6.8–7.7b ≈ 3–4 μm (FULLY OVERSTRETCHED relative to spring model)

This is a **physics-model-level problem**, not a force-distribution fine-tuning problem.

### 8.10.5 次に触るべき点（1つだけ）

**Priority: Motor torque value and/or motor model assumption**

Current value (4e-21 Nm) may need:
1. **Validation against literature**: Is this actually physically reasonable for E. coli motor under canonical conditions?
2. **Alternative: Motor torque distribution in time** (currently constant): gradual ramp-up or impulse structure?
3. **OR: 3-bead stub assumption breach**: Maybe minimal_basal_stub is inherently unstable under actual motor torque loading

**Action NOT to take yet**:
- Further local constraint microadjustments (spring K, preload, weighting coefficients)
- Increased sophistication of force split algorithm
- Changes to basal-bend stiffness or local damping

These are symptom treatments, not root-cause fixes.

### 8.10.6 まとめ

- Motor force split optimized within weighted minimum norm + aggressive spring hardening
- All coefficient variations and split algorithm redesigns yield no significant improvement
- Phase B absolute hard gate remains **3+ orders of magnitude away**
- Indicates problem is not in force-distribution mechanics but in motor model scale or basal structure assumption

---

## 8.11 Split Logic 追加試行（2026-04-09, conservative single-implementation）

### 8.11.1 現行 split の問題点（短要約）

- 現行は `Ta + Tb = Ttot, Ta·ra = 0, Tb·rb = 0` を満たす weighted minimum norm。
- `minimal_basal_stub` では `|ra|` と `|rb|` が同程度に短く、結果的に `Fa` と `Fb` が同時に大きくなりやすい。
- `attach-first` と `first-second` の両方に大きい couple force が入り、局所伸長が同時に悪化する。

### 8.11.2 検討した 2 案（今回は 1 案のみ実装）

1. **Force-space constrained split（実装）**
  - 制約は維持: `r_a×f_a + r_b×f_b = Ttot, f_a·r_a = 0, f_b·r_b = 0`
  - 目的関数: `||f_a||^2 + ||f_b||^2 + λ||f_a-f_b||^2`
  - ねらい: first bead への差分荷重集中を抑える（保守的）
2. **Torque split with explicit physical prior（未実装）**
  - 幾何に応じた split 比を外生で与える方針
  - 拘束との整合性を保つ調整が必要で、今回は見送り

### 8.11.3 実装結果（Phase B 0.1 s, step_summary.csv 実値）

比較対象:
- Before: 現行 weighted minimum norm（現行本線）
- After(trial): force-space constrained split（今回 trial）

| Metric | Before | After(trial) | 判定 |
|---|---:|---:|---|
| local_attach_first_rel_err | 1880.30% | 1926.60% | 悪化 |
| local_first_second_rel_err | 1008.11% | 1038.02% | 悪化 |
| local_second_third_rel_err | 1248.94% | 1245.15% | 微改善（ただし gate 未達） |
| flag_bond_rel_err_max | 1248.94% | 1245.15% | 微改善（ただし gate 未達） |
| local_basal_bend_err_deg | 47.76° | 78.59° | 大幅悪化 |
| motor_degenerate_axis_count | 0 | 0 | 同等 |
| motor_split_rank_deficient_count | 0 | 0 | 同等 |
| motor_bond_length_clipped_count | 0 | 0 | 同等 |

結果:
- numerical finite は維持
- shape-preserving finite は未達のまま
- 特に `local_basal_bend_err_deg` が大幅悪化

このため、trial 実装は **採用しない（revert）**。

### 8.11.4 strict gate に対する未達点

- `local_attach_first_rel_err` 目標 0.5% に対し 1880.30%（3760x over）
- `local_first_second_rel_err` 目標 0.5% に対し 1008.11%（2016x over）
- `local_second_third_rel_err` / `flag_bond_rel_err_max` も 3 桁超過
- `local_basal_bend_err_deg` も閾値 45° を超過

### 8.11.5 次に触るべき点（1つ）

**Motor torque の時系列モデル（定常印加前の ramp 設計）を導入して、初期過渡での局所ピーク荷重を抑えること。**

理由:
- `torque_Nm` の定常値 4e-21 は維持しつつ、禁止事項を回避できる
- preload/midpoint/局所係数の微調整ループに戻らない
- split 係数の静的再配分より、過渡ピークへの直接作用が期待できる

---

## 8.12 Motor Torque 時間構造（linear ramp）試行（2026-04-09）

### 8.12.1 spring hardness の扱い（canonical か一時か）

- `H_over_T_over_b = 20`, `s = 0.2` は、今回以降の Phase B strict-gate 調査における
  **current canonical baseline** とする。
- すなわち、torque 時間構造を検証するときも spring hardness は固定し、
  条件差分は ramp 有無のみとする。

### 8.12.2 ramp モデル仕様（最小実装）

- 対象: motor torque 入力の時間ゲインのみ
- 定常値: `tau_ss = 4e-21`（不変）
- 実効トルク:

  `tau_eff(t) = tau_ss * g(t)`

  `g(t) = min(1, t / t_ramp)` （単調線形, `t_ramp > 0`）

- 今回試行値: `t_ramp = 0.05 s`
- `compute_motor_forces()` / preload / split / local constraint は不変更

### 8.12.3 Phase B 0.1 s 比較（source of truth: step_summary.csv）

| Metric | Before (ramp off) | After (linear ramp 0.05 s) | 判定 |
|---|---:|---:|---|
| local_attach_first_rel_err | 1880.30% | 3757.15% | 悪化 |
| local_first_second_rel_err | 1008.11% | 1250.80% | 悪化 |
| local_second_third_rel_err | 1248.94% | 1858.92% | 悪化 |
| flag_bond_rel_err_max | 1248.94% | 1858.92% | 悪化 |
| local_basal_bend_err_deg | 47.76° | 71.63° | 悪化 |
| motor_degenerate_axis_count | 0 | 0 | 同等 |
| motor_split_rank_deficient_count | 0 | 0 | 同等 |
| motor_bond_length_clipped_count | 0 | 0 | 同等 |

### 8.12.4 strict gate に対する評価

- numerical finite は before/after とも維持
- しかし shape-preserving finite 指標は全般に悪化し、strict gate には遠いまま
- 特に local bond 関連（attach-first/first-second/second-third）が悪化

### 8.12.5 解釈

- この linear ramp（0.05 s）では、初期の低トルク区間で形状が安定化するのではなく、
  後半での追い上げ負荷が局所伸長を増幅した可能性が高い。
- 少なくとも今回の単純 ramp 1種では、Phase B strict-gate 改善は確認できない。

### 8.12.6 次に触るべき点（1つ）

**時間構造は維持しつつ、linear ではなく smoothstep 系 ramp（初期・終端の勾配を0にする）へ切り替えて再評価すること。**

---

## 8.13 Motor Torque 時間構造（smoothstep ramp）試行（2026-04-09）

### 8.13.1 事前方針

- 8.12 の線形 ramp は strict gate 観点で不採用として打ち切る。
- `compute_motor_forces()` / preload / midpoint / local constraint は固定。
- `tau_ss = 4e-21`、`H_over_T_over_b = 20`、`s = 0.2` を固定し、
  時間ゲインのみ smoothstep 1案で評価する。

### 8.13.2 smoothstep ramp 仕様（最小実装）

- 実効トルク: `tau_eff(t) = tau_ss * g(t)`
- ゲイン関数:

  `x = clip(t / t_ramp, 0, 1)`

  `g(t) = 3x^2 - 2x^3`

- 特性: 単調増加、開始・終了勾配が 0
- 試行値: `t_ramp = 0.05 s`

### 8.13.3 Phase B 0.1 s 比較（step_summary.csv 実値）

| Metric | Before (ramp off) | After (linear, 8.12) | After (smoothstep) | 判定 |
|---|---:|---:|---:|---|
| local_attach_first_rel_err | 1880.30% | 3757.15% | 3195.55% | 悪化（linear よりは小） |
| local_first_second_rel_err | 1008.11% | 1250.80% | 1130.47% | 悪化（linear よりは小） |
| local_second_third_rel_err | 1248.94% | 1858.92% | 1674.41% | 悪化（linear よりは小） |
| flag_bond_rel_err_max | 1248.94% | 1858.92% | 1674.41% | 悪化（linear よりは小） |
| local_basal_bend_err_deg | 47.76° | 71.63° | 77.58° | 悪化（linear よりも悪化） |
| motor_degenerate_axis_count | 0 | 0 | 0 | 同等 |
| motor_split_rank_deficient_count | 0 | 0 | 0 | 同等 |
| motor_bond_length_clipped_count | 0 | 0 | 0 | 同等 |

### 8.13.4 判定（継続 or 打ち切り）

- 主要5指標のうち、before を明確改善した項目は 0。
- よって判定ルールに従い、**time-structure 路線（ramp 系）は現時点で一旦打ち切る**。
- strict gate 達成に向けた主戦場としては不十分。

### 8.13.5 次に触るべき点（1つ）

**minimal_basal_stub の幾何モデル妥当性（3-bead stub 前提）を先に検証し、
現在の torque 条件下で shape-preserving finite が理論的に可能かを切り分けること。**

---

## 8.14 motor.torque_Nm = 4e-21 の妥当性確認（2026-04-09）

### 8.14.1 現在の baseline 条件（確定版）

- Phase B canonical: `n_flagella=1`, `stub_mode=minimal_basal_stub`, `duration=0.1 s`
- motor: `torque_Nm=4e-21`, switching off, Brownian off
- spring baseline: `H_over_T_over_b=20`, `s=0.2`
- preload / midpoint / split / local constraint: 現状固定（追加調整しない）
- 判定 source of truth: `step_summary.csv`

### 8.14.2 コード観点の再確認（torque の意味）

`SimulationConfig` と `DynamicsEngine` では、以下のスケーリングを採用している。

- motor torque 入力: `cfg.motor_torque_Nm`
- 力学係数スケール: `torque_for_forces_Nm = |motor_torque_Nm|`
- 代表式:
  - `spring_h = H_over_T_over_b * torque / b`
  - `k_bend = kb_over_T * torque`
  - `k_torsion = kt_over_T * torque`
  - `repulsion_A = A_ss_over_T * torque`

したがって、torque を変更すると motor だけでなく主要剛性側も同時に比例スケールされる。
この実装では、torque 値単独の変更は「絶対力の大きさ」を変える一方、
局所の相対バランスを劇的に変える手段になりにくい。

### 8.14.3 canonical scale との整合（無次元化チェック）

現 baseline (`b=1 um`, `eta=1e-3 Pa·s`) では:

- `eta * b^3 = 1e-21 N·m`
- `motor.torque_Nm = 4e-21 N·m` は `tau* = T/(eta b^3) = 4`
- 単位換算では `4e-21 N·m = 4 pN·nm`

解釈:
- 0037.md に明記の通り、`4e-21` は strict な論文再現値ではなく
  **比較・安定化用の canonical debug 条件**。
- 文献一般で言及される E. coli 系 stall torque はしばしば `1e-18 N·m` オーダー
  （~`10^3 pN·nm`）であり、`4e-21` はそのオーダーより小さい。

### 8.14.4 torque scale vs stub assumption の切り分け

現状観測:
- `4e-21`（既に小さい側）でも Phase B shape 指標は 3桁超過で未達。
- preload / weighting / split / ramp（linear, smoothstep）を触っても未達。

このため、優先仮説は次の順序とする。

1. **先に疑うべき: minimal_basal_stub 仮定（3-bead 近傍幾何）**
   - 現行 torque 条件でも過伸長を抑え切れない。
   - 局所調整の自由度内で改善余地が尽きている。
2. torque scale の再設定は二段目
   - ただし "下げれば通る" は禁止。
   - 文献整合の最終フェーズでは、むしろ大きい torque 側の整合検証が必要。

### 8.14.5 次に実装で触るべき点（1つ）

**minimal_basal_stub の幾何仮定を検証するため、Phase0b/PhaseB 間で
attach-first/first-second の基準長と許容伸長率のモデル整合チェックを追加し、
「現行3-bead stubで strict gate が構造的に到達可能か」を先に判定する。**

---

## 8.15 minimal_basal_stub（3-bead）での strict gate 到達可能性チェック（2026-04-09）

### 8.15.1 比較条件

- source of truth: `outputs/phase_diagnostics_stub_feasibility/*/step_summary.csv`
- 比較対象:
  - Phase0b（motor off, 0.1 s）
  - PhaseB（motor 4e-21, 0.1 s）
- 初期幾何（両者共通）:
  - `stub_mode=minimal_basal_stub`
  - `bead_count=3`
  - `bond_L_over_b=0.58`
  - `initial_attach_first_length_um=0.25`

### 8.15.2 末端 step 実測（0.1 s）

| Metric | Phase0b | PhaseB | PhaseB / Phase0b |
|---|---:|---:|---:|
| local_attach_first_len_over_b | 3.5927 | 4.9508 | 1.38x |
| local_first_second_len_over_b | 4.7418 | 6.4270 | 1.36x |
| local_second_third_len_over_b | 4.5432 | 7.8238 | 1.72x |
| local_attach_first_rel_err | 1337.09% | 1880.30% | 1.41x |
| local_first_second_rel_err | 717.56% | 1008.11% | 1.40x |
| local_second_third_rel_err | 683.30% | 1248.94% | 1.83x |
| local_basal_bend_err_deg | 65.14° | 47.76° | 0.73x |

所見:
- motor off の Phase0b 時点で、local bond 系はすでに strict gate を大幅超過。
- motor on（PhaseB）で悪化倍率は 1.36x〜1.83x だが、崩壊の起点は motor on 前から存在。
- したがって、現状未達は「motor split/ramp のみの問題」ではなく、
  minimal_basal_stub の構造仮定を含む幾何側ボトルネックが主因と判断する。

### 8.15.3 strict gate 項目ごとの到達可能性（現行 3-bead stub 前提）

| 項目 | 現状観測 | 到達可能性判定 |
|---|---|---|
| `pos_all_finite=true` / `any_nan=false` / `any_inf=false` | Phase0b/PhaseB とも達成 | likely reachable（実績あり） |
| motor diagnostics（degenerate/rank_deficient/clipped = 0） | PhaseB で全 0 | likely reachable（実績あり） |
| `local_attach_first_rel_err <= 0.5` | 13.37（Phase0b）〜18.80（PhaseB） | likely structurally impossible（現行前提では困難） |
| `local_first_second_rel_err <= 0.5` | 7.18〜10.08 | likely structurally impossible（現行前提では困難） |
| `local_second_third_rel_err <= 0.5` | 6.83〜12.49 | likely structurally impossible（現行前提では困難） |
| `flag_bond_rel_err_max <= 0.5` | PhaseB 12.49 | likely structurally impossible（現行前提では困難） |
| `local_basal_bend_err_deg <= 45` | 65.14°（Phase0b）, 47.76°（PhaseB） | 境界近傍だが未達（単独改善では不十分） |

結論:
- strict gate は「数値有限性」と「motor 診断」については到達済み。
- 一方、local bond 系 hard gate は 1桁ではなく 1〜2 桁超過しており、
  現行 minimal_basal_stub（3-bead）前提のままでは構造的に到達困難と判定する。

### 8.15.4 次の実装アクション（1つだけ）

**minimal_basal_stub の 3-bead 仮定を置き換える比較実装を 1 本だけ入れる（例: basal 近傍 bead 数を増やす拡張 stub モードを追加）**。

目的は、同一 canonical torque（4e-21）で local bond 系誤差が 1桁以上縮小するかを先に確認し、
strict gate の「構造到達性」を判定すること。

---

## 8.16 比較実装: extended_basal_stub_5（2026-04-09）

### 8.16.1 比較した stub 仮定

- before: `minimal_basal_stub`（bead_count=3）
- after: `extended_basal_stub_5`（bead_count=5）
- 共通条件:
  - `n_flagella=1`
  - `motor.torque_Nm=4e-21`（PhaseB）
  - Phase0b は motor off（`torque_Nm=0`）
  - Brownian off, switching off, `duration_s=0.1`
- source of truth:
  - `outputs/phase_diagnostics_stub_compare_5bead/*/step_summary.csv`

### 8.16.2 Phase0b / PhaseB の before-after 実測（末端 step）

| Metric | Phase0b before (3) | Phase0b after (5) | 改善倍率 (before/after) | PhaseB before (3) | PhaseB after (5) | 改善倍率 (before/after) |
|---|---:|---:|---:|---:|---:|---:|
| local_attach_first_rel_err | 13.3709 | 22.6954 | 0.59x | 18.8030 | 89.9602 | 0.21x |
| local_first_second_rel_err | 7.1756 | 8.7571 | 0.82x | 10.0811 | 31.0942 | 0.32x |
| local_second_third_rel_err | 6.8330 | 8.3120 | 0.82x | 12.4894 | 27.5377 | 0.45x |
| flag_bond_rel_err_max | 7.1756 | 8.7571 | 0.82x | 12.4894 | 31.0942 | 0.40x |
| local_basal_bend_err_deg | 65.1366 | 61.0548 | 1.07x | 47.7553 | 62.9847 | 0.76x |
| motor_degenerate_axis_count | 0 | 0 | 同等 | 0 | 0 | 同等 |
| motor_split_rank_deficient_count | 0 | 0 | 同等 | 0 | 0 | 同等 |
| motor_bond_length_clipped_count | 0 | 0 | 同等 | 0 | 0 | 同等 |

### 8.16.3 1桁以上改善したか

- 判定基準: 改善倍率（before/after）>= 10x を「1桁以上改善」とする。
- 結果: **該当なし（0項目）**。
- 多くの local bond 指標は改善ではなく悪化（改善倍率 < 1）。

### 8.16.4 判定と次の候補（1つだけ）

判定:
- `3 -> 5 bead` の stub 拡張だけでは strict gate 未達の主因を解消できなかった。
- よって、**3-bead 仮定単独を主因とみなすのは不十分**。

次に触るべき点（1つ）:
- **motor torque で剛性・反発も同時スケールしている力学係数設計を分離検証する**
  （`torque_for_forces_Nm` と `motor_torque_Nm` の結合を切り分け、
  canonical torque 固定のまま local bond 指標への寄与を評価）。

---

## 8.17 torque input と剛性/反発スケールの分離検証（2026-04-09）

### 8.17.1 `torque_for_forces_Nm` が支配する係数（実装整理）

`DynamicsEngine.__init__` では `torque = cfg.torque_for_forces_Nm` を基準に、
以下を初期化している。

- spring: `spring_h = H_over_T_over_b * torque / b`
- bend: `k_bend = kb_over_T * torque`
- torsion: `k_torsion = kt_over_T * torque`
- hook: `k_hook = hook.kb_over_T * torque`
- repulsion: `repulsion_A = A_ss_over_T * torque`

一方、motor force には `cfg.motor_torque_Nm` が別経路で使われる。
したがって、従来は `motor_torque_Nm` と `torque_for_forces_Nm` が実質結合していた。

### 8.17.2 分離検証の最小実装案（今回実装）

- 追加: `motor.torque_for_forces_override_Nm`（既定 0.0 = 無効）
- 挙動:
  - `> 0` のとき `torque_for_forces_Nm` をその値で固定
  - `<= 0` のとき既存挙動を維持（後方互換）

これにより、`motor_torque_Nm` は canonical 値のまま、
剛性/反発スケールだけを独立に掃引できる。

### 8.17.3 最小比較実験（PhaseB canonical, 0.1 s）

- 共通条件:
  - `n_flagella=1`
  - `stub_mode=minimal_basal_stub`
  - `motor.torque_Nm=4e-21`（固定）
  - Brownian off, switching off
- before（結合あり）:
  - `motor.torque_for_forces_override_Nm=0.0`（既存結合）
- after（分離）:
  - `motor.torque_for_forces_override_Nm=1e-21`（eta*b^3 相当）
- source of truth:
  - `outputs/phase_diagnostics_torque_force_decouple/*/step_summary.csv`

### 8.17.4 before / after（末端 step）

| Metric | Before (coupled) | After (decoupled) | 改善倍率 (before/after) |
|---|---:|---:|---:|
| local_attach_first_rel_err | 18.8030 | 27.7127 | 0.68x |
| local_first_second_rel_err | 10.0811 | 13.1922 | 0.76x |
| local_second_third_rel_err | 12.4894 | 14.6612 | 0.85x |
| flag_bond_rel_err_max | 12.4894 | 14.6612 | 0.85x |
| local_basal_bend_err_deg | 47.7553° | 69.0063° | 0.69x |
| motor_degenerate_axis_count | 0 | 0 | 同等 |
| motor_split_rank_deficient_count | 0 | 0 | 同等 |
| motor_bond_length_clipped_count | 0 | 0 | 同等 |

### 8.17.5 判定

- strict gate 指標の 1桁改善（>=10x）は 0 項目。
- 今回の分離設定（force scale を 1e-21 へ低下）では、
  local bond 系・basal bend はむしろ悪化。
- よって現時点では、未達の主因を「単純な torque scale」や
  「単純な stiffness/repulsion scale」単独で説明するより、
  **両者の結合設計とその比率設定**が本質課題と判断する。

### 8.17.6 次に触るべき点（1つ）

**`motor_torque_Nm` 固定のまま `torque_for_forces_override_Nm` を上方向（>=4e-21）に1点だけ振る対照実験を行い、
低下側での悪化が「剛性不足」起因かを確認する。**

