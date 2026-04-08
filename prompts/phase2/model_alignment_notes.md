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