# Project Plan

## このドキュメントの役割

このドキュメントは，プロジェクト全体の目的，フェーズ構成，現在地，各フェーズの達成条件を整理する親ドキュメントである。

また，Codex CLIがタスク分解・実装・レビューを行う際に参照する，プロジェクト全体の地図として扱う。

ただし，Codex CLIの通常の入口ではない。Phase 2の現在地はまず `docs/phase2/phase2_current.md` を読み，このファイルは全体地図や phase-level context が必要な場合に参照する。

個別タスクの詳細，Codex実行ログ，ADR，詳細な実験履歴は別ファイルで管理する。

主な正本の優先順位は以下である。

* Repository-wide rules: `AGENTS.md`
* Current Phase 2 entry point: `docs/phase2/phase2_current.md`
* Accepted task status: `docs/phase2/phase2_tasks.md`
* Overall project map: `docs/PROJECT_PLAN.md`
* Historical context: Git history, issue/PR history, and Codex run records

---

## プロジェクトの流れ

本プロジェクトの最終目的は，遊泳細菌の実顕微鏡時系列動画像データから，細菌個体ごとのべん毛本数を推定することである。

入力は時系列動画像データであり，出力は細菌個体ごとのべん毛予測本数である。

実験動画では，細菌ごとのべん毛本数のGround Truthを直接得ることが難しい。そのため，べん毛本数を既知ラベルとして設定できる3D物理シミュレーションを用いて，擬似顕微鏡動画像データを生成する。この擬似データで学習したべん毛本数推定モデルを，最終的に実顕微鏡データへ汎化させる。

プロジェクトは以下の流れで進める。

1. Phase 1: セットアップ・実行基盤整備
2. Phase 2: 3D物理シミュレーション・2D擬似顕微鏡動画生成
3. Phase 3: 動画からの菌体検出・個体クリップ生成
4. Phase 4: べん毛本数推定モデルの学習・評価
5. Phase 5+: 推定結果の可視化・実データ解析支援

---

## フェーズ全体像

* [x] Phase 1: セットアップ・実行基盤整備
* [ ] Phase 2: 3D物理シミュレーション・2D擬似顕微鏡動画生成
* [ ] Phase 3: 動画からの菌体検出・個体クリップ生成
* [ ] Phase 4: べん毛本数推定モデルの学習・評価
* [ ] Phase 5+: 推定結果の可視化・実データ解析支援

チェックは，高レベルなフェーズまたはマイルストーン単位の進捗を示す。

細かいタスクの進捗は，後続で作成する `docs/phase*/phase*_tasks.md`，`docs/planning/`，GitHub Issuesで管理する。

---

## チェック更新条件

このドキュメント内のチェックは，以下を満たした場合のみ更新する。

* 対応する実装または検証が完了している。
* 関連テストまたは定量評価が実行されている。
* 必要に応じてユーザーによる目視・定性評価が実施または記録されている。
* `docs/codex-runs/<run-id>/review_result.json` が `"status": "PASS"` である。
* 完了と判断した根拠が，作業ログ・レビュー結果・関連ドキュメントのいずれかに残っている。

---

## Phase 1: セットアップ・実行基盤整備

Status: Completed

Phase 1では，リポジトリ構成，基本CLI，設定ファイル，出力ディレクトリ，ログ・manifestによる再現性管理の土台を整備した。

### 達成済み項目

* [x] 基本的なCLI入口を作成
* [x] 設定ファイルを読み込む実行形式を整備
* [x] `outputs/YYYY-MM-DD/HHMMSS/` 形式の出力方針を整備
* [x] `run.log` と `manifest.json` による実行条件管理を導入
* [x] `scripts/` と `src/` の責務分離方針を整理

---

## Phase 2: 3D物理シミュレーション・2D擬似顕微鏡動画生成

Status: In progress

Phase 2は，本プロジェクトの中核である。

Phase 2の目的は，単に3D遊泳シミュレーションを動かすことではない。最終的な目的は，べん毛本数を既知ラベルとして持つ擬似顕微鏡動画を生成し，Phase 4のべん毛本数推定モデルの教師データとして利用できるようにすることである。

そのためPhase 2では，以下を同時に満たす必要がある。

* 生物物理モデルとしての妥当性
* 数値計算としての安定性
* 2D動画データとしての利用可能性
* 機械学習ラベルとしての整合性
* 実データへの汎化可能性

### Phase 2 progress

* [x] Phase 2.1: トルク0条件での基本的な形状安定性を検証
* [x] Phase 2.2: 論文モデルに基づく初期形状・幾何条件を検証
* [x] Phase 2.3: body only のトルクあり安定回転を検証
* [x] Phase 2.4: body + hook のトルクあり安定回転を検証
* [x] Phase 2.5: body + hook + single flagellum のトルクあり安定回転を検証
* [x] Phase 2.6: single flagellum の bond / bend / torsion 維持策を探索し，螺旋形状維持と持続的な net 回転をmulti-stepで検証（旧定量 gate の誤検出を修正済み、P2-6-007ではprobe成功まで確認。P2-6-008では `material_twist_local_couple` により局所force coupleでnet 1回転以上とshape gate PASSを確認）
* [x] Phase 2.7: 複数べん毛条件で、後方螺旋中心軸の安定整列を束化成功定義として検証（代表条件は `hook_wrapped_axis_aligned`）
* [ ] Phase 2.8: RUN固定状態で、べん毛本数による菌体推進・姿勢安定性・軸整列の違いを運動指標で検証（Issue #65）
* [ ] Phase 2.9: Tumble状態を段階実装する（設計・診断、motor reversal、polymorph切替、run-and-tumble評価; Issue #69）
* [ ] Phase 2.10: 3D結果から2D擬似顕微鏡動画を生成
* [ ] Phase 2.11: ML学習用データとしてラベル・metadata・出力形式を整備

Phase 2.1は，トルク0条件で，論文モデル条件のべん毛数3および拡張条件のべん毛数9において，ビーズ間距離が発散しないことを確認済みである。ただし，菌体形状がやや膨らむように見える点は，今後の確認対象として残す。

Phase 2.2は，`initial_geometry_summary.json` に初期幾何の target / tolerance / pass-fail を記録し，`step_summary.csv` との整合を pytest で確認できる状態にした。詳細は `docs/phase2/phase2_2_geometry_contract.md` に記録する。

Phase 2.3は，body-only equivalent torque (`body_equiv_load.mode=pure_couple`) の body shape gate を `src/sim_swim/sim/body_shape_gate.py` に集約し，`5.0e-20 N m` を safe representative，`1.0e-19 N m` を `body_spring` first-fail representative として固定した。詳細は `docs/phase2/phase2_3_body_only_baseline.md` に記録する。

Phase 2.4は，`minimal_basal_stub` 条件で `local_hook_scale` sweep を再現可能にし，`body_stiffness_scale=50`, `duration_s=0.02` の診断 baseline において `1.2e-21 N m` を pass，`4.0e-21 N m` を `hook` first-fail representative として固定した。詳細は `docs/phase2/phase2_4_hook_gate.md` に記録する。

Phase 2.5は，`stub_mode=full_flagella`, `n_flagella=1`, `duration_s=0.02` の短時間 motor-on 条件で，`1.2e-21 N m` を safe representative，`4.0e-21 N m` を `flag` first-fail representative として固定した。低トルクでは回転 activity と形状 gate 通過が両立し，高トルクでは body/hook より先に flagellum bond / bend / torsion が破綻する，という仮説を `src/sim_swim/sim/single_flagellum_gate.py` と pytest で検証する。Phase 2.5では維持策そのものの実装は行わず，破綻モードを定量化して次タスクの入力にする。詳細は `docs/phase2/phase2_5_single_flagellum_stability.md` に記録する。

Phase 2.6では，Phase 2.5で特定した flagellum bond / bend / torsion 破綻に対して，維持策を仮説ごとに検証する。候補は，内部時間刻み (`dt_s` / `dt_star`) の縮小，flagellum spring / bend / torsion stiffness の局所 scale 調整，motor torque distribution の見直し，必要最小限の数値安定化補助の導入である。各候補は，物理モデル由来の変更か数値安定化由来の変更かを区別し，multi-step gate で safe/fail representatives を固定してから Phase 2.7 以降へ進める。

Phase 2.6では，従来の `flag_phase_rate_hz` と `median(abs(flag_helix_spin_rate_hz))` が目視の螺旋スピンと一致しないことを確認した。`outputs/phase2_6_dt1e4_spin_review/2026-05-31/230124/` では `median_abs_flag_helix_spin_rate_hz=24.40` だったが，`flag_helix_spin_phase_deg` の累積差は 0.25 s で 0.501 deg (`net_abs_flag_helix_spin_revolutions=0.00139`) しかなく，ユーザー目視でもほぼ回転していなかった。そのため Phase 2.6 gate は，瞬間速度だけでなく `net_abs_flag_helix_spin_revolutions` と `flag_helix_spin_direction_consistency` を hard gate に含める。さらに，既存 topology に存在する `attach-first` hook spring が時間積分で未使用だった実装漏れを修正し，診断用に `motor.force_distribution=distributed_flagellum` を追加した。`time.dt_star` はデフォルトを変えず，必要な実行では CLI override で `time.dt_star=1.0e-4` を明示する。`triplet + hook spring fix` では root 方位の net 回転が螺旋全体の net 回転へ十分に伝わらず，現行 bead-position-only model に material frame / segment twist / axial torque flux がないことが主因候補である。P2-6-007では `axial_torque_flux_probe` と `local_twist_transmission_probe` を追加し，root torque が螺旋全体へ届けば単一べん毛で shape gate PASS と net 1回転以上を両立できることを確認した。P2-6-008では `motor.force_distribution=material_twist_local_couple` を追加し，segment orientation stateを隣接bead対の局所force coupleへ変換した。代表条件 `torque=2.0e-20 N m`, `dt_star=1.0e-4`, `duration_s=0.5`, `local_spring_scale=1.2` では，`net_abs_flag_helix_spin_revolutions=1.11698`, `flag_helix_spin_direction_consistency=0.98175`, `flag_torsion_err_max_deg=6.21047` で shape gate を通過した。P2-6-008後は `material_twist_local_couple` を既定の `motor.force_distribution` とし，`triplet` やprobe系modeは明示的な比較・診断用として残す。torsion force OFF + 新手法 ON では shape gate が fail したため，既存 torsion force は螺旋形状維持，新手法は root torque 伝搬として役割分担する。多本べん毛，束化，遊泳軌跡，2D動画自然さはPhase 2.7以降で扱う。詳細は `docs/phase2/phase2_6_helix_retention_gate.md`, `docs/phase2/phase2_6_triplet_twist_dof_design.md`, `docs/adr/0001_phase2_distributed_flagellar_torque.md`, `docs/adr/0002_phase2_torque_transmission_probes.md`, `docs/adr/0004_phase2_material_frame_twist_transmission.md` に記録する。

Phase 2.7では，旧P2-7（複数べん毛の非崩壊検証）と旧P2-8（後方束化判定）を統合し，後方螺旋中心軸の安定整列を束化成功定義として検証した。代表条件 `n_flagella=3`, `motor.torque_Nm=2.5e-20`, `duration_s=0.5`, `time.dt_star=1.0e-4`, `flagella.initial_helix_axis_from_rear_deg=0` は `hook_wrapped_axis_aligned` となり，平均軸からの最大偏差は `15 deg` 閾値内だった。hook length fail は残るが，flag bond / bend / torsion は大きく破綻しないため，次段階ではRUN状態の本数差評価の前提リスクとして扱う。

Phase 2.8では，RUN固定状態でべん毛本数による菌体挙動の違いを評価する。特に `n_flagella=1,3,6,9` で，body displacement，speed，body axis angle change，angular velocity，wobble RMS，flag helix axis alignment，bundle/body axis angle を比較し，本数差が推進量や姿勢安定性へ与える影響を整理する。このタスクは Issue #65 に対応する。

Phase 2.9では，Tumble状態を段階実装する。まず現行run固定モデルと論文tumble要素の差分を整理し，次に motor reversal，normal / semicoiled / curly1 のpolymorph切替，最後にrun-and-tumble挙動評価へ進む。このタスクは Issue #69 に対応し，各stageで個別にreview_resultを残す。

### Phase 2運用フェーズ（A/B/C/D）と現状

旧 prompt 由来の運用フェーズは，本ドキュメントへ統合して管理する。

現時点（2026-05-30）のコード根拠ベースの整理は以下。

* PhaseA（torque=0, CSV判定）: 実装済み。`shape_pass_nonbody` / `first_fail_category_nonbody` を含む判定が `step_summary.csv` に出力される。
* PhaseB1（first-fail分類）: 実装済み。`run_motor_scale_sweep.py` が `shape_pass` と `first_fail_category` を集計してCSV化する。
* PhaseB2（hard gate）: 実装済みだが，side-attach初期条件の反映により旧 break torque 条件は再同定が必要。
* PhaseC（side-attach, dt_star運用）: 実装済み。`time.dt_star` override で内部刻みを制御できる。
* PhaseD（sweep再開）: D2/D3は実装済み。旧D1（基準トルク再同定）は採択済みタスクへ移管し，P2.3 body-only，P2.4 minimal hook，P2.5 single full flagellum の条件別 representatives として固定した。

### PhaseD相当チェック結果（2026-05-30）

旧 PhaseD 要件に対する確認結果:

1. D1: B2を満たす基準トルク再同定  
旧D1の単一基準トルク運用は採用せず，条件別 representatives として再定義した。P2.3では body-only equivalent torque，P2.4では `minimal_basal_stub` hook gate，P2.5では `full_flagella` single flagellum gate を固定している。P2.6以降では，この短時間基準を multi-step helix retention や複数べん毛条件へ拡張する。

2. D2: 基準トルク固定で scale sweep 再開  
実装済み。`scripts/01_simulate_swimming/run_motor_scale_sweep.py` で `local_hook_scale` / `local_spring_scale` / `local_bend_scale` / `local_torsion_scale` を sweep 可能。

3. D3: pass/fail と first_fail_category の可視化  
実装済み。`scripts/01_simulate_swimming/plot_motor_scale_collapse_heatmap.py` が category heatmap と component pass/fail heatmap（body/hook/flagella）を出力する。

### 参照論文

Phase 2では，以下の論文モデルを主な参照元とする。

* Nobuhiko Watari and Ronald G. Larson, "The Hydrodynamics of a Run-and-Tumble Bacterium Propelled by Polymorphic Helical Flagella", Biophysical Journal, 2010.
* DOI: 10.1016/j.bpj.2009.09.044

この論文では，bead-spring modelにより，菌体・複数べん毛・hydrodynamic/mechanical interaction・モーター回転反転・べん毛の多形変換を扱う。

### 参照論文モデルとの関係

Phase 2は，低レイノルズ数下の細菌遊泳モデルを参照して進める。

基本方針として，本実装は参照論文モデルに基づく。ただし，本プロジェクトは機械学習用データ生成を目的とするため，参照論文モデルをそのまま再現するだけではなく，べん毛本数ラベル付きの2D擬似顕微鏡動画を生成できるように拡張している。

Phase 2では，以下を常に区別する。

* 参照論文モデルに基づく部分
* 本プロジェクト目的のために拡張した部分
* 数値安定化のために導入した近似・制約
* MVPとして一時的にOFFまたは固定している部分
* 将来的に反映すべき部分

重要な差異はADRだけに任せず，冗長にならない範囲で本ドキュメントにも記録する。ADRは，設計判断の背景・理由・影響を詳しく残すために用いる。

### Phase 2モデル要素の整理

| 項目             | 参照論文モデルに基づく部分                        | 本実装・MVPでの扱い                                                                                                      | 注意点                                          |
| -------------- | ------------------------------------ | ---------------------------------------------------------------------------------------------------------------- | -------------------------------------------- |
| 菌体形状           | 菌体をビーズで離散化し，多角柱状の構造として扱う             | MVPでは `body.n_prism = 3` の三角柱として扱う                                                                               | `body.n_prism != 3` はMVP範囲外                  |
| 菌体内部結合         | 隣接ビーズ間をspringで接続する                   | ring方向・長軸方向のspringに加え，層間の斜めbrace springを追加する                                                                     | 論文モデルとの差異。菌体のshear deformationを抑える数値安定化寄りの実装 |
| べん毛本数          | 論文条件では主に3本条件を基準に扱う                   | 本プロジェクトでは教師データ生成のため `0 <= flagella.n_flagella <= 9` を扱う                                                          | ML用ラベル生成のための拡張                               |
| べん毛付着位置        | 菌体表面の点からべん毛を生やす                      | 3本以下では中心層の3点から選ぶ。4本以上では中心層と前後隣接層を含む最大9候補から，seedに基づき重複なしで選ぶ                                                       | 「中心3点 + 隣接する6点」から最大9本を生やす拡張                  |
| べん毛幾何          | normal stateの螺旋形状を使う                 | helix pitch `2.5 b`，helix radius `0.25 b`，隣接距離 `0.58 b` を基本とする                                                   | 初期形状と時間発展後の両方で検証する                           |
| べん毛角度条件        | bending / torsionの平衡角を持つ             | normal stateとして bending angle `142°`，torsion angle `-60°` を扱う                                                    | semicoiled / curly1 は将来的な拡張対象                |
| フック            | 菌体とべん毛の接続部として扱う                      | hook length `0.25 b` を基本とし，attach point, flag0, flag1 の関係で曲げ拘束を扱う                                                | hook lengthのmulti-step維持が重要                  |
| 多形遷移           | normal, semicoiled, curly1などの多形状態を扱う | MVPでは `motor.enable_switching = false` によりnormal stateを固定する                                                      | tumbleやrun-and-tumbleを扱う段階で拡張候補              |
| 壁効果            | 実環境では壁効果が影響し得る                       | 現時点ではwall effectはOFFまたは未反映として扱う                                                                                  | 一時的な簡略化。将来的には検討対象                            |
| Brownian noise | 熱揺らぎを扱い得る                            | 現時点では決定論的条件を基本とし，Brownian noiseは必要に応じて扱う                                                                         | 物理妥当性と数値安定性の切り分けが必要                          |
| 数値安定化制約        | 物理モデルそのものとは区別する                      | 必要に応じてchain-length projection, template projection, repulsion, stiffness scaling, shape constraintsなどを導入する可能性がある | 物理モデル由来か，数値安定化由来かを必ず区別する                     |
| 2D投影           | 論文モデル本体ではなく，本プロジェクト固有の出力処理           | 3Dシミュレーション結果を2D擬似顕微鏡動画へ投影する                                                                                      | ML用データ生成のため必須                                |

### 参照論文モデルに基づく部分

現時点で重要な論文準拠条件は以下である。

#### 菌体

* 菌体はビーズで離散化する。
* MVPでは三角柱状のbodyとして扱う。
* `body.n_prism = 3` をMVP条件とする。
* 菌体の長軸方向と断面方向のspringにより，基本形状を維持する。
* 本実装では，菌体のshear deformationを抑えるため，層間の斜め方向に補助brace springを追加している。この点は参照論文モデルとの差異として扱う。

#### べん毛幾何

* helix pitch: `2.5 b`
* helix diameter: `0.5 b`
* helix radius: `0.25 b`
* adjacent bead distance: `0.58 b`

これらは，初期形状生成時だけでなく，時間発展後にも維持されているかを確認する。

#### べん毛角度条件

参照論文モデルでは，べん毛はnormal，semicoiled，curly1などの多形状態を取り得る。

60-bead coarse-grained modelでの平衡角は以下である。

| State      | Equilibrium bending angle | Equilibrium torsional angle | 備考                 |
| ---------- | ------------------------: | --------------------------: | ------------------ |
| normal     |                    `142°` |                      `-60°` | run状態で使う基本形        |
| semicoiled |                     `90°` |                       `65°` | tumble開始時に切り替わる候補  |
| curly1     |                    `105°` |                      `120°` | tumble中盤以降に切り替わる候補 |

MVPでは，`motor.enable_switching = false` とし，normal stateを固定する。

したがって，semicoiled / curly1 は現時点では実装・検証の必須対象ではない。ただし，tumble，motor reversal，run-and-tumble挙動を扱う段階では，これらの多形状態を導入する候補とする。

#### Polymorphic switching

Polymorphic switchingとは，べん毛フィラメントがnormal，semicoiled，curly1などの多形状態を切り替えることである。

参照論文モデルでは，tumble時にモーター回転が反転すると，対象べん毛の平衡曲げ角・平衡ねじれ角をnormalからsemicoiled，さらにcurly1へ切り替えることで多形変換を表現する。

現行MVPでは，遊泳run状態の再現と2D擬似顕微鏡動画生成を優先するため，polymorphic switchingはOFFとし，normal stateを固定する。

これは最終的な物理モデルとして除外するという意味ではなく，将来的にtumbleやrun-and-tumble挙動を扱う段階で反映する候補である。

#### フック条件

* hook length: `0.25 b`
* attach point, flag0, flag1 の関係に基づき，曲げ拘束や接続条件を扱う。
* hook lengthはmulti-step testで維持を確認する。

#### 低レイノルズ数条件

細菌遊泳は低レイノルズ数条件下の運動として扱う。

慣性支配ではなく，粘性支配の運動として考える。

### 本実装における拡張・差異

#### べん毛本数の拡張

参照論文モデルとは異なり，本プロジェクトでは機械学習用データ生成を目的とするため，べん毛本数を複数条件で扱う。

MVP制約：

* `body.n_prism = 3`
* `0 <= flagella.n_flagella <= 9`

べん毛付着位置は以下の方針で扱う。

* `n_flagella <= 3` の場合，中心層の3点から付着点を選ぶ。
* `n_flagella > 3` の場合，中心層とその前後隣接層を含む最大9点を候補とし，seedに基づいて重複なしで付着点を選ぶ。

この拡張は，べん毛本数推定モデルの教師データを生成するために必要である。

#### 菌体補助spring

本実装では，菌体の隣接層間に斜め方向の補助brace springを追加している。

目的は，菌体のshear deformationを抑え，形状崩壊を避けることである。

これは参照論文モデルそのものではなく，現実装における数値安定化・形状維持のための拡張として扱う。Phase 2でbody deformationを評価する際には，この補助springの影響も考慮する。

#### 2D擬似顕微鏡動画の生成

参照論文モデルではなく，本プロジェクト固有の目的として，3D結果を2Dへ投影し，擬似顕微鏡動画を生成する。

このため，以下が重要になる。

* image size
* pixel scale
* fps
* body centering
* projection metadata
* flagella count label
* output video or frame sequence

2D投影はPhase 3以降の入力データとなるため，3Dシミュレーション側で座標系・長さスケール・scalingを変更した場合は，必ず2D projectionにも反映する。

#### MVPとして一時的に固定またはOFFにするもの

現時点では，数値安定性や検証の単純化のために，一部の物理要素を固定またはOFFにする場合がある。

例：

* run状態の固定
* polymorphic switchingのOFF
* wall effectのOFF
* 一部の流体相互作用の簡略化
* Brownian noiseのOFF
* 数値安定化制約の導入

これらは最終的な物理モデルとして確定したものではなく，将来的に反映・拡張すべき候補として扱う。

#### 数値安定化のための制約

必要に応じて，以下のような数値安定化制約を導入する可能性がある。

* chain-length projection
* template projection
* repulsion
* stiffness scaling
* hook/flagella shape constraints
* time-step control
* torque scaling

これらを導入する場合は，物理モデル由来の条件なのか，数値安定化のための実装上の近似なのかを明示する。

また，scalingを導入または変更した場合は，力学計算，診断値，3D出力，2D projection，metadataへ一貫して反映する。

重要な設計判断であればADRに記録する。ただし，PROJECT_PLANにも冗長にならない範囲で重要な差異を残す。

### Phase 2の成功条件

Phase 2の最終的な成功条件は以下である。

1. 菌体やべん毛が形状崩壊しない。
2. べん毛が束化する。
3. べん毛の回転による菌体の遊泳が再現される。
4. 3D結果から2D擬似顕微鏡動画を生成できる。
5. べん毛本数ラベルと動画が対応している。
6. 下流の学習データとして利用できる。
7. 出力が再現可能である。

ただし，これらは最終目標であり，フェーズやサブフェーズごとに評価可能な定量指標へ落とし込む必要がある。

### 定量評価と定性評価

Phase 2では，pytestや診断量による定量評価と，ユーザーによる目視・定性評価を併用する。

定量評価では，以下のような値を確認する。

* bead distance
* hook length
* helix pitch
* helix radius
* flagellar bond length
* bending angle
* torsion angle
* body deformation
* inter-flagella distance
* body center displacement
* angular velocity
* output file existence
* metadata consistency

一方で，以下の項目は定量指標だけで判断しにくいため，ユーザーによる目視評価も重要である。

* 論文Fig.3 Normalのような後方束化に見えるか
* 側面束化やcollapseではないか
* 菌体やべん毛が不自然に変形していないか
* べん毛回転による遊泳らしさがあるか
* 2D擬似顕微鏡動画として自然か
* ML学習用データとして使えそうか

Codexがユーザー目視評価を必要と判断した場合は，以下を明示する。

* 実行コマンド
* 出力ディレクトリ
* 確認すべき動画・画像・ログ
* 評価すべき観点
* 自動評価で通った点
* 自動評価では判断できない点

ユーザー目視評価が必要で，まだ確認されていない場合，そのタスクは完了扱いにしない。

### Phase 2の既知課題

#### トルク0条件での安定性

トルクなしの静止条件で，ビーズ間距離が発散しないことを確認する。

べん毛数3および拡張条件のべん毛数9でも距離が破綻しないことが重要である。

一方で，菌体形状がやや膨らむように見える場合があり，剛体性・拘束条件・可視化方法のどこに由来するかを区別する必要がある。

#### トルクあり条件での安定回転

現在の主要課題である。

検証は段階的に行う。

1. body only
2. body + hook
3. body + hook + flagella
4. 複数べん毛
5. 束化・遊泳挙動

#### べん毛螺旋の崩壊

トルク増大時に，べん毛が螺旋を維持せず，団子状に潰れる失敗モードがある。

この問題に対して，初期形状が平衡角を満たすこと，bending/torsion条件が維持されること，multi-step testで形状を検証することが重要である。

#### 後方束化

論文Fig.3のNormal条件のように，菌体後方で複数べん毛が束化し，推進する挙動を再現することが目標である。

ただし，現状では側面での束化やcollapse/fly-awayなどの失敗モードが起こり得るため，初期条件，流体相互作用，制約条件，反発，べん毛剛性，投影補正を切り分ける必要がある。

### 失敗モードの扱い

Phase 2では，失敗モードも重要な開発成果として扱う。

失敗モードの例：

* collapse
* fly-away
* hook length divergence
* body deformation
* helix collapse
* side bundling
* posterior bundling failure
* unstable torque response
* unrealistic swimming trajectory

これらは単なる失敗ではなく，モデル差異・数値安定性・制約条件を検証するための診断対象である。

失敗モードを再現できる条件，設定，出力，診断量は保存する。

有用な失敗モードは，review FAILであってもdiagnostic commitとして残してよい。

ただし，失敗モードを残すことと，タスク完了は区別する。

### 出力・ログ・再現性

実行結果は，原則として以下に保存する（JST）。

```text
outputs/YYYY-MM-DD/HHMMSS/
```

各runでは，可能な限り以下を残す。

* `run.log`
* `manifest.json`
* config path
* config overrides
* seed
* input paths
* output paths
* git commit hash
* python version
* environment details

Phase 2では，診断ログとして以下を重視する。

* `step_summary.csv`

`step_summary.csv` は，1 step = 1 row の時系列ログとし，解析に必要な主要診断量を集約する。

`step_summary_full.csv` は再導入しない。ただし，明示的な方針変更がある場合はADRまたはタスクログに記録する。

### Phase 2実行インターフェース（現行コード準拠）

Phase 2実行・sweep・可視化の正本は `scripts/01_simulate_swimming/` の実装とテストとする。

主要コマンド:

* 単一実行: `uv run python -m scripts.01_simulate_swimming [--duration-s ...] [KEY=VALUE ...]`
* sweep: `uv run python -m scripts.01_simulate_swimming.run_motor_scale_sweep --target ... --values ... [--torques ...] [--duration ...]`
* 可視化: `uv run python -m scripts.01_simulate_swimming.plot_motor_scale_collapse_heatmap --summary-csv ... [--output-dir ...]`

sweep集計の主要列:

* `finite_pass`
* `shape_pass_nonbody`
* `body_shape_pass`
* `shape_pass`
* `first_fail_category`
* `motor_degenerate_axis_count`
* `motor_split_rank_deficient_count`
* `motor_bond_length_clipped_count`

run出力で必須とする診断ファイル:

* `sim/step_summary.csv`
* `sim/body_constraint_diagnostics.csv`
* `sim/body_constraint_local_diagnostics.csv`

### テスト方針

手動確認や目視確認だけに依存しない。

pytestによる自動検証を重視する。

基本方針：

* 重要な物理制約はpytestで検証する。
* one-step testだけで完了しない。
* stability-sensitiveな項目はmulti-step testで確認する。
* CLI subprocess testだけに偏らない。
* `src/` のライブラリ層を直接呼ぶテストを優先する。
* 物理量，幾何量，出力形式を分けて検証する。

検証対象の例：

* bead distance
* hook length
* flagellar bond length
* helix pitch
* helix radius
* bending angle
* torsion angle
* body deformation
* flagellar state switching
* body center position in 2D projection
* output file existence
* manifest/log generation

テストやシミュレーションが環境要因で実行できない場合でも，作業内容・未実行理由・代替確認を記録する。

---

## Phase 3: 動画からの菌体検出・個体クリップ生成

Status: Not started

Phase 3では，実験動画または擬似顕微鏡動画から菌体を検出し，個体ごとのクリップ画像列を生成する。

想定手法：

* YOLOによる菌体検出
* トラッキング
* 中心化
* スケール統一
* bbox / track情報の保存

主な対象：

* `scripts/02_detect_bac.py`
* `src/flagella_estimation/`

達成条件：

* 動画から菌体検出結果を出力できる。
* 個体ごとのクリップを生成できる。
* クリップはフレーム間で位置・スケールが揃っている。
* bbox，track ID，frame情報が保存されている。
* Phase 4の学習入力として利用できる。

---

## Phase 4: べん毛本数推定モデルの学習・評価

Status: Not started

Phase 4では，擬似顕微鏡動画または個体クリップを用いて，べん毛本数推定モデルを学習する。

CNNまたはVision Transformerをベースラインとして比較し，最終的には実験動画への汎化を評価する。

主な対象：

* `scripts/03_train_evaluate.py`
* `src/flagella_estimation/`

達成条件：

* 学習・評価パイプラインが動作する。
* 入力形式とラベル形式が定義されている。
* べん毛本数の推定結果を出力できる。
* 評価指標が定義されている。
* ベースライン結果が再現可能である。
* 実顕微鏡動画への適用可能性を評価できる。

---

## Phase 5+: 推定結果の可視化・実データ解析支援

Status: Not started

Phase 5以降では，推定されたべん毛本数を元動画または解析動画へ重畳し，推定結果と遊泳挙動の関係を確認できるようにする。

主な対象：

* `scripts/10_overlay_flagella.py`
* visualization utilities

達成条件：

* 元動画に推定ラベルを重畳できる。
* 推定結果とフレーム・個体IDが対応している。
* 可視化動画や図を再現可能に出力できる。

---

## ディレクトリ方針

### `scripts/`

ユーザーまたは開発者が直接実行するCLI入口を置く。

役割：

* config読み込み
* key=value overrideの処理
* 出力ディレクトリ作成
* logging/manifest初期化
* `src/` のpipeline呼び出し

原則：

* アルゴリズム本体を `scripts/` に増やさない。
* 複雑な処理は `src/` に移す。
* `scripts/` はプロジェクト本体の実行入口であり，Codex運用スクリプトは置かない。

### `src/`

再利用可能な実装本体を置く。

主な責務：

* 物理モデル
* 数値計算
* 投影処理
* データ整形
* モデル定義
* 評価処理

想定パッケージ：

* `src/sim_swim/`: Phase 2
* `src/flagella_estimation/`: Phase 3-4

### `conf/`

実行設定ファイルを置く。

例：

* `conf/sim_swim.yaml`

### `outputs/`

実行成果物を保存する。

原則：

* `outputs/YYYY-MM-DD/HHMMSS/` 形式で保存する。
* 各runに `run.log` と `manifest.json` を残す。
* 出力形式・列名・単位・座標系はフェーズ間で互換性を保つ。

### `tools/`

開発支援ツールを置く。

Codex CLIの自動化スクリプトは，プロジェクト本体の `scripts/` ではなく，以下に配置する。

* `tools/codex/`

---

## 禁止事項

以下は禁止する。

* `main` / `master` への直接作業
* `scripts/` へのアルゴリズム本体の肥大化
* 理由のない大規模リファクタリング
* 理由のない依存関係追加
* `step_summary_full.csv` の再導入
* 出力形式の無断変更
* Ground Truthラベルと出力動画の対応関係を曖昧にすること
* 論文モデルとの差異を記録せずにPhase 2の重要変更を行うこと
* 目視確認のみでPhase 2完了とすること
* review PASSなしでタスク完了扱いにすること
* 秘密情報，API key，token，個人情報のcommit
* 3D座標系，長さスケール，pixel scale，scalingを変更したにもかかわらず，2D projection，metadata，診断出力へ反映しないこと
* 2D projectionだけを変更し，対応する3D座標系・scale・metadataとの整合性を確認しないこと

---

## Codex CLI運用との関係

Codex CLIは，本プロジェクトにおける開発支援エージェントとして使用する。

Codexは以下を行う。

* `AGENTS.md` と本ファイルを読む。
* 必要に応じてタスクを分解する。
* タスク提案を `docs/planning/` に記録する。
* ブランチを作成する。
* 実装する。
* テストを実行する。
* レビューを実施する。
* review PASS後にcommit/pushする。
* 重要な設計判断がある場合のみADRを作成する。
* 作業ログを残す。
* 進捗管理ファイルが導入された後は，review PASS後にチェックを更新する。

Codexは以下を行わない。

* review PASSなしで完了扱いにする。
* main/masterへ直接作業する。
* GitHub issueを自動作成する。
* pull requestを自動作成する。
* mergeする。
* 秘密情報を扱う。
* 論文モデルとの差異を曖昧にしたままPhase 2の重要変更を進める。

---

## Codex用メモ

Codexのトークン使用量を抑えるため，補助的な要約メモを `docs/codex-notes/` に作成してよい。

Codex用メモは正本ではなく，索引・要約・作業補助のための資料である。

以下と矛盾する場合は，正本を優先する。

* `AGENTS.md`
* `docs/PROJECT_PLAN.md`
* `docs/adr/`
* `docs/phase*/phase*_tasks.md`
* accepted GitHub issues
* Git history and accepted issue/PR history

未検証の仮説，観察，推測は，必ずその状態を明示する。検証済みの事実と仮説を混同してはならない。

---

## 進捗管理方針

高レベルのフェーズ進捗は本ファイルで管理する。

細かいタスク進捗は，後続で作成する専用ファイルで管理する。

予定：

* `docs/TASK_MAP.md`: 全フェーズの上位進捗
* `docs/phase2/phase2_tasks.md`: Phase 2の詳細進捗
* `docs/codex-runs/`: Codex作業ログ
* `docs/adr/`: 重要な設計判断
* `docs/planning/`: Codexによるタスク分解案
* `docs/codex-notes/`: Codex用の補助メモ

チェックボックスは，review PASS後にのみ更新する。

---

## 文書の役割分担

* 恒久ルール: `AGENTS.md`
* Codex workflow details: `docs/codex/codex_workflow.md`
* Phase 2の現在地と入口: `docs/phase2/phase2_current.md`
* 採択済み Phase 2 task status: `docs/phase2/phase2_tasks.md`
* 全体方針・フェーズ・重要論点: `docs/PROJECT_PLAN.md`
* フェーズ詳細: `docs/phase*/`
* 重要な設計判断: `docs/adr/`
* Codex run completion records: `docs/codex-runs/`

旧 prompt は現在の source of truth ではない。履歴確認が必要な場合は Git history と accepted issue/PR history を参照する。
