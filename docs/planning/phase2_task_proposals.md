# Phase2 Task Proposals (Codex Draft)

本ファイルは **採択済みタスクではない**。  
`AGENTS.md` / `docs/PROJECT_PLAN.md` / `prompts/00_project_context.md` / `prompts/01_repo_rules.md` を確認した上での、Codexによる提案案である。

## Proposal P2-2-001 (Phase 2.2)
- proposal ID: `P2-2-001`
- title: `論文準拠の初期形状・幾何条件の検証契約を確立する`
- goal: Phase 2.2 の「初期形状・幾何条件」を定量的に判定できる状態にする。
- rationale: Phase 2.2 は後段の安定回転検証の前提条件であり、ここが曖昧だと 2.3 以降の失敗原因切り分けが困難になる。
- scope:
  - 初期形状での `bond/pitch/radius/bending/torsion` の判定ロジックと許容範囲を明文化
  - `step0` と `initial_geometry_summary.json` の整合チェック追加
  - 判定結果を `step_summary.csv` か専用検証出力に残す
- out of scope:
  - 長時間の動的不安定性修正
  - べん毛束化や遊泳の最終再現
- 論文モデル差異フォーカス:
  - 論文の normal state 条件（pitch 2.5b, radius 0.25b, bending 142°, torsion -60°）と現実装の一致度を明示
- 既存実装で確認すべき点:
  - `src/sim_swim/sim/core.py` の初期形状サマリ生成
  - `src/sim_swim/sim/flagella_geometry.py` の幾何生成ロジック
  - `tests/test_simulation.py` / `tests/test_model_builder.py` の既存形状検証
- suggested acceptance criteria:
  1. 初期幾何の主要指標が設定許容範囲を満たすテストが追加される
  2. 許容範囲と根拠（論文値/運用値）が文書化される
  3. 既存の `step_summary.csv` 運用と矛盾しない
- suggested tests:
  - `pytest tests/test_simulation.py -k initial_geometry`
  - `pytest tests/test_model_builder.py -k helix`
- suggested branch name: `feature/phase2-2-geometry-contract`
- suggested GitHub issue title: `Phase2.2: Add paper-aligned initial geometry contract and tests`

## Proposal P2-3-002 (Phase 2.3)
- proposal ID: `P2-3-002`
- title: `body_only トルクあり安定回転の基準帯を再同定する`
- goal: Phase 2.3 の合否判定に使うトルク帯と failure 判定基準を再現可能にする。
- rationale: 現状は side-attach 反映後に旧基準トルクの意味が変化しており、D1未完了の根本原因になっている。
- scope:
  - `body_only` 条件で torque sweep を実施し first-fail を記録
  - break torque と safe band を再同定
  - 結果を `PROJECT_PLAN` の Phase2進捗に反映可能な形で整理
- out of scope:
  - hook/flagella あり条件の最適化
  - 2Dレンダリング改善
- 論文モデル差異フォーカス:
  - 論文の nearly rigid body 前提と現実装（brace spring / stiffness scaling）の挙動差を明示
- 既存実装で確認すべき点:
  - `src/sim_swim/sim/params.py` の `dt_star` / `dt_s` 解釈
  - `tests/test_body_stability.py` の静的安定テスト
  - `src/sim_swim/sim/debug_summary.py` の body diagnostics 出力
- suggested acceptance criteria:
  1. break torque と safe band が再現手順つきで記録される
  2. first-fail category と body指標の優先順位が固定される
  3. 同条件再実行で同じ判定になる
- suggested tests:
  - `pytest tests/test_body_stability.py`
  - `pytest tests/test_simulation.py -k phase0a or phase1`
- suggested branch name: `feature/phase2-3-body-only-baseline`
- suggested GitHub issue title: `Phase2.3: Re-identify body-only torque baseline and first-fail gates`

## Proposal P2-4-003 (Phase 2.4)
- proposal ID: `P2-4-003`
- title: `body+hook トルクあり安定回転と hook 破綻境界を定量化する`
- goal: Phase 2.4 の hook 主導 failure を定量化し、再現可能な gate を確立する。
- rationale: 現在の first-fail で hook が主要カテゴリの一つであり、PhaseD sweep再開にも直結する。
- scope:
  - `minimal_basal_stub` 条件で `local_hook_scale` を含む sweep を標準化
  - `hook_len_rel_err_max` / `local_attach_first_rel_err` の運用閾値候補を定義
  - hook failure と body failure の切り分け基準を文書化
- out of scope:
  - full flagella の螺旋維持改善
  - 多本べん毛束化
- 論文モデル差異フォーカス:
  - hook angle 90°閾値運用と現実装の曲げ拘束挙動の差分整理
- 既存実装で確認すべき点:
  - `scripts/01_simulate_swimming/run_motor_scale_sweep.py`
  - `src/sim_swim/dynamics/engine.py` の motor local scales 適用箇所
  - `tests/test_simulation.py` の PhaseB系テスト
- suggested acceptance criteria:
  1. hook failure 境界を再現できる sweep 条件が確立
  2. `shape_pass_nonbody` と `first_fail_category_nonbody` が期待どおり推移
  3. 判定条件が `PROJECT_PLAN` に反映される
- suggested tests:
  - `pytest tests/test_simulation.py -k phaseb1 or phaseb2`
  - `pytest tests/test_plot_motor_scale_collapse_heatmap.py`
- suggested branch name: `feature/phase2-4-hook-gate`
- suggested GitHub issue title: `Phase2.4: Quantify hook failure boundary under motor-on conditions`

## Proposal P2-5-004 (Phase 2.5)
- proposal ID: `P2-5-004`
- title: `body+hook+single flagellum の安定回転条件を確立する`
- goal: single flagellum 条件でトルクあり短時間安定を再現可能にする。
- rationale: Phase 2.6 以降（螺旋維持・多本化）の前提となる最小実運動条件。
- scope:
  - `stub_mode=full_flagella` の短時間安定条件探索
  - motor split 診断カウンタ (`motor_degenerate_axis_count` 等) を使った fail 判定
  - safe parameter set を1つ確定
- out of scope:
  - 長時間 bundle 再現
  - 2Dデータセット量産
- 論文モデル差異フォーカス:
  - motor torque distribution の忠実度と数値安定化補助項の影響を分離評価
- 既存実装で確認すべき点:
  - `src/sim_swim/dynamics/engine.py` の motor split 実装
  - `src/sim_swim/sim/debug_summary.py` の motor診断列
  - `tests/test_simulation.py` の `phaseb_full_motor_on_first_fail_gate`
- suggested acceptance criteria:
  1. single flagellum で short run の通過条件を1セット以上確立
  2. first-fail が出る条件と出ない条件の差分が説明可能
  3. 診断列の欠損や非有限化がない
- suggested tests:
  - `pytest tests/test_simulation.py -k phaseb_full`
  - `pytest tests/test_motor_forces.py`
- suggested branch name: `feature/phase2-5-single-flagella-stability`
- suggested GitHub issue title: `Phase2.5: Establish stable motor-on single-flagellum regime`

## Proposal P2-6-005 (Phase 2.6)
- proposal ID: `P2-6-005`
- title: `single flagellum の bond / bend / torsion 維持策を探索し multi-step で harden する`
- goal: Phase 2.5 で特定した flagellum-chain dominated failure に対して、維持策を仮説ごとに検証し、1000-2000 step 規模で helix collapse を自動検知できる状態にする。
- rationale: Phase 2.5 では `4.0e-21 N m` 条件で body/hook より先に `flag_bond_rel_err_max` / `flag_bend_err_max_deg` / `flag_torsion_err_max_deg` が破綻することを固定した。次は破綻を観測するだけでなく、どの策で維持できるかを分離評価する必要がある。
- scope:
  - Phase 2.5 の safe representative (`1.2e-21 N m`) と break representative (`4.0e-21 N m`) を入力条件として使う
  - multi-step で `flag_bond_rel_err_max` / `flag_bend_err_max_deg` / `flag_torsion_err_max_deg` を監視
  - `dt_s` / `dt_star`、`local_spring_scale`、`local_bend_scale`、`local_torsion_scale`、motor torque distribution の候補を分けて検証
  - 維持策が回転 activity (`flag_phase_rate_hz`) を潰していないことを確認
  - 時系列劣化の fail 条件（急増・閾値超過）を定義
  - pass case と fail case を再現シナリオとして保存
- out of scope:
  - 多本べん毛の束化最適化
  - 実験動画への適用
- 論文モデル差異フォーカス:
  - “nearly rigid” の実運用定義を、論文記述と現行離散化の間で明文化
  - 物理モデル由来の stiffness 変更と、数値安定化由来の補助拘束を区別する
- 既存実装で確認すべき点:
  - `tests/test_run_state_fixed.py` の長時間テスト設計
  - `src/sim_swim/sim/debug_summary.py` の flag系診断列
  - `conf/sim_swim.yaml` の stiffness・dt_star 設定
- suggested acceptance criteria:
  1. multi-step hard test が追加され CI で再現可能
  2. helix collapse を pass/fail で機械判定できる
  3. bond / bend / torsion 維持策の候補ごとの効果が比較できる
  4. 少なくとも1つの維持策について、回転 activity を保った pass representative または明確な FAIL 診断が記録される
  5. 代表 pass/fail case が文書化される
- suggested tests:
  - `pytest tests/test_run_state_fixed.py`
  - `pytest tests/test_simulation.py -k torsion or flag`
- suggested branch name: `feature/phase2-6-helix-multistep-gate`
- suggested GitHub issue title: `Phase2.6: Add multi-step helix-retention hard gates`

## Proposal P2-7-006 (Phase 2.7)
- proposal ID: `P2-7-006`
- title: `複数べん毛条件での形状崩壊耐性を段階評価する`
- goal: `n_flagella=3..9` での崩壊しない条件帯を把握し、危険域を分類する。
- rationale: 本プロジェクト固有拡張（0〜9本）であり、論文3本条件からの差分管理が必須。
- scope:
  - `n_flagella` を段階的に増やした short run/mid run 検証
  - body/hook/flag の first-fail 分布を整理
  - 不安定化の主因（トルク過大、干渉、剛性不足）仮説を分類
- out of scope:
  - Posterior bundle の定性的最終判断
  - ML学習データ整備
- 論文モデル差異フォーカス:
  - 論文3本条件と、プロジェクト拡張9本条件の運用差異を明示
- 既存実装で確認すべき点:
  - `src/sim_swim/model/builder.py` の付着点選択ロジック
  - `conf/sim_swim.yaml` の `flagella.n_flagella`
  - `scripts/01_simulate_swimming/run_motor_scale_sweep.py` の集計列
- suggested acceptance criteria:
  1. `n_flagella` ごとの安定/不安定帯が再現可能に整理される
  2. first-fail category の分布が出力CSVで確認できる
  3. 既知の catastrophic failure 条件が明示される
- suggested tests:
  - `pytest tests/test_model_builder.py -k flagella`
  - `pytest tests/test_simulation.py -k phase`
- suggested branch name: `feature/phase2-7-multi-flagella-stability`
- suggested GitHub issue title: `Phase2.7: Characterize stability envelope for multi-flagella cases`

## Proposal P2-8-007 (Phase 2.8)
- proposal ID: `P2-8-007`
- title: `後方束化（posterior bundling）判定プロトコルを確立する`
- goal: Phase 2.8 の束化判定を、定量+目視の二層で再現可能にする。
- rationale: 後方束化は数値だけで判定しにくく、ユーザーレビュー運用を組み込む必要がある。
- scope:
  - posterior/side bundling 判定指標の定義（例: inter-flagella距離時系列）
  - 目視レビュー手順（対象動画、観点、FAIL条件）を定義
  - 判定結果を run ログに残す運用案作成
- out of scope:
  - run/tumble の本格導入
  - wall effect 導入
- 論文モデル差異フォーカス:
  - 論文 Fig.3 Normal との見え方差分（side bundling への逸脱）を明文化
- 既存実装で確認すべき点:
  - `scripts/01_simulate_swimming/01_simulate_swimming.py` の出力（render3d/render2d）
  - `AGENTS.md` の User visual review policy
  - `docs/PROJECT_PLAN.md` の定性評価要件
- suggested acceptance criteria:
  1. posterior bundling 判定フローが文書化される
  2. 最低1条件で「定量OK/目視要確認」の判定例を作成
  3. review FAIL/PASS の扱いが `AGENTS.md` と整合
- suggested tests:
  - `pytest tests/test_render_state_and_projection.py`
  - 目視レビュー用の再現コマンド検証（手順テスト）
- suggested branch name: `feature/phase2-8-posterior-bundling-protocol`
- suggested GitHub issue title: `Phase2.8: Define posterior bundling evaluation protocol`

## Proposal P2-9-008 (Phase 2.9)
- proposal ID: `P2-9-008`
- title: `べん毛回転駆動の遊泳挙動を運動指標で検証する`
- goal: 遊泳らしさを位置・速度・角速度指標で定量評価できるようにする。
- rationale: 「回っている」だけでは不十分で、推進・姿勢変化の一貫性検証が必要。
- scope:
  - trajectory と phase 指標から遊泳メトリクスを定義
  - motor torque と推進量の関係を最小限で確認
  - 異常軌跡（fly-away等）の検知ルールを追加
- out of scope:
  - 実データ同化
  - 学習モデル評価
- 論文モデル差異フォーカス:
  - run固定（switching off）でどこまで論文の run 挙動に近づくかを明示
- 既存実装で確認すべき点:
  - `sim/trajectory.csv` 出力項目
  - `step_summary.csv` の `flag_phase_rate_hz` / `flag_body_phase_diff_deg`
  - `tests/test_simulation.py` の phase指標検証
- suggested acceptance criteria:
  1. 遊泳挙動の pass/fail 指標が定義される
  2. 推進不在・異常発散を検知できる
  3. 少なくとも1つの安定条件で指標が基準を満たす
- suggested tests:
  - `pytest tests/test_simulation.py -k phase_rate or trajectory`
  - `pytest tests/test_e2e_cli.py`
- suggested branch name: `feature/phase2-9-swimming-metrics`
- suggested GitHub issue title: `Phase2.9: Add quantitative swimming-behavior validation`

## Proposal P2-10-009 (Phase 2.10)
- proposal ID: `P2-10-009`
- title: `3D→2D擬似顕微鏡出力の整合性契約を確立する`
- goal: 3D状態・2D投影・metadataの整合を機械検証できるようにする。
- rationale: Phase 3入力になるため、座標系・スケール不整合は下流全体の品質劣化に直結する。
- scope:
  - 2D投影の中心化・スケール・fps の契約を定義
  - manifest と出力ファイル群の整合チェックを追加
  - 3D/2Dのタイムサンプル対応を検証
- out of scope:
  - 検出器（YOLO）側の精度改善
  - 学習データ増強設計
- 論文モデル差異フォーカス:
  - 論文非対象の2D投影処理を「プロジェクト拡張」として差分管理対象に固定
- 既存実装で確認すべき点:
  - `src/sim_swim/render/project2d.py`
  - `scripts/01_simulate_swimming/01_simulate_swimming.py` の manifest 更新処理
  - `tests/test_render_state_and_projection.py` / `tests/test_e2e_cli.py`
- suggested acceptance criteria:
  1. 2D出力仕様（fps/pixel_size/centering）が固定化される
  2. manifest と実ファイルの整合が自動テストで担保される
  3. 3D変更時の2D反映漏れを検知できる
- suggested tests:
  - `pytest tests/test_render_state_and_projection.py`
  - `pytest tests/test_e2e_cli.py -k projection`
- suggested branch name: `feature/phase2-10-projection-contract`
- suggested GitHub issue title: `Phase2.10: Define and test 3D-to-2D projection contract`

## Proposal P2-11-010 (Phase 2.11)
- proposal ID: `P2-11-010`
- title: `ML学習用ラベル・metadata出力仕様を確定する`
- goal: Phase 4へ渡すためのラベル対応・メタ情報契約を固定する。
- rationale: べん毛本数ラベルと動画対応が曖昧だと教師データとして使用不能になる。
- scope:
  - `n_flagella` ラベル、seed、設定、run情報の必須項目定義
  - データ分割可能な最小メタ情報スキーマ案作成
  - サンプル出力に対する整合テスト追加
- out of scope:
  - Phase4モデル実装
  - 実験動画への適用評価
- 論文モデル差異フォーカス:
  - 論文再現を超えた「ML教師データ化」の拡張要件を明示
- 既存実装で確認すべき点:
  - `manifest.json` の現行項目
  - `sim/step_summary.csv` / `trajectory.csv` の再利用可能性
  - `scripts/03_train_evaluate.py` との接続前提
- suggested acceptance criteria:
  1. ラベル/metadata必須キーが定義される
  2. 1 run 単位で Phase4入力に必要情報が欠損しない
  3. 出力仕様変更時の互換ポリシーが明記される
- suggested tests:
  - `pytest tests/test_e2e_cli.py -k manifest`
  - `pytest tests/test_run_context.py`
- suggested branch name: `feature/phase2-11-ml-dataset-contract`
- suggested GitHub issue title: `Phase2.11: Finalize ML-ready labels and metadata contract`

## 運用メモ
- 本ファイルは提案集であり、採択はユーザー判断とする。
- 採択後は `docs/phase2/TASKS.md` と GitHub Issue へ移管し、進捗は review PASS ベースで更新する。
