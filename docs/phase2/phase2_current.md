# Phase 2 Current

このファイルは，Phase 2作業を開始するときに読む短い現在地ドキュメントである．
現在のbaseline，Active work，Next queue，blocker，優先参照先のみを記載する．
採択判断と過去結果は`phase2_tasks.md`，詳細仕様はconfig・test・task-specific doc・ADRを参照する．

## Current baseline

- 標準simulation profile: `conf/sim_swim_2010.yaml`
- 2010 project profile: supported baseline
- 2015 project / paper profile: Stage A検証は完了，supported採択まではpending
- Phase 3 / 4 handoff: dataset v1のRUN固定`n_flagella=1,2,3`
- `n_flagella=4`: diagnostic-only
- `n_flagella>=5`: canonical training scope外
- training candidate: 全時間strict passを要求
- CLI例: `KEY=VALUE`形式を第一表記とする
- parameter・閾値・dataset条件: config，test，registryを正本とする

採択判断は`docs/phase2/phase2_tasks.md`のDecision Indexから確認する．

## Active work

### Issue #168 / PR #176: 2015 refined model Stage A

2015 project / paper profileのStage A検証（motor-off `0.1 tau`、motor-on `1 tau`）と、project後方束化の
torque × `dt_star` replay gridを完了した。これは短時間の形状・遊泳安定性の確認であり、supported採択ではない．

完了済み:

- refined 120-bead geometry
- project / paper motor dynamics
- Stage A runner・診断出力・判定閾値
- motor-off pilot
- `dt_star=1e-4` / `1e-5` reference比較（共通defaultは`1e-5`を維持）
- canonical motor-on `1 tau`とreplay
- project後方束化のtorque scale `0.5, 1, 2` × `dt_star=1e-5, 1e-4` の6-cell replay gridとユーザー目視

未完:

- `dt_star`の定量的な有効性・妥当性の説明（Issue #61）
- reference torqueの物理解釈と同一実時間比較のpolicy（Issue #183）
- dataset v2採択向けのtorque・`dt_star`・べん毛数の形状／遊泳安定性評価（Issue #184）
- 2015 profileのsupported採否

参照: Issue #168，PR #176，`phase2_168_2015_stage_a_validation.md`，Decision P2-D18．

## Next queue

1. **Issue #61:** 2010 torque連動time-scaleの`1 tau` campaignは、比較archiveの時刻ずれとmotor反作用不整合を修正して再実行待ち。修正前の2026-08-13出力は採択・比較に使わない。runner/QC契約は`phase2_61_2010_torque_dt_campaign_contract.md`、実行結果の評価後に2015 projectへの反映可否を判断する．
2. **Issue #183:** fixed / tracking reference と同一実時間 / 無次元時間を分離した比較policyを ADR 0016 に固定し、2010先行screenとhook-only diagnosticを比較契約・結果記録へ残した。torque・`dt_star`・dataset v2 の採択は #61・#184 の結果後とする（`phase2_183_reference_torque_comparison_contract.md`）。
3. **Issue #184:** dataset v2採択前に、torque・`dt_star`・べん毛数の安定性と計算効率を検証する．

## Current blockers

- 2015 supported昇格とStage B: Issue #61、#183、#184の定量評価・採択判断待ち
- dataset v2: Issue #158のfailure診断、Issue #184の安定性検証、Issue #157のquality gate確定待ち
- RUN–TUMBLE: dataset v2 RUN core完了後にIssue #69で扱う
- `n_flagella>=4`: 現行training scope外．必要時はIssue #124で安定化する
- flagella-body貫通: 必要時はIssue #93で検証する

## Context routing

1. `AGENTS.md`
2. この`phase2_current.md`
3. 対象Issue / PR
4. `phase2_tasks.md`の関連decision
5. 関連するlive contract / active validation
6. 必要なconfig・schema・test・ADR
7. 必要な場合のみ`review_result.json`，Issue / PR履歴，Git履歴，generated output

長い文書，CSV，run log，`outputs/`配下は，検索して対象範囲を限定してから読む．

## Key references

- Decisions: `docs/phase2/phase2_tasks.md`
- Documentation policy: `docs/codex/phase_document_policy.md`
- 2015 Stage A: `docs/phase2/phase2_168_2015_stage_a_validation.md`
- n=3 diagnostics: `docs/phase2/phase2_158_v1_r1_nf3_proximal_diagnostics.md`
- Phase 3 handoff: `docs/phase2/phase2_8_phase2_to_phase3_handoff.md`
- Dataset registry: `docs/phase2/phase2_8_dataset_version_registry.md`
- Run summary contract: `docs/phase2/phase2_run_summary_contract.md`
- Axis / feature contracts: `phase2_7_flag_helix_axis_diagnostics.md`, `phase2_8_flagella_count_feature_definitions.md`
- Model correspondence: `phase2_163_2010_potential_correspondence.md`, `phase2_167_2015_paper_conditions.md`
- ADRs: `docs/adr/`

## Maintenance rule

- Active workは原則1件，Next queueは最大3件とする
- 完了したIssue・詳細結果・parameter値・output pathを蓄積しない
- 採択判断は`phase2_tasks.md`，詳細仕様は正本へ移す
- 更新時は古くなった記述と参照を削除する
- 本文は目安として80行・4,000文字程度に保つ
