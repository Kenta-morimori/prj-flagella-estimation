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
- 新規runの時間scale: `reference_torque`（`tau_s=eta*b^3/|T|`）。過去再現は`legacy_fixed_tau_s_1`を明示する
- parameter・閾値・dataset条件: config，test，registryを正本とする

採択判断は`docs/phase2/phase2_tasks.md`のDecision Indexから確認する．

## Active work

### Issue #203 / PR #214: uniform torque profile paired campaign

2010 projectの`root_torque_segment_couples + uniform`を、同一seedの既存`diffusive`と
paired comparisonする27-condition campaignがcs10で実行中である。これはprofile候補の
evidence収集であり、2010 default、dataset v2、`dt_star=1e-3`の数値妥当性を変更しない。

次はcanonical aggregateのcompact artifactを転送し、paired解析・composite replayを行い、
User承認後に候補判断だけを#200へ渡す。実行運用は`cs10_runbook.md`、比較・転送契約は
`phase2_203_uniform_torque_profile_runbook.md`を参照する。

## Next queue

1. **Issue #200:** #203のprofile比較結果を受け取るが、`dt_star`収束性の判断責務は#200に残す。
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
- Issue #204 feature-study reference: `docs/phase2/phase2_204_feature_study_reference.md`
- Axis / feature contracts: `phase2_7_flag_helix_axis_diagnostics.md`, `phase2_8_flagella_count_feature_definitions.md`
- Model correspondence: `phase2_163_2010_potential_correspondence.md`, `phase2_167_2015_paper_conditions.md`
- ADRs: `docs/adr/`

## Maintenance rule

- Active workは原則1件，Next queueは最大3件とする
- 完了したIssue・詳細結果・parameter値・output pathを蓄積しない
- 採択判断は`phase2_tasks.md`，詳細仕様は正本へ移す
- 更新時は古くなった記述と参照を削除する
- 本文は目安として80行・4,000文字程度に保つ
