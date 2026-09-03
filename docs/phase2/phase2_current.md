# Phase 2 Current

このファイルは，Phase 2作業を開始するときに読む短い現在地ドキュメントである．
現在のbaseline，Active work，Next queue，blocker，優先参照先のみを記載する．
採択判断と過去結果は`phase2_tasks.md`，詳細仕様はconfig・test・task-specific doc・ADRを参照する．

## Goal

Phase 2は，3D物理simulation modelを構築し，数値・物理妥当性と長時間安定性を検証して，
Phase 3が使用するcanonical physical simulation modelをfreezeする．
2D投影，擬似顕微鏡化，clip生成，べん毛本数の識別可能性評価はPhase 3の責務である．

## Current baseline

- 標準simulation profile: `conf/sim_swim_2010.yaml`
- 2010 project profile: supported baseline
- 2015 project / paper profile: Stage A検証は完了，supported採択まではpending
- Phase 3 handoff: freeze済みcanonical modelのsimulation archive，provenance，physical QC
- `n_flagella=4`: diagnostic-only
- `n_flagella>=5`: canonical training scope外
- training candidate: 全時間strict passを要求
- CLI例: `KEY=VALUE`形式を第一表記とする
- 新規runの時間scale: `reference_torque`（`tau_s=eta*b^3/|T|`）。過去再現は`legacy_fixed_tau_s_1`を明示する
- parameter・閾値・dataset条件: config，test，registryを正本とする

採択判断は`docs/phase2/phase2_tasks.md`のDecision Indexから確認する．

## Active work

### Issue #215: 5.0 s body--flagella axis angle診断

Issue #204の2.0 s feature-study referenceと同じ2010 project tau-linked条件を、5.0 s・
`n=1,2,3,4`・attach / phase seed各3条件の36 independent runsへ延長する。目的は
body--flagella axis angleの2.0--5.0 s推移を3D / 2D時系列とwindow集約で確認することであり、
physical model、dataset、ML policyを採択・変更しない。`execution:cs10`のUser-run campaignで、
runbookは`phase2_215_5s_axis_convergence_runbook.md`を正本とする。

## Next queue

1. **Issue #200:** #203のprofile比較結果を受け取り，2010 projectの`dt_star`収束性を判断する．
2. **Issue #61 / #184:** 2015 projectの`dt_star`，torque，本数条件の安定性と計算効率を検証する．
3. **Issue #205:** 物理・数値妥当性の証拠だけを集約し，canonical modelをfreezeする．

## Current blockers

- 2015 supported昇格とStage B: Issue #61、#183、#184の定量評価・採択判断待ち
- canonical model freeze: #200，#61，#184，#93の必要な物理QC・採択判断待ち（#158のv1 r1 failure診断はv1 r2の5 s stability campaignにより置換済み）
- dataset v2: Phase 3のcontrolled study，観測可能性，robustness，dataset freeze gate待ち
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
- Historical Phase 3 handoff baseline: `docs/phase2/phase2_8_phase2_to_phase3_handoff.md`
- Dataset registry: `docs/phase2/phase2_8_dataset_version_registry.md`
- Run summary contract: `docs/phase2/phase2_run_summary_contract.md`
- Issue #204 feature-study reference: `docs/phase2/phase2_204_feature_study_reference.md`
- Issue #215 5.0 s diagnostic: `docs/phase2/phase2_215_5s_axis_convergence_runbook.md`
- Axis / feature contracts: `phase2_7_flag_helix_axis_diagnostics.md`, `phase2_8_flagella_count_feature_definitions.md`
- Model correspondence: `phase2_163_2010_potential_correspondence.md`, `phase2_167_2015_paper_conditions.md`
- ADRs: `docs/adr/`

## Maintenance rule

- Active workは原則1件，Next queueは最大3件とする
- 完了したIssue・詳細結果・parameter値・output pathを蓄積しない
- 採択判断は`phase2_tasks.md`，詳細仕様は正本へ移す
- 更新時は古くなった記述と参照を削除する
- 本文は目安として80行・4,000文字程度に保つ
