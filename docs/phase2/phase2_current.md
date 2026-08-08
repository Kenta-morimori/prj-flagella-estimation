# Phase 2 Current

このファイルは，Phase 2作業を開始するときに読む短い現在地ドキュメントである．
現在のbaseline，Active work，Next queue，blocker，優先参照先のみを記載する．
採択判断と過去結果は`phase2_tasks.md`，詳細仕様はconfig・test・task-specific doc・ADRを参照する．

## Current baseline

- 標準simulation profile: `conf/sim_swim_2010.yaml`
- 2010 project profile: supported baseline
- 2015 project / paper profile: 実装済み，Issue #168 Stage A完了まではpending
- Phase 3 / 4 handoff: dataset v1のRUN固定`n_flagella=1,2,3`
- `n_flagella=4`: diagnostic-only
- `n_flagella>=5`: canonical training scope外
- training candidate: 全時間strict passを要求
- CLI例: `KEY=VALUE`形式を第一表記とする
- parameter・閾値・dataset条件: config，test，registryを正本とする

採択判断は`docs/phase2/phase2_tasks.md`のDecision Indexから確認する．

## Active work

### Issue #168 / PR #176: 2015 refined model Stage A

2015 project / paper profileのmotor-offおよびmotor-on RUN条件について，短時間安定性とintegration step感度を検証している．

完了済み:

- refined 120-bead geometry
- project / paper motor dynamics
- Stage A runner・診断出力・判定閾値
- motor-off pilot
- `dt_star=1e-4` / `1e-5`比較導線

未完:

- motor-off `0.1 tau` reference run
- motor-on reference / short run
- 刻み感度比較
- replay定性評価
- 2015 profile採否

参照: Issue #168，PR #176，`phase2_168_2015_stage_a_validation.md`，Decision P2-D18．

## Next queue

1. **Issue #154:** Stage A結果から2015 Stage Bの評価条件とproject / paper profileの目的を決める．
2. **Issue #158:** dataset v1 r1の`n_flagella=3`長時間非定常回転とproximal failureを診断する．参照: `phase2_158_v1_r1_nf3_proximal_diagnostics.md`，P2-D19．
3. **Issue #157:** 修正後モデルでRUN固定`n_flagella=1,2,3`のdataset v2 coreをfreezeする．5 s source，attach × phase独立run，全時間strict passを候補とする．

## Current blockers

- 2015 supported昇格とStage B: Issue #168の刻み感度・安定性・定性評価待ち
- dataset v2: Issue #158のfailure診断とIssue #157のquality gate確定待ち
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
