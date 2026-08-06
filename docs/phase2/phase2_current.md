# Phase 2 Current

このファイルは，CodexがPhase 2作業を開始するときに読む短い現在地ドキュメントである．

現在のbaseline，Active work，Next queue，blocker，優先参照先のみを記載する．  
採択判断と過去結果は`phase2_tasks.md`，詳細仕様はconfig・test・task-specific doc・ADRを参照する．

## Current baseline

Phase 2は，3D物理simulationと2D擬似顕微鏡動画を，機械学習用データとして利用可能な形へ整えるフェーズである．

現在の基準:

- 標準simulation profile: `conf/sim_swim_2010.yaml`
- 2010 project profile: supported baseline
- 2015 project / paper profile: 実装済み，Issue #168のStage A完了まではpending
- Phase 3 / 4 handoff baseline: dataset v1のRUN固定`n_flagella=1,2,3`
- `n_flagella=4`: diagnostic-only
- `n_flagella>=5`: canonical training scope外
- training candidate: 全時間strict passを要求
- 新しいCLI例: `KEY=VALUE`形式を第一表記とする
- parameter値・閾値・dataset条件: config，test，registryを正本とする

採択判断は`docs/phase2/phase2_tasks.md`のDecision Indexから確認する．

## Active work

### Issue #168 / PR #176: 2015 refined model Stage A

2015 project / paper profileについて，motor-offおよびmotor-on RUN条件の短時間安定性を検証している．

完了済み:

- refined 120-bead geometry
- project / paper motor dynamics
- Stage A runnerと診断出力
- motor-off pilot
- Stage A判定閾値
- `dt_star=1e-4` / `1e-5`比較導線

未完:

- `dt_star=1e-4` motor-off `0.1 tau` reference run
- `dt_star=1e-4` motor-on `1 tau` reference run
- `dt_star=1e-5` motor-on short runとの刻み感度比較
- replayによる定性評価
- 2015 profileの採否判断

参照先:

- Issue #168
- PR #176
- `docs/phase2/phase2_168_2015_stage_a_validation.md`
- Decision P2-D18

## Next queue

### 1. Issue #154: 2015 refined modelの次段階

Issue #168の結果に基づき，Stage Bの評価条件，project / paper profileの目的，Kong et al. (2015)との定量比較範囲を決める．

### 2. Issue #158: n=3長時間failure診断

dataset v1 r1の`n_flagella=3`で残った非定常回転とproximal flag bond failureを診断し，dataset v2前のモデル修正要否を判断する．

参照先:

- `docs/phase2/phase2_158_v1_r1_nf3_proximal_diagnostics.md`
- Decision P2-D19

### 3. Issue #157: dataset v2 core

修正後の物理モデルを用いて，RUN固定`n_flagella=1,2,3`のdataset v2 coreをfreezeする．

候補条件:

- source duration: 5 s
- attach seed × phase seedを独立runとして評価
- 全時間strict passを要求
- RUN–TUMBLEは別scope

## Current blockers

- 2015 profileのsupported昇格: Issue #168 Stage A完了待ち
- 2015 Stage B: Stage Aの刻み感度・安定性・定性評価待ち
- dataset v2: Issue #158のfailure診断とIssue #157のquality gate確定待ち
- RUN–TUMBLE: dataset v2 RUN core完了後にIssue #69で扱う
- `n_flagella>=4`: 現行training scope外．必要時はIssue #124で安定化する
- flagella-body貫通: 必要時はIssue #93で検証する

## Context routing

Phase 2作業では，以下の順に必要な範囲だけ読む．

1. `AGENTS.md`
2. この`phase2_current.md`
3. 対象Issue / PR
4. `phase2_tasks.md`の関連decision
5. currentまたはdecisionから参照されるlive contract / active validation
6. 必要なconfig・schema・test・ADR
7. 必要な場合のみ`review_result.json`，Issue / PR履歴，Git履歴，generated output

長い文書，CSV，run log，`outputs/`配下は，検索して対象範囲を限定してから読む．

## Key references

- Decisions: `docs/phase2/phase2_tasks.md`
- Documentation policy: `docs/codex/phase_document_policy.md`
- 2015 Stage A: `docs/phase2/phase2_168_2015_stage_a_validation.md`
- n=3 long-run diagnostics: `docs/phase2/phase2_158_v1_r1_nf3_proximal_diagnostics.md`
- Phase 2 → Phase 3 handoff: `docs/phase2/phase2_8_phase2_to_phase3_handoff.md`
- Dataset registry: `docs/phase2/phase2_8_dataset_version_registry.md`
- Flagella-axis contract: `docs/phase2/phase2_7_flag_helix_axis_diagnostics.md`
- Feature definitions: `docs/phase2/phase2_8_flagella_count_feature_definitions.md`
- 2010 potential mapping: `docs/phase2/phase2_163_2010_potential_correspondence.md`
- 2015 paper conditions: `docs/phase2/phase2_167_2015_paper_conditions.md`
- ADRs: `docs/adr/`

## Maintenance rule

- Active workは原則1件とし，意図的な並列作業では複数件を許容する
- Next queueは最大3件とする
- 完了したIssueと詳細結果は削除する
- parameter値，閾値，output path，実行日時，pass件数を蓄積しない
- 採択判断は`phase2_tasks.md`へ移す
- Phase固有の詳細仕様はlive contract，config，test，ADRへ置く
- 更新時は古くなった記述と参照も削除する
- 本文は目安として80行・4,000文字程度に保つ
