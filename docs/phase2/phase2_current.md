# Phase 2 Current

このファイルは、CodexがPhase 2作業を開始するときに読む短い現在地ドキュメントである。

現在の採用条件、進行中タスク、次候補、blocker、参照先のみを記載する。
詳細な仕様、数値結果、実験履歴、実行ログは、個別文書、ADR、Issue、PRへ分離する。

## Current baseline

Phase 2は、3D物理シミュレーションと2D擬似顕微鏡動画を、機械学習用データとして利用できる形へ整えるフェーズである。

現在の基準は以下とする。

* 標準simulation profileは `conf/sim_swim_2010.yaml`
* 2010 project profileは現行のsupported baseline
* 2015 project / paper profileは実装済みだが、Issue #168のStage A完了まではpending
* Phase 2からPhase 3 / 4へ渡すbaselineはdataset v1のRUN固定 `n_flagella=1,2,3`
* `n_flagella=4`はdiagnostic-only
* `n_flagella>=5`は現行のcanonical training scope外
* 新しいPhase 2 CLI例は `KEY=VALUE` 形式を第一表記とする
* 詳細なparameter値は、このファイルへ複製せずconfigを正本とする

現在の判断と根拠の一覧は、`docs/phase2/phase2_decision_index.md` を参照する。

## Active work

### Issue #168 / PR #176: 2015 refined model Stage A

2015 project / paper profileについて、motor-offおよびmotor-on RUN条件の短時間安定性を検証している。

完了済み:

* 2015 refined geometryの実装
* project / paper motor dynamicsの実装
* Stage A専用runnerと診断出力の実装
* motor-off pilotの完走
* Stage A判定閾値の固定
* `dt_star=1e-4` / `1e-5`比較用profileと解析導線の実装

未完:

* `dt_star=1e-4` motor-off `0.1 tau`の参考run
* `dt_star=1e-4` motor-on `1 tau`の参考run
* `dt_star=1e-5` motor-on `0.1 tau`との刻み感度比較
* replayによる定性評価
* 必要に応じたcanonical motor-on `1 tau`の実行
* 2015 profileの採否判断

詳細:

* Issue #168
* PR #176
* `docs/phase2/phase2_168_2015_stage_a_validation.md`

## Next queue

Phase 2では、以下を次候補とする。必要に応じて並列に進める。

### 1. Issue #154: 2015 refined modelの次段階

Issue #168のStage A結果を受けて、2015 modelの次段階へ進む。

主な候補:

* Stage Bの評価条件を定義する
* 遊泳速度、body rotation、flagellar rotation、wobbleを評価する
* Kong et al. (2015)の定量値との比較範囲を決める
* project profileとpaper profileの採用目的を分ける

Stage B用の子Issueを作成した場合は、Issue #154ではなく子Issueをここへ記載する。

### 2. Issue #158: n=3長時間非定常回転とproximal failure

dataset v1 r1の `n_flagella=3` 長時間runで残った、非定常回転とproximal flag bond failureの原因を分離する。

結果はdataset v2の物理モデルとquality gateへ反映する。

### 3. Issue #157: dataset v2生成・品質改善

修正後の物理モデルを用いて、RUN固定 `n_flagella=1,2,3` のdataset v2 coreを先にfreezeする。

主な条件:

* source durationは5秒を候補とする
* attach × phaseの独立runを評価する
* training candidateは全時間strict-passを要求する
* RUN–TUMBLEはv2 core完了後の別scopeとして扱う

## Current blockers

* 2015 profileの採用は、Issue #168のStage A完了待ち
* 2015 Stage Bの開始条件は、Stage Aの刻み感度、安定性、定性評価結果待ち
* dataset v2は、Issue #158の長時間failure診断とIssue #157のquality gate確定待ち
* RUN–TUMBLEは、dataset v2 RUN core完了後にIssue #69で扱う
* `n_flagella>=4`は現行training scope外とし、必要時はIssue #124で物理安定化を扱う
* flagella-body貫通の検証が必要な場合はIssue #93で扱う

open Issueの一覧はこのファイルへ列挙せず、Phase 2 roadmapまたは各親Issueを参照する。

## Context routing

Phase 2作業では、以下の順に必要な範囲だけ読む。

1. `AGENTS.md`
2. この `docs/phase2/phase2_current.md`
3. 対象Issue / PR
4. 対象Issueに対応する `docs/phase2/phase2_*.md`
5. 採択済みtask statusが必要な場合のみ `docs/phase2/phase2_tasks.md`
6. 設計判断の根拠が必要な場合のみ `docs/phase2/phase2_decision_index.md` と `docs/adr/`
7. 過去の作業結果が必要な場合のみ、対応する `review_result.json`
8. 上記で不足する場合のみ、Git履歴、PR履歴、generated outputを確認する

原則として最初に全文を読まない:

* `docs/PROJECT_PLAN.md`
* `work_log.md`
* 大規模CSV
* 長いrun log
* `outputs/` 配下のgenerated artifacts

長い文書や出力を確認する場合は、先に検索して対象箇所、列、step範囲を限定する。

## Key references

現在のPhase 2作業で優先する参照先:

* 採択済みtask status: `docs/phase2/phase2_tasks.md`
* 判断根拠一覧: `docs/phase2/phase2_decision_index.md`
* 2015 Stage A: `docs/phase2/phase2_168_2015_stage_a_validation.md`
* Phase 2からPhase 3へのhandoff: `docs/phase2/phase2_8_phase2_to_phase3_handoff.md`
* Dataset version registry: `docs/phase2/phase2_8_dataset_version_registry.md`
* Phase 2 roadmap: Issue #134
* 2015 refined model parent: Issue #154
* Dataset v2 roadmap: Issue #157
* n=3 long-run diagnostics: Issue #158
* 設計判断: `docs/adr/`

個別の実装仕様や過去の定量結果は、対応するIssue、PR、ADR、個別文書を参照する。

## Maintenance rule

このファイルは以下の規則で維持する。

* Active workは原則1件とする
* 並列作業中はActive workを複数記載してよい
* Next queueは最大3件とする
* 完了したIssueはActive workから削除する
* 完了したIssueの詳細結果を追記し続けない
* 完了後も必要な判断は `phase2_decision_index.md` から根拠へ辿れるようにする
* parameter値、閾値、output path、実行日時、pass件数、failure条件一覧を記載しない
* config、ADR、個別文書と同じ説明を重複して記載しない
* 更新時は、新しい情報の追加と同時に古くなった情報を削除する
* 本文は目安として80行・4,000文字以内に保つ
* 上限を超える場合は、情報を個別文書またはdecision indexへ移す

Phase 2 taskの完了条件、commit、push、PR、review resultの一般規則は `AGENTS.md` と `docs/codex/` を正本とし、このファイルへ重複記載しない。
