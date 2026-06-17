# Phase 2 Current

このファイルは，Codex CLI が Phase 2 作業の最初に読む短い現在地ドキュメントである。

原則として 100〜150 行以内に保つ。詳細な実験履歴，長い失敗条件，過去 run の詳細はここに書かず，必要な参照先だけを置く。

## Source of truth priority

* Repository-wide rules: `AGENTS.md`
* Current Phase 2 entry point: `docs/phase2/phase2_current.md`
* Accepted task status: `docs/phase2/phase2_tasks.md`
* Overall project map: `docs/PROJECT_PLAN.md`
* Historical prompts: source of truth ではない

## Context reading policy

* 長い文書を開く前に `rg -n` で関連箇所を絞る。
* 長い Markdown，logs，CSV，`work_log.md` はデフォルトで全文表示しない。
* 過去 run は `review_result.json` を先に読み，必要な場合だけ `work_log.md` を読む。
* `outputs/` の大きな成果物は，文書・summary・review_result で不足するときだけ読む。

## Current Phase 2 status

Phase 2 は，3D物理シミュレーションと2D擬似顕微鏡動画生成を，ML教師データに使える形へ整えるフェーズである。

完了済みの大枠:

* Phase 2.1〜2.5: 初期幾何，body-only，body+hook，single flagellum の短時間安定条件を整理済み。
* Phase 2.6: `material_twist_local_couple` により，単一べん毛の螺旋形状維持と net 回転を両立する代表条件を確認済み。
* 3D出力の間引き用 `output_sampling.fps_out_3d` は導入済み。

現在の主対象:

* Phase 2.7: 複数べん毛で，崩壊せず後方側へ束化または軸整列する条件を探索中。
* Source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/58`
* Main branch for ongoing work: `feature/phase2-58-posterior-bundling-swim`

## Current blockers

* 後方螺旋中心軸条件で，`duration_s=0.5` の代表 sweep とユーザー目視レビューが未完了。
* posterior bundle / partial bundle の成功条件はまだ確定していない。
* 既存診断では，代表トルクで hook drift / hook length fail が先行し，低トルクでは no_bundle のまま残る条件がある。
* Phase 2.8 以降の遊泳挙動・2D動画自然さ・MLデータ整備は，Phase 2.7 の代表条件確定後に進める。

## Default read path

Phase 2 task の開始時は，次の順で必要な範囲だけ読む。

1. `AGENTS.md`
2. この `docs/phase2/phase2_current.md`
3. 対象 issue / PR 本文とコメント
4. `docs/phase2/phase2_tasks.md` の該当 task 周辺
5. 該当する `docs/phase2/phase2_*.md`
6. 関連する `docs/codex-runs/*/review_result.json`

`docs/PROJECT_PLAN.md` は，全体地図や phase-level context が必要な場合だけ読む。

## Key references

* Accepted Phase 2 tasks: `docs/phase2/phase2_tasks.md`
* Overall project map: `docs/PROJECT_PLAN.md`
* Phase 2.6 single flagellum gate: `docs/phase2/phase2_6_helix_retention_gate.md`
* Phase 2.6 torque model evaluation: `docs/phase2/phase2_6_torque_transmission_model_evaluation.md`
* Phase 2.7 bundling plan: `docs/phase2/phase2_7_bundling_stability_plan.md`
* Phase 2.7 helix axis diagnostics: `docs/phase2/phase2_7_flag_helix_axis_diagnostics.md`
* Material twist ADR: `docs/adr/0004_phase2_material_frame_twist_transmission.md`

## Recent run anchors

* `docs/codex-runs/20260616_121312_phase2_7_issue58_helix_axis_integration/review_result.json`
* `docs/codex-runs/20260609_150211_phase2_issue58_posterior_bundling_diagnostic/review_result.json`
* `docs/codex-runs/20260609_014607_phase2_fps_out_3d_sampling/review_result.json`

## Usually do not read first

* `prompts/`: historical prompt files, not source of truth.
* `work_log.md`: read only after the matching `review_result.json` is insufficient.
* Large CSVs and generated outputs under `outputs/`: inspect only when summary docs cannot answer the question.

## Completion update rule

Phase 2 task completion should update only the documents needed by the result:

* Always: `docs/codex-runs/<run-id>/review_result.json`
* Usually: `docs/phase2/phase2_current.md`
* When accepted task status changes: `docs/phase2/phase2_tasks.md`
* When project-level phase status changes: `docs/PROJECT_PLAN.md`
* When a modeling or workflow decision is significant: `docs/adr/`

Do not mark a task complete unless the relevant review result is `PASS`.
