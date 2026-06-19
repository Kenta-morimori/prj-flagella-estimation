# Phase 2 Current

このファイルは，Codex CLI が Phase 2 作業の最初に読む短い現在地ドキュメントである。

原則として 100〜150 行以内に保つ。詳細な実験履歴，長い失敗条件，過去 run の詳細はここに書かず，必要な参照先だけを置く。

## Source of truth priority

* Repository-wide rules: `AGENTS.md`
* Current Phase 2 entry point: `docs/phase2/phase2_current.md`
* Accepted task status: `docs/phase2/phase2_tasks.md`
* Overall project map: `docs/PROJECT_PLAN.md`
* Historical context: Git history, issue/PR history, and Codex run records

## Context reading policy

* 長い文書を開く前に `rg -n` で関連箇所を絞る。
* 長い Markdown，logs，CSV，`work_log.md` はデフォルトで全文表示しない。
* 過去 run は `review_result.json` を先に読み，必要な場合だけ `work_log.md` を読む。
* `outputs/` の大きな成果物は，文書・summary・review_result で不足するときだけ読む。

## Current Phase 2 status

Phase 2 は，3D物理シミュレーションと2D擬似顕微鏡動画生成を，ML教師データに使える形へ整えるフェーズである。

完了済みの大枠:

* Phase 2.1〜2.5: 初期幾何，body-only，body+hook，single flagellum の短時間安定条件を整理済み。
* Phase 2.6: `material_twist_local_couple` により，単一べん毛の螺旋形状維持と net 回転を両立する代表条件を確認済み。Issue #53 / #54 は完了。
* Phase 2.7: 複数べん毛の束化検証は，近接束ではなく螺旋中心軸の安定整列を成功定義として完了。代表条件は `hook_wrapped_axis_aligned`。
* 3D出力の間引き用 `output_sampling.fps_out_3d` は導入済み。

現在の主対象:

* Phase 2.8 / Issue #65: 親Issue #71 のRUN固定べん毛数差分析に向けて，特徴量カテゴリとdataset出力方針を定義済み。
* Phase 2.8 / Issue #72: 親Issue #71 のRUN固定べん毛数差分析に向けて，複数条件実行スクリプトとdataset builderを追加済み。標準datasetは `seed=0` 固定の `fc_nf1_2_3_6_seed1_dur0p5`。
  raw sample には `trajectory.csv` と `state_archive.npz` を残し，後から 3D / 2D render を再生成できる。
  CLI は `scripts/02_phase2_analysis/`，設定は `conf/phase2_analysis/` に置く。
* Phase 2.8 / Issue #71: 次は生成datasetの分布可視化・plot・解釈を進める。
* Phase 2.9 / Issue #69: Tumble状態を段階実装する。まず設計・診断，次にmotor reversal，polymorph切替，run-and-tumble評価へ分ける。

## Current blockers

* P2-7代表条件は `hook_wrapped_axis_aligned` であり，hook length fail は残る。#71 のRUN本数差評価では前提リスクとして扱う。
* #71 では #72 のdataset出力を使い，body displacement, speed, body axis angle change, angular velocity, wobble RMS, flag helix axis alignment を本数別に比較する必要がある。標準datasetは `seed=0` 固定で，本数差の見え方をまず確認する。
* #69 は一括実装にせず，設計・診断から段階化する必要がある。

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
* Phase 2.8 flagella count feature definitions: `docs/phase2/phase2_8_flagella_count_feature_definitions.md`
* Material twist ADR: `docs/adr/0004_phase2_material_frame_twist_transmission.md`
* Codex workflow details: `docs/codex/codex_workflow.md`

## Recent run anchors

* `docs/codex-runs/20260616_144322_phase2_7_axis_alignment_definition/review_result.json`
* `docs/codex-runs/20260616_211555_phase2_7_axis_plot_readability/review_result.json`
* `docs/codex-runs/20260616_121312_phase2_7_issue58_helix_axis_integration/review_result.json`

## Usually do not read first

* `work_log.md`: read only after the matching `review_result.json` is insufficient.
* Large CSVs and generated outputs under `outputs/`: inspect only when summary docs cannot answer the question.

## Completion update rule

Phase 2 task completion should update only the documents needed by the result:

* Always: `docs/codex-runs/<run-id>/review_result.json`
* Usually: `docs/phase2/phase2_current.md`
* When accepted task status changes: `docs/phase2/phase2_tasks.md`
* When project-level phase status changes: `docs/PROJECT_PLAN.md`
* When a modeling or workflow decision is significant: `docs/adr/`
* When Codex workflow details change: `docs/codex/codex_workflow.md`

Do not mark a task complete unless the relevant review result is `PASS`.
