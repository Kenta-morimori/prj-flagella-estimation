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
* Phase 2.8 / Issue #72: 親Issue #71 のRUN固定べん毛数差分析に向けて，複数条件実行スクリプトとdataset builderを追加済み。
  raw sample には `trajectory.csv` と `state_archive.npz` を残し，後から 3D / 2D render を再生成できる。
  Issue #78 により，後出し render CLI はデフォルトで3D全step描画を避け，`--fps-out-3d` / `--fps-out-2d` でfpsを指定できる。
  CLI は `scripts/02_phase2_analysis/`，設定は `conf/phase2_analysis/` に置く。
* Phase 2.8 / Issue #76: seed依存の初期条件ばらつきとして，付着点選択 `attach_seed` と初期helix phase `phase_seed` を分離した。標準datasetは `attach_seeds=[0,1,2]` と `phase_seeds=[0,1,2]` の直積を使う `fc_nf1_2_3_6_as3_ps3_dur0p5`。
* Phase 2.8 / Issue #73: `summary.csv` から本数別の特徴量分布plot，QC plot，要約統計，NaN集計，簡易スクリーニングを生成するCLIを追加済み。標準datasetの分析出力は `plots/distributions/`，`plots/qc/`，`analysis/` に生成できる。
* Phase 2.8 / Issue #81: `run_flagella_count_behavior_sweep.py` に `runner.step_summary_stride` / `runner.state_stride` / `runner.sample_order` を追加し，診断・archiveを間引く fast config `conf/phase2_analysis/flagella_count_behavior_dataset_fast.yaml` を追加済み。標準configは従来どおり全step保存を維持する。
* Phase 2 / Issue #84: 実験簡略化のため，render frame保存defaultをOFFにした。`conf/sim_swim.yaml` の `motor.torque_Nm` default は，Phase 2.6 torque評価の第一候補，Phase 2.7代表条件，Phase 2.8 dataset条件と揃え，`2.5e-20` とする。`1.0e-4` は `time.dt_star` の実行時overrideとして使う値であり，モータートルクdefaultにはしない。3D render は RUN/TUMBLE・時刻・実効トルク・`follow_camera_3d` を併記し，dataset directory 指定で全raw sampleを `replays/<sample_id>/` へ一括再描画できる。残タスクの `seeded_surface` 付着点配置では，`attach_seed` 前半を `n_flagella` ごとの center-triangle priority 配置にし，`n_flagella=1..9` の前半seed数を `3,3,1,6,15,20,15,6,1` として明文化した。前半seedのみのdataset用に `sweep.attach_seed_mode=center_priority_prefix` と `conf/phase2_analysis/flagella_count_behavior_dataset_center_prefix.yaml` を追加し，`phase_seeds=[0]` で27 sampleを生成できる。
* Phase 2.8 / Issue #71: 次は #73 の可視化結果を読み，seed数を増やすか，評価時間を伸ばすか，3D特徴量の初期分離性評価へ進むかを判断する。
* Phase 2.9 / Issue #69: Tumble状態を段階実装する。まず設計・診断，次にmotor reversal，polymorph切替，run-and-tumble評価へ分ける。

## Current blockers

* P2-7代表条件は `hook_wrapped_axis_aligned` であり，hook length fail は残る。#71 のRUN本数差評価では前提リスクとして扱う。
* #71 では #72/#73/#76 のdataset・plot出力を使い，body displacement, speed, body axis angle change, angular velocity, wobble RMS, flag helix axis alignment を本数別に比較する必要がある。標準datasetは付着点seedとphase seedを分けるため，seed要因を層別して解釈する。
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
