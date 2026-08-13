# Phase 2 analysis CLI

`01_simulate_swimming/` が作成した既存の simulation / campaign 出力を、再simulationせずに
plan、集計、plot、replay、診断するための entrypoint を置く。各 script は user-facing な薄い
CLI であり、実装本体は原則 `src/sim_swim/analysis/` にある。

新しい作業では、まず下の **共通運用** を使う。Issue 固有 CLI は、対応する active contract または
結果記録の再現が必要な場合だけ使う。物理パラメータ、config schema、出力契約はこの directory では
定義しない。

## 共通運用

| CLI | 用途 | 実装 | 正本 |
| --- | --- | --- | --- |
| `build_run_summary.py` | 既存 diagnostics から first-read 用 `run_summary.json` を生成 | `sim_swim.analysis.run_summary` | [run summary contract](../../docs/phase2/phase2_run_summary_contract.md) |
| `inspect_step_summary.py` | `run_summary.json` を入口に `step_summary.csv` の限定区間だけを読む | `sim_swim.analysis.step_summary_inspection` | [run summary contract](../../docs/phase2/phase2_run_summary_contract.md) |
| `plot_heatmap.py` | sweep / multi-run summary の profile-driven plot を作る | `sim_swim.analysis.phase2_heatmap` | `scripts/README.md` の Heatmap 節 |
| `render_phase2_replay.py` | 保存済み archive から comparison plot / 3D replay を作る | `sim_swim.analysis.phase2_replay` | `scripts/README.md` の Replay Render 節 |

`step_summary.csv` は通常読み込まず、まず `run_summary.json` を確認する。詳細調査が必要なときだけ
`inspect_step_summary.py` に時間窓と列を明示して渡す。

## Active diagnostic / campaign

| CLI | owner | 用途 | 実装 | 詳細手順 |
| --- | --- | --- | --- | --- |
| `plan_2010_torque_dt_stability.py` | #61（実験基盤 #190） | 2010 torque-linked `dt_star` campaign を simulation なしで展開 | `sim_swim.analysis.torque_dt_stability` | [#61 contract](../../docs/phase2/phase2_61_2010_torque_dt_campaign_contract.md) |
| `analyze_2010_torque_dt_stability.py` | #61 | 既存 torque × `dt_star` campaign の QC / 比較集計を再生成 | `sim_swim.analysis.torque_dt_stability_campaign` | [#61 contract](../../docs/phase2/phase2_61_2010_torque_dt_campaign_contract.md) |
| `visualize_2010_torque_dt_stability.py` | #61 | campaign の torque × `dt_star` diagnostic heatmap を作る | `sim_swim.analysis.torque_dt_stability_visuals` | `scripts/README.md` の Issue #190 / #61 節 |
| `diagnose_v1_r1_nf3_failures.py` | #158 | 既存 v1 r1 n=3 raw output の proximal failure を診断 | `sim_swim.analysis.phase2_158_diagnostics` | [#158 diagnostic](../../docs/phase2/phase2_158_v1_r1_nf3_proximal_diagnostics.md) |

これらは現行の検証・診断に必要であり、削除・統合の対象ではない。#193 の fixed-real-time
performance analyzer は PR #194 merge後に、この表へ同じ分類で追加する。

## 結果再現 / policy 参照

| CLI | owner | 用途 | 実装 | 詳細手順 |
| --- | --- | --- | --- | --- |
| `plan_reference_torque_comparison.py` | #183 | fixed / tracking reference torque 比較を simulation なしで計画 | `sim_swim.analysis.reference_torque_comparison` | [#183 contract](../../docs/phase2/phase2_183_reference_torque_comparison_contract.md) |
| `analyze_spring_formulations.py` | #163 | 2010 `legacy` / `fene_fraenkel` 比較を再集計 | `sim_swim.analysis.spring_formulations_workflow` | [#163 correspondence](../../docs/phase2/phase2_163_2010_potential_correspondence.md) |
| `analyze_2015_stage_a.py` | #168 | 2015 refined model Stage A の結果を再評価 | `sim_swim.analysis.stage_a_2015_workflow` | [#168 validation](../../docs/phase2/phase2_168_2015_stage_a_validation.md) |
| `plot_distributions.py` | Phase 2.8 / #118 | 既存 behavior dataset の分布を plot | `sim_swim.analysis.behavior_dataset_distributions` | [dataset registry](../../docs/phase2/phase2_8_dataset_version_registry.md) |
| `analyze_2d_separability.py` | #125 | behavior dataset v1 の 2D separability を再解析 | `sim_swim.analysis.behavior_dataset_separability` | [Phase 3 handoff](../../docs/phase2/phase2_8_phase2_to_phase3_handoff.md) |

closed Issue の CLI は、採択済み結果・比較policyを再現する入口として保持する。通常の新規campaignでは
共通運用 CLI と active diagnostic の対応contractを優先する。

## 実行規約

- 公開commandは `uv run python scripts/02_phase2_analysis/<script>.py ...` とする。
- 実行条件と閾値は config / live contract を正本とし、ここには複製しない。
- 長時間 simulation、sweep、training、render を起動しない。ここで扱うのは既存出力の plan / 集計 / 描画である。
