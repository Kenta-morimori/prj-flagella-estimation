# Phase 2 Issue #183: 2010 reference torque 先行screen結果

## 目的と範囲

これはreference torqueを採択するための結果ではない。[比較契約](phase2_183_reference_torque_comparison_contract.md)とADR 0016に従い、#61 が比較条件を固定して `dt_star` と実行性能を評価し、#184 が候補を検証する前の、2010 projectにおける短い安全性screenである。この実行は明示した `legacy_fixed_tau_s_1` controlであるため、`tracking-reference` は物性scale連動を確認する条件であって時間相似の主張ではない。新規runの既定値はADR 0019の `reference_torque` である。

全runは `same-real-time=1 s`、`dt_star=1e-4`、10,000 steps、default seedで実行した。確認には各runの `sim/run_summary.json` を用い、巨大な `step_summary.csv` は読み込んでいない。body gateは全runでunavailableである。

## 定量結果

| policy | motor scale | finite | non-body shape | 最初のFAIL [s] | hook最大相対誤差 | flag bond最大相対誤差 | compact evidence |
| --- | ---: | --- | --- | ---: | ---: | ---: | --- |
| fixed-reference | 0.5 | PASS | FAIL | 0.9918 | 25.60 | 18.22 | `.../fixed_reference_same_real_time_motor_scale_0.5/2026-08-11/125458/sim/run_summary.json` |
| fixed-reference | 1.0 | PASS | PASS | - | 0.0912 | 0.0949 | `.../fixed_reference_same_real_time_motor_scale_1/2026-08-11/130553/sim/run_summary.json` |
| fixed-reference | 2.0 | PASS | FAIL | 0.2698 | 32.81 | 25.61 | `.../fixed_reference_same_real_time_motor_scale_2/2026-08-11/131650/sim/run_summary.json` |
| tracking-reference | 0.5 | PASS | PASS | - | 0.0779 | 0.0756 | `.../tracking_reference_same_real_time_motor_scale_0.5/2026-08-11/132750/sim/run_summary.json` |
| tracking-reference | 2.0 | PASS | FAIL | 0.0001 | 43.70 | 20.15 | `.../tracking_reference_same_real_time_motor_scale_2/2026-08-11/133847/sim/run_summary.json` |

`tracking-reference` scale 1.0は、fixed-reference scale 1.0と数値条件が同一であるため、重複実行していない。

fixed-reference scale 0.5は、0.9918 sでflagが4 step続けてFAILし、その後0.9999 sまでhook FAILが78 step継続した。従って、renderで0.96–1.00 sに見える形状崩壊は、定量gate上は最後の8.1 msのpersistent failureである。有限値PASSは形状・遊泳安定性PASSを意味しない。

## hook-only diagnostic

tracking-reference scale 2.0のfailureをhookだけで改善できるかを調べるため、`local_hook_scale` だけを `1.0 / 1.25 / 1.5 / 2.0` に変えた0.02 s（200 step）の診断を行った。入力configは `conf/phase2_multi_run/reference_torque_2010_tracking_scale2_hook_diagnostic.yaml`、集計は `outputs/phase2_reference_torque/2010_project/tracking_reference_same_real_time_motor_scale_2_hook_diagnostic/2026-08-12/223908/summary.csv` である。

| local hook scale | non-body shape | 最初のFAIL [s] | hook最大相対誤差 | flag bond最大相対誤差 |
| ---: | --- | ---: | ---: | ---: |
| 1.0 | FAIL | 0.0001 | 43.70 | 11.87 |
| 1.25 | FAIL | 0.0001 | 51.51 | 14.23 |
| 1.5 | FAIL | 0.0001 | 46.48 | 8.34 |
| 2.0 | FAIL | 0.0001 | 49.57 | 12.43 |

全条件が最初の観測stepからhook FAILであり、局所hook scale単独ではこの条件を安全な候補へ戻せなかった。そのため、この組合せの1 s confirmationは開始しない。

## 定性確認と次の比較

既存の `fixed_reference_same_real_time_motor_scale_0.5/.../render/swim3d.mp4` はfailure時刻の視覚的な補助証拠として保持する。ただし、FAIL条件の動画を並べた比較replayは、policy優劣・遊泳品質・採択を示すものではないため、このIssueでは作成しない。

次に比較replayを作るのは、body gateもPASSした候補に限る。その際は同一camera、同一frame rate、同一実時間窓、同一seedを固定し、動画とともにspeed、body roll、axis-center rotation、形状gateの時系列を提示する。現時点でnon-body shapeがPASSしたのはfixed 1.0とtracking 0.5のみであり、遊泳安定性・reference torqueの採択については結論を出さない。
