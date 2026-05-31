# Phase 2.2 初期幾何検証契約

## 目的

Phase 2.2 では、参照論文 normal state に基づく初期形状を、動画目視ではなく数値契約として検証する。

対象は `initial_geometry_summary.json` と `step_summary.csv` である。`initial_geometry_summary.json` は初期配置そのものの正本とし、`step_summary.csv` は短時間実行後に同じ診断軸が破綻していないことを確認する。

## 参照目標値

| 指標 | 目標値 | 根拠 |
| --- | ---: | --- |
| べん毛内部 bond 長 | `0.58 b` | Watari and Larson Table 1 / project MVP |
| 螺旋半径 | `0.25 b` | Watari and Larson normal state |
| 螺旋 pitch | `2.5 b` | Watari and Larson normal state |
| bending 角 | `142 deg` | Watari and Larson normal state |
| torsion 角 | `-60 deg` | Watari and Larson normal state |

## 運用許容範囲

| 指標 | 許容範囲 | 理由 |
| --- | ---: | --- |
| bond 平均 | `±2%` | project MVP requirement |
| bond min/max | `±5%` | project MVP requirement |
| 螺旋半径 | 絶対誤差 `±0.035 b` | 現行 `paper_table1` 生成は bend/torsion を優先し、半径は約 `0.279 b` として導出されるため |
| 螺旋 pitch | `±5%` | 現行の離散 chain が `2.5 b` 近傍に収まるため |
| bending 最大誤差 | `<= 5 deg` | paper normal state 互換 gate |
| torsion 最大誤差 | `<= 5 deg` | paper normal state 互換 gate |
| 基部接線と後方方向 | `90 ± 10 deg` | 現行 side-attach 初期条件では first bead を半径方向、最初の chain segment を接線方向に置くため |

## 現行モデルとの差分

`paper_table1` mode では、`bond_L_over_b`, `potentials.bend.theta0_deg.normal`, `potentials.torsion.phi0_deg.normal` を正本として扱う。そのため、この mode における radius と pitch は直接の生成パラメータではなく、導出診断値である。

現行 default full-flagellum 条件で観測される導出値は、おおよそ以下である。

- radius: `0.279 b`
- pitch: `2.40 b`

Phase 2.2 では、現行 paper-table 初期化の主契約を bend/torsion 一致に置くため、この差分を許容する。将来、radius/pitch と bend/torsion を同時に厳密に満たす生成器へ変更する場合は、この許容範囲を狭め、モデル差分として記録する。

## 検証

自動検証:

- `tests/test_model_builder.py::test_paper_table1_mode_matches_target_bend_and_torsion`
- `tests/test_simulation.py::test_phase2_initial_geometry_summary_contract_matches_step0`

想定成果物:

- `sim/initial_geometry_summary.json`
- `sim/step_summary.csv`
