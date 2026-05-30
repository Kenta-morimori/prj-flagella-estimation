# Phase 2.2 Initial Geometry Contract

## Purpose

Phase 2.2 では、参照論文 normal state に基づく初期形状を、動画目視ではなく数値契約として検証する。

対象は `initial_geometry_summary.json` と `step_summary.csv` である。`initial_geometry_summary.json` は初期配置そのものの source of truth とし、`step_summary.csv` は短時間実行後に同じ診断軸が破綻していないことを確認する。

## Reference Targets

| metric | target | source |
| --- | ---: | --- |
| flag intra bond length | `0.58 b` | Watari and Larson Table 1 / project MVP |
| helix radius | `0.25 b` | Watari and Larson normal state |
| helix pitch | `2.5 b` | Watari and Larson normal state |
| bending angle | `142 deg` | Watari and Larson normal state |
| torsion angle | `-60 deg` | Watari and Larson normal state |

## Operational Tolerances

| metric | tolerance | reason |
| --- | ---: | --- |
| bond mean | `±2%` | project MVP requirement |
| bond min/max | `±5%` | project MVP requirement |
| helix radius | `±0.035 b` absolute | current `paper_table1` construction prioritizes bend/torsion and derives radius around `0.279 b` |
| helix pitch | `±5%` | current discrete chain remains close to `2.5 b` |
| bending max error | `<= 5 deg` | paper-normal compatibility gate |
| torsion max error | `<= 5 deg` | paper-normal compatibility gate |
| base tangent vs rear direction | `90 ± 10 deg` | current side-attach initial condition places the first bead radially and the first chain segment tangentially |

## Current Model Difference

`paper_table1` mode uses `bond_L_over_b`, `potentials.bend.theta0_deg.normal`, and `potentials.torsion.phi0_deg.normal` as the source of truth. Therefore radius and pitch are derived diagnostics rather than direct generation parameters in this mode.

Current observed derived values for the default full-flagellum condition are approximately:

- radius: `0.279 b`
- pitch: `2.40 b`

This is accepted for Phase 2.2 because bend/torsion agreement is the primary contract for the current paper-table initialization. If a future task changes the generator to satisfy radius/pitch exactly, this tolerance should be tightened and the modeling difference should be recorded.

## Verification

Automated checks:

- `tests/test_model_builder.py::test_paper_table1_mode_matches_target_bend_and_torsion`
- `tests/test_simulation.py::test_phase2_initial_geometry_summary_contract_matches_step0`

Expected artifacts:

- `sim/initial_geometry_summary.json`
- `sim/step_summary.csv`
