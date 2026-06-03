# Phase 2.3 body-only トルク基準

## 目的

Phase 2.3 では、flagella / hook を含めない body-only 条件で、トルク負荷に対する body 形状 gate を再現可能にする。

現実装では body-only 条件に実モーターは存在しないため、`body_equiv_load.mode=pure_couple` を使って同等トルクを body に加える。

## Gate 指標

`body_constraint_diagnostics.csv` を正本とする。判定は `src/sim_swim/sim/body_shape_gate.py` に集約する。

| 優先順位 | category | fail 条件 |
| ---: | --- | --- |
| 1 | `body_nonfinite` | いずれかの gate 指標が非有限 |
| 2 | `body_spring` | `body_spring_max_stretch_ratio > 1.0` |
| 3 | `body_bend` | `body_bend_max_error_deg > 60.0` |
| 4 | `body_centerline` | `body_centerline_max_deviation_um > 2.0` |
| 5 | `body_area` | `body_triangle_area_ratio_min < 0.5` |
| pass | `none` | すべての判定を通過 |

## 再同定した基準

条件:

- `n_flagella=0`
- `body_equiv_load.enabled=true`
- `body_equiv_load.mode=pure_couple`
- `duration_s=0.05`
- `dt_star=1.0e-3`
- 決定論的条件、Brownian off

代表結果:

| target torque | 期待結果 | first-fail category |
| ---: | --- | --- |
| `5.0e-20 N m` | pass | `none` |
| `1.0e-19 N m` | fail | `body_spring` |

したがって、現行 gate 下での再現可能な safe representative は `5.0e-20 N m`、最初の break representative は `1.0e-19 N m` である。

## 検証

自動検証:

- `tests/test_simulation.py::test_phase23_body_only_torque_baseline_safe_and_first_fail`

想定成果物:

- `sim/body_constraint_diagnostics.csv`
- `sim/body_constraint_local_diagnostics.csv`
