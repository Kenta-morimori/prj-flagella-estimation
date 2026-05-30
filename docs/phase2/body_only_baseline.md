# Phase 2.3 Body-Only Torque Baseline

## Purpose

Phase 2.3 では、flagella / hook を含めない body-only 条件で、トルク負荷に対する body 形状 gate を再現可能にする。

現実装では body-only 条件に実モーターは存在しないため、`body_equiv_load.mode=pure_couple` を使って同等トルクを body に加える。

## Gate Metrics

`body_constraint_diagnostics.csv` を source of truth とする。判定は `src/sim_swim/sim/body_shape_gate.py` に集約する。

| priority | category | fail condition |
| ---: | --- | --- |
| 1 | `body_nonfinite` | any gate metric is non-finite |
| 2 | `body_spring` | `body_spring_max_stretch_ratio > 1.0` |
| 3 | `body_bend` | `body_bend_max_error_deg > 60.0` |
| 4 | `body_centerline` | `body_centerline_max_deviation_um > 2.0` |
| 5 | `body_area` | `body_triangle_area_ratio_min < 0.5` |
| pass | `none` | all checks pass |

## Re-Identified Baseline

Condition:

- `n_flagella=0`
- `body_equiv_load.enabled=true`
- `body_equiv_load.mode=pure_couple`
- `duration_s=0.05`
- `dt_star=1.0e-3`
- deterministic, Brownian off

Representative results:

| target torque | expected result | first-fail category |
| ---: | --- | --- |
| `5.0e-20 N m` | pass | `none` |
| `1.0e-19 N m` | fail | `body_spring` |

Thus the current reproducible safe representative is `5.0e-20 N m`, and the first observed break representative is `1.0e-19 N m` under the gate above.

## Verification

Automated check:

- `tests/test_simulation.py::test_phase23_body_only_torque_baseline_safe_and_first_fail`

Expected artifacts:

- `sim/body_constraint_diagnostics.csv`
- `sim/body_constraint_local_diagnostics.csv`
