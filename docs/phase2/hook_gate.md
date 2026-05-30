# Phase 2.4 Hook Gate Baseline

## Purpose

Phase 2.4 では、`body + hook` 相当の `minimal_basal_stub` 条件で、motor-on 時の hook 主導 failure を再現可能にする。

この段階では full flagellum の螺旋維持や多本べん毛の束化は対象外とする。対象は basal 近傍、特に attach-first hook link と first-second link の診断である。

## Standard Sweep Command

```bash
uv run python -m scripts.01_simulate_swimming.run_motor_scale_sweep \
  --target local_hook_scale \
  --values 1,8 \
  --torques 1.2e-21,4e-21 \
  --duration 0.02 \
  --body-stiffness-scale 50 \
  --output-dir outputs/phase2_4_hook_gate
```

This writes:

- `outputs/phase2_4_hook_gate/local_hook_scale_sweep_summary.csv`
- per-condition `step_summary.csv`
- per-condition `body_constraint_diagnostics.csv`

## Gate Metrics

The first source of truth is `step_summary.csv`.

| metric | role |
| --- | --- |
| `shape_pass_nonbody` | hook / flagella-side shape pass |
| `first_fail_category_nonbody` | first non-body failure category |
| `local_attach_first_rel_err` | attach-first hook link relative error |
| `hook_len_rel_err_max` | maximum hook length relative error |
| `hook_angle_err_max_deg` | hook bend angle error |
| `local_first_second_rel_err` | basal first-second link relative error |

`body_constraint_diagnostics.csv` is also recorded. In the current diagnostic baseline, the failing high-torque condition can also deform body shape by the last step; therefore the hook/body split is interpreted by priority:

1. `first_fail_category_nonbody == "hook"` means the basal hook gate detects the non-body failure.
2. `body_shape_pass` and body metrics are retained in the sweep summary to show whether body deformation is also present by the final step.
3. If `first_fail_category_nonbody == "none"` but body gate fails, the combined `first_fail_category` is a body category.

## Re-Identified Representatives

Condition:

- `stub_mode=minimal_basal_stub`
- `body_stiffness_scale=50`
- `duration_s=0.02`
- `target=local_hook_scale`
- `values=1,8`
- deterministic, Brownian off

| torque | local_hook_scale | expected result | first-fail category |
| ---: | ---: | --- | --- |
| `1.2e-21 N m` | `1` | pass | `none` |
| `1.2e-21 N m` | `8` | pass | `none` |
| `4.0e-21 N m` | `1` | fail | `hook` |
| `4.0e-21 N m` | `8` | fail | `hook` |

The current boundary is torque-dominated for this minimal sweep: `1.2e-21 N m` is the pass representative and `4.0e-21 N m` is the hook first-fail representative for both tested hook scales.

## Modeling Note

The hook potential threshold remains `90 deg`. Current failures in this diagnostic baseline are dominated by hook length / attach-first relative error rather than hook angle error alone. This means `local_hook_scale` is not a simple monotonic control of failure in the current model; the sweep is used to keep the failure mode traceable, not to claim a final physical optimum.

## Verification

Automated checks:

- `tests/test_motor_scale_sweep.py::test_phase24_local_hook_scale_sweep_reproduces_hook_gate`
- `tests/test_simulation.py::test_phaseb1_known_failure_case_detected_minimal_stub_005s`
- `tests/test_simulation.py::test_phaseb2_break_torque_observed_at_1_2e21_minimal_stub_005s`
- `tests/test_plot_motor_scale_collapse_heatmap.py`
