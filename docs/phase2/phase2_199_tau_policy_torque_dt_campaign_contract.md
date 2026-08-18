# Phase 2 Issue #199 タスクA: τ policy × torque × `dt_star`

2010 project の experiment-only campaign として、固定実時間 `0.05 s` における
`tau_fixed_control` と `torque_linked_tau` を比較する。baseline、2015 project、dataset採択を変更しない。

## Conditions

- torque: `1e-21`, `2.5e-20`, `1e-19 N m` per flagellum
- `dt_star`: `1e-3`, `1e-4`、`phase_seed=0`、`stiffness_scales.body=1.0`
- `tau_fixed_control`: 2010 project の `profile_default`（`tau_s=1 s`）
- `torque_linked_tau`: `tau_s=eta*b^3/T`、`dt_internal_s=dt_star*tau_s`
- motor/reference/force torque は同一、force override は使用しない

合計12条件である。`T=1e-21 N m` は両policyの同値controlであり、解析で一致度を記録する。

## Execution

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_tau_policy_torque_dt_0p05s.yaml \
  dry_run=true
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_tau_policy_torque_dt_0p05s.yaml
uv run python scripts/02_phase2_analysis/analyze_2010_torque_dt_stability.py \
  --config conf/phase2_multi_run/2010_project_tau_policy_torque_dt_0p05s.yaml \
  --run-dir <campaign-root>
uv run python scripts/02_phase2_analysis/visualize_2010_torque_dt_stability.py \
  --run-dir <campaign-root>
uv run python scripts/02_phase2_analysis/render_phase2_replay.py \
  config=conf/phase2_multi_run/2010_project_tau_policy_torque_dt_0p05s.yaml \
  run_dir=<campaign-root>
```

解析はsimulationを再起動しない。campaign root の `qc_summary.json`、`dt_comparison.csv`、
`tau_policy_comparison.csv`、policy別heatmap、`replay/` を確認する。これらは採択判断ではなく、
#200 の収束・複数seed評価に先立つscreenである。
