# Phase 2 Issue #199 タスクA: τ policy × torque × `dt_star`

2010 project の experiment-only campaign として、固定実時間 `0.05 s` における
`tau_fixed_control` と `torque_linked_tau` を比較する。baseline、2015 project、dataset採択を変更しない。

## Conditions

- torque: `1e-21`, `2.5e-20`, `1e-19 N m` per flagellum
- `dt_star`: `1e-3`, `1e-4`、`phase_seed=0`、`stiffness_scales.body=1.0`
- `tau_fixed_control`: 明示的な`legacy_fixed_tau_s_1`（`tau_s=1 s`）
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
uv run python scripts/03_dataset_building/replay_dataset.py \
  --run-dir <campaign-root> --view 3d
```

解析はsimulationを再起動しない。campaign root の `qc_summary.json`、`dt_comparison.csv`、
`tau_policy_comparison.csv`、policy別heatmap、`replay/` を確認する。これらは採択判断ではなく、
#200 の収束・複数seed評価に先立つscreenである。

## Task D: 2015 project τ-linked torque × `dt_star` safety screen

2015 project profile を対象に、Task Aと同じ物理トルク `1e-21`, `2.5e-20`, `1e-19 N m / flagellum`
を使う。`motor.torque_Nm` と `motor.reference_torque_Nm` は各条件で同じ物理値とし、
`time.scale_policy=reference_torque` を使う。これにより material force scale と `tau_s` も条件ごとに
連動する。`duration=1 tau`、`n_flagella=3`、seed 0、`dt_star=1e-3` / `1e-4` の6条件を実行する。

これは locked Stage A threshold による safety screen と replay/heatmap 用の診断であり、
`dt_star=1e-3`の採用、2015 profile のsupported化、dataset条件の変更を意味しない。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_task_d_physical_torque_tau_linked_dt1e3.yaml
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_task_d_physical_torque_tau_linked_dt1e4.yaml

uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind task-d-2015 \
  --reference-run <dt1e-4-run-root> \
  --candidate-run <dt1e-3-run-root> \
  --threshold-contract conf/phase2_validation/2015_stage_a_thresholds.yaml \
  --output-dir <analysis-root>

uv run python scripts/03_dataset_building/replay_dataset.py --run-dir \
  --input-dir <analysis-root>/replay_input \
  --output-dir <analysis-root>/replay \
  --mode render-only --target-frame-count 201 --max-panels-per-grid 6
```
