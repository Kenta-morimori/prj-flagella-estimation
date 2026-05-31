# Phase 2.6 dt_star=1.0e-4 spin search

- ユーザー指示に従い、デフォルト設定を変更せず、実行時 override の `time.dt_star=1.0e-4` 固定で条件探索した。
- 前回の `3.0e-20 N m`, `dt_star=6.25e-5` は 0.25 s で約1.1回転相当のため、定性レビューには弱いと判断した。
- `dt_star=1.0e-4`, `duration_s=0.10` でトルクと local scale を探索した。
- `torque=3.0e-20`, local scale `(4,2,2,2)` は `flag_bond_rel_err_max=0.314 > 0.25` で fail した。
- `torque=2.5e-20`, local scale `(4,2,2,2)` は 0.25 s / 2500 steps で pass した。
- 新しい定性レビュー候補として `outputs/phase2_6_dt1e4_spin_review/2026-05-31/230124/` を生成した。
- `AGENTS.md`, `docs/PROJECT_PLAN.md`, `docs/phase2/phase2_tasks.md`, `docs/phase2/phase2_6_helix_retention_gate.md` に `dt_star=1.0e-4` override 方針と新代表条件を記録した。
