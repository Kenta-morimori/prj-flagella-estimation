# P2-6-005 single flagellum 螺旋維持 gate

## 作業内容

- `src/sim_swim/sim/helix_retention_gate.py` を追加し、`step_summary.csv` から multi-step helix retention を判定する gate を実装した。
- `tests/test_run_state_fixed.py` に Phase 2.6 の代表条件を追加した。
  - P2-5 break representative (`4.0e-21 N m`, `dt_star=1.0e-3`) は `flag` fail として再現。
  - `dt_star=2.5e-4` のみでは bond / bend / torsion は維持できるが `motor_no_rotation` になることを固定。
  - `dt_star=2.5e-4`, `local_bend_scale=8` では 200-step CI representative で helix retention と回転 activity が両立することを固定。
- `scripts.01_simulate_swimming.run_motor_scale_sweep` に `--dt-star` を追加し、summary CSV に `dt_star` と `dt_internal_s` を記録するようにした。
- `docs/phase2/phase2_6_helix_retention_gate.md`, `docs/phase2/phase2_tasks.md`, `docs/PROJECT_PLAN.md`, `scripts/README.md` を更新した。

## 仮説と結果

仮説:

- `dt_star` 縮小だけでは形状維持と回転 activity が両立しない。
- `dt_star=2.5e-4` と `local_bend_scale=8` の組み合わせで、`4.0e-21 N m` 条件でも bond / bend / torsion と回転 activity を両立できる。

結果:

- default (`dt_star=1.0e-3`): `flag` fail。P2-5 の break representative を再現。
- `dt_star=2.5e-4` only: `motor_no_rotation`。形状は維持するが安定回転条件としては不十分。
- `dt_star=2.5e-4`, `local_bend_scale=8`: PASS。200-step CI representative で hard gate を通過。
- ローカル probe: 同じ pass 条件で `duration_s=0.25` の 1000-step run も PASS。これは単一べん毛・決定論的条件での定量 gate 確認であり，多本べん毛，束化，遊泳軌跡，2D動画自然さは未確認。

## 検証

- `uv run pytest tests/test_run_state_fixed.py -k phase26`
- `uv run pytest tests/test_run_state_fixed.py`
- `uv run pytest tests/test_motor_scale_sweep.py`
- `uv run pytest tests/test_simulation.py -k phase2_execution_and_diagnostics_guidance_exists_in_project_plan`
- `uv run ruff format --check .`
- `uv run ruff check .`
- `uv run pytest -q`

## 結果

- Review result: PASS
- User visual review: 不要
- ADR: 不要
