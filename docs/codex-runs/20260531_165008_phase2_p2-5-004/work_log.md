# P2-5-004 single flagellum 短時間安定回転 gate

## 作業内容

- `stub_mode=full_flagella`, `n_flagella=1` の短時間 motor-on 条件について、safe/fail representatives を固定した。
- `src/sim_swim/sim/single_flagellum_gate.py` を追加し、以下を同時に判定する gate を作成した。
  - finite pass
  - non-body shape pass
  - body shape pass
  - motor split counter
  - `flag_phase_rate_hz` による回転 activity
  - hook / flagellum の主要形状診断
- `1.2e-21 N m` を safe representative、`4.0e-21 N m` を `flag` first-fail representative として pytest で固定した。
- `scripts.01_simulate_swimming.run_motor_scale_sweep` に `--stub-mode` と `--n-flagella` を追加し、full flagellum sweep を標準コマンドで実行できるようにした。
- `scripts/plot_motor_scale_collapse_heatmap.py` は phase-specific 実体への wrapper であり、Phase2 実行系の責務境界を明確にするため削除した。テストは `scripts/01_simulate_swimming/plot_motor_scale_collapse_heatmap.py` を直接読むように更新した。
- `docs/phase2/single_flagellum_stability.md`, `docs/phase2/TASKS.md`, `docs/PROJECT_PLAN.md`, `scripts/README.md` を更新した。

## 仮説と結果

仮説:

- 低トルクでは full flagellum が短時間の回転 activity を持ちながら形状 gate を通過する。
- 高トルクでは hook/body ではなく flagellum bond / bend / torsion が first-fail になる。

結果:

- `1.2e-21 N m`: pass。body / non-body gate を通過し、`median(abs(flag_phase_rate_hz)) > 1.0`。
- `4.0e-21 N m`: fail。body gate は pass、hook relative errors は hook threshold 未満のまま、`flag_bond_rel_err_max`, `flag_bend_err_max_deg`, `flag_torsion_err_max_deg` が gate を超える。

## 検証

- `uv run pytest tests/test_simulation.py -k "phase25 or phaseb_full or phase3_full"`
- `uv run pytest tests/test_motor_scale_sweep.py`
- `uv run pytest tests/test_plot_motor_scale_collapse_heatmap.py`
- `uv run pytest tests/test_motor_forces.py`
- `uv run pytest tests/test_simulation.py -k phase2_execution_and_diagnostics_guidance_exists_in_project_plan`
- `uv run ruff format --check .`
- `uv run ruff check .`
- `uv run pytest -q`

## 結果

- Review result: PASS
- User visual review: 不要
- ADR: 不要
