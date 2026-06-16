# Phase 2.7 issue #58 helix-axis integration

## 実施内容

- `origin/main` を `feature/phase2-58-posterior-bundling-swim` へmergeし、skill関連の変更を取り込んだ。
- PR #63 の螺旋軸診断差分を PR #60 ブランチへ統合した。
- 旧 `initial_flagellum_axis_from_rear_deg` / `posterior_aligned` / `initial_orientation_mode` は現行実装から外し、未merge互換aliasを残さない方針にした。
- 現行の後方初期条件は `flagella.initial_helix_axis_from_rear_deg` に統一した。

## 検証

- `uv run ruff check ...`
- `uv run ruff format --check ...`
- `uv run pytest tests/test_helix_axis.py tests/test_model_builder.py tests/test_params.py tests/test_phase2_7_bundling_sweep.py tests/test_render_state_and_projection.py tests/test_simulation.py -q`
- `git diff --check`

## 残課題

- `duration_s=0.5`, `time.dt_star=1.0e-4`, `flagella.initial_helix_axis_from_rear_deg=0` の代表sweepと目視レビューは未完了。
- Phase 2.7 全体はまだ `FAIL / diagnostic` 状態で、後方束化成功条件または失敗原因の確定が必要。
