# Phase 2 Tasks

このファイルは Phase 2 の採択済みタスクを管理する。提案段階の案は `docs/planning/phase2_task_proposals.md` に残し、採択後に本ファイルへ移す。

チェックボックスは review PASS 後にのみ更新する。

## Phase 2.2: 初期形状・幾何条件

### P2-2-001: 論文準拠の初期形状・幾何条件の検証契約を確立する

- status: completed
- source proposal: `docs/planning/phase2_task_proposals.md#proposal-p2-2-001-phase-22`
- branch: `feature/phase2-geometry-body-baseline`
- goal: 初期形状の bond / pitch / radius / bending / torsion を定量判定できる状態にする。
- acceptance criteria:
  - [x] `initial_geometry_summary.json` に target / tolerance / pass/fail が記録される。
  - [x] step summary と初期幾何 summary の整合を pytest で確認する。
  - [x] 許容範囲と論文値・現実装差分を文書化する。
- verification:
  - `uv run pytest tests/test_simulation.py -k phase2_initial_geometry`
  - `uv run pytest tests/test_model_builder.py -k paper_table1`
- docs:
  - `docs/phase2/geometry_contract.md`

## Phase 2.3: body only トルクあり安定回転

### P2-3-002: body_only トルクあり安定回転の基準帯を再同定する

- status: completed
- source proposal: `docs/planning/phase2_task_proposals.md#proposal-p2-3-002-phase-23`
- branch: `feature/phase2-geometry-body-baseline`
- goal: body-only equivalent torque 条件の safe band と first-fail 境界を再現可能にする。
- acceptance criteria:
  - [x] body-only safe torque と break torque の代表条件を pytest で固定する。
  - [x] body diagnostics の pass/fail 閾値と first-fail priority を文書化する。
  - [x] 同条件再実行で同じ判定になる。
- verification:
  - `uv run pytest tests/test_simulation.py -k phase23_body_only_torque_baseline`
- docs:
  - `docs/phase2/body_only_baseline.md`

## Phase 2.4: body + hook トルクあり安定回転

### P2-4-003: body+hook トルクあり安定回転と hook 破綻境界を定量化する

- status: completed
- source proposal: `docs/planning/phase2_task_proposals.md#proposal-p2-4-003-phase-24`
- branch: `feature/phase2-4-hook-gate`
- goal: `minimal_basal_stub` 条件で hook first-fail boundary を再現可能にする。
- acceptance criteria:
  - [x] `local_hook_scale` を含む sweep 条件が標準化される。
  - [x] `shape_pass_nonbody` と `first_fail_category_nonbody` の期待推移が pytest で固定される。
  - [x] hook failure と body failure の切り分け方針が文書化される。
- verification:
  - `uv run pytest tests/test_motor_scale_sweep.py`
  - `uv run pytest tests/test_simulation.py -k "phaseb1 or phaseb2"`
  - `uv run pytest tests/test_plot_motor_scale_collapse_heatmap.py`
- docs:
  - `docs/phase2/hook_gate.md`
