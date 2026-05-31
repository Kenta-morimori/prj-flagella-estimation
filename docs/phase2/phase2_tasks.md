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
  - `docs/phase2/phase2_2_geometry_contract.md`

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
  - `docs/phase2/phase2_3_body_only_baseline.md`

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
  - `docs/phase2/phase2_4_hook_gate.md`

## Phase 2.5: body + hook + single flagellum トルクあり安定回転

### P2-5-004: body+hook+single flagellum の安定回転条件を確立する

- status: completed
- source proposal: `docs/planning/phase2_task_proposals.md#proposal-p2-5-004-phase-25`
- branch: `feature/phase2-5-single-flagella-stability`
- goal: `stub_mode=full_flagella`, `n_flagella=1` 条件で短時間 motor-on の safe/fail representatives を再現可能にする。
- hypothesis:
  - 低トルクでは single full flagellum が短時間の回転 activity を持ちながら形状 gate を通過する。
  - 高トルクでは hook/body ではなく flagellum bond / bend / torsion が first-fail になる。
- acceptance criteria:
  - [x] single flagellum の short run 通過条件を1セット以上確立する。
  - [x] first-fail が出る条件と出ない条件の差分を文書化する。
  - [x] 診断列の欠損、非有限値、motor split counter の異常を pytest で検出できる。
- verification:
  - `uv run pytest tests/test_simulation.py -k "phase25 or phaseb_full or phase3_full"`
  - `uv run pytest tests/test_motor_scale_sweep.py`
  - `uv run pytest tests/test_motor_forces.py`
- docs:
  - `docs/phase2/phase2_5_single_flagellum_stability.md`

## Phase 2.6: single flagellum 螺旋形状維持

### P2-6-005: single flagellum の bond / bend / torsion 維持策を探索し multi-step で harden する

- status: user_review_required
- source proposal: `docs/planning/phase2_task_proposals.md#proposal-p2-6-005-phase-26`
- branch: `feature/phase2-6-helix-retention-gate`
- goal: Phase 2.5 の flagellum-chain dominated failure に対して、回転 activity を保った multi-step 螺旋維持条件を固定する。
- hypothesis:
  - `flag_phase_rate_hz` は root azimuth 由来の proxy であり、目視の螺旋スピンと一致しない場合がある。
  - `flag_helix_spin_rate_hz` を用いると、形状維持だけでなく螺旋そのものの安定回転を判定できる。
  - デフォルト設定を変えず、実行時 override の `time.dt_star=1.0e-4` で条件を探索する。
  - `torque=2.5e-20`, `dt_star=1.0e-4`, `local_hook_scale=4`, `local_spring_scale=2`, `local_bend_scale=2`, `local_torsion_scale=2` の組み合わせで、0.25 s / 2500 steps の螺旋維持と螺旋スピンを両立できる。
- acceptance criteria:
  - [x] Phase 2.5 break representative が multi-step gate で `flag` fail として再現される。
  - [x] `dt_star` 縮小や旧代表条件では形状を保っても螺旋スピンが出ず `motor_no_rotation` になることを pytest で固定する。
  - [x] `torque=8.0e-21`, `dt_star=1.25e-4`, `local_hook_scale=8`, `local_spring_scale=5`, `local_bend_scale=8`, `local_torsion_scale=4` が 400-step CI representative で pass する。
  - [x] `torque=2.5e-20`, `dt_star=1.0e-4`, `local_hook_scale=4`, `local_spring_scale=2`, `local_bend_scale=2`, `local_torsion_scale=2` が 2500-step local representative で pass する。
  - [x] hard gate の指標と閾値を文書化する。
  - [ ] 生成動画をユーザーが目視し、単一べん毛の定性的な安定回転を確認する。
- verification:
  - `uv run pytest tests/test_run_state_fixed.py -k phase26`
  - `uv run pytest tests/test_run_state_fixed.py`
  - `uv run pytest tests/test_motor_scale_sweep.py`
- docs:
  - `docs/phase2/phase2_6_helix_retention_gate.md`
