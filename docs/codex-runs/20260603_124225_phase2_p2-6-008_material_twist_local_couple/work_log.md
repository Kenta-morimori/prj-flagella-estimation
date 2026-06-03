# P2-6-008 material twist local couple

## Summary

P2-6-008として、`motor.force_distribution=material_twist_local_couple` を追加した。

このmodeは、segment軸まわりのscalar `orientation` を内部状態として使い、root側のtorque activityを先端側へ拡散・緩和させる。そのactivityを各segmentの重みに変換し、隣接bead対へ作用反作用のforce coupleを置く。P2-6-007のprobeと異なり、flagellum beads全体へ一括で接線forceを入れない。

## Verification

- `uv run ruff format --check src/sim_swim/dynamics/forces.py src/sim_swim/dynamics/engine.py src/sim_swim/sim/helix_retention_gate.py tests/test_motor_forces.py tests/test_run_state_fixed.py tests/test_helix_retention_gate.py`
- `uv run ruff check .`
- `uv run pytest tests/test_motor_forces.py tests/test_helix_retention_gate.py tests/test_run_state_fixed.py -q`
- `uv run pytest -q`
- commit hook: `uv run ruff format --check .`, `uv run ruff check .`, `uv run pytest -q`

## Representative Result

条件:

- `motor.force_distribution=material_twist_local_couple`
- `motor.torque_Nm=2.0e-20`
- `time.dt_star=1.0e-4`
- `duration_s=0.5`
- `local_hook_scale=1.0`
- `local_spring_scale=1.2`
- `local_bend_scale=1.0`
- `local_torsion_scale=1.0`

結果:

- `helix_retention_pass=True`
- `first_fail_category=none`
- `net_abs_flag_helix_spin_revolutions=1.1169834466127073`
- `flag_helix_spin_direction_consistency=0.981747947248682`
- `helix_to_root_net_rotation_ratio=8.771178783249558`
- `max_hook_len_rel_err=0.40113477213486026`
- `max_flag_bond_rel_err=0.17372724278455898`
- `max_flag_bend_err_deg=4.3491832837972915`
- `max_flag_torsion_err_deg=6.210473374082372`
- `local_twist_tip_activity_ratio=0.3483590883311593`

## Torsion OFF Counterexample

`local_torsion_scale=0.0` かつ `stiffness_scales.flag_torsion=0.0` では、同じ0.5 s条件で `helix_retention_pass=False` になった。

- `first_fail_category=flag`
- `first_fail_step=146`
- `max_hook_len_rel_err=0.7655563853323256`
- `max_flag_bond_rel_err=0.4783030996100029`
- `max_flag_torsion_err_deg=179.9956070777144`

したがって、現時点では既存torsion forceを置き換えない。既存torsion forceは螺旋形状維持、`material_twist_local_couple` はroot torque伝搬として役割分担する。

## Decision

P2-6-008はPASSとする。新modeは完全な3軸material frame / Cosserat rodではないが、非局所force injectionを避け、0.5 s代表条件でnet 1回転以上とshape gate PASSを満たす最小物理実装である。
