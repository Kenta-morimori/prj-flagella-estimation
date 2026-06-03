# Phase 2.6 helix spin gate rework

- ユーザー目視で `outputs/phase2_6_qualitative_review/2026-05-31/214109/render/swim3d.mp4` が「全く回転していない」と判定されたため、旧 `flag_phase_rate_hz` 指標を見直した。
- `flag_phase_rate_hz` は root azimuth 由来で、螺旋そのもののスピンとは一致しないことを確認した。
- `src/sim_swim/sim/debug_summary.py` に `flag_helix_spin_phase_deg`, `flag_helix_spin_rate_hz`, `flag_helix_spin_fit_r2` を追加した。
- `src/sim_swim/sim/helix_retention_gate.py` の回転判定を `flag_helix_spin_rate_hz` に変更した。
- `4.0e-21 N m`, `dt_star=2.5e-4`, local scale `(8,5,8,4)` は形状を保つが `motor_no_rotation` として fail することを pytest で固定した。
- `8.0e-21 N m`, `dt_star=1.25e-4`, local scale `(8,5,8,4)` は 400-step CI representative で pass することを pytest で固定した。
- ローカル probe で `3.0e-20 N m`, `dt_star=6.25e-5`, local scale `(4,2,2,2)`, `duration_s=0.25` が 4000-step gate を通過することを確認した。
- ユーザー定性レビュー用に `outputs/phase2_6_spin_review/2026-05-31/222332/` を生成した。
- `uv run pytest tests/test_run_state_fixed.py -k phase26` と `uv run pytest` を実行し、どちらも通過した。
