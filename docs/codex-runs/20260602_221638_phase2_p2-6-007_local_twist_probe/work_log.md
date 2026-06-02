# P2-6-007 local twist probe work log

## 目的

Phase 2.6の主課題は、単一べん毛で螺旋形状を崩さず、`dt_star=1.0e-4` 条件で長時間のnet回転を確認することである。

既存の `triplet + hook spring fix` では、root方位proxyは動くが、螺旋全体のnet回転へほぼ伝わらなかった。原因仮説は、現行bead-position-only modelにsegment軸まわりのorientation/twist状態がなく、root torqueを内部状態として保存・輸送できないことである。

## 実装

- `motor.force_distribution=local_twist_transmission_probe` を追加した。
- flagellum segmentごとに内部orientation stateを持ち、root motor torqueでroot側orientationを進める。
- orientation activityをdiffusion/relaxationで先端側へ伝える。
- 伝搬したactivityをsegment重みとして、root軸まわりのbead接線forceへ変換する。
- 反作用torqueはbody側へ与え、全体のnet force / root軸net torqueを釣り合わせる。
- `step_summary.csv` に以下の診断列を追加した。
  - `local_twist_root_orientation_deg`
  - `local_twist_tip_orientation_deg`
  - `local_twist_abs_mean_deg`
  - `local_twist_abs_max_deg`
  - `local_twist_tip_activity_ratio`

## 実験結果

短時間比較 (`duration_s=0.05`, `dt_star=1.0e-4`, `torque=2.0e-20`, `local_scale=(1,1.2,1,1)`):

- `helix_retention_pass=True`
- `net_abs_flag_root_revolutions=0.02072`
- `net_abs_flag_helix_spin_revolutions=0.12133`
- `helix_to_root_net_rotation_ratio=5.85492`
- `flag_helix_spin_direction_consistency=1.00000`
- `local_twist_root_orientation_deg=17.731`
- `local_twist_tip_orientation_deg=0.0134`
- `local_twist_tip_activity_ratio=0.0007566`

長時間代表条件 (`duration_s=0.5`, `dt_star=1.0e-4`, `torque=2.0e-20`, `local_scale=(1,1.2,1,1)`):

- `helix_retention_pass=True`
- `net_abs_flag_helix_spin_revolutions=1.04745`
- `flag_helix_spin_direction_consistency=0.81507`
- `hook_len_rel_err_max=0.31502`
- `flag_bond_rel_err_max=0.17314`
- `flag_bend_err_max_deg=2.52590`
- `flag_torsion_err_max_deg=7.77964`
- `local_twist_root_orientation_deg=66.809`
- `local_twist_tip_orientation_deg=23.273`
- `local_twist_abs_max_deg=8.855`
- `local_twist_tip_activity_ratio=0.34836`

比較条件:

- `torque=1.5e-20`, `duration_s=0.5` は `net_abs_flag_helix_spin_revolutions=0.85138` で `motor_no_rotation` fail。
- `torque=2.0e-20`, `local_scale=(1,1,1,1)` でも `net_abs_flag_helix_spin_revolutions=1.04745` でshape gate PASS。

## 考察

形状崩壊せずに回転できた理由は、root近傍3点の局所変形だけにtorqueを閉じ込めず、orientation activityとして先端側へ伝える内部状態を持たせたためである。これにより、螺旋bead群へ一方向の接線速度が生じ、net回転として観測できる。

ただし、このprobeは完全なmaterial frame / segment twist物理モデルではない。local twist potentialから厳密にbead forceを導出しているわけではなく、orientation activityをbead接線forceの重みに使う近似である。論文モデルへの重大な違反として既定値を変えることは避け、`triplet`を既定値のまま残し、実験用modeとして明示的に選択する。

Phase 2.6の方針としては、`distributed_flagellum`を上限診断、`axial_torque_flux_probe`を短期比較用、`local_twist_transmission_probe`を内部twist状態を持つ妥協案候補として整理する。

## 確認

- `uv run ruff format --check .`: PASS
- `uv run ruff check .`: PASS
- `uv run pytest tests/test_motor_forces.py tests/test_run_state_fixed.py::test_phase26_local_twist_transmission_probe_tracks_twist_state -q`: PASS
- `uv run pytest tests/test_run_state_fixed.py -q`: PASS
- `uv run pytest -q`: PASS, 141 passed
- 代表動画生成: PASS
  - `outputs/phase2_6_local_twist_review/2026-06-02/221954/render/swim3d.mp4`
  - `outputs/phase2_6_local_twist_review/2026-06-02/221954/render2d/projection.mp4`
  - `outputs/phase2_6_local_twist_review/2026-06-02/221954/sim/step_summary.csv`

代表動画生成runのhelix retention summary:

- `helix_retention_pass=True`
- `first_fail_category=none`
- `net_abs_flag_helix_spin_revolutions=1.0474479497573188`
- `flag_helix_spin_direction_consistency=0.8150711113472757`
- `max_hook_len_rel_err=0.31502368424398897`
- `max_flag_bond_rel_err=0.17313621315973723`
- `max_flag_bend_err_deg=2.5258957111137565`
- `max_flag_torsion_err_deg=7.779637196442446`
