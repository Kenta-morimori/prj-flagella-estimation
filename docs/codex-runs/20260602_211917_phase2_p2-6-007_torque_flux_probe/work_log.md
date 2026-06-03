# Work log

## 対象

- Phase 2.6 P2-6-007 継続
- material frame / segment twist 予定内容をP2-6-007へ吸収
- torque transmission 実装検証

## 実施内容

- `motor.force_distribution=axial_torque_flux_probe` を追加した。
- `axial_torque_flux_probe` は root bead を原点、flagellum の root-to-tip 主軸を torque flux 軸として、root から先端へ緩やかに減衰する接線力を flagellum beads に分布させる。
- 反作用 torque は body 側へ与え、全体の合力と軸まわり torque balance を保つ。
- `tests/test_motor_forces.py` に torque balance の単体テストを追加した。
- `tests/test_run_state_fixed.py` に短時間 integration gate を追加した。
- P2-6-008として切り出していた内容をP2-6-007へ吸収し、`docs/phase2/phase2_tasks.md` からP2-6-008を削除した。
- `docs/adr/0002_phase2_axial_torque_flux_probe.md` を作成した。
- `docs/phase2/phase2_6_triplet_twist_dof_design.md` に実装検証結果と考察、有効アプローチを追記した。

## 結果

短時間比較 (`duration_s=0.05`, `dt_star=1.0e-4`, `torque=2.0e-20`, `local_scale=(1,1.2,1,1)`):

- `triplet`: fail, `net_abs_flag_helix_spin_revolutions=0.00023`, `helix_to_root_net_rotation_ratio=0.00618`
- `axial_torque_flux_probe`: pass, `net_abs_flag_helix_spin_revolutions=0.11566`, `helix_to_root_net_rotation_ratio=2.81489`
- `distributed_flagellum`: pass, `net_abs_flag_helix_spin_revolutions=0.10966`, `helix_to_root_net_rotation_ratio=1.83882`

長時間確認 (`duration_s=0.5`, `dt_star=1.0e-4`):

- `axial_torque_flux_probe`, `torque=2.0e-20`, `local_scale=(1,1.2,1,1)`: pass
- `net_abs_flag_helix_spin_revolutions=1.08546`
- `flag_helix_spin_direction_consistency=1.0`
- `hook_len_rel_err_max=0.39986`
- `flag_bond_rel_err_max=0.17206`
- `flag_bend_err_max_deg=1.12490`
- `flag_torsion_err_max_deg=3.97803`

## 考察

- `triplet + hook spring fix` だけでは、root方位の運動が螺旋全体のnet回転へ伝わりにくい。
- `axial_torque_flux_probe` は、material frameを持たない現行モデルでも、軸方向torque flux近似が有効であることを示した。
- ただし、これは離れたbeadへ接線力を直接分布させる近似であり、segment twistを状態として保存・輸送する物理モデルではない。
- 短期的には `axial_torque_flux_probe` を有効アプローチとして比較対象に残す。
- 論文モデルへの忠実性を高める本命候補は、material frame / segment twist の明示導入である。

## 検証

- `uv run pytest tests/test_motor_forces.py -q`
- `uv run pytest tests/test_motor_forces.py tests/test_run_state_fixed.py::test_phase26_axial_torque_flux_probe_drives_net_helix_spin -q`
- `uv run pytest tests/test_run_state_fixed.py -q`
