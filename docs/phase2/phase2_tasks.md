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

- status: complete
- source proposal: `docs/planning/phase2_task_proposals.md#proposal-p2-6-005-phase-26`
- branch: `feature/phase2-6-helix-retention-gate`
- goal: Phase 2.5 の flagellum-chain dominated failure に対して、回転 activity を保った multi-step 螺旋維持条件を固定する。
- hypothesis:
  - `flag_phase_rate_hz` は root azimuth 由来の proxy であり、目視の螺旋スピンと一致しない場合がある。
  - `median(abs(flag_helix_spin_rate_hz))` だけでは、フィット jitter や往復揺れを持続回転として誤判定する。
  - `flag_helix_spin_phase_deg` の累積差と方向一貫性を用いると、目視で見える螺旋全体の net 回転を判定できる。
  - デフォルト設定を変えず、実行時 override の `time.dt_star=1.0e-4` で条件を探索する。
  - 従来の `triplet` motor では root 方位や螺旋フィット位相は揺れるが、螺旋全体の累積回転へ十分に伝達されない。
  - `attach-first` hook spring force を復元し、診断用 `motor.force_distribution=distributed_flagellum` を使うと、単一べん毛の螺旋形状維持と net 回転を両立できる。
- acceptance criteria:
  - [x] Phase 2.5 break representative が multi-step gate で `flag` fail として再現される。
  - [x] `dt_star` 縮小や旧代表条件では形状を保っても螺旋スピンが出ず `motor_no_rotation` になることを pytest で固定する。
  - [x] `median(abs(flag_helix_spin_rate_hz))` が高くても net 回転がなければ fail することを pytest で固定する。
  - [x] hard gate の指標と閾値を文書化する。
  - [x] motor torque が螺旋全体の累積回転へ伝達されない原因を実装レベルで切り分ける。
  - [x] `dt_star=1.0e-4` で、螺旋形状維持と net 1回転以上を両立する代表条件を固定する。
  - [x] 生成動画をユーザーが目視し、単一べん毛の定性的な安定回転を確認する。
- verification:
  - `uv run pytest tests/test_helix_retention_gate.py`
  - `uv run pytest tests/test_run_state_fixed.py -k phase26`
  - `uv run pytest tests/test_run_state_fixed.py`
  - `uv run pytest tests/test_motor_scale_sweep.py`
- docs:
  - `docs/phase2/phase2_6_helix_retention_gate.md`

### P2-6-006: motor torque の螺旋 net 回転伝達を診断・修正する

- status: complete
- branch: `feature/phase2-6-helix-retention-gate`
- goal: motor torque が root 方位の揺れや局所変形に消えず、単一 full flagellum の螺旋全体を継続回転させる条件または実装修正を確立する。
- background:
  - `torque=2.5e-20`, `dt_star=1.0e-4`, `local scale=(4,2,2,2)` は形状 gate と瞬間スピン rate では pass したが、0.25 s の累積螺旋回転は 0.00139 回転であり、目視でもほぼ回転しなかった。
  - `dt_star=1.0e-4`, 0.05 s の短時間スクリーニングでは、net 回転が大きい条件はすべて `shape_pass_nonbody=False` となり、形状維持と持続回転の両立条件は見つかっていない。
  - 追加診断で、`attach-first` hook spring が topology には存在するが時間積分で未使用だったことを確認した。
- tasks:
  - [x] motor force couple の作用軸、body/hook/root/flagellum への torque 分配、拘束力との打ち消しを診断する。
  - [x] root azimuth, helix phase, body phase の関係を比較し、root の揺れと螺旋 net 回転を分離する。
  - [x] `attach-first` hook spring force を時間積分へ復元する。
  - [x] `motor.force_distribution=distributed_flagellum` を追加し、参照論文モデルとの差分を ADR に記録する。
  - [x] 修正後に `time.dt_star=1.0e-4` で形状維持と net 回転 gate を再探索する。
  - [x] 代表動画を生成し、ユーザー目視レビューを受ける。

### P2-6-007: triplet motor でのねじれ・回転自由度不足を明確化しprobeで検証する

- status: complete
- branch: `feature/phase2-6-helix-retention-gate`
- goal: `distributed_flagellum` を最終解とせず、root motor torque が hook/root から flagellum chain へ伝搬するために必要な回転自由度・ねじれ自由度を明確化し、実験用probeで仮説を検証する。
- background:
  - `distributed_flagellum` では、トルクが螺旋全体へ伝われば単一べん毛が綺麗に回転することを確認した。
  - しかし `triplet + hook spring fix` では、root 方位の net 回転と螺旋全体の net 回転が分離している。
  - `torque=2.5e-20`, `dt_star=1.0e-4`, `local_scale=(4,2,2,2)` の `triplet` 条件では、0.25 s で `net_abs_flag_root_revolutions=0.8796` に対して `net_abs_flag_helix_spin_revolutions=0.00139`、`helix_to_root_net_rotation_ratio=0.00158` だった。
- terms:
  - material frame: 各 flagellum segment に付随する局所座標系で、segment 接線方向だけでなく断面の向きを表す。bead 位置だけでは表せない「軸まわりの向き」を状態として持つための概念。
  - segment twist: 隣接 segment の material frame 同士が、segment 軸まわりにどれだけ相対回転しているかを表すねじれ量。root torque を弾性的な torsional deformation として chain に蓄える・伝えるために必要。
  - axial torque flux: root から flagellum 先端側へ、segment 軸方向に伝わる torque の流れ。現行の bead-position-only model ではこれを明示的に保存・輸送する状態量がない。
- tasks:
  - [x] 現行 bead-position-only モデルで、material frame / segment twist / axial torque flux が明示されていないことを整理する。
  - [x] root torque を segment chain へ伝える候補として、material frame 導入、segment torsional torque flux、quasi-rigid helical body approximation を比較する。
  - [x] `distributed_flagellum` を diagnostic upper-bound とし、triplet 系モデルの改善量を `helix_to_root_net_rotation_ratio` で評価する。
  - [x] `axial_torque_flux_probe` を実装し、現行 bead-position-only model で torque flux 近似が有効か検証する。
  - [x] `local_twist_transmission_probe` を実装し、segment orientation state が root から先端側へ伝わるか検証する。
  - [x] `time.dt_star=1.0e-4`, 0.5 s, net 1回転以上、shape gate PASS を満たす条件を確認する。
  - [x] 実装検証の結果と考察を記録し、有効なアプローチを提案する。
  - [x] 物理モデル変更を行う場合は ADR を作成し、`time.dt_star=1.0e-4`, 0.5 s, net 1回転以上、shape gate PASS を受入基準にする。
- decision:
  - 本PRはprobe成功までをscopeとし、完全な material frame / segment twist 物理実装は次タスクへ分離する。
  - `axial_torque_flux_probe` は短期比較用の有効アプローチとして残す。
  - `local_twist_transmission_probe` は、orientation/local_twist状態を持つため、主probeとして `axial_torque_flux_probe` より優先する。
  - ただし bead force 変換はまだ近似であり、完全な material frame / segment twist 物理モデルとしては扱わない。
  - `distributed_flagellum` は、螺旋全体へ torque が届いた場合の上限診断として残す。
  - 論文モデルへの忠実性を高める本命候補は material frame / segment twist の導入である。
  - local twist transmission は、べん毛全体に1つのねじれ量を与えるのではなく、隣接segment同士の軸まわり向きの差を局所twistとして扱う。
  - 完全物理実装では、orientation/material frame診断、root torque入力、twist potential、局所force couple変換の順に段階化する。
- verification:
  - `uv run pytest tests/test_helix_retention_gate.py -q`
  - `uv run pytest tests/test_motor_forces.py -q`
  - `uv run pytest tests/test_run_state_fixed.py -q`
- docs:
  - `docs/phase2/phase2_6_triplet_twist_dof_design.md`
  - `docs/adr/0002_phase2_torque_transmission_probes.md`
  - `docs/adr/0003_phase2_local_twist_transmission.md`
  - `docs/adr/0004_phase2_material_frame_twist_transmission.md`

### P2-6-008: material frame / segment twist による局所 torque transmission を実装検証する

- status: complete
- branch: `feature/phase2-6-008-local-twist-physics`
- goal: `triplet` root motor torque を、非局所force injectionではなく material frame / segment twist / 局所force couple により flagellum chain へ伝搬し、単一べん毛の螺旋形状維持とnet回転を両立できるか検証する。
- background:
  - P2-6-007では、`axial_torque_flux_probe` と `local_twist_transmission_probe` により、root torque が螺旋全体へ届けば安定回転できることを確認した。
  - ただし両probeは、bead force への最終変換で離れたflagellum beadへ直接接線力を入れる近似を含む。
  - 現行 torsion force は4 bead dihedral angle による螺旋形状復元力であり、root motor torque を時間発展する内部ねじれ状態として保存・伝搬する仕組みではない。
- tasks:
  - [x] `orientation` と `material frame` の状態量・初期化・出力列を定義する。
  - [x] `local_twist_i = wrap_angle(orientation_{i+1} - orientation_i - rest_twist_i)` を診断する。
  - [x] `U_twist_i = 0.5 * k_twist * local_twist_i^2` を第一候補として扱う。
  - [x] root motor torque を root側orientation/material frameへ入力する。
  - [x] twist state を隣接segmentの局所bead対への force couple に変換し、非局所force injectionを使わない。
  - [x] 既存 torsion force との役割分担を評価し、形状維持と torque transmission の役割分担を記録する。
  - [x] `time.dt_star=1.0e-4`, `duration_s>=0.5` で net 1回転以上、shape gate PASS を確認する。
  - [x] 物理モデル拡張の結果を ADR 0004 に記録する。
- decision:
  - `motor.force_distribution=material_twist_local_couple` を追加し、P2-6-008の主実装とする。
  - 代表条件では `duration_s=0.5`, `time.dt_star=1.0e-4`, `motor.torque_Nm=2.0e-20`, `local_spring_scale=1.2` で `net_abs_flag_helix_spin_revolutions=1.11698` と shape gate PASS を確認した。
  - torsion force OFF + 新手法 ON では shape gate が fail したため、現時点では既存 torsion force を置き換えない。
  - 既存 torsion force は螺旋形状維持、新手法は root torque の保存・伝搬として役割分担する。
- acceptance criteria:
  - default `triplet` と paper-compatible geometry 契約を壊さない。
  - 新しい material-frame系挙動は明示的なmodeまたは設定で有効化する。
  - `net_abs_flag_helix_spin_revolutions >= 1.0`
  - `flag_helix_spin_direction_consistency >= 0.5`
  - `hook_len_rel_err_max <= 0.5`
  - `flag_bond_rel_err_max <= 0.25`
  - `flag_bend_err_max_deg <= 30`
  - `flag_torsion_err_max_deg <= 60`
- docs:
  - `docs/adr/0004_phase2_material_frame_twist_transmission.md`
  - `docs/phase2/phase2_6_triplet_twist_dof_design.md`

## Phase 2.7: multi flagella 非崩壊検証

### P2-7-006: 複数べん毛条件での形状崩壊耐性を段階評価する

- status: proposed
- source proposal: `docs/planning/phase2_task_proposals.md#proposal-p2-7-006-phase-27`
- goal: `n_flagella=3..9` での非崩壊条件帯を把握し、特に `motor.force_distribution=distributed_flagellum` で複数べん毛が破綻しないかを検証する。
- background:
  - Phase 2.6 では `distributed_flagellum` により単一べん毛の螺旋形状維持と net 回転を確認した。
  - 複数べん毛に拡張した場合、distributed torque が flagellum 間干渉、hook 変形、body 変形を増幅して collapse/fly-away を起こす可能性がある。
  - posterior bundling の成立可否は別 issue として扱い、本タスクでは束化の良否ではなく非崩壊性を先に評価する。
- tasks:
  - [ ] `motor.force_distribution=distributed_flagellum`, `time.dt_star=1.0e-4` を明示した `n_flagella=3` representative を作る。
  - [ ] `n_flagella` を段階的に増やし、body/hook/flag の first-fail 分布を整理する。
  - [ ] collapse/fly-away が出る条件を再現可能な diagnostic output として保存する。
  - [ ] posterior bundling 判定は Phase 2.8 または別 issue の scope として分離する。
