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
  - `time.dt_star` のデフォルトは変えず、実行時 override の `time.dt_star=1.0e-4` で条件を探索する。
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
  - paper-compatible geometry 契約を壊さない。
  - P2-6-008後のdefault `motor.force_distribution` は `material_twist_local_couple` とし、`triplet` は比較・診断用modeとして明示指定で残す。
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

### P2-6-009: トルク伝搬機構を修正した拡張モデルの詳細評価

- status: in_progress
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/54`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/53`
- branch: `feature/phase2-54-torque-transmission-eval`
- goal: `material_twist_local_couple` 導入後の単一べん毛モデルについて、高トルク条件での形状安定性と local stiffness scaling の必要性を定量評価する。
- result:
  - `local_*_scale=1.0`, `duration_s=0.5`, `time.dt_star=1.0e-4` で、`2.0e-20` から `3.0e-20 N m` は PASS、`3.5e-20 N m` 以上は flag 破綻、`1.0e-19 N m` も FAIL。
  - `3.5e-20 N m` で local spring / bend / torsion / hook を個別に `2.0` へ上げても PASS しないため、local scaling は高トルク安定化の主手段として採用しない。
  - `torque x local-scale mode` heatmap でも、`3.5e-20 N m` 以上は `all=1`, `spring=2`, `bend=2`, `torsion=2`, `hook=2` の全 mode で `flag` fail となることを確認した。
  - `local_*_scale=1.0` を baseline とし、非 1.0 は論文モデルとの差分を伴う診断条件として warning を出す。
  - Issue #58 へ渡す代表候補は `2.5e-20 N m`、上限側候補は `3.0e-20 N m`、失敗境界代表は `3.5e-20 N m` とする。
- background:
  - P2-6-008 では、`motor.torque_Nm=2.0e-20`, `time.dt_star=1.0e-4`, `duration_s=0.5`, `motor.local_spring_scale=1.2` で shape gate PASS と net 1回転以上を確認した。
  - ただし `local_spring_scale=1.2` が物理的に必要な拡張なのか、数値安定化なのか、あるいは不要なのかは未評価である。
  - 後方束化・遊泳検証へ進む前に、単一べん毛で安定トルク帯、scaling の必要性、既存 torsion force の役割を整理する。
- tasks:
  - [x] `local_*_scale=1.0` で `1.0e-19 N m` までの torque sweep を行い、安定境界を確認する。
  - [x] `local_spring_scale`, `local_bend_scale`, `local_torsion_scale`, `local_hook_scale` の one-factor sweep を行う。
  - [x] local-scale mode heatmap で、local scale 変更が高トルク collapse を救済しないことを可視化する。
  - [x] shape gate, first-fail category, net 回転数, 方向一貫性, helix/root 回転比を集計する。
  - [x] 代表 PASS 条件で、菌体重心変位・平均速度・body axis 角度変化を確認する。
  - [x] 既存 torsion force は残す前提で評価し、torsion force OFF は必要時の追加診断に留める。
  - [x] Issue #58 へ渡す多べん毛評価用の代表条件を提示する。
- acceptance criteria:
  - [x] `local_*_scale=1.0` で `1.0e-19 N m` までの torque 安定境界が報告されている。
  - [x] `motor.local_spring_scale=1.2` が必要か、不要か、条件付きで必要かが説明されている。
  - [x] scaling が必要な場合、どの局所項が効いているかが one-factor sweep または heatmap で示されている。
  - [x] 代表 PASS / FAIL 条件が `duration_s>=0.5`, `time.dt_star=1.0e-4` で再現可能である。
  - [x] 代表 PASS 条件で、菌体重心変位または平均速度が報告されている。
  - [x] 既存 torsion force を残す前提が評価結果と矛盾しない。
- docs:
  - `docs/phase2/phase2_6_torque_transmission_model_evaluation.md`

## Phase 2.7: multi flagella 後方束化・非崩壊条件探索

### P2-7-006: 複数べん毛で崩壊せず後方束化する条件を探索する

- status: in_progress
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/58`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/53`
- source proposal: `docs/planning/phase2_task_proposals.md#proposal-p2-7-006-phase-27`
- supersedes: `P2-8-007`
- branch: `feature/phase2-58-posterior-bundling-swim`
- goal: `n_flagella=3..9` で、螺旋形状・hook・bodyが崩壊せず、かつ複数べん毛の螺旋中心軸方向が安定的に揃う条件帯を把握する。特にトルク、初期螺旋軸角度、べん毛本数の関係を整理する。
- background:
  - Phase 2.6 では `material_twist_local_couple` により、単一べん毛の螺旋形状維持とnet回転を確認した。
  - 次に必要なのは、複数べん毛で collapse/fly-away せず、複数べん毛軸が安定的に揃う条件を探索することである。
  - 旧P2-7の「非崩壊性」と旧P2-8の「後方束化判定」は分離すると同じ実験を二度行うため、本タスクで統合する。ただし本PRでは「束化」を近接ではなく軸方向の安定整列として定義する。
  - 短時間で後方束化候補を観察するには、第1ビーズを菌体長軸に対して垂直外向きに保ちつつ、初期べん毛軸を菌体後方へ向けた代表条件を作るのが有効である。hook角度は目的指標ではなく、破綻してはいけない制約として扱う。
  - 近接束化は完全な二値判定が難しいため、本PRでは複数べん毛の螺旋中心軸方向が後半80%で平均軸から `15 deg` 以内に揃うことを主判定にする。1本だけ軸方向が外れる場合は、時系列plot上でも外れとして確認する。
  - Issue #54 / PR #59 で、単一べん毛の代表条件として `motor.torque_Nm=2.5e-20`, `time.dt_star=1.0e-4`, `local_*_scale=1.0` を多べん毛評価へ渡すことにした。
  - PR #55 は診断用WIPであり、最新 `main` から派生した本ブランチに必要な実装だけを移植する。
- tasks:
  - [ ] `motor.force_distribution=material_twist_local_couple`, `time.dt_star=1.0e-4`, `motor.torque_Nm=2.5e-20` を基本条件とした `n_flagella=3` representative を作る。
  - [ ] `distributed_flagellum` は診断用比較条件として残し、必要に応じて同じ `n_flagella` / torque で比較する。
  - [ ] `flagella.initial_helix_axis_from_rear_deg=0` を主条件として、第2ビーズ以降の螺旋中心軸を同一の菌体後方方向へ揃えた診断条件を作る。
  - [ ] `motor.torque_Nm` をsweepし、軸整列する条件、軸整列しない条件、hook巻き付き候補、崩壊する条件を分類する。
  - [ ] `n_flagella=3,6,9` を段階評価し、body/hook/flag の first-fail 分布を整理する。
  - [ ] 軸整列候補指標を実装・記録する。候補は、べん毛軸同士のpair angle、平均軸からの偏差、alignment order、後方角である。
  - [ ] 第1ビーズ外向き維持の診断指標として、`local_attach_first_vs_body_axis_angle_deg` と `local_attach_first_vs_body_axis_err_deg` を記録する。
  - [ ] 各べん毛の第2ビーズ以降から螺旋中心軸を推定し、`flag_helix_axis_vs_rear_angle_deg` として後方向きかを記録する。
  - [ ] hook length drift は自動判定だけで採否を決めず、3D螺旋軸overlayと併せて定性評価する。
  - [ ] 1本のみ軸方向が外れるケースを、`axis_not_aligned` として時系列plot上でも確認できるようにする。
  - [ ] collapse/fly-away が出る条件を再現可能な diagnostic output として保存する。
  - [ ] 代表PASS/FAIL/PARTIAL条件について、定量結果と目視レビュー対象動画を記録する。
- acceptance criteria:
  - [ ] `n_flagella=3` で、0.5 s以上、shape gate PASS の代表条件が1つ以上ある。
  - [ ] `n_flagella` と torque ごとの collapse / axis_not_aligned / hook_wrapped_axis_aligned / axis_aligned_stable の分類表がある。
  - [ ] 軸整列の定量指標が `step_summary.csv` または別CSVへ記録される。
  - [ ] 代表条件でべん毛軸角度の時系列plotを報告できる。
  - [ ] 目視レビューが必要な条件では、対象動画・確認観点・自動判定の限界をreview_resultに記録する。
- current diagnostic notes:
  - 2026-06-09時点では、旧 `initial_flagellum_axis_from_rear_deg=10`, `local_*_scale=1.0`, `duration_s=0.5` の条件で、`motor.torque_Nm=2.5e-20` は `n_flagella=3,6,9` すべて `collapse / hook` となった。この旧条件はhook近傍の接線制御であり、以後の主条件にはしない。
  - `0.5e-20..2.0e-20` へ下げても、全条件で `collapse / hook` となり、`bundle_participation_ratio=0.0`, `flag_flag_close_pair_count=0` だった。
  - `2.5e-21` まで下げると `n_flagella=3,9` は0.5秒のshape gateを通過したが、どちらも `no_bundle` であり、束化・接触・反発は確認されなかった。
  - 現時点の未達理由は、代表トルク帯では `hook_drift` が先に出ること、形状を保てる低トルクでは `no_bundle_drive` になることである。
  - 2026-06-15の追加診断では、初期geometryは「第1ビーズ外向き90度、べん毛軸後方10度」を満たしたが、トルク入力後に `attach -> first` が110度前後へ傾き、hook角度ではなくhook長誤差が1.0を越えてfailした。
  - hook巻き付きは実在挙動として許容し得るため、既存strict判定は残しつつ、hook長のみ `2.0` まで許容する `shape_pass_nonbody_hook_len_relaxed` と `phase27_class_hook_len_relaxed` を併記する。
  - hook長relaxed判定では、`n_flagella=3`, `torque=2.5e-20`, `time.dt_star=1.0e-4` のstrict停止時点は `no_bundle` として扱えるが、0.5秒まで継続すると `hook_len_rel_err_max=2.3783` まで伸びて relaxed 判定でも `collapse / hook` になる。
  - `torque=5.0e-21` では0.5秒まで hook長relaxedで `no_bundle` を維持したが、`bundle_participation_ratio=0.0`, `flag_flag_close_pair_count=0` であり、束化駆動は確認できていない。
  - 2026-06-16時点では、`initial_helix_axis_from_rear_deg=0`, `duration_s=0.5`, `time.dt_star=1.0e-4` により、初期螺旋中心軸を菌体後方へ揃える診断条件を追加した。`n_flagella=3` では束化候補なし、`n_flagella=6`, `motor.torque_Nm=5.0e-21` では `shape_pass_nonbody=True` のまま close pair が出た。
  - `n_flagella=6`, `motor.torque_Nm=1.0e-20` と `2.5e-20` では hook fail 後も close pair が残るため、fail条件でも3D動画で束化候補として確認する。
  - 2026-06-16のユーザー目視評価により、明示的な近接束ではなくても、複数べん毛軸が揃い菌体軸の揺れが小さい条件は成功候補とみなせることを確認した。本PRでは菌体軸の揺れ・移動距離は扱わず、軸整列の定量定義へ更新する。
  - 新定義では、`n_flagella=3`, `motor.torque_Nm=2.5e-20`, `duration_s=0.5`, `time.dt_star=1.0e-4`, `initial_helix_axis_from_rear_deg=0` が `hook_wrapped_axis_aligned` になった。hook長failは残るが、平均軸からの最大偏差は `11.8841 deg` で `15 deg` 閾値を満たす。
  - 今後は、hook drift の原因解明、束化駆動の有無評価、多本数拡張、遊泳挙動評価、2D/ML用データ妥当性評価を目的別に分ける。
  - 詳細は `docs/phase2/phase2_7_bundling_stability_plan.md` に記録する。
- docs:
  - `docs/phase2/phase2_7_bundling_stability_plan.md`
  - `docs/phase2/phase2_7_flag_helix_axis_diagnostics.md`

## Phase 2.8: 遊泳挙動の運動指標検証

### P2-8-008: 束化条件で菌体の推進・姿勢変化を定量化する

- status: accepted
- source proposal: `docs/planning/phase2_task_proposals.md#proposal-p2-9-008-phase-29`
- renumbered from proposal: `P2-9-008`
- branch: `feature/phase2-8-swimming-motion-metrics`
- goal: 複数べん毛が回転する条件で、菌体がどの程度推進し、どの程度姿勢を変えるかを定量化する。
- background:
  - べん毛が回転し束化しても、遊泳らしい推進や姿勢安定性が出るとは限らない。
  - べん毛束軸と菌体軸がずれると、菌体が揺れる・首振りする・角度を変える可能性がある。
  - 本数差により、推進速度、姿勢揺らぎ、角度変化が変わる可能性がある。
- tasks:
  - [ ] body重心trajectory、body axis、body angular velocityを時系列出力・集計する。
  - [ ] べん毛束軸を定義し、bundle axis と body axis の角度を出力する。
  - [ ] 遊泳速度、直進性、姿勢角の累積変化、姿勢揺らぎRMS、角速度RMSを指標化する。
  - [ ] `n_flagella=1,3,6,9` で、推進量と姿勢安定性の差を比較する。
  - [ ] fly-awayや異常回転を運動指標で検出する。
- acceptance criteria:
  - [ ] 遊泳指標CSVが再現可能に出力される。
  - [ ] 少なくとも1条件で、推進速度・body axis角度変化・bundle-body軸角度を報告できる。
  - [ ] 本数差による挙動比較の最小表がある。

## Phase 2.9: 動画出力・サンプリング整備

### P2-9-009: 3D/2D動画出力のレビュー向けサンプリングを整備する

- status: accepted
- branch: `feature/phase2-9-output-sampling`
- background:
  - `dt_star=1.0e-4` では内部step数が多く、3D出力を全step保存するとファイル数・動画生成時間が大きくなる。
  - 2D側には `fps_out_2d` があるが、3D側にも同様の `fps_out_3d` が必要である。
- tasks:
  - [ ] `output_sampling.fps_out_3d` を追加する。
  - [ ] 3D render / frame保存 / manifest に `fps_out_3d` を反映する。
  - [ ] 既存の `output_sampling.out_all_steps_3d` との優先順位を定義する。
  - [ ] CLI overrideで `output_sampling.fps_out_3d=100` のように指定できることを確認する。
- acceptance criteria:
  - [ ] `dt_star=1.0e-4`, `duration_s>=0.5` でも3D動画フレーム数を制御できる。
  - [ ] manifestに3D出力サンプリング条件が記録される。
  - [ ] 既存の全step保存挙動を必要時に維持できる。
