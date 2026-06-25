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
  - `motor.force_distribution=root_torque_segment_couples` を追加し、P2-6-008の主実装とする。旧名 `material_twist_local_couple` は Issue #88 以降 deprecated alias として扱う。
  - 代表条件では `duration_s=0.5`, `time.dt_star=1.0e-4`, `motor.torque_Nm=2.0e-20`, `local_spring_scale=1.2` で `net_abs_flag_helix_spin_revolutions=1.11698` と shape gate PASS を確認した。
  - torsion force OFF + 新手法 ON では shape gate が fail したため、現時点では既存 torsion force を置き換えない。
  - 既存 torsion force は螺旋形状維持、新手法は root torque の保存・伝搬として役割分担する。
- acceptance criteria:
  - paper-compatible geometry 契約を壊さない。
  - P2-6-008後のdefault `motor.force_distribution` は `root_torque_segment_couples` とし、`triplet` は比較・診断用modeとして明示指定で残す。
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

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/54`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/53`
- branch: `feature/phase2-54-torque-transmission-eval`
- goal: `root_torque_segment_couples` 導入後の単一べん毛モデルについて、高トルク条件での形状安定性と local stiffness scaling の必要性を定量評価する。
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

### P2-6-010: `motor.force_distribution` の正式名を整理する

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/88`
- branch: `feature/phase2-88-force-distribution-naming`
- goal: root torque を flagellum 全体へ分配する方式を一貫した命名で扱い、正式運用modeと過去の診断probeを分離する。
- result:
  - 正式名を `triplet`, `root_torque_segment_couples`, `root_torque_axis_projection` に整理した。
  - 旧 `material_twist_local_couple` は `root_torque_segment_couples`、旧 `distributed_flagellum` は `root_torque_axis_projection` の deprecated alias として正規化する。
  - `axial_torque_flux_probe` と `local_twist_transmission_probe` はコード・設定上の実行modeから削除した。
  - `conf/sim_swim.yaml` と Phase 2.8 dataset config は新名 `root_torque_segment_couples` を使う。
- acceptance criteria:
  - [x] 新名指定で既存代表条件が動作する。
  - [x] 旧名 alias は warning 付きで新名へ正規化される。
  - [x] 削除済み probe mode は設定読み込み時に拒否される。
  - [x] 関連ドキュメントに旧名と新名の対応が残っている。
- docs:
  - `docs/phase2/phase2_current.md`
  - `docs/phase2/phase2_6_helix_retention_gate.md`
  - `docs/adr/0004_phase2_material_frame_twist_transmission.md`

## Phase 2.7: multi flagella 後方束化・非崩壊条件探索

### P2-7-006: 複数べん毛で崩壊せず後方束化する条件を探索する

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/58`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/53`
- source proposal: `docs/planning/phase2_task_proposals.md#proposal-p2-7-006-phase-27`
- supersedes: `P2-8-007`
- branch: `feature/phase2-58-posterior-bundling-swim`
- goal: `n_flagella=3..9` で、螺旋形状・hook・bodyが崩壊せず、かつ複数べん毛の螺旋中心軸方向が安定的に揃う条件帯を把握する。特にトルク、初期螺旋軸角度、べん毛本数の関係を整理する。
- background:
  - Phase 2.6 では `root_torque_segment_couples` により、単一べん毛の螺旋形状維持とnet回転を確認した。
  - 次に必要なのは、複数べん毛で collapse/fly-away せず、複数べん毛軸が安定的に揃う条件を探索することである。
  - 旧P2-7の「非崩壊性」と旧P2-8の「後方束化判定」は分離すると同じ実験を二度行うため、本タスクで統合する。ただし本PRでは「束化」を近接ではなく軸方向の安定整列として定義する。
  - 短時間で後方束化候補を観察するには、第1ビーズを菌体長軸に対して垂直外向きに保ちつつ、初期べん毛軸を菌体後方へ向けた代表条件を作るのが有効である。hook角度は目的指標ではなく、破綻してはいけない制約として扱う。
  - 近接束化は完全な二値判定が難しいため、本PRでは複数べん毛の螺旋中心軸方向が後半80%で平均軸から `15 deg` 以内に揃うことを主判定にする。1本だけ軸方向が外れる場合は、時系列plot上でも外れとして確認する。
  - Issue #54 / PR #59 で、単一べん毛の代表条件として `motor.torque_Nm=2.5e-20`, `time.dt_star=1.0e-4`, `local_*_scale=1.0` を多べん毛評価へ渡すことにした。
  - PR #55 は診断用WIPであり、最新 `main` から派生した本ブランチに必要な実装だけを移植する。
- tasks:
  - [x] `motor.force_distribution=root_torque_segment_couples`, `time.dt_star=1.0e-4`, `motor.torque_Nm=2.5e-20` を基本条件とした `n_flagella=3` representative を作る。
  - [x] `root_torque_axis_projection` は軸投影比較条件として残し、必要に応じて同じ `n_flagella` / torque で比較する。
  - [x] `flagella.initial_helix_axis_from_rear_deg=0` を主条件として、第2ビーズ以降の螺旋中心軸を同一の菌体後方方向へ揃えた診断条件を作る。
  - [x] `motor.torque_Nm` をsweepし、軸整列する条件、軸整列しない条件、hook巻き付き候補、崩壊する条件を分類する。
  - [x] `n_flagella=3,6,9` を段階評価し、body/hook/flag の first-fail 分布を整理する。
  - [x] 軸整列候補指標を実装・記録する。候補は、べん毛軸同士のpair angle、平均軸からの偏差、alignment order、後方角である。
  - [x] 第1ビーズ外向き維持の診断指標として、`local_attach_first_vs_body_axis_angle_deg` と `local_attach_first_vs_body_axis_err_deg` を記録する。
  - [x] 各べん毛の第2ビーズ以降から螺旋中心軸を推定し、`flag_helix_axis_vs_rear_angle_deg` として後方向きかを記録する。
  - [x] hook length drift は自動判定だけで採否を決めず、3D螺旋軸overlayと併せて定性評価する。
  - [x] 1本のみ軸方向が外れるケースを、`axis_not_aligned` として時系列plot上でも確認できるようにする。
  - [x] collapse/fly-away が出る条件を再現可能な diagnostic output として保存する。
  - [x] 代表条件について、定量結果と任意確認対象plotを記録する。
- acceptance criteria:
  - [x] `n_flagella=3` で、0.5 s以上、`hook_wrapped_axis_aligned` の代表条件が1つ以上ある。
  - [x] `n_flagella` と torque ごとの collapse / axis_not_aligned / hook_wrapped_axis_aligned / axis_aligned_stable の分類表がある。
  - [x] 軸整列の定量指標が `step_summary.csv` または別CSVへ記録される。
  - [x] 代表条件でべん毛軸角度の時系列plotを報告できる。
  - [x] 目視レビューは必須ではなく任意確認とし、対象plot・確認観点・自動判定の限界をreview_resultに記録する。
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

## Phase 2.8: べん毛数の違いによるRUN遊泳挙動分析

### P2-8-008: べん毛数分析用の特徴量を定義する

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/65`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/71`
- branch: `feature/phase2-65-feature-definitions`
- goal: 親Issue #71 のRUN固定データセット分析に向けて、特徴量カテゴリ、代表変数、NaN方針、分析用特徴量とML候補特徴量の分離方針を定義する。
- background:
  - #71 は、RUN固定条件の3Dシミュレーション結果を使い、べん毛数差が推進、直進性、姿勢安定性、wobble、べん毛軸と菌体軸の関係に与える影響をデータセットとして評価する親Issueである。
  - #65 はその小タスクとして、実行スクリプトやdataset構築の前に特徴量カテゴリと出力方針を固定する。
  - RUN本数差評価は #71 に含まれる後続実装であり、本タスクでは特徴量定義までを完了範囲とする。
- default conditions:
  - `n_flagella = 1, 2, 3, 6`
  - `motor.enable_switching=false`
  - `motor.force_distribution=material_twist_local_couple`
  - `time.dt_star=1.0e-4`
  - `flagella.initial_helix_axis_from_rear_deg=0`
  - `duration_s=0.5`
- tasks:
  - [x] 特徴量カテゴリ `metadata`, `quality`, `cell_translation`, `cell_orientation`, `flagella_axis`, `cell_flagella_relation`, `diagnostics` を定義する。
  - [x] Feature Registry の正本YAML `conf/phase2_analysis/flagella_count_behavior_features.yaml` を作成し、カテゴリごとの代表変数名を記録する。
  - [x] 各カテゴリが cell / flagella / relation / QC のどれを対象にするか明記する。
  - [x] 定義不能な特徴量を `NaN` として扱い、plot / summary で除外数を記録する方針を明記する。
  - [x] 分析用特徴量とML候補特徴量を分け、`quality` と `diagnostics` を原則ML入力候補に含めない方針を明記する。
  - [x] `dataset_id` と `sample_id` を持つ出力成果物要件を定義する。
- acceptance criteria:
  - [x] 特徴量カテゴリ名が定義されている。
  - [x] 各カテゴリが cell / flagella / relation / QC のどれを対象とするか明確である。
  - [x] 簡易 YAML 形式で、カテゴリと代表変数名が整理されている。
  - [x] 定義不能な特徴量は `NaN` として扱う方針が明記されている。
  - [x] 分析用特徴量と ML 候補特徴量を分ける方針が明記されている。
  - [x] 後続のデータセット構築・分布可視化Issueの入力仕様として利用できる。
- docs:
  - `conf/phase2_analysis/flagella_count_behavior_features.yaml`
  - `docs/phase2/phase2_8_flagella_count_feature_definitions.md`
- follow-up:
  - #71 配下で、複数条件実行、raw output保存、`summary.csv` / `timeseries/<sample_id>.csv` 生成、べん毛数ごとの分布可視化を小タスクに分けて進める。
  - 後続の初期datasetでは `n_flagella=1,2,3,6` と seed 0固定を基本にする。Issue #76 で付着点seedとphase seedの直積条件へ拡張済み。

### P2-8-009: 複数条件の実行スクリプトとdataset構築を追加する

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/72`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/71`
- branch: `feature/phase2-72-flagella-count-dataset`
- goal: RUN固定のべん毛本数差分析に向けて、`n_flagella = 1, 2, 3, 6` を `seed=0` 固定でまとめて実行し、`summary.csv` と sample別 `timeseries` を持つdatasetへ変換できるようにする。
- default conditions:
  - `n_flagella = 1, 2, 3, 6`
  - `seed = 0`
  - `duration_s = 0.5`
  - `time.dt_star = 1.0e-4`
  - `motor.torque_Nm = 2.5e-20`
  - `motor.enable_switching = false`
  - `motor.force_distribution = material_twist_local_couple`
  - `flagella.initial_helix_axis_from_rear_deg = 0`
- tasks:
  - [x] Phase2 analysis用dataset config `conf/phase2_analysis/flagella_count_behavior_dataset.yaml` を追加する。
  - [x] `scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py` で条件表、sample config、sample raw output、`run_manifest.json` を生成する。
  - [x] `scripts/02_phase2_analysis/build_flagella_count_behavior_dataset.py` で `run_manifest.json` から `summary.csv`、`qc_summary.csv`、`timeseries/<sample_id>.csv`、`dataset_manifest.json`、`feature_schema_used.yaml` を生成する。
  - [x] 既存 sample はデフォルトで上書きせず、`--overwrite` 指定時のみ再実行・再生成する。
  - [x] plot、動画出力、ML分類は対象外のまま維持する。
- acceptance criteria:
  - [x] analysis用configが追加されている。
  - [x] `n_flagella = 1, 2, 3, 6`、seed 0固定の条件表を生成できる。
  - [x] 各sampleに一意な `sample_id` が付与される。
  - [x] 各シミュレーション結果が `runs/<run_batch_id>/samples/<sample_id>/raw/` に保存される。
  - [x] `run_manifest.json` が作成される。
  - [x] dataset builder が `run_manifest.json` を入力として dataset を作成できる。
  - [x] `summary.csv` が `datasets/<dataset_id>/summary.csv` に出力される。
  - [x] `timeseries/<sample_id>.csv` が sampleごとに出力される。
  - [x] `dataset_manifest.json` が作成される。
  - [x] `feature_schema_used.yaml` が dataset directory に保存される。
  - [x] `NaN`、fail、relaxed sample を壊さず処理できる。
  - [x] plot、動画出力、ML分類は本Issueの対象外として維持されている。
- tests/checks:
  - `uv run pytest tests/test_flagella_count_behavior_dataset.py tests/test_phase2_7_bundling_sweep.py`
  - `uv run ruff check scripts/02_phase2_analysis tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff format --check scripts/02_phase2_analysis tests/test_flagella_count_behavior_dataset.py`

### P2-8-012: seed依存の初期配置ばらつきを導入する

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/76`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/10`
- branch: `feature/phase2-76-seeded-initial-conditions`
- goal: dataset増強に向けて，seed値により別サンプルと見なせる初期条件ばらつきを導入する。
- default conditions:
  - `flagella.placement_mode = seeded_surface`
  - `flagella.initial_phase_mode = seeded`
  - `seed.attach_seed` はべん毛付着点選択に使う。
  - `seed.phase_seed` は初期helix phase選択に使う。
  - 未指定時はどちらも `seed.global_seed` にfallbackする。
- tasks:
  - [x] `n_flagella=1,2,3,6` の全条件で，seedにより初期付着点またはphaseが変わるようにする。
  - [x] 旧互換の `uniform` 配置と等間隔phaseを明示modeとして残す。
  - [x] Phase 2.8 sweepで `attach_seeds` と `phase_seeds` の直積条件を生成できるようにする。
  - [x] manifest，dataset summary，initial geometry summaryに分離seedを記録する。
- acceptance criteria:
  - [x] 同じ `attach_seed` / `phase_seed` では初期条件が再現する。
  - [x] `attach_seed` だけを変えると付着点が変わる。
  - [x] `phase_seed` だけを変えると付着点は同じで初期helix phaseが変わる。
  - [x] 旧 `sweep.seeds` 形式は `attach_seed=seed`, `phase_seed=seed` として互換動作する。
- tests/checks:
  - `uv run pytest tests/test_params.py tests/test_model_builder.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff check src/sim_swim/model src/sim_swim/sim src/sim_swim/analysis scripts/02_phase2_analysis tests/test_model_builder.py tests/test_params.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff format --check src/sim_swim/model/builder.py src/sim_swim/model/types.py src/sim_swim/sim/params.py src/sim_swim/sim/core.py src/sim_swim/analysis/flagella_count_behavior.py scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py scripts/02_phase2_analysis/build_flagella_count_behavior_dataset.py tests/test_model_builder.py tests/test_params.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py --sample-limit 1 --overwrite --progress-interval 5000`
  - `uv run python scripts/02_phase2_analysis/build_flagella_count_behavior_dataset.py --overwrite`
- docs:
  - `conf/phase2_analysis/flagella_count_behavior_dataset.yaml`
  - `docs/codex-runs/20260618_173743_phase2_72_flagella_count_dataset/review_result.json`

### P2-8-010: 特徴量分布の可視化を追加する

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/73`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/71`
- branch: `feature/phase2-73-feature-distributions`
- goal: #72 の `summary.csv` dataset から，べん毛本数ごとの特徴量分布，QC偏り，NaN，簡易スクリーニングを確認できる分析出力を生成する。
- tasks:
  - [x] `scripts/02_phase2_analysis/plot_flagella_count_behavior_distributions.py` を追加する。
  - [x] `--dataset-dir` または `--dataset-id` で dataset directory を指定できるようにする。
  - [x] category別の sample点付き分布図を `plots/distributions/` に出力する。
  - [x] `quality_class` と `use_for_analysis` の違いを確認できる図を出力する。
  - [x] `quality_summary.csv`，`nan_summary.csv`，`feature_summary_by_n_flagella.csv`，`feature_screening_summary.csv` を出力する。
  - [x] 既存分析出力はデフォルトで上書きせず，`--overwrite` 指定時のみ再生成する。
- acceptance criteria:
  - [x] dataset directory を指定して `summary.csv` を読み込める。
  - [x] `n_flagella` ごとの特徴量分布図が出力される。
  - [x] sample ごとの値が図上で確認できる。
  - [x] `quality_class` を区別した可視化ができる。
  - [x] `use_for_analysis == true` のみの可視化ができる。
  - [x] quality / diagnostics の本数別集計が出力される。
  - [x] feature ごとの NaN 数が `nan_summary.csv` に出力される。
  - [x] feature ごとの要約統計が `feature_summary_by_n_flagella.csv` に出力される。
  - [x] 差がありそうな特徴量候補が `feature_screening_summary.csv` に出力される。
- tests/checks:
  - `uv run pytest tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff check scripts/02_phase2_analysis tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff format --check scripts/02_phase2_analysis tests/test_flagella_count_behavior_dataset.py`
  - `uv run python scripts/02_phase2_analysis/plot_flagella_count_behavior_distributions.py --dataset-id fc_nf1_2_3_6_seed1_dur0p5 --overwrite`
- docs:
  - `scripts/README.md`
  - `docs/codex-runs/20260619_191237_phase2_73_feature_distributions/review_result.json`

## Completed support task: 動画出力・サンプリング整備

### P2-9-009: 3D/2D動画出力のレビュー向けサンプリングを整備する

- status: complete
- branch: `feature/phase2-9-output-sampling`
- background:
  - `dt_star=1.0e-4` では内部step数が多く、3D出力を全step保存するとファイル数・動画生成時間が大きくなる。
  - 2D側には `fps_out_2d` があるが、3D側にも同様の `fps_out_3d` が必要である。
- tasks:
  - [x] `output_sampling.fps_out_3d` を追加する。
  - [x] 3D render / frame保存 / manifest に `fps_out_3d` を反映する。
  - [x] 既存の `output_sampling.out_all_steps_3d` との優先順位を定義する。
  - [x] CLI overrideで `output_sampling.fps_out_3d=100` のように指定できることを確認する。
- acceptance criteria:
  - [x] `dt_star=1.0e-4`, `duration_s>=0.5` でも3D動画フレーム数を制御できる。
  - [x] manifestに3D出力サンプリング条件が記録される。
  - [x] 既存の全step保存挙動を必要時に維持できる。
- docs:
  - `docs/codex-runs/20260609_014607_phase2_fps_out_3d_sampling/review_result.json`

### P2-8-078: dataset raw sample replay render のfps指定を追加する

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/78`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/71`
- branch: `feature/phase2-78-render-fps-sampling`
- background:
  - `dt_star=1.0e-4`, `duration_s=0.5` 条件では archive state が約5000件になり，後出し3D描画を全stepで行うとファイル数と処理時間が大きくなる。
  - 既存の `output_sampling.out_all_steps_3d` / `fps_out_3d` / `fps_out_2d` を raw sample replay CLI から指定できる必要がある。
- tasks:
  - [x] `render_flagella_count_behavior_sample.py` に `--fps-out-3d` / `--fps-out-2d` を追加する。
  - [x] replay render のデフォルトを `out_all_steps_3d=false` とし，必要時のみ `--out-all-steps-3d` で全step描画できるようにする。
  - [x] `manifest.json` と `run.log` に effective な render sampling 条件を記録する。
  - [x] `conf/phase2_analysis/flagella_count_behavior_dataset.yaml` に標準 sampling 条件を明記する。
  - [x] 3D/2D mp4 は H.264 (`avc1` / `H264`) を優先し，encoder がない環境では `mp4v` にfallbackする。
- acceptance criteria:
  - [x] 既存 raw archive を再シミュレーションせず，3D/2D replay fps を指定できる。
  - [x] デフォルト実行で3D全step描画を避ける。
  - [x] sampling 条件が manifest から追跡できる。
  - [x] 実際に使われた動画codecが manifest から追跡できる。
- tests/checks:
  - `uv run pytest tests/test_flagella_count_behavior_dataset.py`
  - `uv run pytest tests/test_render_state_and_projection.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff check scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff check src/sim_swim/render scripts/01_simulate_swimming/01_simulate_swimming.py scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py tests/test_render_state_and_projection.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff format --check scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff format --check src/sim_swim/render scripts/01_simulate_swimming/01_simulate_swimming.py scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py tests/test_render_state_and_projection.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run python scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py --help`
  - `uv run python scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py --sample-dir outputs/phase2_analysis/flagella_count_behavior/runs/fc_nf1_2_3_6_seed1_dur0p5/samples/nf03_seed000 --output-dir /private/tmp/phase2_codec_probe_nf03 --fps-out-3d 25 --fps-out-2d 25`
  - `ffprobe -hide_banner /private/tmp/phase2_codec_probe_nf03/render/swim3d.mp4`
  - `ffprobe -hide_banner /private/tmp/phase2_codec_probe_nf03/render2d/projection.mp4`
- docs:
  - `scripts/README.md`
  - `docs/phase2/phase2_current.md`
  - `docs/codex-runs/20260619_220923_phase2_78_render_fps_sampling/review_result.json`

### P2-8-081: sweep実行時間の短縮

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/81`
- branch: `feature/phase2-76-seeded-initial-conditions`
- goal: `run_flagella_count_behavior_sweep.py` の後半 sample ほど実行時間が長く見える問題に対し，標準互換を維持しつつ軽量実行手段を追加する。
- result:
  - 実測ログでは `n_flagella` 増加に伴い 1 sample の計算時間が増え，標準順序が `1,2,3,6` のため後半ほど遅く見えていた。
  - `runner.sample_order=interleave_n_flagella` で seed 条件ごとに `n_flagella` を混ぜて実行できるようにした。
  - 2026-06-22後続判断により，保存段階で情報を消す `runner.step_summary_stride` / `runner.state_stride` と `conf/phase2_analysis/flagella_count_behavior_dataset_fast.yaml` は廃止した。Phase 2.8 raw sample は全step保存し，軽量化は可視化側samplingで行う。
- acceptance criteria:
  - [x] 既定挙動は既存 runner / dataset builder と互換である。
  - [x] manifest に runner 設定と sample timing が記録される。
  - [x] `interleave_n_flagella` で sample 順序を制御できる。
  - [x] `runner.step_summary_stride` / `runner.state_stride` は指定すると error になる。
- verification:
  - `uv run pytest tests/test_flagella_count_behavior_dataset.py tests/test_simulation.py tests/test_model_builder.py`
  - `uv run ruff check scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py src/sim_swim/sim/core.py tests/test_flagella_count_behavior_dataset.py tests/test_simulation.py`
  - `uv run ruff format --check scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py src/sim_swim/sim/core.py tests/test_flagella_count_behavior_dataset.py tests/test_simulation.py`
- docs:
  - `scripts/README.md`
  - `docs/phase2/phase2_current.md`
  - `docs/codex-runs/20260620_230903_phase2_81_sweep_runtime_shortening/review_result.json`
  - `docs/codex-runs/20260622_183135_phase2_remove_raw_stride_sampling/review_result.json`

### P2-8-084: 実験簡略化のための設定・再描画導線整理

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/84`
- branch: `feature/phase2-84-experiment-simplification`
- goal: 実験時の設定編集と動画確認を簡略化するため，標準config default，3D render表示，dataset一括再描画CLIを整理する。
- result:
  - `conf/sim_swim.yaml` の `motor.torque_Nm` default は，Phase 2.6 torque評価の第一候補，Phase 2.7代表条件，Phase 2.8 dataset条件と揃え，`2.5e-20` とした。`-1` sentinel の意味はコメントに残した。
  - `1.0e-4` は `time.dt_star` の実行時overrideとして使う値であり，モータートルクdefaultにはしない。
  - `flagella.placement_mode` と `flagella.initial_phase_mode` の取りうる値をconfigコメントへ明記した。
  - `render.save_frames_3d` / `render.save_frames_2d` の default を `false` にし，mp4 と final image は維持した。
  - 3D render に RUN/TUMBLE，時刻，実効トルク，`follow_camera_3d` を併記するようにした。
  - `render_flagella_count_behavior_sample.py --dataset-dir` で dataset 内の全raw sampleを `replays/<sample_id>/` へ一括再描画できるようにした。
- acceptance criteria:
  - [x] `conf/sim_swim.yaml` と parser fallback default が一致する。
  - [x] `-1` sentinel と明示的な `motor.torque_Nm` override は維持される。
  - [x] 3D render の表示内容を単体テストで確認できる。
  - [x] dataset directory から `dataset_manifest.json` / `run_manifest.json` を辿って全sampleを再描画できる。
  - [x] dataset一括再描画は既存 `replays/<sample_id>/` を置き換える。
- verification:
  - `uv run pytest tests/test_params.py tests/test_render_state_and_projection.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run python -c "import yaml; yaml.safe_load(open('conf/sim_swim.yaml', encoding='utf-8'))"`
  - `uv run ruff check src/sim_swim/sim/params.py src/sim_swim/render/render3d.py scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py tests/test_params.py tests/test_render_state_and_projection.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff format --check src/sim_swim/sim/params.py src/sim_swim/render/render3d.py scripts/02_phase2_analysis/render_flagella_count_behavior_sample.py tests/test_params.py tests/test_render_state_and_projection.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run python -m scripts.01_simulate_swimming time.duration_s=0.001 time.dt_star=1.0e-4 motor.torque_Nm=0 output.base_dir=/private/tmp/phase2_issue84_smoke`
- docs:
  - `scripts/README.md`
  - `docs/phase2/phase2_current.md`
  - `docs/codex-runs/20260622_124330_phase2_84_experiment_simplification/review_result.json`

### P2-8-084b: seed値によるべん毛初期配置のcenter-priority化

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/84#issuecomment-4764990562`
- branch: `feature/phase2-84-seeded-attach-placement`
- goal: `seeded_surface` の小さい `attach_seed` で付着点が1辺へ寄る問題を解消し，`n_flagella` ごとに中心三角形を優先的に埋める seed 範囲を明文化する。
- result:
  - `placement_mode=seeded_surface` では，`attach_seed` を付着点組合せ列のindexとして扱う。
  - 各 `n_flagella` の前半 seed は center layer の三角形を優先する。`n_flagella<=3` は center から必要数を選び，`n_flagella>=4` は center 3点を必ず含めて残りを前後隣接層から選ぶ。
  - center-priority seed 数は `n_flagella=0..9` で `1,3,3,1,6,15,20,15,6,1` とする。つまり `n_flagella=1..9` の center-priority `attach_seed` 範囲は `0..2`, `0..2`, `0`, `0..5`, `0..14`, `0..19`, `0..14`, `0..5`, `0` である。
  - center-priority 範囲を超えた `attach_seed` は，9候補全体から center-priority と重複しない制約なし配置へ進み，候補列全体を modulo で循環する。
- acceptance criteria:
  - [x] `n_flagella=0..9` の center-priority seed 数がテストで固定される。
  - [x] center-priority 範囲内では，`n_flagella<=3` はcenter layerのみ，`n_flagella>=4` はcenter三角形3点を必ず含む。
  - [x] center-priority 範囲の直後ではcenter優先制約が外れる。
  - [x] `phase_seed` だけを変えると付着点は変わらない既存保証を維持する。
- verification:
  - `uv run pytest tests/test_model_builder.py tests/test_params.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff check src/sim_swim/model/builder.py tests/test_model_builder.py`
  - `uv run ruff format --check src/sim_swim/model/builder.py tests/test_model_builder.py`
  - `uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py --dry-run --config conf/phase2_analysis/flagella_count_behavior_dataset_center_prefix.yaml --sample-limit 6`
- docs:
  - `docs/phase2/phase2_current.md`
  - `docs/phase2/phase2_tasks.md`
  - `docs/codex-runs/20260622_142503_phase2_84_seeded_attach_placement/review_result.json`

### P2-8-084c: center-priority前半seedのみdataset option

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/84#issuecomment-4764990562`
- branch: `feature/phase2-84-seeded-attach-placement`
- goal: `n_flagella` ごとに異なる center-priority 前半 `attach_seed` だけを使うdatasetを，手動列挙なしで生成できるようにする。
- result:
  - `sweep.attach_seed_mode=center_priority_prefix` を追加し，`build_conditions()` で `n_flagella` ごとの前半seedを自動展開する。
  - `attach_seed_mode` と `attach_seeds` の同時指定は曖昧なため error にする。
  - 専用config `conf/phase2_analysis/flagella_count_behavior_dataset_center_prefix.yaml` を追加した。`n_flagella=[1,2,3,6]`, `phase_seeds=[0]` で27 sampleを生成する。
  - 既存の明示 `attach_seeds` / `phase_seeds` と legacy `seeds` は互換維持する。
- acceptance criteria:
  - [x] `n_flagella=[1,2,3,6]` で前半 `attach_seed` が `0..2`, `0..2`, `0`, `0..19` に展開される。
  - [x] 専用configから27 sampleが生成され，sample idは既存 split seed 形式を使う。
  - [x] `phase_seed` は専用configでは `0` のみにできる。
  - [x] `attach_seed_mode` と `attach_seeds` の同時指定は error になる。
- verification:
  - `uv run pytest tests/test_flagella_count_behavior_dataset.py tests/test_model_builder.py`
  - `uv run ruff check scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py src/sim_swim/analysis/flagella_count_behavior.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run ruff format --check scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py src/sim_swim/analysis/flagella_count_behavior.py tests/test_flagella_count_behavior_dataset.py`
  - `uv run python scripts/02_phase2_analysis/run_flagella_count_behavior_sweep.py --dry-run --config conf/phase2_analysis/flagella_count_behavior_dataset_center_prefix.yaml --sample-limit 6`
- docs:
  - `conf/phase2_analysis/flagella_count_behavior_dataset_center_prefix.yaml`
  - `docs/phase2/phase2_current.md`
  - `docs/phase2/phase2_tasks.md`
  - `docs/codex-runs/20260622_175854_phase2_84_center_prefix_dataset/review_result.json`

### P2-8-082: hook過伸長対策の局所補強候補を比較する

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/82`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/10`
- branch: `feature/phase2-82-hook-overstretch`
- goal: `root_torque_segment_couples` 条件で body とべん毛根元が離れる hook過伸長に対し，論文モデルdefaultを変えずに局所補強候補を比較できるようにする。
- implementation notes:
  - `motor.local_attach_first_spring_scale` は body attach bead とべん毛第1ビーズの距離補強である。
  - `motor.local_attach_first_body_axis_angle_scale` は，`attach -> first` ベクトルが body長軸に対して90度を保つための補強である。これは `hook_triplets=(attach, first, second)` のなす角ではない。
  - `motor.local_first_second_spring_scale` はべん毛第1ビーズと第2ビーズの距離補強である。
  - 3つの新scaleはいずれもdefault `1.0` とし，標準configの挙動は変えない。非defaultは診断用の paper model extension として扱う。
- result:
  - `n_flagella=1`, `duration_s=0.02`, `time.dt_star=1.0e-4` の軽量代表 sweep で，body軸90度補強は `local_attach_first_vs_body_axis_err_deg` を baseline の `28.8650 deg` から `0.4656 deg` へ改善した。
  - 同条件で `hook_len_rel_err_max` は baseline の `0.4676` から `0.2210` へ下がり，`shape_pass_nonbody=True` を維持した。
  - 強い combined scale `4.0` は hook fail したため，補強scaleは標準defaultへ昇格せず，診断用 extension として残す。
  - `n_flagella=3`, `duration_s=0.5` のフル代表 sweep は既存 segment repulsion 計算が重く，ローカル実行では中断した。
  - 2026-06-25追補として，実行条件を `body-first-grid` と `first-second-grid` に分け，body-first 距離・body軸90度補正を先に評価し，必要な場合だけ第1-第2ビーズ距離補正を追加評価できるようにした。
  - `plot_phase2_82_hook_overstretch_heatmap.py` を追加し，破綻カテゴリ，shape pass/fail，hook距離誤差，body-first距離誤差，body軸90度誤差，第1-第2ビーズ距離誤差を heatmap 化できるようにした。
- acceptance criteria:
  - [x] 3つの新scaleが config override で個別に指定できる。
  - [x] body長軸90度補強が既存 hook三点角度拘束と混同されない形で実装・テストされる。
  - [x] Issue #82 用 sweep CLI が baseline と局所補強候補を比較し，`step_summary.csv` の hook距離・body軸角度・螺旋回転指標を summary に集約できる。
  - [x] 短時間代表条件で hook過伸長の改善有無を記録する。
  - [x] body-first 距離・body軸90度補正の grid sweep と heatmap 出力ができる。
  - [x] 第1-第2ビーズ距離補正追加時の sweep と heatmap 出力ができる。
- verification:
  - `uv run pytest tests/test_params.py tests/test_motor_forces.py tests/test_simulation.py`
  - `uv run ruff check src/sim_swim/dynamics src/sim_swim/sim scripts/01_simulate_swimming/run_phase2_82_hook_overstretch_sweep.py tests/test_params.py tests/test_motor_forces.py tests/test_simulation.py`
  - `uv run ruff format --check src/sim_swim/dynamics src/sim_swim/sim scripts/01_simulate_swimming/run_phase2_82_hook_overstretch_sweep.py tests/test_params.py tests/test_motor_forces.py tests/test_simulation.py`
  - `uv run python scripts/01_simulate_swimming/run_phase2_82_hook_overstretch_sweep.py --dry-run --sample-limit 3`
  - `uv run python scripts/01_simulate_swimming/run_phase2_82_hook_overstretch_sweep.py --duration-s 0.02 --n-flagella 1 --overwrite --progress-interval 5000`
  - `uv run pytest tests/test_phase2_82_hook_overstretch_sweep.py tests/test_phase2_82_hook_overstretch_heatmap.py`
- user-run commands:
  - Stage 1 body-first grid:
    `uv run python scripts/01_simulate_swimming/run_phase2_82_hook_overstretch_sweep.py --mode body-first-grid --duration-s 0.5 --dt-star 1.0e-4 --torque-nm 2.5e-20 --n-flagella 3 --attach-seed 0 --phase-seed 0 --attach-first-spring-scales 1,1.25,1.5,2,3 --body-axis-angle-scales 1,1.25,1.5,2,3 --output-dir outputs/phase2_82/body_first_grid --overwrite --progress-interval 5000`
  - Stage 1 heatmap:
    `uv run python scripts/01_simulate_swimming/plot_phase2_82_hook_overstretch_heatmap.py --summary-csv outputs/phase2_82/body_first_grid/phase2_82_hook_scale_sweep_summary.csv --mode body-first-grid --output-dir outputs/phase2_82/body_first_grid/plots`
  - Stage 2 first-second grid:
    `uv run python scripts/01_simulate_swimming/run_phase2_82_hook_overstretch_sweep.py --mode first-second-grid --duration-s 0.5 --dt-star 1.0e-4 --torque-nm 2.5e-20 --n-flagella 3 --attach-seed 0 --phase-seed 0 --fixed-attach-first-spring-scale 2 --fixed-body-axis-angle-scale 2 --first-second-spring-scales 1,1.25,1.5,2,3 --output-dir outputs/phase2_82/first_second_grid --overwrite --progress-interval 5000`
  - Stage 2 heatmap:
    `uv run python scripts/01_simulate_swimming/plot_phase2_82_hook_overstretch_heatmap.py --summary-csv outputs/phase2_82/first_second_grid/phase2_82_hook_scale_sweep_summary.csv --mode first-second-grid --output-dir outputs/phase2_82/first_second_grid/plots`
- docs:
  - `docs/phase2/phase2_current.md`
  - `docs/phase2/phase2_tasks.md`
  - `docs/codex-runs/20260625_205543_phase2_82_hook_overstretch/review_result.json`

## Phase 2.9: Tumble状態の段階実装

### P2-9-010: Tumble状態を段階実装する

- status: accepted
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/69`
- parent issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/10`
- branch: `feature/phase2-69-tumble-state`
- goal: 現行のRUN固定モデルを前提に、Tumble状態を設計・診断・motor reversal・polymorph切替・run-and-tumble評価へ段階分割して実装する。
- background:
  - 参照論文では、tumble時にモーター回転反転と normal / semicoiled / curly1 などのpolymorph切替が関係する。
  - 現行MVPでは `motor.enable_switching=false` によりRUN固定を基本にしている。
  - #65 でRUN固定条件の本数差と姿勢安定性を評価した後、Tumbleを導入することで方向転換やrun-and-tumble挙動を扱えるようにする。
- stages:
  - [ ] Stage 1: 現行run固定実装と論文tumble要素の差分を整理し、必要なstate, config, diagnostics, testsを設計する。
  - [ ] Stage 2: motor reversal の最小実装を追加し、短時間で形状・トルク伝搬・ログが壊れないことを検証する。
  - [ ] Stage 3: normal / semicoiled / curly1 のpolymorph state切替を導入し、平衡 bending / torsion 値と state transition log を検証する。
  - [ ] Stage 4: run-and-tumble挙動として、姿勢変化、方向転換、動画自然さを評価する。
- acceptance criteria:
  - [ ] 各stageでreview_resultを分け、Tumble全体の完了を一括で主張しない。
  - [ ] `motor.enable_switching=false` のRUN固定挙動を壊さない。
  - [ ] motor reversal とpolymorph切替の状態遷移が `step_summary.csv` または別CSVに記録される。
  - [ ] normal / semicoiled / curly1 の平衡角と論文モデルとの差分が文書化される。
  - [ ] 代表run-and-tumble条件について、姿勢変化または方向転換指標と目視レビュー要否が記録される。
