# Phase 2 Current

このファイルは，Codex CLI が Phase 2 作業の最初に読む短い現在地ドキュメントである。

原則として 100〜150 行以内に保つ。詳細な実験履歴，長い失敗条件，過去 run の詳細はここに書かず，必要な参照先だけを置く。

## Source of truth priority

* Repository-wide rules: `AGENTS.md`
* Current Phase 2 entry point: `docs/phase2/phase2_current.md`
* Accepted task status: `docs/phase2/phase2_tasks.md`
* Overall project map: `docs/PROJECT_PLAN.md`
* Historical context: Git history, issue/PR history, and Codex run records

## Context reading policy

* 長い文書を開く前に `rg -n` で関連箇所を絞る。
* 長い Markdown，logs，CSV，`work_log.md` はデフォルトで全文表示しない。
* 過去 run は `review_result.json` を先に読み，必要な場合だけ `work_log.md` を読む。
* `outputs/` の大きな成果物は，文書・summary・review_result で不足するときだけ読む。

## Current Phase 2 status

Phase 2 は，3D物理シミュレーションと2D擬似顕微鏡動画生成を，ML教師データに使える形へ整えるフェーズである。

完了済みの大枠:

* Phase 2.1〜2.5: 初期幾何，body-only，body+hook，single flagellum の短時間安定条件を整理済み。
* Phase 2.6: `root_torque_segment_couples` により，単一べん毛の螺旋形状維持と net 回転を両立する代表条件を確認済み。Issue #53 / #54 は完了。旧名 `material_twist_local_couple` は deprecated alias として受け付ける。
* Phase 2.7: 複数べん毛の束化検証は，近接束ではなく螺旋中心軸の安定整列を成功定義として完了。代表条件は `hook_wrapped_axis_aligned`。
* 3D出力の間引き用 `output_sampling.fps_out_3d` は導入済み。
* Phase 2 CLI の新規コマンド例は `config=...`，`time.duration_s=...`，`time.dt_star=...` などの `KEY=VALUE` 形式を第一表記にする。`--config` / `--duration-s` / `--fps-out` などの option 形式は legacy compatibility としてのみ残す。

現在の主対象:

* Phase 2.8 / Issue #65: 親Issue #71 のRUN固定べん毛数差分析に向けて，特徴量カテゴリとdataset出力方針を定義済み。
* Phase 2.8 / Issue #72: 親Issue #71 のRUN固定べん毛数差分析に向けて，複数条件実行スクリプトとdataset builderを追加済み。
  raw sample には `trajectory.csv` と `state_archive.npz` を残し，後から 3D / 2D render を再生成できる。
  Issue #78 により，後出し render CLI はデフォルトで3D全step描画を避け，`--fps-out-3d` / `--fps-out-2d` でfpsを指定できる。
  CLI は `scripts/02_phase2_analysis/`，設定は `conf/phase2_analysis/` に置く。
* Phase 2.8 / Issue #76: seed依存の初期条件ばらつきとして，付着点選択 `attach_seed` と初期helix phase `phase_seed` を分離した。標準datasetは `attach_seeds=[0,1,2]` と `phase_seeds=[0,1,2]` の直積を使う `fc_nf1_2_3_6_as3_ps3_dur0p5`。
* Phase 2.8 / Issue #73: `summary.csv` から本数別の特徴量分布plot，QC plot，要約統計，NaN集計，簡易スクリーニングを生成するCLIを追加済み。標準datasetの分析出力は `plots/distributions/`，`plots/qc/`，`analysis/` に生成できる。
* Phase 2.8 / Issue #81: `run_flagella_count_behavior_sweep.py` の保存間引き軽量化は廃止した。Phase 2.8 raw sample は `step_summary.csv`，`trajectory.csv`，`state_archive.npz` を全step保存し，軽量化は replay render / visualization 側の sampling で行う。情報を消さない `runner.sample_order=interleave_n_flagella` は維持する。
* Phase 2 / Issue #84: 実験簡略化のため，render frame保存defaultをOFFにした。`conf/sim_swim.yaml` の `motor.torque_Nm` default は，Phase 2.6 torque評価の第一候補，Phase 2.7代表条件，Phase 2.8 dataset条件と揃え，`2.5e-20` とする。`time.dt_star` default は Phase 2 標準の内部積分刻みとして `1.0e-4` とし，`time.dt_s=1.0e-3` は出力・記録間隔として扱う。`time.dt_star: null` を明示した場合だけ，内部積分刻みは従来互換として `time.dt_s/tau_s` に戻る。通常の3D render default は `output_sampling.out_all_steps_3d=false` とし，`fps_out_3d` で間引く。全内部stepを3D描画したい診断runだけ `output_sampling.out_all_steps_3d=true` を明示する。3D render は RUN/TUMBLE・時刻・実効トルク・`follow_camera_3d` を併記し，dataset directory 指定で全raw sampleを `replays/<sample_id>/` へ一括再描画できる。残タスクの `seeded_surface` 付着点配置では，`attach_seed` 前半を `n_flagella` ごとの center-triangle priority 配置にし，`n_flagella=1..9` の前半seed数を `3,3,1,6,15,20,15,6,1` として明文化した。前半seedのみのdataset用に `sweep.attach_seed_mode=center_priority_prefix` と `conf/phase2_analysis/flagella_count_behavior_dataset_center_prefix.yaml` を追加し，`phase_seeds=[0]` で27 sampleを生成できる。
* Phase 2 / Issue #88: `motor.force_distribution` の正式名を `triplet`, `root_torque_segment_couples`, `root_torque_axis_projection` に整理した。旧名 `material_twist_local_couple` と `distributed_flagellum` は deprecated alias として新名へ正規化する。診断用 probe mode `axial_torque_flux_probe` / `local_twist_transmission_probe` は削除した。
* Phase 2 / Issue #82: `root_torque_segment_couples` 条件での hook過伸長対策として，局所補強候補を比較できるようにした。`motor.local_attach_first_body_axis_angle_scale` の90度補強は，`hook_triplets=(attach, first, second)` の角度ではなく，`attach -> first` ベクトルが body長軸に対して垂直に生えることを指す。default `1.0` では標準挙動を変えず，非defaultを診断用 extension として扱う。body-first距離・body軸90度補正，第1-第2ビーズ距離補正，per-flag hook診断，first fail / 全期間最大hook event指標，`motor.local_attach_frame_position_scale` / `motor.local_attach_frame_tangent_scale` と `attach-frame-grid` sweep / heatmap に対応済み。2026-06-27時点の代表評価では，Stage A の `fp=3, ft=1.5` が baseline `hook_len_rel_err_max=3.2865` を `0.0145` まで抑え，後方条件0.5 sの定性評価でも hook根元挙動は問題なし。ただし first fail は `flag` へ移り，Stage B で `local_first_second_spring_scale=1..3` を足しても `flag_bond_rel_err_max` は `1.1875..1.1961` に残った。dt sweep では `dt_star=1.0e-4, 5.0e-5, 2.5e-5` で破綻種別が解消しないため，時間刻みではなく拘束・力の釣り合い不足として扱う。長時間3D定性評価 `frame_fp3_ft1p5_fs1p5` では first fail `flag`, `hook_len_rel_err_max=0.0157`, `flag_bond_rel_err_max=1.0006`, final `flag_bond_rel_err_max=2.0500` となり，hookではなく flagellum bond 過伸長が残課題である。反作用トルクは消しておらず body 側へ入っているが，attach-frame補強で body-root 間の相対運動が抑えられ，一体回転に近く見える。Issue #94 では flag bond 最大伸長の `flag_id` / bead pair / per-flag 最大値を `step_summary.csv` と #82 sweep summary に出し，sweep summary では local bead pair も出す。さらに proximal local bond `0-1..4-5` の rel err を per-flag で `step_summary.csv` と sweep summary の final / first fail / max event に出す。`fp=3, ft=1.5` 固定の短時間 `fs=1..3` 比較では `max_flag_bond_rel_err` はほぼ同等で明確な改善は見えず，0.6 s `fs=1..3` 比較でも全条件の first fail は `t=0.4363 s`, `flag_id=1`, global bead `29-30`, local bead `3-4` に固定された。`fs=1.25` は max event が flag 0 local `0-1` へ飛び，悪化候補として扱う。first-second や second-third ではないため，`fs` 増強は採用しない。2026-06-29 の torque切り分けでは `fs=1.0, fp=3, ft=1.5` で `1.0e-20..2.0e-20` は0.6 s first failなし，`2.5e-20` は非有限座標に起因すると見られる SVD crash で summary化できなかった。attach-frame強度切り分けでは `ft=1` が local `0-1` 破綻を誘発しやすく不採用，`ft=1.5` は維持候補とする。`fp=3, ft=1.5` は first fail が `0.4363 s` と最遅だが final/max `flag_bond_rel_err_max=1.4964`，`fp=2, ft=1.5` は first fail `0.4217 s` で final/max `1.3315` まで下がるため，attach-frame診断代表は `fs=1.0, ft=1.5` 固定で `fp=2` と `fp=3` を併記候補にする。SVD crash guard 実装後は `fp=3, ft=1.5, torque=2.5e-20, duration_s=0.6` でも Python 例外ではなく summary CSV まで出力できる。ただし first fail は `t=0.0009 s`, `hook`，max `flag_bond_rel_err` は極大値となるため，物理条件としては破綻扱いを継続する。body-flagella 貫通検証は Issue #93 側へ分離する。Issue #82 配下の残タスクは，`#94=flag bond 過伸長の診断・抑制`, `#97=torque distribution 比較と採用判断`, `#100=user-facing sweep/heatmap 導線の generic 化` として分担する。#100 は `hook_overstretch` 起源の実験基盤整理であり，#97 の比較再開前に先行実施する sub-issue とする。
* Phase 2 / Issue #94 follow-up: 2026-06-30時点の追加切り分けでは，`run_phase2_82_hook_overstretch_sweep.py` が `flagella.initial_helix_axis_from_rear_deg=0` を内部overrideするため，数値sweepは後方束化初期条件で評価している。`af=1, axis=1, fs=1, fp=2, ft=1.5` 固定では `torque=1.9e-20` だけが 0.6 s の自動 shape gate を通った。`fp=3, ft=1.5` では `torque=1.9e-20` と `2.0e-20` が shape gate を通り，`2.0e-20` は `max_flag_bond_rel_err=0.9063` で比較候補中もっとも余裕があった。merge判定用の後方動画 `outputs/phase2_94/posterior_merge_review/qual_videos/fp3_ft1p5_torque2p0_dur0p6/2026-06-30/170248` は user定性評価で概ね良好と判断されたため，PR #95 の採用候補は `fp=3, ft=1.5, fs=1, torque=2.0e-20` とする。自動指標上は後方軸整列が時間とともに弱まるため，べん毛が螺旋軸中心に回っているかの厳密評価も含め，トルク分散・回転伝達の深掘りは Issue #97 `[Phase2] hookとべん毛形状を保つトルク分散方法を見直す` で扱う。
* Phase 2 / Issue #97: 後方条件での `body / hook / flagella / attach-bead` 回転安定性を主目的に，トルク分散比較を整理した。`motor.torque_distribution_profile` を導入し，旧 `motor.torque_segment_weight_profile` は deprecated alias とした。旧 `local_twist_activity` 系 profile は `diffusive`, `diffusive_sqrt`, `diffusive_floor_0p2`, `diffusive_floor_0p4` に改名し，`root_torque_axis_projection` でも同じ profile 名を使えるようにした。後方主評価 Stage E では `root_torque_segment_couples + diffusive` だけが 0.6 s を shape pass したため，比較面を `root_torque_segment_couples / root_torque_axis_projection` × `diffusive / uniform` の 2x2 に縮小した。Issue #100 で現役導線は `shape_stability_grid` sweep / heatmap へ移し，replay utility `render_shape_stability_grid_replay.py` は `state_archive.npz` / `trajectory.csv` を使った再シミュレーションなしの定性確認に対応した。2026-07-07 の最終評価では，後方 1.0 s `outputs/phase2_97/stage_g_posterior_distribution_grid_fp3_ft1p5_torque2p0_dur1p0/summary.csv` は4条件すべて `flag` fail，側方 0.6 s `outputs/phase2_97/stage_h_lateral_distribution_grid_fp3_ft1p5_torque2p0_dur0p6/summary.csv` は4条件すべて shape pass した。`root_torque_segment_couples` では `diffusive` が `uniform` より first fail を遅らせる一方，側方 user定性評価 `stage_h_lateral_grid_qual_fp3_ft1p5_torque2p0_dur0p6` ではどの2x2条件も螺旋軸中心の回転ではなく菌体を含む剛体回転に見えた。従って #97 では現行2x2候補を採用せず，剛体回転の主因は torque profile より attach-frame補強による body-root 一体化 / basal自由度不足として扱い，診断を Issue #103 として #82 配下へ分離した。
* Phase 2 / Issue #103: attach-frame補強による剛体回転を診断するため，`body_roll_phase_deg` と `axis_center_body_relative_phase_deg` を追加し，`shape_stability_grid` summary で body roll と螺旋軸中心 spin を分離できるようにした。`motor.local_attach_frame_tangent_mode=basal_bearing` は，`attach -> first` 軸まわり方位角を拘束せず，`first -> second` の軸方向成分と垂直半径だけを保つ #103 診断用 extension である。#103 の比較とユーザー定性評価により，現行 `local_attach_frame_tangent_scale` は「接線方向だけ」を保つのではなく，body attach frame で見た `first -> second` ベクトル全体を初期 target に戻すため，root の軸まわり相対 spin も抑えやすく，側方・後方の両条件で剛体回転様の見え方を招くと整理した。言い換えると，attach-frame tangent を強くすると root tangent が body 座標系へ縛られ，べん毛だけが軸中心に回るより先に菌体ごと一緒に回る見え方になりやすい。続く position-only 比較 `conf/phase2_sweeps/basal_freedom_position_only_sweep.yaml` では，側方 0.6 s は `no_frame` がもっとも自然な軸中心回転に見え，`fp>=1.25` は hook を大きく改善する一方で body-relative spin / body-roll ratio を `103.5 -> 約48.6` まで落とした。後方 1.0 s では `no_frame` が `final_shape_pass_nonbody=True` でも `t=0.0449 s` に `hook` first fail しており，無破綻ではない。これに対し `fp1p25` は `t=0.3295 s` まで `flag` fail を遅らせ，`max_hook_len_rel_err=0.1207`, `body_roll_net_abs_revolutions=0.0135`, `axis_center_to_body_roll_ratio_mean=151.2` で最良の折衷となった。`fp1p5` は `first_fail_t_s=0.3498 s` とわずかに遅いが body roll が増え，`fp2+` は `flag_bond_rel_err_max` 悪化で不利だった。したがって default は tangent を外し，`motor.local_attach_frame_position_scale=1.25` の最小補強だけを採る。採用判断は `final_shape_pass_nonbody` 単独ではなく `first_fail_category_nonbody` / `first_fail_t_s` を優先して読む。さらに #94 既存結果どおり `local_first_second_spring_scale` は主要 failure である proximal flag bond 過伸長を改善しておらず，現時点では優先する補強候補ではない。`basal_freedom_diagnostic.yaml` は #103 の診断再現用として保持し，追加の長時間 sweep / replay render はユーザー実行待ちである。
* Phase 2 / Issue #96: `scripts/01_simulate_swimming` の sweep / heatmap 導線を整理した。現役CLIは `run_sweep.py config=conf/phase2_sweeps/<profile>.yaml` と `plot_heatmap.py config=conf/phase2_sweeps/<profile>.yaml` に統一し，issue番号付きの個別scriptは用途別 module under `src/sim_swim/analysis/` へ移した。新しい sweep summary の標準名は `summary.csv`。`shape_stability_grid.yaml` / `shape_stability_heatmap.yaml` を canonical profile とし，`hook_overstretch.yaml` / `hook_overstretch_heatmap.yaml` は historical alias として残す。`conf/phase2_sweeps/*.yaml` は `metadata` に role / canonical / 推奨対応先を持ち，wrapper から `list_canonical_profiles=true` と `describe_profile=true` で確認できる。`scripts/README.md` は初見ユーザー向けの現行操作説明だけに更新した。未使用だった `src/flagella_estimation` package は削除し，wheel 配布対象は `sim_swim` に切り替えた。
* Phase 2 / Issue #100: `sweep` と `multi-run / replay` の導線を整理した。`run_sweep.py` は `conf/phase2_sweeps/` の task-specific diagnostic sweep 専用とし，user-facing な複数条件実行は `run_multi_run.py config=conf/phase2_multi_run/<profile>.yaml` に寄せる。`conf/phase2_multi_run/` の profile は run / plot / replay で 1 枚共用し，`output.timestamp_subdir=false` の profile では `output.base_dir` を固定 run root として `plot_heatmap.py config=...` が `plots/` を，`render_shape_stability_grid_replay.py config=...` が `replay/` を作れる。さらに `shape_stability_grid` sweep 側も `condition_order` / `base_config` / `axes` / `condition_label` を含む共通 campaign 契約へ寄せ，`render_shape_stability_grid_replay.py` は sweep / multi-run の両方を同じ reader で読めるようにした。標準 profile は `latest_model_torque_shape_stability` とし，最新モデルの torque 複数条件比較を 3 コマンドで回せる。tracked local wrapper は不要になったため削除した。`plot_heatmap.py` / replay の log / manifest 強化は別PRで扱う。なお `shape_stability_grid` の heatmap 実装は sweep 診断用 mode ごとの集計ロジックを残しており，generic 化しているのは wrapper / config / summary contract と multi-run 側の plot 導線である。
* Phase 2.8 / Issue #71: 次は #73 の可視化結果を読み，seed数を増やすか，評価時間を伸ばすか，3D特徴量の初期分離性評価へ進むかを判断する。
* Phase 2.9 / Issue #69: Tumble状態を段階実装する。まず設計・診断，次にmotor reversal，polymorph切替，run-and-tumble評価へ分ける。

## Current blockers

* P2-7代表条件は `hook_wrapped_axis_aligned` であり，hook length fail は残る。#71 のRUN本数差評価では前提リスクとして扱う。
* #82/#94 の attach-frame 補強は hook過伸長を 0.5 s 代表条件で抑えるが，長時間では flagellum bond 過伸長が支配的になる。#103 の整理後は `ft=1.5` を維持候補から外し，default は `motor.local_attach_frame_position_scale=1.25`, `motor.local_attach_frame_tangent_scale=1.0`, `motor.local_first_second_spring_scale=1.0` とする。Issue #94 後は破綻 flag / bead pair を CSV で追え，非有限座標・SVD失敗も fail-safe に記録できるが，根本的な flag bond 安定化は未解決である。出力整理では `step_summary.csv`，`trajectory.csv`，`state_archive.npz`，`manifest.json` を再現用として保持し，途中停止run・重複run・報告に使わない再生成可能な動画/frameだけを削除対象にする。
* #97 の 2x2 torque distribution 候補は，後方 1.0 s で全条件 `flag` fail，側方 0.6 s で剛体回転様の user定性評価となったため採用しない。#103 の初期診断も `ft` / `basal_bearing` の継続採用には至らなかった。position-only 比較の結果，現時点の best posterior compromise は `fp1p25`，次点は `fp1p5`，lateral の自然回転参照点は `no_frame` である。次の blocker は `fp1p25` / `fp1p5` のどちらを #82 へ返す代表候補にするかの最終判断であり，採用判断では `final_shape_pass_nonbody` 単独ではなく `first_fail_category_nonbody` / `first_fail_t_s` を優先して読む。
* #71 では #72/#73/#76 のdataset・plot出力を使い，body displacement, speed, body axis angle change, angular velocity, wobble RMS, flag helix axis alignment を本数別に比較する必要がある。標準datasetは付着点seedとphase seedを分けるため，seed要因を層別して解釈する。
* #69 は一括実装にせず，設計・診断から段階化する必要がある。

## Default read path

Phase 2 task の開始時は，次の順で必要な範囲だけ読む。

1. `AGENTS.md`
2. この `docs/phase2/phase2_current.md`
3. 対象 issue / PR 本文とコメント
4. `docs/phase2/phase2_tasks.md` の該当 task 周辺
5. 該当する `docs/phase2/phase2_*.md`
6. 関連する `docs/codex-runs/*/review_result.json`

`docs/PROJECT_PLAN.md` は，全体地図や phase-level context が必要な場合だけ読む。

## Key references

* Accepted Phase 2 tasks: `docs/phase2/phase2_tasks.md`
* Overall project map: `docs/PROJECT_PLAN.md`
* Phase 2.6 single flagellum gate: `docs/phase2/phase2_6_helix_retention_gate.md`
* Phase 2.6 torque model evaluation: `docs/phase2/phase2_6_torque_transmission_model_evaluation.md`
* Phase 2.7 bundling plan: `docs/phase2/phase2_7_bundling_stability_plan.md`
* Phase 2.7 helix axis diagnostics: `docs/phase2/phase2_7_flag_helix_axis_diagnostics.md`
* Phase 2.8 flagella count feature definitions: `docs/phase2/phase2_8_flagella_count_feature_definitions.md`
* Root torque segment-couple ADR: `docs/adr/0004_phase2_material_frame_twist_transmission.md`
* Codex workflow details: `docs/codex/codex_workflow.md`

## Recent run anchors

* `docs/codex-runs/20260616_144322_phase2_7_axis_alignment_definition/review_result.json`
* `docs/codex-runs/20260616_211555_phase2_7_axis_plot_readability/review_result.json`
* `docs/codex-runs/20260616_121312_phase2_7_issue58_helix_axis_integration/review_result.json`

## Usually do not read first

* `work_log.md`: read only after the matching `review_result.json` is insufficient.
* Large CSVs and generated outputs under `outputs/`: inspect only when summary docs cannot answer the question.

## Completion update rule

Phase 2 task completion should update only the documents needed by the result:

* Always: `docs/codex-runs/<run-id>/review_result.json`
* Usually: `docs/phase2/phase2_current.md`
* When accepted task status changes: `docs/phase2/phase2_tasks.md`
* When project-level phase status changes: `docs/PROJECT_PLAN.md`
* When a modeling or workflow decision is significant: `docs/adr/`
* When Codex workflow details change: `docs/codex/codex_workflow.md`

Do not mark a task complete unless the relevant review result is `PASS`.
