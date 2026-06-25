# Phase 2.7: 複数べん毛の後方束化・非崩壊条件探索

## 目的

Phase 2.7 では、単一べん毛で確認した `root_torque_segment_couples` を複数べん毛へ拡張し、螺旋形状・hook・body が崩壊しないまま複数べん毛の螺旋中心軸が安定的に揃う条件を探す。

本PRでは「束化」を、べん毛同士が近接して1本の束に見えることではなく、複数べん毛軸方向が安定的に揃うこととして定義する。菌体軸の揺れ、移動距離、推進速度などの遊泳特徴量は別PRで評価する。

## 2026-06-16 現在の作業方針

本作業は issue #58 `[Phase2] 後方束化による安定的な形状維持・遊泳検証` として、最新 `main` から派生した新ブランチ `feature/phase2-58-posterior-bundling-swim` で進める。

PR #55 は古い `main` から派生した診断用WIPとして扱い、必要な実装だけをこのブランチへ移植する。#55 の内容をそのまま継続しない理由は、PR #59 で確定した Phase 2.6 の `local_*_scale=1.0` baseline、`fps_out_3d`、トルク境界評価と重複する差分があるためである。

PR #63 で進めていた螺旋軸診断は、PRを分けずに本ブランチへ統合する。今後の主条件は、第1ビーズを含む `attach-first-second` の接線方向ではなく、第2ビーズ以降から推定した螺旋中心軸を菌体後方へ向ける `flagella.initial_helix_axis_from_rear_deg` とする。

今回の初期代表条件は次のように固定する。

- `motor.force_distribution=root_torque_segment_couples`
- `motor.torque_Nm=2.5e-20`
- `time.dt_star=1.0e-4`
- `duration_s>=0.5`
- `flagella.initial_helix_axis_from_rear_deg=0`
- `motor.local_*_scale=1.0`
- `n_flagella=3,6,9` を主対象とし、必要に応じて `n_flagella=1` を比較基準にする。

分類は二値判定にしない。`axis_aligned_stable`, `hook_wrapped_axis_aligned`, `axis_not_aligned`, `collapse` を区別する。1本だけ軸方向が外れる状態は plot 上でも外れとして確認できるようにし、`axis_not_aligned` として記録する。

今回の最優先は、後方または側方の初期条件で多べん毛が形状崩壊せず動き、複数べん毛軸方向が平均軸から `15 deg` 以内で安定的に揃うかを、`step_summary.csv`、sweep summary、時系列plotで再現可能にすることである。hook角度そのものは目的指標ではなく、hook長failは菌体巻き付き由来の候補として strict collapse とは分けて記録する。

## 基本方針

検証は次の順で行う。

1. 自然初期条件で多べん毛の形状安定トルク帯を求める。
2. 第2ビーズ以降の螺旋中心軸を菌体後方または側方へ向けた初期条件を作る。
3. 複数べん毛軸方向が安定的に揃うトルク帯を探索し、軸整列候補を分類する。
4. 軸整列候補条件の3D動画と時系列plotを確認する。

この順序にする理由は、最初から後方へ揃えた条件だけを見ると、自然初期条件での多べん毛耐性と、人工的に後方へ向けた初期条件の影響を切り分けにくいためである。

## 実行条件

- `motor.force_distribution=root_torque_segment_couples` を主条件とする。
- `root_torque_axis_projection` は、torque がべん毛全体へ届いた場合の軸投影比較modeとして残す。旧名 `distributed_flagellum` は deprecated alias。
- `time.dt_star=1.0e-4` はCLI overrideで指定し、デフォルト設定は変更しない。
- 最低検証時間は `duration_s=0.5` とする。
- 代表条件が見えたら `duration_s=1.0` 以上も検討する。
- 3D動画は `output_sampling.fps_out_3d` で間引き、長時間レビュー可能なフレーム数にする。

## torque sweep

最初の探索は `n_flagella=3` に固定し、single flagellum の代表条件 `2.0e-20 N m` を中心に粗く見る。

初期候補:

- `0.5e-20 N m`
- `1.0e-20 N m`
- `1.5e-20 N m`
- `2.0e-20 N m`
- `2.5e-20 N m`
- `3.0e-20 N m`

多べん毛では flagella 間干渉で不安定化する可能性があるため、上限は結果を見てから広げる。

## 後方螺旋中心軸条件

`flagella.initial_helix_axis_from_rear_deg` は、自然な初期条件ではなく、短時間で後方束化を観察しやすくするための診断条件として扱う。

目標は、各べん毛の螺旋本体の中心軸が菌体長軸の後方成分を持つことである。ここでいう螺旋中心軸は、各べん毛の第2ビーズ以降からPCAで推定する。第1ビーズはhook側の接続ビーズであり、菌体への巻き付きやhook length driftの影響を受けやすいため、螺旋本体の軸推定には含めない。

幾何の意図は、`attach -> first` を菌体長軸に対してほぼ垂直な外向き方向に保ち、螺旋本体の中心軸だけを後方へ向けることである。したがって、`hook.threshold_deg=90` は後方束化の目的変数ではなく、hookが潰れすぎないための制約である。

2026-06-16時点では、後方軸からの角度を連続値として扱う方針にする。

- `flagella.initial_helix_axis_from_rear_deg=0`: 第2ビーズ以降の螺旋中心軸を菌体後方軸と平行にする主診断条件。
- `flagella.initial_helix_axis_from_rear_deg=10`: ほぼ後方向きだが、完全平行ではなく少し外側へ開く比較条件。
- `flagella.initial_helix_axis_from_rear_deg=null`: 従来の `side_attach` 条件。

旧PR #60では `initial_tangent_vs_rear_deg` や `initial_flagellum_axis_from_rear_deg` という名前を検討したが、これらはhook近傍の局所接線と螺旋本体の中心軸を混同しやすい。そのため、未mergeの互換aliasは残さず、今後のドキュメント・実行コマンドでは `initial_helix_axis_from_rear_deg` を使う。

角度sweepの初期候補は `0, 10, 30, 60 deg` とする。まずは `0 deg` と代表トルク `2.5e-20 N m` を優先し、螺旋軸が本当に後方へ揃った状態で束化・部分束化が見られるかを確認する。

重要な確認点は、初期からべん毛同士が近すぎて repulsion dominated になっていないかである。そのため、初期stepから flagella-flagella minimum distance と close-contact pair count を記録する。

## 分類

各条件は次のカテゴリに分類する。

- `axis_aligned_stable`: 後半80%の全stepで、各べん毛軸の平均軸からの最大偏差が `15 deg` 以内で、hook以外のshape gateが保たれる。
- `hook_wrapped_axis_aligned`: 軸整列は `axis_aligned_stable` と同等だが、hook長failが出る。菌体への巻き付き候補として、strict collapseとは分ける。
- `axis_not_aligned`: hook以外のshapeは保つが、後半80%で軸整列が `15 deg` 以内に収まらない。1本だけ外れる場合もここに含める。
- `collapse`: finite、flag bond、flag bend、flag torsion など、hook長fail以外の非body shape gateが破綻する。

後方か側方かは主分類ではなく、副指標として `flag_helix_axis_vs_rear_angle_deg_*` と時系列plotで確認する。

## 指標

形状安定性:

- `shape_pass_nonbody`
- `first_fail_category_nonbody`
- `shape_pass_nonbody_strict`
- `first_fail_category_nonbody_strict`
- `shape_pass_nonbody_hook_len_relaxed`
- `first_fail_category_nonbody_hook_len_relaxed`
- `hook_len_strict_limit`
- `hook_len_relaxed_limit`
- `hook_len_rel_err_max`
- `flag_bond_rel_err_max`
- `flag_bend_err_max_deg`
- `flag_torsion_err_max_deg`
- `net_abs_flag_helix_spin_revolutions`
- `local_attach_first_vs_body_axis_angle_deg`
- `local_attach_first_vs_body_axis_err_deg`

軸整列:

- `flag_helix_axis_pair_angle_deg_mean`
- `flag_helix_axis_pair_angle_deg_max`
- `flag_helix_axis_mean_deviation_deg_max`
- `flag_helix_axis_alignment_order`
- `flag_helix_axis_vs_rear_angle_deg_*`
- `flag_helix_axis_angles_timeseries.png`

べん毛間反発:

- flagella-flagella bead-pair minimum distance
- flagella-flagella close-contact pair count
- flagella-flagella repulsion mean / max
- root/hook近傍のrepulsion mean / max
- 束中心軸距離とrepulsionの関係

近接・反発は補助指標として残すが、本PRの成功条件にはしない。

## 完了条件

- `n_flagella=3`, `duration_s>=0.5`, `time.dt_star=1.0e-4` で shape gate PASS の代表条件がある。
- 自然初期条件と後方螺旋中心軸条件のそれぞれについて、torque と軸整列分類結果の表がある。
- `axis_aligned_stable` または `hook_wrapped_axis_aligned` の候補条件が1つ以上ある。
- 軸整列候補条件で、べん毛軸角度の時系列plotが出力される。
- 代表動画と確認観点を review_result に記録する。

## 2026-06-16 現行実装

PR #60 には、PR #63で分離していた螺旋軸診断を統合する。

- `src/sim_swim/sim/helix_axis.py` を追加し、菌体後方軸と各べん毛の螺旋中心軸を推定する。
- 螺旋中心軸は、各べん毛の第2ビーズ以降からPCAで推定する。第1ビーズはhook接続部として扱い、軸推定には含めない。
- `flagella.initial_helix_axis_from_rear_deg` を追加する。`0` は螺旋中心軸を菌体後方へ揃え、`null` は従来のside-attach初期条件を使う。
- `render.show_flagella_helix_axis_3d=true` で、3D動画に推定螺旋軸を重ねる。
- `step_summary.csv` に `flag_helix_axis_vs_rear_angle_deg_*`, `flag_helix_axis_rearward_projection_min`, `flag_helix_axis_fit_r2_min`, `flag_helix_axis_degenerate_count` を追加する。
- `step_summary.csv` に `flag_helix_axis_pair_angle_deg_*`, `flag_helix_axis_mean_deviation_deg_max`, `flag_helix_axis_alignment_order` を追加する。
- `flag_helix_axis_diagnostics.csv` に、各べん毛ごとの軸方向・角度・fit品質を出力する。
- `run_phase2_7_bundling_sweep.py` は `--helix-axis-angles-deg` を受け取り、旧 `--tangent-angles-deg` は採用しない。
- `run_phase2_7_bundling_sweep.py` は各run directoryに `flag_helix_axis_angles_timeseries.png` を自動生成する。

この実装で確認する主な問いは、「複数べん毛の螺旋本体の軸が、後半80%で平均軸から `15 deg` 以内に安定して揃うか」である。`flag_helix_axis_pair_angle_deg_max` は補助指標として残す。後方か側方かは `flag_helix_axis_vs_rear_angle_deg_*` とplot上の `0 deg` / `90 deg` 参照線で区別する。

代表コマンド:

```bash
uv run python scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py \
  --n-flagella 3,6,9 \
  --torques 2.5e-20 \
  --duration-s 0.5 \
  --dt-star 1.0e-4 \
  --helix-axis-angles-deg 0 \
  --output-dir outputs/phase2_7_issue58_helix_axis_bundling_v1 \
  output_sampling.out_all_steps_3d=false \
  output_sampling.fps_out_3d=25 \
  render.show_flagella_helix_axis_3d=true \
  render.save_frames_3d=false
```

現時点では、後方螺旋中心軸条件での定量sweepとユーザー定性レビューは未完了である。そのため、Phase 2.7 全体のreview statusはまだ `FAIL / diagnostic` とする。

### 2026-06-16 軸整列定義での代表条件

新しい軸整列定義では、`n_flagella=3`, `motor.torque_Nm=2.5e-20`, `duration_s=0.5`, `time.dt_star=1.0e-4`, `flagella.initial_helix_axis_from_rear_deg=0` を代表条件として評価した。

結果は `hook_wrapped_axis_aligned` だった。strict / hook-relaxed のどちらでも hook長failは残るが、flag bond / bend / torsion は大きく破綻せず、後半80%の軸整列は `15 deg` 閾値を満たした。

| 指標 | 値 |
| --- | ---: |
| `phase27_class` | `hook_wrapped_axis_aligned` |
| `phase27_axis_alignment_stable` | `True` |
| `phase27_axis_alignment_stable_fraction` | `1.0` |
| `flag_helix_axis_mean_deviation_deg_max` | `11.8841` |
| `flag_helix_axis_pair_angle_deg_max` | `20.5403` |
| `flag_helix_axis_vs_rear_angle_deg_max` | `11.9216` |
| `hook_len_rel_err_max` | `2.6763` |
| `flag_bond_rel_err_max` | `0.1869` |
| `flag_bend_err_max_deg` | `3.7965` |
| `flag_torsion_err_max_deg` | `9.9368` |

出力:

- `outputs/phase2_7_axis_alignment_v1/phase2_7_bundling_sweep_summary.csv`
- `outputs/phase2_7_axis_alignment_v1/helix_axis_angle_0deg/n_3/torque_2.50e-20/flag_helix_axis_angles_timeseries.png`
- `outputs/phase2_7_axis_alignment_v1/helix_axis_angle_0deg/n_3/torque_2.50e-20/step_summary.csv`
- `outputs/phase2_7_axis_alignment_v1/helix_axis_angle_0deg/n_3/torque_2.50e-20/flag_helix_axis_diagnostics.csv`

解釈:

- 複数べん毛軸は後方へ揃い、平均軸からの最大偏差は `15 deg` 以内に収まった。
- 明示的な近接束ではないが、本PRの新定義では軸整列成功候補である。
- hook長failは菌体への巻き付き候補として別カテゴリにしたため、flag形状崩壊とは分けて扱う。
- 菌体軸の揺れ、移動距離、推進速度は別PRで評価する。

## 2026-06-09 main派生ブランチでの旧接線制御診断結果

この節は、旧 `initial_flagellum_axis_from_rear_deg` / `initial_orientation_mode=posterior_aligned` に基づく診断結果である。hook近傍の局所接線制御と螺旋本体の中心軸制御が混ざるため、以後の主条件にはしない。ただし、hook length drift と `no_bundle_drive` の失敗モードを確認した履歴として残す。

### 追加・移植した実装

本ブランチでは、PR #55 のWIP実装をそのまま継続せず、PR #59 merge後の最新 `main` から必要部分だけを移植した。評価条件は Phase 2.6 の結論に合わせ、`motor.local_*_scale=1.0` を維持する。

- 旧実装では `flagella.initial_orientation_mode` を追加した。
  - `side_attach`: 従来どおり、側面付着方向へ初期べん毛を伸ばす。
  - `posterior_aligned`: 診断用。hook付着は側面のまま、hook近傍の初期軸を菌体後方へ向ける。
- 旧実装では `flagella.initial_flagellum_axis_from_rear_deg` を追加した。指定時は `initial_orientation_mode` より優先し、0度を菌体後方軸、90度を従来の側面付着相当として扱った。
- `step_summary.csv` に束化・遊泳指標を追加した。
- `scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py` を追加し、`n_flagella`, `torque_Nm`, 初期軸条件のsweep結果を `phase2_7_bundling_sweep_summary.csv` へ集約できるようにした。
- sweep用途として `Simulator.run(..., stop_on_shape_fail=True)` を追加し、shape gate fail後の不要な計算を止められるようにした。通常の `scripts.01_simulate_swimming` 実行では既定値 `False` のままである。
- `local_attach_first_vs_body_axis_angle_deg` と `local_attach_first_vs_body_axis_err_deg` を追加し、第1ビーズが菌体長軸に対して垂直外向きに保たれているかを診断できるようにした。

べん毛間距離・反発指標は、長時間の `n_flagella=9` sweepを現実的な時間で回すため、厳密なsegment-segment最近接ではなく、異なるべん毛に属するビーズ間pairの高速近似として記録する。列名は実装内容に合わせて `flag_flag_bead_pair_dist_*` とする。

### 今回の代表条件

- `motor.force_distribution=material_twist_local_couple`
- `motor.local_*_scale=1.0`
- 旧 `flagella.initial_flagellum_axis_from_rear_deg=10`
- `time.dt_star=1.0e-4`
- `duration_s=0.5`
- `n_flagella=3,6,9`

### `torque_Nm=2.5e-20` の結果

ユーザー指定の初期代表トルクでは、全本数で0.5秒まで到達する前にhook gateで停止した。

| n_flagella | final_t_s | 判定 | first_fail_category_nonbody | bundle_participation_ratio | flag_flag_close_pair_count | body_displacement_um |
| ---: | ---: | --- | --- | ---: | ---: | ---: |
| 3 | 0.0713 | collapse | hook | 0.0 | 0 | 0.0565 |
| 6 | 0.0331 | collapse | hook | 0.0 | 0 | 0.0542 |
| 9 | 0.0513 | collapse | hook | 0.0 | 0 | 0.1308 |

出力:

- `outputs/phase2_7_issue58_posterior_bundling_sweep_v3/phase2_7_bundling_sweep_summary.csv`

### 低トルク sweep の結果

`0.5e-20..2.0e-20` でも、全条件で `collapse / hook` となった。トルクを上げるほど fail時刻は早くなり、束化指標はすべて0のままだった。

| n_flagella | torque_Nm | final_t_s | 判定 | first_fail_category_nonbody | bundle_participation_ratio | flag_flag_close_pair_count |
| ---: | ---: | ---: | --- | --- | ---: | ---: |
| 3 | 0.5e-20 | 0.4430 | collapse | hook | 0.0 | 0 |
| 3 | 1.0e-20 | 0.1974 | collapse | hook | 0.0 | 0 |
| 3 | 1.5e-20 | 0.1207 | collapse | hook | 0.0 | 0 |
| 3 | 2.0e-20 | 0.0847 | collapse | hook | 0.0 | 0 |
| 6 | 0.5e-20 | 0.2112 | collapse | hook | 0.0 | 0 |
| 6 | 1.0e-20 | 0.0919 | collapse | hook | 0.0 | 0 |
| 6 | 1.5e-20 | 0.0570 | collapse | hook | 0.0 | 0 |
| 6 | 2.0e-20 | 0.0409 | collapse | hook | 0.0 | 0 |
| 9 | 0.5e-20 | 0.2739 | collapse | hook | 0.0 | 0 |
| 9 | 1.0e-20 | 0.1360 | collapse | hook | 0.0 | 0 |
| 9 | 1.5e-20 | 0.0873 | collapse | hook | 0.0 | 0 |
| 9 | 2.0e-20 | 0.0625 | collapse | hook | 0.0 | 0 |

出力:

- `outputs/phase2_7_issue58_posterior_bundling_low_torque_v1/phase2_7_bundling_sweep_summary.csv`

### さらに低い `torque_Nm=2.5e-21` の結果

トルクを1桁下げると、`n_flagella=3` と `n_flagella=9` は0.5秒のshape gateを通過した。ただし、どちらも `no_bundle` であり、束化参加率・close-contact・repulsionは0のままだった。`n_flagella=6` は0.478秒でhook gateをわずかに超えた。

| n_flagella | final_t_s | 判定 | shape_pass_nonbody | first_fail_category_nonbody | net_abs_flag_helix_spin_revolutions | bundle_participation_ratio | flag_flag_close_pair_count | body_displacement_um |
| ---: | ---: | --- | --- | --- | ---: | ---: | ---: | ---: |
| 3 | 0.4999 | no_bundle | True | none | 0.1489 | 0.0 | 0 | 0.0217 |
| 6 | 0.4781 | collapse | False | hook | 0.2259 | 0.0 | 0 | 0.0321 |
| 9 | 0.4999 | no_bundle | True | none | 0.2386 | 0.0 | 0 | 0.0627 |

出力:

- `outputs/phase2_7_issue58_posterior_bundling_very_low_torque_v1/phase2_7_bundling_sweep_summary.csv`

### 現時点の原因整理

今回の条件では、後方束化成功条件は見つからなかった。主因は次の2つである。

1. 代表トルク `2.5e-20` と低トルク `0.5e-20..2.0e-20` では、束化や接触の前にhook length gateが破綻する。これは `hook_drift` である。
2. `2.5e-21` まで下げると一部条件で形状は保てるが、回転量・推進量が小さく、べん毛間の接近も起きない。これは `no_bundle_drive` である。

また、すべての実行で `flag_flag_close_pair_count=0` かつ repulsion proxy が0であった。したがって、今回の失敗原因は「複数べん毛が絡まって大きな反発で壊れた」ではなく、「代表トルクではhookが先に崩れ、形状を保てる低トルクでは束化を駆動できない」と解釈する。

### 次に詰めるべき点

- 後方10度初期条件でhookが伸びる原因を、hook長・hook曲げ・attach-first springの整合として調査する。単純な `local_*_scale` 補強は Phase 2.6 のbaseline方針と矛盾するため、まずは局所形状と力の内訳を評価する。
- 現行のspring-spring repulsionは反発であり、べん毛同士を束へ引き寄せる相互作用ではない。回転する螺旋同士の流体相互作用や、束化を促す相互作用が現行モデルに不足していないかを別途評価する。
- `2.5e-21` は形状維持の診断条件としては有用だが、束化・遊泳代表条件としてはトルク不足の可能性が高い。

## 2026-06-15 以降の目的別タスク整理

今後は、後方束化・非崩壊・遊泳を1つの条件探索に混ぜず、次の目的ごとに分けて進める。

1. hook drift の原因解明と安定化
   - 後方10度初期条件で、hook長・hook角度・attach-first springがなぜ崩れるかを局所診断する。
   - まずは `local_*_scale` 補強ではなく、幾何と力の内訳から原因を確定する。
   - 最小成功条件は `torque=2.5e-20`, `time.dt_star=1.0e-4`, `duration_s=0.5`, `n_flagella=3` で hook gate PASS とする。
2. 束化駆動の有無を評価する
   - hookが保てる条件で、べん毛同士が近づかない原因を調べる。
   - 現行のspring-spring repulsionは反発であり、束へ引き寄せる駆動ではないため、流体相互作用または束化を促す相互作用の不足を評価する。
   - `posterior_bundle` または `partial_bundle` 候補が1条件以上あることを目標にする。
3. 多本数へ拡張する
   - `n_flagella=3` の代表条件を確定してから `n_flagella=6,9` へ広げる。
   - 本数ごとの `collapse / no_bundle / partial_bundle / posterior_bundle` 表を作る。
4. 遊泳挙動を評価する
   - 束化候補条件に限って、body displacement、body speed、body axis wobble、bundle-body angleを評価する。
   - `duration_s>=0.5` の定量結果を最小条件とし、代表条件は必要に応じて `1-5s` の動画レビューへ進める。
5. 2D/ML用データとしての妥当性を評価する
   - 3Dで安定・束化・遊泳の代表条件が見えた後に、2D疑似顕微鏡動画として使えるかを別途確認する。

## 2026-06-15 第1ビーズ外向き診断結果

`hook.threshold_deg=90` は目的指標ではなく、後方束化を見るためにhookが潰れすぎないようにする制約である。そこで、第1ビーズが菌体長軸に対して垂直外向きに保たれているかを `local_attach_first_vs_body_axis_angle_deg` と `local_attach_first_vs_body_axis_err_deg` で確認した。

初期geometryでは、旧 `initial_flagellum_axis_from_rear_deg=10` の条件で、hook近傍の局所配置は期待通りになっていた。

| 指標 | 値 |
| --- | ---: |
| `initial_tangent_vs_rear_direction_angle_deg` | `10.0 deg` |
| `attach_first_vs_body_axis_angle_deg` | `90.0 deg` |
| `initial_hook_angle_deg` | `100.0 deg` |
| `initial_geometry_pass` | `True` |

短時間 `duration_s=0.02`, `torque_Nm=2.5e-20`, `n_flagella=3` ではshape gateは通ったが、すでにhook driftが進んでいた。

| t_s | shape_pass_nonbody | local_attach_first_vs_body_axis_angle_deg | local_attach_first_vs_body_axis_err_deg | hook_len_rel_err_max | bundle_participation_ratio | flag_flag_close_pair_count |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| `0.0000` | `True` | `90.0001` | `0.0746` | `0.0063` | `0.0` | `0` |
| `0.0199` | `True` | `98.9621` | `11.0711` | `0.5009` | `0.0` | `0` |

0.5秒代表条件 `duration_s=0.5`, `torque_Nm=2.5e-20`, `n_flagella=3` は、従来通り `0.0713 s` で `collapse / hook` となった。fail時点では `hook_angle_err_max_deg=0.5446` と小さいため、hook角度が鋭く潰れたことではなく、`local_attach_first_rel_err` / `hook_len_rel_err_max` が `1.0` を越えたことが直接のfail理由である。

| torque_Nm | final_t_s | 判定 | hook_len_rel_err_max | local_attach_first_vs_body_axis_angle_deg | local_attach_first_vs_body_axis_err_deg | hook_angle_err_max_deg | bundle_participation_ratio | flag_flag_close_pair_count |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `0.5e-20` | `0.4430` | collapse / hook | `1.0002` | `111.9687` | `31.8109` | `0.5400` | `0.0` | `0` |
| `1.0e-20` | `0.1974` | collapse / hook | `1.0002` | `111.9874` | `31.5932` | `0.5602` | `0.0` | `0` |
| `1.5e-20` | `0.1207` | collapse / hook | `1.0000` | `111.9549` | `31.0070` | `0.5753` | `0.0` | `0` |
| `2.0e-20` | `0.0847` | collapse / hook | `1.0008` | `111.9587` | `30.4827` | `0.5802` | `0.0` | `0` |
| `2.5e-20` | `0.0713` | collapse / hook | `1.0017` | `110.7503` | `27.9341` | `0.5446` | `0.0` | `0` |

出力:

- `outputs/phase2_7_issue58_posterior_axis_metrics_short_v1/phase2_7_bundling_sweep_summary.csv`
- `outputs/phase2_7_issue58_posterior_axis_metrics_representative_v1/phase2_7_bundling_sweep_summary.csv`
- `outputs/phase2_7_issue58_posterior_axis_metrics_torque_sweep_v1/phase2_7_bundling_sweep_summary.csv`

現時点の結論:

- 初期配置は「第1ビーズは菌体長軸に垂直外向き、べん毛軸は後方10度」を満たしている。
- しかし、トルク入力後に `attach -> first` が外向き90度から110度前後へ傾き、同時にhook長が基準長から大きく伸びる。
- `hook_angle_err_max_deg` は小さいため、今回のhook failは「hook角度が90度から崩れる」問題ではなく、「第1ビーズ外向き配置を保てず、hook長が伸びる」問題である。
- 全条件で `bundle_participation_ratio=0.0`, `flag_flag_close_pair_count=0` であり、後方束化や接触の前にhook長が破綻している。

## 2026-06-15 hook長relaxed gate方針

ユーザー目視確認では、今回のhook挙動は菌体に巻き付いている状態として許容し得る。そのため、現行のstrict gateは残したまま、hook長だけを緩めた評価列を追加する。

- `shape_pass_nonbody` / `shape_pass_nonbody_strict`: 既存互換のstrict判定。
- `first_fail_category_nonbody` / `first_fail_category_nonbody_strict`: strict判定でのfail理由。
- `shape_pass_nonbody_hook_len_relaxed`: hook長のみ `2.0` まで許容した判定。hook角度、flag bond、flag bend、flag torsionはstrictと同じ。
- `first_fail_category_nonbody_hook_len_relaxed`: hook長relaxed判定でのfail理由。
- `hook_len_strict_limit`: `1.0`。
- `hook_len_relaxed_limit`: `2.0`。
- `phase27_class`: strict判定による分類。
- `phase27_class_hook_len_relaxed`: hook長relaxed判定による分類。

この分け方により、「厳密にはhook長が伸びているが、hook巻き付きとして許容した場合に束化・遊泳評価へ進めるか」を区別して記録する。

### hook長relaxed gateの評価結果

`n_flagella=3`, 旧 `flagella.initial_flagellum_axis_from_rear_deg=10`, `time.dt_star=1.0e-4` で、strict gate と hook長relaxed gate の差分を評価した。

| 条件 | final_t_s | strict分類 | relaxed分類 | hook_len_rel_err_max | bundle_participation_ratio | flag_flag_close_pair_count | body_displacement_um |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: |
| `torque=2.5e-20`, strict停止 | `0.0713` | `collapse / hook` | `no_bundle` | `1.0017` | `0.0` | `0` | `0.0565` |
| `torque=2.5e-20`, 0.5秒継続 | `0.4999` | `collapse / hook` | `collapse / hook` | `2.3783` | `0.0` | `0` | `0.2682` |
| `torque=5.0e-21`, 0.5秒継続 | `0.4999` | `collapse / hook` | `no_bundle` | `1.1296` | `0.0` | `0` | `0.0469` |

代表トルク `2.5e-20 N m` では、strict停止時点 `t=0.0713 s` で `hook_len_rel_err_max=1.0017` となり、strictでは `collapse / hook` だが hook長relaxedでは `no_bundle` だった。この時点では hook角度、flagellum bond/bend/torsion は大きく崩れておらず、strict gate が hook長のわずかな超過を拾っていた。

一方、`--no-stop-on-shape-fail` で同条件を `0.5 s` まで継続すると、`t=0.3007 s` 付近で `hook_len_rel_err_max=2.0658` となり、hook長relaxed gate も fail した。最終時刻 `t=0.4999 s` では `hook_len_rel_err_max=2.3783` で、`phase27_class_hook_len_relaxed=collapse` だった。つまり、代表トルクでは「hook wrappingとして許容できる軽微な伸び」を越えて、hook長が継続的に伸びる。

低トルク `5.0e-21 N m` では、`0.5 s` まで継続しても `hook_len_rel_err_max=1.1296` に留まり、hook長relaxedでは `no_bundle` を維持した。ただし `bundle_participation_ratio=0.0`, `flag_flag_close_pair_count=0`, `body_displacement_um=0.0469` であり、後方束化や十分な遊泳駆動は確認できなかった。

strict停止のtorque sweepでは、`0.5e-20..2.5e-20 N m` の全条件が strictでは `collapse / hook`、hook長relaxedでは `no_bundle` になった。これは、fail直後の分類としてはhook長relaxedが過剰な早期停止を分離できることを示すが、代表トルクで0.5秒の安定条件を満たすことは示さない。

出力:

- `outputs/phase2_7_issue58_hook_len_relaxed_representative_v1/phase2_7_bundling_sweep_summary.csv`
- `outputs/phase2_7_issue58_hook_len_relaxed_representative_full_v1/phase2_7_bundling_sweep_summary.csv`
- `outputs/phase2_7_issue58_hook_len_relaxed_torque_sweep_v1/phase2_7_bundling_sweep_summary.csv`
- `outputs/phase2_7_issue58_hook_len_relaxed_low_torque_full_v1/phase2_7_bundling_sweep_summary.csv`

現時点では、hook長relaxed gate は診断として有効だが、後方束化成功条件はまだ見つかっていない。代表トルクでは relaxed gate でも長時間維持できず、形状を維持できる低トルクでは束化駆動が不足している。

## 2026-06-03 旧orientation mode実装・診断結果

この節は、旧 `initial_orientation_mode` に基づく診断履歴である。現行PRではこのconfig keyを採用せず、螺旋中心軸制御 `flagella.initial_helix_axis_from_rear_deg` に置き換える。

### 当時追加した実装

- `output_sampling.fps_out_3d` を追加し、`output_sampling.out_all_steps_3d=false` のときに3D動画を指定fpsで間引けるようにした。
- 旧 `flagella.initial_orientation_mode` を追加した。
  - `side_attach`: 既定値。既存の側面付着条件で、初期べん毛軸は菌体後方に対して約90度。
  - `posterior_aligned`: 診断用。hookは側面付着のまま、べん毛初期軸を菌体後方へ向ける。
- `step_summary.csv` に束化・遊泳候補指標を追加した。
- `scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py` を追加し、旧orientation mode / n_flagella / torque のsweep結果を `phase2_7_bundling_sweep_summary.csv` に集約できるようにした。

### 短時間スクリーニング

条件:

- `n_flagella=3`
- `time.dt_star=1.0e-4`
- `duration_s=0.02`
- `motor.force_distribution=material_twist_local_couple`
- `torque_Nm = 1.0e-20, 2.0e-20, 3.0e-20`
- 旧 `initial_orientation_mode = side_attach, posterior_aligned`

結果:

- 6条件すべてで `shape_pass_nonbody=True`。
- `side_attach` は束軸が菌体後方に対して約90度で、後方束化候補ではない。
- 旧 `posterior_aligned` は束軸自体は後方を向くが、`bundle_participation_ratio=0.0` のままで、短時間では束化しない。
- 出力: `outputs/phase2_7_bundling_sweep/2026-06-03_n3_dt1e4_dur0p02/phase2_7_bundling_sweep_summary.csv`

### 0.5秒代表条件

条件:

- `n_flagella=3`
- 旧 `initial_orientation_mode=posterior_aligned`
- `time.dt_star=1.0e-4`
- `duration_s=0.5`

結果:

| torque_Nm | local_spring_scale | shape_pass_nonbody | first_fail_category_nonbody | net_abs_flag_helix_spin_revolutions | bundle_axis_vs_rear_angle_deg | bundle_participation_ratio | body_displacement_um | 判定 |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | --- |
| `2.0e-20` | `1.0` | `False` | `hook` | `1.4087` | `18.18` | `0.0` | `0.1952` | collapse |
| `1.0e-20` | `1.0` | `False` | `hook` | `0.6685` | `7.77` | `0.0` | `0.0932` | collapse |
| `0.5e-20` | `1.0` | `False` | `hook` | `0.2573` | `2.31` | `0.0` | `0.0784` | collapse |
| `0.5e-20` | `1.2` | `False` | `hook` | `0.2628` | `2.24` | `0.0` | `0.0789` | collapse |

`motor.local_spring_scale=1.2` はhook破綻を遅らせたが、0.5秒最終stepでは `hook_len_rel_err_max=1.0713` となり、現行gateの上限 `1.0` を超えた。

### 現時点の解釈

旧実装では、`posterior_aligned` によって全べん毛の束軸を菌体後方へ向けることはできた。しかし、0.5秒の長時間条件では低トルクでもhook長が徐々に崩れ、後方束化以前にhook gateで落ちる。

また、後方へ向けた条件でも `bundle_participation_ratio=0.0` であり、少なくとも今回の `n_flagella=3`, `0.5e-20..3.0e-20` 範囲では、べん毛同士が束中心軸へ近づく挙動は確認できていない。つまり、短時間で「後方を向いている」状態は作れるが、「後方束化している」状態までは到達していない。

### 残課題

- hook破綻を抑える条件探索が必要である。候補は `motor.local_spring_scale`, `motor.local_hook_scale`, hook spring / bend の扱い、posterior aligned 時のhook初期形状である。
- 束化しない原因を切り分ける必要がある。候補は、べん毛間の接触・反発・流体相互作用の強さ、初期配置の距離、posterior aligned が平行配置に近く束化駆動を作らないこと。
- Phase 2.7 は現時点では未完了であり、今回の結果は診断進捗として扱う。

## 2026-06-04 追加方針

次の旧調査では、`posterior_aligned` の0度固定ではなく、菌体後方軸からの角度を振った。

優先順:

1. 旧 `initial_flagellum_axis_from_rear_deg=10`
2. 旧 `initial_flagellum_axis_from_rear_deg=30`
3. 旧 `initial_flagellum_axis_from_rear_deg=0`
4. 旧 `initial_flagellum_axis_from_rear_deg=60`
5. 旧 `initial_flagellum_axis_from_rear_deg=90`

実行は全grid長時間ではなく、短時間screeningと代表長時間検証に分ける。

1. `duration_s=0.02` で `0, 10, 30, 60, 90 deg` と低トルク候補をscreeningする。
2. `shape_pass_nonbody=True` かつhook指標が良い角度を選ぶ。
3. 選んだ角度だけ `duration_s=0.5` へ伸ばす。
4. 0.5 sで shape gate PASSした条件だけ、束化・遊泳指標を解釈する。

失敗時は、失敗原因を次のどれかに分類する。

- `hook_drift`: hook length / local attach-first が時間とともに伸びる。
- `repulsion_dominated`: flagella-flagella距離が小さく、repulsion force が上がる。
- `no_bundle_drive`: 後方を向くがflagella間距離が縮まらず、repulsionも小さい。
- `flag_shape_fail`: flagellum bond / bend / torsion が先に崩れる。
- `body_motion_unstable`: shapeは通るがbody姿勢変化が大きすぎる。

## 2026-06-04 旧接線角度sweep・反発診断結果

### 当時追加した実装

- 旧 `flagella.initial_flagellum_axis_from_rear_deg` を追加した。指定時は `initial_orientation_mode` より優先し、0度を菌体後方軸、90度を従来 `side_attach` 相当として扱った。
- `scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py` に旧 `--tangent-angles-deg` を追加した。
- `step_summary.csv` と sweep summary に、べん毛同士のbead-pair距離、close-contact数、flagella-flagella repulsion force、basal近傍repulsion forceを追加した。

### 短時間screening

条件:

- `n_flagella=3`
- `time.dt_star=1.0e-4`
- `duration_s=0.02`
- `motor.force_distribution=material_twist_local_couple`
- `motor.local_spring_scale=1.2`
- 旧 `initial_flagellum_axis_from_rear_deg = 0, 10, 30, 60, 90`
- `torque_Nm = 0.5e-20, 1.0e-20`

結果:

- 10条件すべてで `shape_pass_nonbody=True`。
- すべて `phase27_class=no_bundle` で、`bundle_participation_ratio=0.0`。
- すべて `flag_flag_close_pair_count=0`, `flag_flag_repulsion_force_max_N=0.0`。
- つまり、短時間では「初期から近すぎて反発で壊れる」状態ではなく、「べん毛同士が接触・束化する距離まで近づいていない」状態だった。
- 出力: `outputs/phase2_7_bundling_sweep/2026-06-04_n3_angle_screen_dur0p02/phase2_7_bundling_sweep_summary.csv`

### 0.5秒代表条件

条件:

- `n_flagella=3`
- `time.dt_star=1.0e-4`
- `duration_s=0.5`
- `motor.force_distribution=material_twist_local_couple`
- `motor.local_spring_scale=1.2`
- `torque_Nm=0.5e-20`

結果:

| 旧 initial_flagellum_axis_from_rear_deg | shape_pass_nonbody | first_fail_category_nonbody | net_abs_flag_helix_spin_revolutions | hook_len_rel_err_max | bundle_axis_vs_rear_angle_deg | bundle_participation_ratio | flag_flag_close_pair_count | 判定 |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| `0` | `False` | `hook` | `0.2628` | `1.0713` | `2.24` | `0.0` | `0` | collapse |
| `10` | `True` | `none` | `0.2927` | `0.9982` | `2.66` | `0.0` | `0` | no_bundle |
| `30` | `False` | `hook` | `0.3748` | `1.0948` | `39.04` | `0.0` | `0` | collapse |

10度条件では0.5秒のshape gateをぎりぎり通過した。しかし `bundle_participation_ratio=0.0` で、flagella間距離も縮まらず、flagella-flagella repulsionは発生しなかった。

出力:

- `outputs/phase2_7_bundling_sweep/2026-06-04_n3_angle_representative_dur0p5/phase2_7_bundling_sweep_summary.csv`

### 10度条件でのtorque上昇

旧 `initial_flagellum_axis_from_rear_deg=10`, `duration_s=0.5`, `motor.local_spring_scale=1.2` で torque を上げると、束化より先にhookが破綻した。

| torque_Nm | shape_pass_nonbody | first_fail_category_nonbody | first fail t_s | net_abs_flag_helix_spin_revolutions | hook_len_rel_err_max | bundle_participation_ratio | flag_flag_close_pair_count |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `1.0e-20` | `False` | `hook` | `0.2247` | `0.7206` | `1.5328` | `0.0` | `0` |
| `2.0e-20` | `False` | `hook` | `0.0967` | `1.4328` | `1.8936` | `0.0` | `0` |

`motor.local_spring_scale=2.0` まで補強すると `1.0e-20` はshape gateを通過したが、束化は起きなかった。

| torque_Nm | local_spring_scale | shape_pass_nonbody | first_fail_category_nonbody | net_abs_flag_helix_spin_revolutions | hook_len_rel_err_max | bundle_participation_ratio | flag_flag_close_pair_count | 判定 |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | --- |
| `1.0e-20` | `2.0` | `True` | `none` | `0.6894` | `0.9908` | `0.0` | `0` | no_bundle |
| `2.0e-20` | `2.0` | `False` | `hook` | `1.3843` | `1.4057` | `0.0` | `0` | collapse |

出力:

- `outputs/phase2_7_bundling_sweep/2026-06-04_n3_angle10_torque_escalation_dur0p5/phase2_7_bundling_sweep_summary.csv`
- `outputs/phase2_7_bundling_sweep/2026-06-04_n3_angle10_spring2_dur0p5/phase2_7_bundling_sweep_summary.csv`

### 現時点の原因整理

今回の範囲では、後方束化成功条件は見つからなかった。主因は次の2つに分けられる。

1. 形状が保てる低トルク条件では、べん毛同士が近づかず、flagella-flagella repulsionも発生しない。これは `no_bundle_drive` であり、初期べん毛軸を後方へ向けるだけでは束化を駆動できていない。
2. トルクを上げると螺旋回転量は増えるが、束化や接触の前にhook length gateが破綻する。これは `hook_drift` である。

したがって、現行の `material_twist_local_couple` + 後方寄り初期べん毛軸 + 現行spring-spring repulsion だけでは、`duration_s=0.5` で「崩壊せず後方束化する」条件は未達である。

### 次に詰めるべき点

- 後方束化を作る物理要素として、回転する螺旋同士の流体相互作用が現行実装で十分に表現されているかを確認する。
- hook drift を抑える方法を、単なる局所spring補強ではなく、hookの平衡長・曲げ拘束・付着直後の初期形状との整合として見直す。
- `distributed_flagellum` または診断用modeで多本べん毛を回した場合に束化傾向が出るかを比較し、問題がtorque伝達不足か、相互作用モデル不足かを切り分ける。
- 束化判定は `bundle_participation_ratio` だけでなく、flagella間距離の時間変化とclose-contact発生率を併用する。
