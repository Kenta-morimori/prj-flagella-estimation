# Phase 2.7: 複数べん毛の後方束化・非崩壊条件探索

## 目的

Phase 2.7 では、単一べん毛で確認した `material_twist_local_couple` を複数べん毛へ拡張し、螺旋形状・hook・body が崩壊しないまま後方束化する条件を探す。

最終的には、後方束化したべん毛束によって菌体がどの程度推進し、姿勢がどの程度揺れるかも確認する。単に「束化して見える」だけでなく、遊泳挙動として利用できるかを評価する。

## 基本方針

検証は次の順で行う。

1. 自然初期条件で多べん毛の形状安定トルク帯を求める。
2. hook角度または初期接線方向を調整し、全べん毛を菌体後方へ向けた `posterior_aligned` 条件を作る。
3. `posterior_aligned` 条件で安定トルク帯を再探索し、後方束化候補を分類する。
4. 後方束化候補条件で、菌体の推進量と姿勢変化を定量化する。

この順序にする理由は、最初から後方へ揃えた条件だけを見ると、自然初期条件での多べん毛耐性と、人工的に後方へ向けた初期条件の影響を切り分けにくいためである。

## 実行条件

- `motor.force_distribution=material_twist_local_couple` を主条件とする。
- `distributed_flagellum` は、torque がべん毛全体へ届いた場合の診断用比較modeとして残す。
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

## posterior aligned 条件

`posterior_aligned` は、自然な初期条件ではなく、短時間で後方束化を観察しやすくするための診断条件として扱う。

目標は、各べん毛の初期接線方向が菌体長軸の後方成分を持つことである。完全に同じ軸へ重ねるのではなく、付着点と螺旋位相の差は残す。これにより、後方束化しやすい条件を作りつつ、べん毛間干渉や部分束化も観察できるようにする。

## 分類

各条件は次のカテゴリに分類する。

- `posterior_bundle`: 多数のべん毛が後方側で近接し、束軸が菌体後方を向く。
- `partial_bundle`: 2本以上は束化するが、1本以上が独立する。
- `no_bundle`: 崩壊しないが束化しない。
- `collapse`: shape gate fail、fly-away、hook破綻、flagellum形状崩壊のいずれか。

1本のみ独立する状態は即FAILとはせず、`partial_bundle` として記録する。

## 指標

形状安定性:

- `shape_pass_nonbody`
- `first_fail_category_nonbody`
- `hook_len_rel_err_max`
- `flag_bond_rel_err_max`
- `flag_bend_err_max_deg`
- `flag_torsion_err_max_deg`
- `net_abs_flag_helix_spin_revolutions`

束化:

- 束中心軸への距離
- 束参加率
- 独立べん毛数
- flagella間距離
- bundle axis と body axis の角度
- bundle axis の後方成分

遊泳:

- body重心移動量
- body平均速度
- body axis の累積角度変化
- body angular velocity RMS
- body姿勢揺らぎRMS
- bundle axis と body axis の角度

## 完了条件

- `n_flagella=3`, `duration_s>=0.5`, `time.dt_star=1.0e-4` で shape gate PASS の代表条件がある。
- 自然初期条件と `posterior_aligned` 条件のそれぞれについて、torque と分類結果の表がある。
- `posterior_bundle` または `partial_bundle` の候補条件が1つ以上ある。
- 束化候補条件で、菌体推進量と姿勢変化の指標が出力される。
- 代表動画と確認観点を review_result に記録する。

## 2026-06-03 実装・診断結果

### 追加した実装

- `output_sampling.fps_out_3d` を追加し、`output_sampling.out_all_steps_3d=false` のときに3D動画を指定fpsで間引けるようにした。
- `flagella.initial_orientation_mode` を追加した。
  - `side_attach`: 既定値。既存の側面付着条件で、初期接線は菌体後方に対して約90度。
  - `posterior_aligned`: 診断用。hookは側面付着のまま、べん毛初期接線を菌体後方へ向ける。
- `step_summary.csv` に束化・遊泳候補指標を追加した。
- `scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py` を追加し、orientation mode / n_flagella / torque のsweep結果を `phase2_7_bundling_sweep_summary.csv` に集約できるようにした。

### 短時間スクリーニング

条件:

- `n_flagella=3`
- `time.dt_star=1.0e-4`
- `duration_s=0.02`
- `motor.force_distribution=material_twist_local_couple`
- `torque_Nm = 1.0e-20, 2.0e-20, 3.0e-20`
- `initial_orientation_mode = side_attach, posterior_aligned`

結果:

- 6条件すべてで `shape_pass_nonbody=True`。
- `side_attach` は束軸が菌体後方に対して約90度で、後方束化候補ではない。
- `posterior_aligned` は束軸自体は後方を向くが、`bundle_participation_ratio=0.0` のままで、短時間では束化しない。
- 出力: `outputs/phase2_7_bundling_sweep/2026-06-03_n3_dt1e4_dur0p02/phase2_7_bundling_sweep_summary.csv`

### 0.5秒代表条件

条件:

- `n_flagella=3`
- `initial_orientation_mode=posterior_aligned`
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

現行実装では、`posterior_aligned` によって全べん毛の束軸を菌体後方へ向けることはできた。しかし、0.5秒の長時間条件では低トルクでもhook長が徐々に崩れ、後方束化以前にhook gateで落ちる。

また、後方へ向けた条件でも `bundle_participation_ratio=0.0` であり、少なくとも今回の `n_flagella=3`, `0.5e-20..3.0e-20` 範囲では、べん毛同士が束中心軸へ近づく挙動は確認できていない。つまり、短時間で「後方を向いている」状態は作れるが、「後方束化している」状態までは到達していない。

### 残課題

- hook破綻を抑える条件探索が必要である。候補は `motor.local_spring_scale`, `motor.local_hook_scale`, hook spring / bend の扱い、posterior aligned 時のhook初期形状である。
- 束化しない原因を切り分ける必要がある。候補は、べん毛間の接触・反発・流体相互作用の強さ、初期配置の距離、posterior aligned が平行配置に近く束化駆動を作らないこと。
- Phase 2.7 は現時点では未完了であり、今回の結果は診断進捗として扱う。
