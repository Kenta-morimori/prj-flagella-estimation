# Phase 2.7: べん毛螺旋軸の後方向き評価

## 目的

後方束化条件を探索する前提として、各べん毛の螺旋軸が本当に菌体後方へ向いているかを定量・定性の両方で確認できるようにする。

## 螺旋軸の定義

本タスクでは、べん毛軸を `flagella_indices[f][1:]`、つまり各べん毛の第2ビーズ以降から推定する螺旋中心軸と定義する。

`flagella_indices[f][0]` は hook側の第1ビーズであり、菌体との接続や hook length drift の影響を強く受ける。そのため、螺旋本体の中心軸評価には含めない。

軸推定は第2ビーズ以降のPCA第1主成分で行い、`第2ビーズ -> 先端ビーズ` の方向と同じ向きになるよう符号を揃える。この軸と菌体後方方向の角度を `flag_helix_axis_vs_rear_angle_deg` とする。

## 出力

`step_summary.csv` には、各stepの全べん毛を集約した列を出力する。

- `flag_helix_axis_vs_rear_angle_deg_min`
- `flag_helix_axis_vs_rear_angle_deg_mean`
- `flag_helix_axis_vs_rear_angle_deg_max`
- `flag_helix_axis_rearward_projection_min`
- `flag_helix_axis_fit_r2_min`
- `flag_helix_axis_degenerate_count`
- `flag_flag_helix_bead_dist_min_um`
- `flag_flag_helix_close_pair_count`
- `flag_helix_bundle_radius_mean_um`
- `flag_helix_bundle_radius_max_um`

各べん毛ごとの詳細は `flag_helix_axis_diagnostics.csv` に出力する。

- `step`, `t_s`, `flag_id`
- `flag_helix_axis_vs_rear_angle_deg`
- `flag_helix_axis_rearward_projection`
- `flag_helix_axis_fit_r2`
- `axis_origin_x_um`, `axis_origin_y_um`, `axis_origin_z_um`
- `axis_dir_x`, `axis_dir_y`, `axis_dir_z`

## 3D可視化

`render.show_flagella_helix_axis_3d=true` のとき、3D動画に第2ビーズ以降から推定した螺旋中心軸を重ねる。デフォルトは `false` とし、通常の動画出力の見た目は変えない。

hookは従来どおり黄色で描画する。hook length drift は数値だけでは許容可否を決めず、螺旋軸 overlay と併せて、菌体への巻き付きとして許容できるか、非物理的な伸びとして扱うべきかを目視レビューで確認する。

## 後方螺旋軸初期条件

`flagella.initial_helix_axis_from_rear_deg` を指定すると、各べん毛の第2ビーズ以降から推定される螺旋中心軸を菌体後方から指定角度だけ傾ける。`0` は全べん毛の螺旋軸を同一の菌体後方方向へ揃える条件である。`null` は従来の side-attach 条件を維持する。

この設定は第1ビーズを菌体側hookとして残したまま、螺旋本体の軸だけを後方へ揃える診断条件である。そのため、hook長やhook近傍の形状破綻は別途確認する必要がある。

## 確認観点

- `flag_helix_axis_vs_rear_angle_deg_*` が意図した後方向き条件を反映しているか。
- 第1ビーズの位置ずれや hook伸長が、螺旋軸推定に混入していないか。
- `flag_helix_axis_fit_r2_min` が極端に低い条件では、螺旋軸の解釈が不安定でないか。
- `flag_flag_helix_bead_dist_min_um` と `flag_helix_bundle_radius_*` が減少するかを見て、fail判定時でも束化らしい接近が起きていないかを確認する。
- `hook_len_rel_err_max` が大きい条件では、3D動画で hook の巻き付き・伸長を定性評価する。

## 2026-06-16 後方螺旋軸条件での定量probe

`flagella.initial_helix_axis_from_rear_deg=0`, `time.dt_star=1.0e-4`, `duration_s=0.5` で、全べん毛の第2ビーズ以降から推定した螺旋軸を初期状態で菌体後方へ揃えた。

判定では `shape_pass_nonbody` が `False` になっても打ち切らず、最終stepまで `flag_flag_helix_bead_dist_min_um`, `flag_flag_helix_close_pair_count`, `flag_helix_bundle_radius_mean_um` を確認した。これは、hook fail が出ても目視では束化している場合があり得るためである。

出力:

- `outputs/2026-06-16/012334_phase2_7_posterior_axis_bundling_probe/summary_recomputed.csv`
- `outputs/2026-06-16/014710_phase2_7_posterior_axis_bundling_probe_n6/summary_recomputed.csv`

結果:

| n_flagella | torque_Nm | shape_pass_nonbody_final | first_fail | min flag-flag dist [um] | max close pairs | classification |
|---:|---:|---|---|---:|---:|---|
| 3 | 0 | True | none | 0.736 | 0 | no_bundle |
| 3 | 5.0e-21 | True | none | 0.655 | 0 | no_bundle |
| 3 | 1.0e-20 | False | hook at 0.430 s | 0.665 | 0 | fail_no_bundle |
| 3 | 2.5e-20 | False | hook at 0.177 s | 0.417 | 0 | fail_no_bundle |
| 6 | 5.0e-21 | True | none | 0.268 | 9 | bundle_candidate |
| 6 | 1.0e-20 | False | hook at 0.248 s | 0.224 | 9 | fail_bundle_candidate |
| 6 | 2.5e-20 | False | hook at 0.104 s | 0.200 | 6 | fail_bundle_candidate |

暫定結論:

- 後方螺旋軸初期条件は機能している。初期の `flag_helix_axis_vs_rear_angle_deg_max` は n=3/n=6 の全条件で 0.03 deg 未満だった。
- n=3 では 0.5 s の範囲で束化候補は出なかった。高トルクでは hook fail が先行し、fail 後も close pair は 0 のままだった。
- n=6 では `5.0e-21 N m` が shape pass のまま close pair を持つ束化候補になった。
- n=6 の `1.0e-20`, `2.5e-20 N m` では hook fail が出るが、fail 後も close pair が残る。これは「failなので不採用」と即断せず、3D動画で hook巻き付きとして許容できるかを確認する対象である。
- `flag_helix_bundle_radius_mean_um` は大きく縮まっていないため、現在の proxy は「完全な束化」ではなく「部分的な近接・束化候補」を示す。定性評価では、1本のみ独立する部分束化、hookの菌体巻き付き、非物理的なhook伸長を分けて見る。

次の定性評価対象:

- 第一候補: `n_flagella=6`, `motor.torque_Nm=5.0e-21`, `flagella.initial_helix_axis_from_rear_deg=0`
- 比較候補: `n_flagella=6`, `motor.torque_Nm=1.0e-20`, `flagella.initial_helix_axis_from_rear_deg=0`
- 高トルクfail候補: `n_flagella=6`, `motor.torque_Nm=2.5e-20`, `flagella.initial_helix_axis_from_rear_deg=0`
