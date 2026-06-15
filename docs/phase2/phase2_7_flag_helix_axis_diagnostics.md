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

## 確認観点

- `flag_helix_axis_vs_rear_angle_deg_*` が意図した後方向き条件を反映しているか。
- 第1ビーズの位置ずれや hook伸長が、螺旋軸推定に混入していないか。
- `flag_helix_axis_fit_r2_min` が極端に低い条件では、螺旋軸の解釈が不安定でないか。
- `hook_len_rel_err_max` が大きい条件では、3D動画で hook の巻き付き・伸長を定性評価する。

