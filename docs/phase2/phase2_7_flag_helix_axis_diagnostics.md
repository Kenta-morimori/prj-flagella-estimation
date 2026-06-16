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
- `flag_helix_axis_pair_angle_deg_mean`
- `flag_helix_axis_pair_angle_deg_max`
- `flag_helix_axis_mean_deviation_deg_max`
- `flag_helix_axis_alignment_order`
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

## 時系列plot

`run_phase2_7_bundling_sweep.py` は、各run directoryに `flag_helix_axis_angles_timeseries.png` を出力する。

plotは上下2段で確認する。

- 上段: 各べん毛の `flag_helix_axis_vs_rear_angle_deg`。`0 deg` は菌体後方、`90 deg` は側方を示す参照線である。
- 下段: 各べん毛軸の平均軸からの偏差。`15 deg` を軸整列の閾値線として表示する。

後半80%の評価区間は背景で示す。1本だけ軸方向が外れる場合は、そのflagの線だけが他と離れるため、plot上でも外れとして確認できる。

## 軸整列の定量定義

本PRでは、束化を「近接した1本の束に見えること」ではなく、「複数べん毛の螺旋中心軸方向が安定的に揃うこと」と定義する。

- 主指標は `flag_helix_axis_mean_deviation_deg_max`。`flag_helix_axis_pair_angle_deg_max` は軸同士の最大開き幅を示す補助指標として使う。
- 成功閾値は `15 deg`。
- `duration_s=0.5` の後半80%の全stepで、平均軸からの最大偏差が閾値以内なら、軸整列を安定とみなす。
- 後方か側方かは主分類ではなく、`flag_helix_axis_vs_rear_angle_deg_*` と時系列plotで副次的に確認する。
- hook長failのみが出る場合は、菌体への巻き付き候補として strict collapse とは分けて記録する。

## 後方螺旋軸初期条件

`flagella.initial_helix_axis_from_rear_deg` を指定すると、各べん毛の第2ビーズ以降から推定される螺旋中心軸を菌体後方から指定角度だけ傾ける。`0` は全べん毛の螺旋軸を同一の菌体後方方向へ揃える条件である。`null` は従来の side-attach 条件を維持する。

この設定は第1ビーズを菌体側hookとして残したまま、螺旋本体の軸だけを後方へ揃える診断条件である。そのため、hook長やhook近傍の形状破綻は別途確認する必要がある。

## 確認観点

- `flag_helix_axis_vs_rear_angle_deg_*` が意図した後方向き条件を反映しているか。
- `flag_helix_axis_mean_deviation_deg_max` が後半80%で `15 deg` 以内に収まるか。
- 第1ビーズの位置ずれや hook伸長が、螺旋軸推定に混入していないか。
- `flag_helix_axis_fit_r2_min` が極端に低い条件では、螺旋軸の解釈が不安定でないか。
- `flag_helix_axis_angles_timeseries.png` で、後方/側方の区別と1本外れの有無を確認する。
- `flag_flag_helix_bead_dist_min_um` と `flag_helix_bundle_radius_*` が減少するかを見て、fail判定時でも束化らしい接近が起きていないかを確認する。
- `hook_len_rel_err_max` が大きい条件では、3D動画で hook の巻き付き・伸長を定性評価する。
