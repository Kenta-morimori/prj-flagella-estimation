# Phase 2.6: トルク伝搬拡張モデルの詳細評価計画

## 目的

Issue #54 では、`material_twist_local_couple` 導入後の単一べん毛モデルについて、高トルク条件での形状安定性と local stiffness scaling の必要性を定量評価する。

特に、論文モデル相当の `local_*_scale=1.0` でどこまで安定するか、scaling が必要な場合にどの局所項が破綻防止に効くかを明確にする。

この結果を確定してから、Issue #58 の多べん毛・後方束化・遊泳検証へ進む。

## 背景

Phase 2.6 では、従来の `triplet` motor で root torque が螺旋全体の net 回転へ十分に伝わらないことを確認した。

その後、probe 系 mode で「トルクが螺旋全体へ届けば綺麗に回る」ことを診断し、P2-6-008 では `material_twist_local_couple` を採用した。代表条件では、`motor.torque_Nm=2.0e-20`, `time.dt_star=1.0e-4`, `duration_s=0.5`, `motor.local_spring_scale=1.2` で shape gate PASS と net 1回転以上を確認した。

一方で、代表条件には `local_spring_scale=1.2` が含まれている。これが物理的に必要な拡張なのか、数値安定化なのか、あるいは新しいトルク伝搬機構の導入後は不要なのかは未整理である。

## 評価対象

主対象:

- `n_flagella=1`
- `flagella.stub_mode=full_flagella`
- `motor.force_distribution=material_twist_local_couple`
- `time.dt_star=1.0e-4`
- `duration_s=0.5`
- Brownian off の決定論的条件

比較対象:

- 原則として必須比較は置かず、`material_twist_local_couple` を単独で評価する。
- `distributed_flagellum` や probe 系 mode は、評価結果の解釈に詰まった場合の追加診断に限定する。

理由は、Issue #54 の主目的が「採用済み mode の安定境界と scaling 必要性」を明確にすることであり、診断用 mode との再比較を主軸にすると評価対象が広がりすぎるためである。

## 評価軸

### 1. torque sweep

まず `local_*_scale=1.0` を基準に、torque に対する破綻境界を確認する。

初期候補:

- `0.5e-20 N m`
- `1.0e-20 N m`
- `2.0e-20 N m`
- `4.0e-20 N m`
- `6.0e-20 N m`
- `8.0e-20 N m`
- `1.0e-19 N m`

`2.0e-20 N m` は P2-6-008 の代表条件であるため、必ず含める。

上限は `1.0e-19 N m` とする。`material_twist_local_couple` 未反映の単一べん毛条件でも `time.dt_star=1.0e-4` でこのトルク帯を扱えた実績があるため、採用済み拡張モデルでも同じ上限まで形状安定性と駆動力を評価する。

全 torque 条件は `duration_s=0.5` で評価する。短時間 screening は使わず、比較対象を同じ時間条件に揃える。

### 2. local scaling sweep

torque sweep で破綻する条件が見えたら、以下を個別に sweep する。

- `motor.local_spring_scale`
- `motor.local_bend_scale`
- `motor.local_torsion_scale`
- `motor.local_hook_scale`

初期候補:

- `1.0`
- `1.1`
- `1.2`
- `1.5`
- `2.0`

まず one-factor sweep で支配的な項を探し、その後、支配的な local scale が見えた場合だけ必要な2軸を heatmap 化する。

### 3. 代表条件での遊泳駆動確認

shape gate と net 回転を満たす代表条件について、菌体の重心変位と速度を確認する。

Issue #54 の主対象は単一べん毛モデル評価であり、多べん毛の後方束化・遊泳検証は Issue #58 に分離する。ただし、単一べん毛代表条件で菌体位置が長時間でほとんど変化しない場合、その torque 条件は後続の遊泳検証に渡す候補として弱い。

そのため、代表 PASS 条件では以下も記録する。

- body center displacement
- body mean speed
- body axis 角度変化
- torque と重心変位の関係

この指標は「多べん毛で遊泳する」ことの判定ではなく、Issue #58 へ渡す torque 条件が駆動力として十分かを判断するための前段評価である。

### 4. torsion force の役割確認

既存 torsion force は、root torque を伝搬するものではなく、螺旋形状を normal state へ戻す形状復元力である。

ただし、P2-6-008 ですでに torsion force OFF + `material_twist_local_couple` ON は shape gate fail となることを確認している。そのため、Issue #54 の初回評価では torsion force OFF を必須項目にしない。

今回の基本方針は以下とする。

- 既存 torsion force は螺旋形状維持のために残す。
- `material_twist_local_couple` は root torque の保存・伝搬を担う。
- torsion force OFF は、評価結果で置き換え可能性を再検討する必要が出た場合の追加診断に限る。

追加診断を行う場合は、以下を確認する。

- `material_twist_local_couple` が torsion force の代替になっていないこと。
- torsion force がないと、どの形状指標が最初に破綻するか。
- torsion force を残す必要性が shape gate と first-fail category で説明できるか。

torsion force OFF + 新手法 ON で一部条件が通った場合でも、Phase 2.2 の paper normal state 幾何契約、長時間条件、多べん毛条件を再評価するまでは、既存 torsion force の置き換えとは扱わない。

## 指標

最低限記録する指標:

- `helix_retention_pass`
- `shape_pass_nonbody`
- `first_fail_category_nonbody`
- `net_abs_flag_helix_spin_revolutions`
- `flag_helix_spin_direction_consistency`
- `helix_to_root_net_rotation_ratio`
- `hook_len_rel_err_max`
- `local_attach_first_rel_err`
- `flag_bond_rel_err_max`
- `flag_bend_err_max_deg`
- `flag_torsion_err_max_deg`
- `local_twist_root_orientation_deg`
- `local_twist_tip_orientation_deg`
- `local_twist_tip_activity_ratio`
- `body_displacement_um`
- `body_mean_speed_um_s`
- `body_axis_angle_change_deg`

解釈では、shape gate PASS だけでなく、net 回転と方向一貫性を同時に見る。

## 成果物

- `local_*_scale=1.0` の torque 安定境界表
- one-factor sweep による local scaling rescue 可否の表
- 支配的な local scale が見えた場合のみ、torque x scaling の破綻 heatmap
- 支配的な local scale が見えた場合のみ、torque x scaling の first-fail category heatmap
- 代表 PASS / FAIL 条件表
- `local_*_scale=1.0` の限界条件
- scaling が必要な場合の最小条件
- scaling の役割に関する考察
- 代表条件での菌体重心変位・速度
- Issue #58 へ渡す代表条件

## 受け入れ条件

- `local_*_scale=1.0` での torque 安定境界が報告されている。
- `motor.local_spring_scale=1.2` が必要か、不要か、または条件付きで必要かが説明されている。
- scaling が必要な場合、どの局所項が効いているかが one-factor sweep または heatmap で示されている。
- 代表 PASS 条件と代表 FAIL 条件が `duration_s>=0.5`, `time.dt_star=1.0e-4` で再現可能である。
- 既存 torsion force の役割が明確化されている。
- 代表 PASS 条件で、菌体重心変位または平均速度が報告されている。
- 多べん毛・後方束化に進むための代表条件が1つ以上提示されている。

## 評価結果（2026-06-06）

実行条件:

- `motor.force_distribution=material_twist_local_couple`
- `flagella.n_flagella=1`
- `flagella.stub_mode=full_flagella`
- `time.dt_star=1.0e-4`
- `duration_s=0.5`
- `motor.local_hook_scale=1.0`
- `motor.local_spring_scale=1.0`
- `motor.local_bend_scale=1.0`
- `motor.local_torsion_scale=1.0`
- Brownian off

実行コマンド:

```bash
uv run python -m scripts.01_simulate_swimming.run_phase2_6_torque_model_evaluation --scale-torques none
uv run python -m scripts.01_simulate_swimming.run_phase2_6_torque_model_evaluation --torques 2.5e-20,3.0e-20,3.5e-20 --scale-torques none
uv run python -m scripts.01_simulate_swimming.run_phase2_6_torque_model_evaluation --torques 3.5e-20 --scale-torques all --scale-values 2.0
```

出力:

- `outputs/2026-06-06/195006/phase2_6_torque_model_evaluation/phase2_6_torque_model_evaluation_summary.csv`
- `outputs/2026-06-06/201058/phase2_6_torque_model_evaluation/phase2_6_torque_model_evaluation_summary.csv`
- `outputs/2026-06-06/201631/phase2_6_torque_model_evaluation/phase2_6_torque_model_evaluation_summary.csv`

### local scale 1.0 の torque 境界

| torque [N m] | 判定 | 主な理由 | net helix rev | direction consistency | flag bond max | bend max [deg] | torsion max [deg] | body displacement [um] | mean speed [um/s] |
|---:|---|---|---:|---:|---:|---:|---:|---:|---:|
| `0.5e-20` | FAIL | `motor_no_rotation` | 0.288 | 1.000 | 0.172 | 3.22 | 3.21 | 0.030 | 0.060 |
| `1.0e-20` | FAIL | `motor_no_rotation` | 0.564 | 1.000 | 0.172 | 3.82 | 3.93 | 0.073 | 0.146 |
| `2.0e-20` | PASS | none | 1.120 | 0.985 | 0.174 | 4.35 | 6.19 | 0.057 | 0.114 |
| `2.5e-20` | PASS | none | 1.242 | 0.682 | 0.176 | 6.33 | 27.07 | 0.075 | 0.151 |
| `3.0e-20` | PASS | none | 1.297 | 0.704 | 0.249 | 18.61 | 51.95 | 0.119 | 0.239 |
| `3.5e-20` | FAIL | `flag` | 1.328 | 0.486 | 1.648 | 96.48 | 105.45 | 0.190 | 0.379 |
| `4.0e-20` | FAIL | `flag` | 1.359 | 0.304 | 0.627 | 62.97 | 157.39 | 0.262 | 0.524 |
| `6.0e-20` | FAIL | `flag` | 1.488 | 0.152 | 2.025 | 102.31 | 151.59 | 0.619 | 1.239 |
| `8.0e-20` | FAIL | `flag` | 0.612 | 0.043 | 2.712 | 101.34 | 159.34 | 0.969 | 1.938 |
| `1.0e-19` | FAIL | `flag` | 1.834 | 0.070 | 6.646 | 131.47 | 175.22 | 1.254 | 2.507 |

解釈:

- `local_*_scale=1.0` のまま、`2.0e-20` から `3.0e-20 N m` は 0.5 s で shape gate と net 1回転以上を満たした。
- `0.5e-20` と `1.0e-20 N m` は形状破綻ではなく、0.5 s 内に net 1回転へ届かないため FAIL とした。
- `3.5e-20 N m` 以上は net 回転自体は出るが、flag bond / bend / torsion が gate を超え、螺旋形状を保てない。
- 上限 `1.0e-19 N m` は、菌体変位は大きいが flag 形状破綻が支配的であり、Issue #58 へ渡す代表条件にはしない。

### local scaling の必要性

`3.5e-20 N m` で各 local scale を個別に `2.0` へ上げた one-factor rescue を行った。

| target | value | 判定 | 主な理由 | net helix rev | direction consistency | flag bond max | bend max [deg] | torsion max [deg] |
|---|---:|---|---|---:|---:|---:|---:|---:|
| baseline | 1.0 | FAIL | `flag` | 1.328 | 0.486 | 1.648 | 96.48 | 105.45 |
| `local_spring_scale` | 2.0 | FAIL | `flag` | 1.323 | 0.466 | 1.648 | 96.48 | 105.45 |
| `local_bend_scale` | 2.0 | FAIL | `flag` | 1.315 | 0.438 | 1.821 | 99.41 | 107.06 |
| `local_torsion_scale` | 2.0 | FAIL | `flag` | 1.113 | 0.083 | 1.623 | 89.09 | 157.12 |
| `local_hook_scale` | 2.0 | FAIL | `flag` | 1.327 | 0.487 | 1.648 | 96.48 | 105.45 |

結論:

- `motor.local_spring_scale=1.2` は、少なくとも `motor.torque_Nm=2.0e-20`, `duration_s=0.5`, `time.dt_star=1.0e-4` では不要である。scale 1.0 のまま PASS した。
- `3.5e-20 N m` の flag 破綻は、単純に local spring / bend / torsion / hook を `2.0` へ上げても救済できなかった。
- one-factor rescue で支配的な local scale が見えなかったため、今回の評価では torque x scaling heatmap は作成していない。
- したがって、今回の結果では local scaling を高トルク安定化の主手段として採用しない。高トルク側をさらに安定させたい場合は、局所剛性倍率ではなく、flag geometry、repulsion、時間刻み、または motor force の分配式そのものを別タスクで検討する。

### 代表条件

Issue #58 の多べん毛・後方束化検証へ渡す候補は以下とする。

- 第一候補: `motor.torque_Nm=2.5e-20`
  - shape gate PASS、net helix rev 1.242、body mean speed 0.151 um/s。
  - `3.0e-20` より shape margin が広く、初期代表条件として安全。
- 上限側候補: `motor.torque_Nm=3.0e-20`
  - shape gate PASS、net helix rev 1.297、body mean speed 0.239 um/s。
  - torsion max 51.95 deg で gate 上限 60 deg に近いため、後続では注意して扱う。
- 失敗境界代表: `motor.torque_Nm=3.5e-20`
  - flag 破綻が再現するため、崩壊診断・比較用に残す。

単一べん毛では菌体重心変位は確認できるが、body axis 角度変化も大きい。これは「単一べん毛で安定遊泳する」確認ではなく、Issue #58 で多べん毛・後方束化条件へ進むための駆動力評価として扱う。

## 今回の実装範囲

最初のPRでは、評価のための sweep tooling、集計、可視化、ドキュメント整理を中心にする。

物理モデル本体の変更は、評価で明確な不足が示された場合に限定する。新しい物理実装を追加する場合は、論文モデルとの差分、数値安定化との区別、ADR要否を review_result に記録する。

## 実装前確認事項と方針

1. torque sweep の上限をどこまで含めるか。
2. `duration_s=0.5` を全条件に適用するか、短時間 screening 後に代表条件だけ 0.5 s へ伸ばすか。
3. heatmap の主軸を `local_spring_scale` にするか、one-factor sweep の結果で決めるか。
4. `distributed_flagellum` と probe 系 mode をどこまで比較対象に含めるか。
5. torsion force OFF 診断を今回の必須項目にするか、追加診断に分けるか。

2026-06-06時点の方針:

- torque sweep は上限を `1.0e-19 N m` とする。
- 全条件を `duration_s=0.5` で評価する。
- heatmap の主軸は one-factor sweep の結果で決める。今回の結果では支配的な local scale が見えないため、heatmap は作成しない。
- `distributed_flagellum` と probe 系 mode は必須比較に含めない。
- torsion force OFF 診断は必須にしない。
