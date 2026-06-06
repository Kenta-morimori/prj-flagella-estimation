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
- `duration_s>=0.5`
- Brownian off の決定論的条件

比較対象:

- `motor.force_distribution=triplet`
- `motor.force_distribution=distributed_flagellum`
- 必要に応じて probe 系 mode

比較対象は最終候補ではなく、`material_twist_local_couple` の安定性・回転性を解釈するための reference として扱う。

## 評価軸

### 1. torque sweep

まず `local_*_scale=1.0` を基準に、torque に対する破綻境界を確認する。

初期候補:

- `1.0e-20 N m`
- `1.5e-20 N m`
- `2.0e-20 N m`
- `2.5e-20 N m`
- `3.0e-20 N m`

`2.0e-20 N m` は P2-6-008 の代表条件であるため、必ず含める。

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

まず one-factor sweep で支配的な項を探し、その後必要な2軸だけ heatmap 化する。

### 3. torsion force の役割確認

既存 torsion force は、root torque を伝搬するものではなく、螺旋形状を normal state へ戻す形状復元力である。

そのため、torsion force OFF は置き換え可否の即時判断ではなく、以下を確認する診断条件として扱う。

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

解釈では、shape gate PASS だけでなく、net 回転と方向一貫性を同時に見る。

## 成果物

- torque x scaling の破綻 heatmap
- torque x scaling の first-fail category heatmap
- 代表 PASS / FAIL 条件表
- `local_*_scale=1.0` の限界条件
- scaling が必要な場合の最小条件
- scaling の役割に関する考察
- Issue #58 へ渡す代表条件

## 受け入れ条件

- `local_*_scale=1.0` での torque 安定境界が報告されている。
- `motor.local_spring_scale=1.2` が必要か、不要か、または条件付きで必要かが説明されている。
- scaling が必要な場合、どの局所項が効いているかが one-factor sweep または heatmap で示されている。
- 代表 PASS 条件と代表 FAIL 条件が `duration_s>=0.5`, `time.dt_star=1.0e-4` で再現可能である。
- 既存 torsion force の役割が明確化されている。
- 多べん毛・後方束化に進むための代表条件が1つ以上提示されている。

## 今回の実装範囲

最初のPRでは、評価のための sweep tooling、集計、可視化、ドキュメント整理を中心にする。

物理モデル本体の変更は、評価で明確な不足が示された場合に限定する。新しい物理実装を追加する場合は、論文モデルとの差分、数値安定化との区別、ADR要否を review_result に記録する。

## 事前に詰めるべき点

現時点では、以下を実装開始前の確認事項として残す。

1. torque sweep の上限を `3.0e-20 N m` で十分とするか、`4.0e-20 N m` 以上も含めるか。
2. `duration_s=0.5` を全条件に適用するか、短時間 screening 後に代表条件だけ 0.5 s へ伸ばすか。
3. heatmap の主軸を `local_spring_scale` にするか、one-factor sweep の結果で決めるか。
4. `distributed_flagellum` と probe 系 mode をどこまで比較対象に含めるか。
5. torsion force OFF 診断を今回の必須項目にするか、追加診断に分けるか。
