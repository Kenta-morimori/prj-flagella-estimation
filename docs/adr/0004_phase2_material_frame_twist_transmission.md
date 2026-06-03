# ADR 0004: Phase 2 material frame twist transmission

- status: accepted
- date: 2026-06-03
- scope: Phase 2.6 follow-up

## Context

P2-6-007では、`axial_torque_flux_probe` と `local_twist_transmission_probe` により、root torque が螺旋全体へ届くと単一べん毛が螺旋形状を保ちつつnet回転できることを確認した。

一方で、これらのprobeは完全物理実装ではない。特に、bead force への最終変換では離れた flagellum beads へ接線力を直接入れる近似が残る。

現行 `triplet` motor が失敗する主因は、現行 bead-position-only model に segment の軸まわり姿勢、つまり material frame / segment twist を保存・輸送する状態量がないことである。現行 torsion force は4点dihedral angleの形状復元力であり、root motor torque を先端側へ伝える内部状態ではない。

## Terms

- `axial_torque_flux_probe`: root bead を起点に、flagellum のroot-to-tip主軸まわりへ接線力を分布させる診断用mode。root側を強く、先端側を弱くすることで「rootから先端へtorqueが流れたらどう見えるか」を近似する。ただし、離れたbeadへ直接forceを入れるため、最終モデルではなくprobeとして扱う。
- `local_twist_transmission_probe`: segmentごとのscalar `orientation` を内部状態として持ち、root側のorientation activityを先端側へ拡散・緩和させる診断用mode。bead forceへの変換では、activityを重みとしてflagellum beads全体へ接線forceを分布させるため、非局所force injectionが残る。
- `material_twist_local_couple`: P2-6-008で採用したmode。`local_twist_transmission_probe` と同じくscalar `orientation` を使うが、bead forceへの変換を隣接bead対の局所force coupleに限定する。これにより、flagellum全体へ一括でforceを入れる非局所force injectionを避ける。
- `orientation`: segment軸まわりの向きを表す1自由度の角度。完全な3軸material frameではなく、material frameの軸まわり成分だけを取り出した最小状態である。
- `material frame`: 各segmentに付随する局所座標系。segment接線方向だけでなく、断面が軸まわりにどちらを向いているかを表すための概念。
- `segment twist`: 隣接segmentのmaterial frame同士が、segment軸まわりにどれだけ相対回転しているかを表す局所ねじれ量。

## Decision

`motor.force_distribution=material_twist_local_couple` を追加する。

この mode は、material frame の軸まわり成分を最小化した scalar `orientation` を segment ごとに持ち、root torque を局所的な twist state として保存・伝搬させる。

`orientation` は、segment 軸まわりの向きを表す1自由度の角度であり、material frame の軸まわり成分の簡略表現として扱う。本実装では、完全な3軸material frameではなく、この scalar orientation を最小物理状態として採用する。

今回確定した実装は、explicitなtwist potentialからforceを導出するenergy-based実装ではない。root torque を `orientation_0` へ入力し、orientation activityをdiffusion/relaxationで先端側へ伝え、そのactivityを隣接bead対の局所force coupleへ変換する最小実装である。

## Confirmed Implementation

P2-6-008で採用した実装は以下である。

1. root torque を root側 `orientation_0` に入力する。
2. orientation activity を diffusion/relaxation で先端側へ伝える。
3. segmentごとの activity を torque weight として正規化する。
4. 各segmentの隣接bead対に作用反作用のforce coupleを置く。
5. 各force coupleはflagellum主軸まわりの局所torqueを作る。
6. flagellum側の合計torqueと反対向きのreaction torqueをbody側へ与える。

つまり、**P2-6-008で確定したのは harmonic potential 型でも cosine potential 型でもなく、diffusion/relaxationで伝搬したorientation activityを局所force coupleへ変換する実装**である。

## Future Candidate: Energy-Based Twist Potential

将来、完全なmaterial frame / segment twistモデルへ寄せる場合は、twist potentialからtorqueとbead forceを導出するenergy-based実装を比較候補にする。

第一候補は harmonic potential である。

```text
local_twist_i = wrap_angle(orientation_{i+1} - orientation_i - rest_twist_i)
U_twist_i = 0.5 * k_twist * local_twist_i^2
torque_i = -dU_twist / d orientation_i
```

この式は、小変形の弾性エネルギーがずれの2乗に比例すること、現行 spring / bend / torsion potential の設計思想と整合すること、rod/filament のねじれ弾性を離散化した形として自然であることを根拠にする。

ただし、これはP2-6-008の確定実装ではない。角度は周期量であるため、大回転で harmonic potential が不安定または不自然な場合は、将来タスクで `k_twist * (1 - cos(local_twist_i))` 型を比較候補にする。

## Torsion force との役割分担

現行 torsion force は、4つのbead位置からdihedral angleを計算し、螺旋形状が基準から崩れたときに戻す力である。これは形状維持には効くが、root motor torque を時間発展する内部ねじれ状態として保存・伝搬する仕組みではない。

本実装では、まず以下の役割分担で検証する。

- 現行 torsion force: paper normal state の螺旋形状維持
- material frame / segment twist: root torque の保存・伝搬

torsion force OFF + `material_twist_local_couple` ON の診断では、0.5 s 条件で `flag_torsion_err_max_deg=179.996`, `flag_bond_rel_err_max=0.478`, `hook_len_rel_err_max=0.766` となり、shape gate が fail した。したがって、現時点では既存 torsion force を置き換えない。

この結果から、既存 torsion force は螺旋形状維持を担い、新しい material twist local couple は root torque の伝搬を担う、という役割分担を採用する。将来置き換えを検討する場合は、Phase 2.2で固定した paper normal state の幾何契約、長時間条件、複数べん毛条件を再検証する。

## Difference From Previous Probes

`material_twist_local_couple` は `local_twist_transmission_probe` と同じく内部orientation状態を持つが、force変換が異なる。

- `local_twist_transmission_probe`: activityをflagellum beads全体の接線force重みに使う。root torque伝搬の診断として有効だが、離れたbeadへ直接forceを入れる。
- `material_twist_local_couple`: activityを隣接bead対ごとの局所force coupleに変換する。各segment近傍でtorqueを作るため、非局所force injectionを避ける。

`axial_torque_flux_probe` はmaterial frameやorientation状態を持たず、root-to-tip主軸に沿って接線forceを分布させるだけの診断modeである。P2-6-008の採用対象ではない。

## Acceptance criteria

P2-6-008の受入条件は、少なくとも以下とする。

- default `triplet` と既存 paper-compatible geometry は変更しない。
- 新しい material-frame系挙動は明示的な mode または設定で有効化する。
- 非局所 force injection を使わず、隣接segmentまたは局所bead群への force couple として torque を伝える。
- `time.dt_star=1.0e-4` を明示する。
- `duration_s >= 0.5`、可能なら `1.0` で評価する。
- `net_abs_flag_helix_spin_revolutions >= 1.0` を満たす。
- `flag_helix_spin_direction_consistency >= 0.5` を満たす。
- `hook_len_rel_err_max <= 0.5` を満たす。
- `flag_bond_rel_err_max <= 0.25` を満たす。
- `flag_bend_err_max_deg <= 30` を満たす。
- `flag_torsion_err_max_deg <= 60` を満たす。

## Verification result

代表条件:

- `motor.force_distribution=material_twist_local_couple`
- `motor.torque_Nm=2.0e-20`
- `time.dt_star=1.0e-4`
- `duration_s=0.5`
- `local_hook_scale=1.0`
- `local_spring_scale=1.2`
- `local_bend_scale=1.0`
- `local_torsion_scale=1.0`

結果:

- `helix_retention_pass=True`
- `first_fail_category=none`
- `net_abs_flag_helix_spin_revolutions=1.11698`
- `flag_helix_spin_direction_consistency=0.98175`
- `helix_to_root_net_rotation_ratio=8.77118`
- `hook_len_rel_err_max=0.40113`
- `flag_bond_rel_err_max=0.17373`
- `flag_bend_err_max_deg=4.34918`
- `flag_torsion_err_max_deg=6.21047`
- `local_twist_root_orientation_deg=66.809`
- `local_twist_tip_orientation_deg=23.273`
- `local_twist_tip_activity_ratio=0.34836`

この条件は、`duration_s=0.5`, `time.dt_star=1.0e-4`, net 1回転以上、shape gate PASS、非局所force injectionなし、というP2-6-008の最低条件を満たす。

## Consequences

- material frame / segment twist は参照論文モデルそのものの完全再現ではなく、現行 bead-spring 実装への拡張である。
- scalar orientation は完全な3軸material frameより簡略だが、root torque を内部ねじれ状態として保存・伝搬し、局所force coupleへ変換するため、P2-6-007のprobeより物理的に筋が通る。
- 現時点では既存 torsion force を置き換えず、螺旋形状維持として残す。
