# ADR 0004: Phase 2 material frame twist transmission

- status: proposed
- date: 2026-06-03
- scope: Phase 2.6 follow-up

## Context

P2-6-007では、`axial_torque_flux_probe` と `local_twist_transmission_probe` により、root torque が螺旋全体へ届くと単一べん毛が螺旋形状を保ちつつnet回転できることを確認した。

一方で、これらのprobeは完全物理実装ではない。特に、bead force への最終変換では離れた flagellum beads へ接線力を直接入れる近似が残る。

現行 `triplet` motor が失敗する主因は、現行 bead-position-only model に segment の軸まわり姿勢、つまり material frame / segment twist を保存・輸送する状態量がないことである。現行 torsion force は4点dihedral angleの形状復元力であり、root motor torque を先端側へ伝える内部状態ではない。

## Decision

完全物理実装の次タスクでは、material frame / segment twist を導入し、root torque を局所的な twist state として保存・伝搬させる。

`orientation` は、segment 軸まわりの向きを表す1自由度の角度であり、material frame の軸まわり成分の簡略表現として扱う。最小検証では orientation を使い、最終採用候補は material frame + segment twist とする。

`material frame` は、各 segment に付随する局所座標系であり、segment 接線方向だけでなく断面の向きを表す。隣接する material frame 同士の軸まわり相対回転を `segment twist` として定義する。

twist potential の第一候補は harmonic potential とする。

```text
local_twist_i = wrap_angle(orientation_{i+1} - orientation_i - rest_twist_i)
U_twist_i = 0.5 * k_twist * local_twist_i^2
torque_i = -dU_twist / d orientation_i
```

この式は、小変形の弾性エネルギーがずれの2乗に比例すること、現行 spring / bend / torsion potential の設計思想と整合すること、rod/filament のねじれ弾性を離散化した形として自然であることを根拠にする。

角度は周期量であるため、実装では `wrap_angle` により `[-pi, pi]` へ畳む。大回転で harmonic potential が不安定または不自然な場合は、`k_twist * (1 - cos(local_twist_i))` 型を比較候補にする。

## Torsion force との役割分担

現行 torsion force は、4つのbead位置からdihedral angleを計算し、螺旋形状が基準から崩れたときに戻す力である。これは形状維持には効くが、root motor torque を時間発展する内部ねじれ状態として保存・伝搬する仕組みではない。

次タスクでは、まず以下の役割分担で検証する。

- 現行 torsion force: paper normal state の螺旋形状維持
- material frame / segment twist: root torque の保存・伝搬

二重カウントが見えた場合は、既存 torsion force と新しい twist potential の役割を再定義する。いきなり既存 torsion force を置き換えない。Phase 2.2で固定した paper normal state の幾何契約を壊さないためである。

## Acceptance criteria

次タスクの受入条件は、少なくとも以下とする。

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

## Consequences

- material frame / segment twist は参照論文モデルそのものの完全再現ではなく、現行 bead-spring 実装への拡張である。
- ただし、root torque を内部ねじれ状態として保存・伝搬し、非局所force injectionを避けるため、probeより物理的に筋が通る。
- 実装時はADR更新または追加ADRで、既存 torsion force との二重カウント、force couple の局所性、安定性gateを記録する。
