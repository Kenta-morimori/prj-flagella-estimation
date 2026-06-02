# Phase 2.6 triplet motor のねじれ・回転自由度設計メモ

## 目的

`distributed_flagellum` を最終的な物理モデルとして扱わず、従来の `triplet` motor と hook spring fix を中心に、root torque が flagellum chain へ伝わらない原因を明確にする。

本メモでは、現行 bead-position-only model に不足している状態量を整理し、torque transmission model の候補を実装検証する。

## 現行モデルで起きていること

`triplet` motor は `attach-first-second` の3点に対して力カップルを与える。`attach-first` hook spring force を復元しても、root 方位の net 回転と螺旋全体の net 回転は分離している。

代表的な `triplet` 条件:

- `motor.force_distribution=triplet`
- `torque=2.5e-20 N m`
- `time.dt_star=1.0e-4`
- `local_scale=(4,2,2,2)`
- `duration_s=0.25`

この条件では、`net_abs_flag_root_revolutions=0.8796` に対して `net_abs_flag_helix_spin_revolutions=0.00139`、`helix_to_root_net_rotation_ratio=0.00158` だった。root 近傍は動くが、螺旋全体の一方向回転へほぼ伝わっていない。

## 不足している状態量

噛み砕くと、現状の問題は「根本をねじっているが、べん毛内部にそのねじれを隣へ渡す部品がない」ということである。

時間積分は順次進むため、bead位置の変化は次step以降に周囲へ影響する。しかし、それは位置の押し引きであり、segmentが自分の軸まわりにどれだけ回ったかを保存して隣へ渡す仕組みとは別である。

そのため、根本の3点に入った torque はroot近傍の小さな動きとしては見えるが、「segment0がねじれたのでsegment1へねじれを渡す」という内部伝達にはならない。

現時点で明確になっている失敗機序は以下である。

1. `triplet` motor は `attach-first-second` の局所3点だけに力カップルを入れる。
2. その torque は root 近傍の bead positions を動かすが、現行モデルには segment の軸まわり姿勢を保持する状態量がない。
3. 現行 torsion force は4点 dihedral angle の形状復元力であり、root から先端へ twist angular momentum や torque flux を輸送する状態ではない。
4. そのため root 方位 proxy は動いても、螺旋全体の各beadに一方向の接線速度が作られず、局所変形・弾性復元・フィットjitterとして消える。

つまり、`triplet + hook spring fix` が失敗している主因は、hook spring の有無だけではなく、root torque を flagellum chain の軸まわり回転として伝える自由度と伝達則が現行 bead-position-only model にないことである。

### material frame

各 flagellum segment に付随する局所座標系で、segment 接線方向だけでなく断面の向きを表す。現行モデルは bead positions だけを状態として持つため、segment の「軸まわり姿勢」を直接保持しない。

### segment twist

隣接 segment の material frame 同士が、segment 軸まわりにどれだけ相対回転しているかを表すねじれ量。root torque を弾性的な torsional deformation として chain に蓄える・伝えるために必要である。

これは、べん毛全体に1つのねじれ量を与えるものではない。隣接するsegment同士の向きの差から、局所的に判断する。

```text
segment0 と segment1 の向きの差 = local_twist_0
segment1 と segment2 の向きの差 = local_twist_1
```

現行の torsion force は4点 dihedral angle の位置エネルギー勾配であり、形状のねじれを保つ力としては働く。しかし、segment 自身が軸まわりにどの向きを持つかを独立状態として持たないため、root からの motor torque を material twist として保存・輸送する経路にはなっていない。

### axial torque flux

root から先端側へ、segment 軸方向に伝わる torque の流れ。現行 `triplet` は root 近傍に局所的な力カップルを入れるが、各 segment 間で torque flux を保存・減衰・伝搬させる明示的な状態や計算はない。

## 設計候補

| 候補 | 内容 | 論文モデルとの距離 | 利点 | リスク | P2-6での扱い |
| --- | --- | --- | --- | --- | --- |
| material frame + segment twist | 各 segment に局所 frame または twist angle を持たせ、root torque を torsional spring/damping として隣接 segment へ伝える | 中。bead-spring model を拡張するが、軸まわり回転自由度として物理的に明確 | torque の保存・伝搬を説明しやすい。`triplet` の root 駆動に近い | 状態変数、積分器、安定性 gate、ADR が必要 | 次の本命候補 |
| segment torsional torque flux | bead positions は維持しつつ、隣接 segment に軸方向 torque flux を近似分配する | 中-大。状態量を省略した有効モデル | 実装量は material frame より小さい可能性がある | 何を保存しているかが曖昧になりやすい。force couple への変換で数値不安定化しやすい | material frame 案の比較対象 |
| axial_torque_flux_probe | root bead を原点、root-to-tip 主軸を torque flux 軸として、root から先端へ弱く減衰する接線力を分布させる | 中-大。material frame なしの実験用 probe | 現行モデルのまま torque flux 近似の有効性を検証できる | 離れた bead に直接力を入れるため、真の segment twist 伝搬ではない | P2-6-007で実装検証 |
| quasi-rigid helical body approximation | 単一 flagellum 螺旋をほぼ剛体として root torque で全体回転させる | 大。ML用生成には実用的だが論文モデルから離れる | 形状維持と回転は作りやすい | collapse 診断や elastic failure の意味が薄くなる | Phase 2.6の主案にはしない |
| distributed_flagellum | flagellum 全体に主軸まわりのゼロ合力 torque drive を分布させる | 大。root motor ではなく診断用 upper-bound | torque が螺旋全体へ伝わった場合の見え方を確認できる | 離れた bead に直接 torque を入れるため、現実の root torque 伝搬ではない | mode として残すが最終解にしない |

## `axial_torque_flux_probe` と `distributed_flagellum` の違い

両者はどちらも、離れた flagellum beads に直接接線力を入れる近似であり、root motor torque が material twist として物理的に伝わる実装ではない。この点ではどちらも論文モデルから離れている。

違いは、先端側で減衰させているかだけではない。

| 観点 | `distributed_flagellum` | `axial_torque_flux_probe` |
| --- | --- | --- |
| 回転中心 | flagellum 全体の重心 | root bead |
| 回転軸 | flagellum 点群の主軸 | root-to-tip 向きに揃えた主軸 |
| 力分布 | 全flagellum beadへほぼ一様な接線drive | root側を強く、先端側を弱くする減衰drive |
| 物理的解釈 | 「螺旋全体へtorqueが届いたらどう見えるか」のupper-bound | 「rootから先端へtorque fluxが流れたらどう見えるか」の妥協近似 |
| 非現実性 | 離れたbeadへ同時に直接torqueを入れる | 離れたbeadへ直接力を入れるが、root起点・減衰で伝搬らしく寄せる |
| 本PRでの位置づけ | 診断用baseline | 最終妥協案の候補 |

したがって、最終的な妥協案としてどちらかを残すなら、`axial_torque_flux_probe` の方が妥当である。理由は、root beadを起点にし、先端側へ弱くなる力分布にしているため、`distributed_flagellum` より「root motor torque が軸方向に伝わる」という解釈に近いからである。

## 実装検証

P2-6-007では、`motor.force_distribution=axial_torque_flux_probe` を追加した。これは material frame を明示しないが、root bead を原点とする torque flux 軸に沿って flagellum beads へ接線力を分布させる probe である。詳細なモデル判断は `docs/adr/0002_phase2_axial_torque_flux_probe.md` に記録する。

短時間比較 (`duration_s=0.05`, `dt_star=1.0e-4`, `torque=2.0e-20`, `local_scale=(1,1.2,1,1)`):

| mode | pass | first fail | root net rev | helix net rev | transfer ratio | direction consistency | 備考 |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `triplet` | false | `motor_no_rotation` | `0.03789` | `0.00023` | `0.00618` | `0.01690` | root 近傍の揺れに留まる |
| `axial_torque_flux_probe` | true | `none` | `0.04109` | `0.11566` | `2.81489` | `1.00000` | torque flux 近似で一方向 net 回転が出る |
| `distributed_flagellum` | true | `none` | `0.05964` | `0.10966` | `1.83882` | `1.00000` | diagnostic upper-bound |

長時間確認 (`duration_s=0.5`, `dt_star=1.0e-4`):

| mode | torque | local scales `(hook,spring,bend,torsion)` | pass | helix net rev | direction consistency | hook err | bond err | bend err deg | torsion err deg |
| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `axial_torque_flux_probe` | `2.0e-20` | `(1,1.2,1,1)` | true | `1.08546` | `1.00000` | `0.39986` | `0.17206` | `1.12490` | `3.97803` |
| `axial_torque_flux_probe` | `1.5e-20` | `(1,1.2,1,1)` | false | `0.85969` | `1.00000` | `0.39961` | `0.17199` | `0.79433` | `0.89821` |
| `axial_torque_flux_probe` | `2.0e-20` | `(1,1,1,1)` | true | `1.08547` | `1.00000` | `0.40020` | `0.17206` | `1.12497` | `3.97803` |

## 判断

P2-6-007では、`distributed_flagellum` は diagnostic upper-bound として残し、最終的な triplet 系改善の評価基準に使う。現時点で `distributed_flagellum` を物理モデルの本解として採用しない。

`axial_torque_flux_probe` は、現行 bead-position-only model のままでも torque flux 近似が単一べん毛の螺旋 net 回転に有効であることを示した。0.5 s で net 1回転以上と shape gate PASS を満たすため、短期的な有効アプローチである。

ただし、`axial_torque_flux_probe` は material frame / segment twist を明示しない。離れた bead へ直接接線力を分布させるため、root motor torque が segment twist として保存・輸送される物理モデルではない。したがって、論文モデルへの忠実性を高める本命候補は引き続き material frame / segment twist の導入である。

本PRでの優先判断は以下とする。

- なぜ伝わらないか: 現行 `triplet` は局所3点の力カップルであり、bead-position-only model には root torque を軸まわり twist として保存・輸送する状態量がないため。
- どうすれば伝わるか: flagellum beads に螺旋主軸まわりの接線速度を作る必要がある。物理的には material frame / segment twist を導入する。妥協案では root起点の `axial_torque_flux_probe` で軸方向torque fluxを近似する。
- 最終妥協案候補: `distributed_flagellum` ではなく `axial_torque_flux_probe` を優先する。`distributed_flagellum` は上限診断として残す。

quasi-rigid helical body approximation は、将来 ML 用データ生成で必要になった場合の実用近似としては残すが、Phase 2.6 の triplet torque transmission 修正の主案にはしない。

## 評価方法

`triplet` 系の改善は、少なくとも以下で評価する。

- `time.dt_star=1.0e-4`
- `duration_s >= 0.5`
- `net_abs_flag_helix_spin_revolutions >= 1.0`
- `flag_helix_spin_direction_consistency >= 0.5`
- `helix_to_root_net_rotation_ratio` が `triplet + hook spring fix` baseline より改善する
- `hook_len_rel_err_max <= 0.5`
- `flag_bond_rel_err_max <= 0.25`
- `flag_bend_err_max_deg <= 30`
- `flag_torsion_err_max_deg <= 60`

`helix_to_root_net_rotation_ratio` を比較に使うため、helix retention gate は `flag_phase_deg` と `flag_phase_rate_hz` を必須列として扱う。root phase が欠損している出力は、torque transmission の評価対象にしない。

## 有効なアプローチ

1. 短期: `axial_torque_flux_probe` を Phase 2.6 の実験用 candidate として残し、`triplet` baseline と `distributed_flagellum` upper-bound の間に置いて比較する。
2. 中期: material frame / segment twist を明示する本実装へ進む。`axial_torque_flux_probe` の pass 条件を、material frame 実装の最低限の定量目標にする。
3. 採用判断: 複数べん毛や束化へ進む前に、単一べん毛で `dt_star=1.0e-4`, 0.5 s, net 1回転以上、shape gate PASS を維持できることを条件にする。

local twist transmission の実装方針は `docs/adr/0003_phase2_local_twist_transmission.md` に記録する。まずはsegmentごとのorientation/local_twistを状態として持ち、root torqueでroot側のlocal_twistが増え、twist potentialで先端側へ伝わることを確認する。その後、local twist torqueをbead forceへ変換して、螺旋全体のnet回転に接続する。
