# Phase 2.6 triplet motor のねじれ・回転自由度設計メモ

## 目的

`distributed_flagellum` を最終的な物理モデルとして扱わず、従来の `triplet` motor と hook spring fix を中心に、root torque が flagellum chain へ伝わらない原因を明確にする。

本メモでは、現行 bead-position-only model に不足している状態量を整理し、次に実装すべき torque transmission model の候補を比較する。

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

### material frame

各 flagellum segment に付随する局所座標系で、segment 接線方向だけでなく断面の向きを表す。現行モデルは bead positions だけを状態として持つため、segment の「軸まわり姿勢」を直接保持しない。

### segment twist

隣接 segment の material frame 同士が、segment 軸まわりにどれだけ相対回転しているかを表すねじれ量。root torque を弾性的な torsional deformation として chain に蓄える・伝えるために必要である。

現行の torsion force は4点 dihedral angle の位置エネルギー勾配であり、形状のねじれを保つ力としては働く。しかし、segment 自身が軸まわりにどの向きを持つかを独立状態として持たないため、root からの motor torque を material twist として保存・輸送する経路にはなっていない。

### axial torque flux

root から先端側へ、segment 軸方向に伝わる torque の流れ。現行 `triplet` は root 近傍に局所的な力カップルを入れるが、各 segment 間で torque flux を保存・減衰・伝搬させる明示的な状態や計算はない。

## 設計候補

| 候補 | 内容 | 論文モデルとの距離 | 利点 | リスク | P2-6での扱い |
| --- | --- | --- | --- | --- | --- |
| material frame + segment twist | 各 segment に局所 frame または twist angle を持たせ、root torque を torsional spring/damping として隣接 segment へ伝える | 中。bead-spring model を拡張するが、軸まわり回転自由度として物理的に明確 | torque の保存・伝搬を説明しやすい。`triplet` の root 駆動に近い | 状態変数、積分器、安定性 gate、ADR が必要 | 次の本命候補 |
| segment torsional torque flux | bead positions は維持しつつ、隣接 segment に軸方向 torque flux を近似分配する | 中-大。状態量を省略した有効モデル | 実装量は material frame より小さい可能性がある | 何を保存しているかが曖昧になりやすい。force couple への変換で数値不安定化しやすい | material frame 案の比較対象 |
| quasi-rigid helical body approximation | 単一 flagellum 螺旋をほぼ剛体として root torque で全体回転させる | 大。ML用生成には実用的だが論文モデルから離れる | 形状維持と回転は作りやすい | collapse 診断や elastic failure の意味が薄くなる | Phase 2.6の主案にはしない |
| distributed_flagellum | flagellum 全体に主軸まわりのゼロ合力 torque drive を分布させる | 大。root motor ではなく診断用 upper-bound | torque が螺旋全体へ伝わった場合の見え方を確認できる | 離れた bead に直接 torque を入れるため、現実の root torque 伝搬ではない | mode として残すが最終解にしない |

## 判断

P2-6-007では、`distributed_flagellum` は diagnostic upper-bound として残し、最終的な triplet 系改善の評価基準に使う。現時点で `distributed_flagellum` を物理モデルの本解として採用しない。

次に実装する場合の第一候補は、material frame / segment twist を明示する方法である。これは現行 bead-position-only model への物理モデル拡張であり、実装前にADRを作成する。

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

## 次アクション

1. material frame / segment twist を入れる最小設計ADRを作成する。
2. 単一 flagellum の root segment から先端側へ twist state を伝える最小実装を検討する。
3. `distributed_flagellum` 代表条件を upper-bound とし、`triplet` 系改善案の `helix_to_root_net_rotation_ratio` と shape gate を比較する。
