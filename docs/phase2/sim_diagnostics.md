# Phase2 Simulation Diagnostics Reference

## 概要

このドキュメントは、シミュレーション診断出力の列説明を管理するための参照資料です。

各フェーズで重要な診断指標をまとめ、破綻検出・原因分析・段階的ゲート設定に用いられます。

**出力対象**:
- `step_summary.csv` (main diagnostics)
- `body_constraint_diagnostics.csv` (body-only mode)
- `body_constraint_local_diagnostics.csv` (body spring/bend details)

---

## フェーズ別重要指標マトリクス

| 指標 | Phase0a | Phase0b | Phase1 | Phase2 | Phase3 |
|------|---------|---------|--------|--------|--------|
| `pos_all_finite` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `any_nan`, `any_inf` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `F_total_mean_all` | ✓ | ✓ | ✓ | ✓ | ✓ |
| `local_attach_first_rel_err` | - | ✓ | - | ✓ | ✓ |
| `local_first_second_rel_err` | - | ✓ | - | ✓ | ✓ |
| `motor_*` | - | - | - | ✓ | ✓ |
| `flag_bond_rel_err_max` | - | - | - | ✓ | ✓ |
| `flag_bend_err_max_deg` | - | - | - | △ | ✓ |
| `flag_torsion_err_max_deg` | - | - | - | - | ✓ |

✓: 主要監視対象, △: 補助監視, -: 不適用 (phase該当条件で出力されない可能性)

---

## 主要診断列の詳細

### A. 基本有限性チェック

#### `pos_all_finite`
- **意味**: 全 bead の座標が有限値（±∞ でない）か
- **値**: boolean (0/1)
- **用途**: NaN/Inf 発生の最初の兆候; すべてのフェーズで必須チェック
- **破綻時の見方**: 1 → 0 に遷移した step が first fail

#### `any_nan`
- **意味**: いずれかの計算中間変数に NaN が発生したか
- **値**: boolean (0/1)
- **用途**: 数値計算の発散（0割など）を検出
- **破綻時の見方**: 1 → 1 のままなら計算崩壊；diagnostics が以後出力されない可能性

#### `any_inf`
- **意味**: いずれかの計算中間変数に ±∞ が発生したか
- **値**: boolean (0/1)
- **用途**: 力が発散したことの指標（spring stretch など）
- **破綻時の見方**: 1 が立つと通常は `pos_all_finite=0` に引き続く

---

### B. 全体力統計

#### `F_total_mean_body`
- **意味**: body 上の全 bead に作用する全力の平均ノルム (N)
- **値**: float (N)
- **用途**: body 側の力バランス診断; Phase0a での基準値設定に重要
- **正常範囲**: Phase0a 静的では μN オーダー; motor ON で數 pN に増加

#### `F_total_mean_flag`
- **意味**: flagellum 上の全 bead に作用する全力の平均ノルム (N)
- **値**: float (N)
- **用途**: flagellum 側の力バランス診断; Phase2/3 で basal vs distal の差を見る
- **正常範囲**: Phase1 で μN 以下; Phase2/3 では flagellum 長に応じて増加可能

#### `F_total_mean_all`
- **意味**: 全 bead の全力平均ノルム (N)
- **値**: float (N)
- **用途**: 全体的な力バランス監視; hard gate の一次指標候補
- **破綻時の見方**: 値が指数関数的に増加する場合は spring stretch または motor 破綻の兆候

---

### C. Basal / Hook 局所診断

これらは `n_flagella≥1` で、flagellum が少なくとも 3 bead を含む場合に出力されます。

#### `local_attach_first_rel_err`
- **意味**: body attach bead と flagellum first bead 間の spring 相対誤差 (dimensionless)
- **定義**: `(actual_dist - rest_length) / rest_length`
- **値**: float
- **用途**: basal link の stretch 度; **Phase0b/2 での一次 failure 指標**
- **正常範囲**: Phase0b 静的では <10% (typical: <5%); Phase2 motor ON では up to 100%+ 可能
- **破綻時の見方**: >100% は basal link の大幅な拉伸; >200% は近い first-fail

#### `local_first_second_rel_err`
- **意味**: flagellum first と second bead 間の spring 相対誤差
- **定義**: `(actual_dist - rest_length) / rest_length`
- **値**: float
- **用途**: basal 直後の intra-flagellar link; attach-first と組で basal region 診断
- **正常範囲**: Phase0b では <10%; Phase2 では first と並行監視
- **破綻時の見方**: attach-first より先に悪化する場合は basal bend 強度不足の可能性

#### `local_second_third_rel_err`
- **意味**: flagellum second と third bead 間の spring 相対誤差
- **定義**: `(actual_dist - rest_length) / rest_length`
- **値**: float
- **用途**: basal region より少し下流の link; trend 監視用
- **正常範囲**: Phase0b では <10%; Phase2 では attach-first・first-second と三点監視
- **破綻時の見方**: attach/first/second より後発して悪化する場合は flagellum 内の連鎖崩れの指標

---

### D. Basal Bending & Torsion

#### `local_basal_bend_err_deg`
- **意味**: body attachment 点での flagellum bend angle の誤差 (degrees)
- **定義**: `actual_angle_deg - theta0_target_deg`
- **値**: float (degrees)
- **用途**: basal region の bending 状態; 90° 超過で hook repulsion が off
- **正常範囲**: Phase0a で 0; Phase0b-2 では ±30° 以内が目安
- **破綻時の見方**: >90° で bend potential off (flagellar torsion dominant); >120° は基本破綻

#### `local_first_torsion_err_deg`
- **意味**: flagellum first bead の torsion angle 誤差 (degrees)
- **定義**: `actual_torsion_deg - phi0_target_deg`
- **値**: float (degrees) or `nan`
- **用途**: basal torsion の捻れ保持; Phase2/3 重要
- **正常範囲**: Phase0b では出力されない可能性（2 bead stub で torsion 未定義）; Phase2 では ±30° 以内
- **破綻時の見方**: NaN は torsion 定義不可（bead <4); 値が大きく累積する場合は torsion stiffness 不足

---

### E. Flagellum 全体統計

#### `flag_intra_count`
- **意味**: flagellum 内部の spring 数 (body attach 除外)
- **値**: int
- **用途**: flagellum 構造が正しく構築されたかの確認; phase 定義確認
- **典型値**: Phase0b `minimal_basal_stub` で 1 (first-second link のみ); Phase3 `full_flagella` で 10

#### `flag_bond_rel_err_max`
- **意味**: 全 intra-flagellar bond の相対誤差の最大値
- **定義**: `max((actual_dist - rest_length) / rest_length) over all flag bonds`
- **値**: float
- **用途**: flagellum 全体での最悪ケース spring stretch; Phase2/3 で early warning
- **正常範囲**: Phase0b では <50%; Phase2 では <150% 目安; >300% は重度破綻

#### `flag_bend_err_max_deg`
- **意味**: 全 bend angle の誤差の最大値 (degrees)
- **定義**: `max(|angle_error_deg|) over all triplets within flagellum`
- **値**: float (degrees)
- **用途**: flagellum 全体での helix 保持度; Phase3 重要
- **正常範囲**: Phase3 では <45° が目安
- **破綻時の見方**: >90° は helix collapse の兆候

#### `flag_torsion_err_max_deg`
- **意味**: 全 torsion angle の誤差の最大値 (degrees)
- **定義**: `max(|torsion_error_deg|) over all quads within flagellum`
- **値**: float (degrees) or `nan`
- **用途**: flagellum 全体での torsion 保持度; Phase3 重要
- **正常範囲**: Phase3 では <45° が目安
- **破綻時の見方**: NaN は torsion 定義不可; >90° は counter-rotating や twist-up の兆候

---

### F. Motor 関連診断

これらは `motor.torque_Nm != 0` (Phase2 以降) で出力されます。

#### `motor_ta_dot_ra_abs`
- **意味**: motor torque が attach bead の lever arm と平行か（内積の絶対値）
- **値**: float (0-1 scale, 理想的には接近 1)
- **用途**: motor torque decomposition の数値安定性確認; rank deficiency detector
- **正常範囲**: >0.9 が目安
- **破綻時の見方**: <0.5 は motor axis が lever arm に垂直に近い（特異性）

#### `motor_tb_dot_rb_abs`
- **意味**: motor torque が first-second lever arm と平行か
- **値**: float (0-1 scale)
- **用途**: motor torque の二次 decomposition 安定性
- **正常範囲**: >0.8 が目安

#### `motor_split_residual_norm`
- **意味**: motor torque を 3 point forces に分解したときの residual ノルム
- **値**: float (Nm)
- **用途**: decomposition の成功度; ゼロに近いほど正確
- **正常範囲**: <1e-20
- **破綻時の見方**: >1e-15 は motor axis が特異に近い；>1e-10 は実装エラーの可能性

#### `motor_degenerate_axis_count`
- **意味**: motor axis が特異（degenerate）と判定された回数（累積）
- **値**: int
- **用途**: motor decomposition 失敗の統計
- **正常値**: 0

#### `motor_split_rank_deficient_count`
- **意味**: 3-point decomposition が rank deficient と判定された回数（累積）
- **値**: int
- **用途**: motor force split の数値安定性監視
- **正常値**: 0

---

### G. Body Equivalent Load (Phase1 のための Surrogate Torque)

#### `body_equiv_load_mode`
- **意味**: body-only surrogate torque の mode
- **値**: string ("off" / "pure_couple" / ...)
- **用途**: Phase1 実行条件を diagnostics から確認
- **正常値**: Phase1 では "pure_couple"; 他は "off"

#### `body_equiv_load_target_torque_Nm`
- **意味**: 設定された surrogate torque の target 値 (Nm)
- **値**: float (Nm)
- **用途**: surrogate torque 入力条件の記録
- **Phase1 典型値**: 1e-20 オーダー

#### `body_equiv_load_target_force_N`
- **意味**: 設定された surrogate force の target 値 (N)
- **値**: float (N)
- **用途**: Phase1 で torque 以外の等価荷重を使う場合の記録
- **Phase1 典型値**: 0.0

---

### H. Projection State (Diagnostics 用)

#### `body_spring_max_stretch_ratio` (in body_constraint_diagnostics.csv)
- **意味**: body 内 spring の最大 stretch ratio (projected / non-projected)
- **値**: float (ratio)
- **用途**: constraint projection の効果を diagnostics（非修正条件）で観察
- **正常範囲**: 1.0 に近い（stretch なし）

---

## 破綻の読み取り方

### Scenario 1: Immediate Divergence (Step 1-2)
```
pos_all_finite: 1 → 0
any_nan: 0 → 1 or 0 → 0
any_inf: 0 → 1
```
- **原因候補**: 初期化エラー、極端な parameter 設定
- **対応**: 初期 config 検査

### Scenario 2: Gradual Spring Stretch (Phase2 Motor ON 典型)
```
local_attach_first_rel_err: 0% → 50% → 100% → (diverges)
flag_bond_rel_err_max: 0% → 30% → 100% → 200% → (diverges)
F_total_mean_all: μN × 1 → μN × 10 → μN × 100 → (inf)
pos_all_finite: 1 → 1 → 1 → 0
```
- **原因候補**: spring stiffness 不足、motor torque 過大
- **対応**: k_spring / stiffness_scale 調整、motor torque 減尐

### Scenario 3: Helix Collapse (Phase3 典型)
```
flag_bend_err_max_deg: 0° → 20° → 50° → 90° → (diverges)
flag_torsion_err_max_deg: 0° → 15° → 40° → nan
F_total_mean_flag: μN × 1 → μN × 5 → (inf)
```
- **原因候補**: helix bending/torsion stiffness 不足
- **対応**: θ, φ stiffness 調整

### Scenario 4: Motor Decomposition Failure
```
motor_split_residual_norm: 1e-20 → 1e-15 → 1e-10 → (inf)
motor_degenerate_axis_count: 0 → 1 → N (cumulative)
motor_split_rank_deficient_count: 0 → 1 → N
```
- **原因候補**: motor axis geometry 特異、motor torque 入力エラー
- **対応**: motor attach point 位置確認、torque direction 確認

---

## 診断 CSV の使い方

### パイソン例：First-fail step を見つける
```python
import pandas as pd

df = pd.read_csv("step_summary.csv")
first_fail = df[df['pos_all_finite'] == 0].iloc[0] if (~df['pos_all_finite']).any() else None
if first_fail is not None:
    print(f"First fail at step {first_fail['step']}, t_s={first_fail['t_s']}")
    print(f"  local_attach_first_rel_err: {first_fail['local_attach_first_rel_err']}")
    print(f"  flag_bond_rel_err_max: {first_fail['flag_bond_rel_err_max']}")
    print(f"  F_total_mean_all: {first_fail['F_total_mean_all']}")
```

### トレンド監視
```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(3, 1, figsize=(10, 8))
axes[0].semilogy(df['t_s'], df['F_total_mean_all'], label='F_total_mean_all')
axes[0].set_ylabel('Force (N)')
axes[0].legend()

axes[1].plot(df['t_s'], df['local_attach_first_rel_err'], label='attach-first');
axes[1].plot(df['t_s'], df['local_first_second_rel_err'], label='first-second');
axes[1].set_ylabel('Rel Err (ratio)')
axes[1].legend()

axes[2].plot(df['t_s'], df['flag_bend_err_max_deg'], label='bend');
axes[2].plot(df['t_s'], df['flag_torsion_err_max_deg'], label='torsion');
axes[2].set_ylabel('Angle Error (deg)')
axes[2].set_xlabel('Time (s)')
axes[2].legend()

plt.tight_layout()
plt.savefig("diagnostics_trend.pdf")
```

---

## 更新履歴

| Date | Author | Change |
|------|--------|--------|
| 2026-04-08 | Codex | Issue #37: Initial phase diagnostics reference |

---

## 補足: Diagnostics 列追加時のチェックリスト

新規列を追加した場合は、以下のセクションに追記してください。

- [ ] 列名
- [ ] 意味（数学的定義）
- [ ] 出力値のデータ型・単位
- [ ] 用途（どのフェーズで重要か）
- [ ] 正常範囲（典型値）
- [ ] 破綻時の見方（異常値の解釈）

