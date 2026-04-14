# 2026-04-14 Motor Scale Sweep Report

## 目的
- body + hook 条件を基点に、`motor_local_hook_scale` を先に弱めたときの回転角と形状崩壊の出方を観測する。
- torque ramp は使わず、固定トルク条件でどこが先に崩れるかを step_summary.csv ベースで確認する。
- レポートは観測結果の整理を主目的とし、閾値ベースの合否判定は行わない。

## 実行条件
- Script: `scripts/run_motor_scale_sweep.py`
- Base motor torque: `4.0e-21 N·m`
- Flagella: `n_flagella=1`
- Stub: `minimal_basal_stub`
- torque ramp: disabled
- Sweep target: `motor.local_hook_scale`
- Sweep values:
  - Smoke test: `8, 0` at `duration_s=0.001`
  - Main observation: `8, 2, 0` at `duration_s=0.05`

## 要約
- どの sweep でも `pos_all_finite=True`, `any_nan=False`, `any_inf=False` を維持した。
- `motor_degenerate_axis_count`, `motor_split_rank_deficient_count`, `motor_bond_length_clipped_count` は全ケースで 0 のままだった。
- hook scale を 8 から 0 へ下げても、短時間では first-fail は発生しなかった。
- ただし、形状誤差は徐々に悪化し、特に `local_first_second_rel_err` と `flag_bond_rel_err_max` が先に増えた。
- 回転系では `flag_phase_deg` の unwrap 表現は振れたが、`flag_phase_rate_hz` は scale 低下とともに絶対値が大きくなる傾向を示した。

## Main observation results

### local_hook_scale = 8.0
- `pos_all_finite=True`
- `local_attach_first_rel_err = 6.592%`
- `local_first_second_rel_err = 0.899%`
- `hook_angle_err_max_deg = 12.693`
- `hook_len_rel_err_max = 6.592%`
- `flag_bond_rel_err_max = 2.216e-02`
- `flag_bend_err_max_deg = 1.344`
- `flag_phase_rate_hz = -43.144`

### local_hook_scale = 2.0
- `pos_all_finite=True`
- `local_attach_first_rel_err = 6.348%`
- `local_first_second_rel_err = 14.756%`
- `hook_angle_err_max_deg = 11.110`
- `hook_len_rel_err_max = 6.348%`
- `flag_bond_rel_err_max = 4.968e-01`
- `flag_bend_err_max_deg = 1.476`
- `flag_phase_rate_hz = -42.890`

### local_hook_scale = 0.0
- `pos_all_finite=True`
- `local_attach_first_rel_err = 9.007%`
- `local_first_second_rel_err = 153.899%`
- `hook_angle_err_max_deg = 13.975`
- `hook_len_rel_err_max = 9.007%`
- `flag_bond_rel_err_max = 1.539`
- `flag_bend_err_max_deg = 13.229`
- `flag_phase_rate_hz = 48.839`

## 分析

### 1. 破綻の先行点
- 今回の短時間 sweep では、明確な NaN/Inf や退化カウンタ増加は出なかった。
- ただし、崩壊に向かう最初の兆候は `local_first_second_rel_err` と `flag_bond_rel_err_max` に現れた。
- `local_attach_first_rel_err` も悪化するが、先に大きく振れるのは first-second link 側だった。

### 2. 回転角の見え方
- `flag_root_azimuth_deg` と `flag_phase_deg` は unwrap の影響を強く受けるため、単独の瞬時値よりも `flag_phase_rate_hz` と合わせて読むのが安全。
- `local_hook_scale` を下げると、回転速度はゼロ近傍に潰れるのではなく、むしろ符号や大きさが変化した。
- したがって、hook 補助を弱めることは「単純に回転を止める」よりも、「回転の向きと位相の安定性を乱す」方向に効いている。

### 3. 形状崩壊の解釈
- hook の補助を外しても、まず body / hook 接続が即座に破綻するわけではなかった。
- 一方で `flag_bond_rel_err_max` と `flag_bend_err_max_deg` が先に伸びており、べん毛側の維持が弱くなっている。
- このため、hook scale は主に basal 近傍の保形を支えるが、下流の flagella 保持にも間接的に効いていると考えられる。

### 4. 実行上の補足
- 10ms の smoke sweep でも 50ms の main observation でも、現時点では first-fail に到達しなかった。
- つまり、今回の設定では「発散境界」まではまだ十分に近づいていない。
- 次に掘るなら `local_spring_scale` または `local_bend_scale` / `local_torsion_scale` を順に落とすのが自然。

## 気づいた点
- `flag_phase_deg` は 360 度折り返しの影響を受けるので、今後の比較では `flag_phase_rate_hz` を主指標にした方が解釈しやすい。
- hook weakening 単独では、今回の時間長では「即発散」ではなく「徐々に形状が崩れる」タイプの挙動だった。
- body + hook のみの sweep としては、今回の結果は「観測開始点」として十分で、次は長手方向の spring / bend / torsion の順に弱める方がよい。

## 次の候補
1. `motor.local_spring_scale` を同じ条件で sweep する。
2. `motor.local_bend_scale` と `motor.local_torsion_scale` を分けて sweep し、flagella 側の崩れ方を比較する。
3. 既存 sweep helper を使って、`flag_phase_rate_hz` と `local_first_second_rel_err` の関係を表形式で並べる。
