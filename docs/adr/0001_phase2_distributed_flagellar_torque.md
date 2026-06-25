# ADR 0001: Phase 2 distributed flagellar torque drive

## Status

Accepted for Phase 2 diagnostic use.

## Naming Update

Issue #88 で旧 `distributed_flagellum` は `root_torque_axis_projection` の deprecated alias になった。このADR内の `distributed_flagellum` は当時の実装名として読む。新規設定・コマンドでは `motor.force_distribution=root_torque_axis_projection` を使う。

## Context

Phase 2.6 では、単一 full flagellum が螺旋形状を保ちつつ、べん毛全体として継続的に回転する条件を探している。

従来の motor force は `attach-first-second` の3点に力カップルを与える `triplet` 方式だった。この方式では root 方位や螺旋位相フィットの瞬間値は揺れるが、螺旋全体の累積 net 回転がほとんど進まない条件があった。`outputs/phase2_6_dt1e4_spin_review/2026-05-31/230124/` では `median_abs_flag_helix_spin_rate_hz=24.40` だったが、0.25 s の net 回転は 0.00139 回転しかなかった。

また、実装上 `attach-first` の hook spring 行を分類していたにもかかわらず、時間積分で hook spring force を計算していなかった。これにより、分布型 torque で net 回転を作ると長時間で hook link が伸びやすかった。

## Decision

Phase 2.6 の診断用に `motor.force_distribution=distributed_flagellum` を追加する。

この方式では、各 flagellum 点群の主軸まわりにゼロ合力の接線力を分布させ、同じ軸まわりの反作用トルクを body 側へ与える。これにより、root 近傍だけではなく、螺旋全体へ motor torque を伝える。

このADR作成時点では既定値を `triplet` のままとし、`distributed_flagellum` は明示的な実行時 override でのみ使う判断とした。

P2-6-008後の既定値は `material_twist_local_couple` である。`distributed_flagellum` は引き続き診断用modeであり、使用時は明示的に指定する。

あわせて、hook spring force を時間積分に復元する。これは新しい物理モデルではなく、既存 topology に存在する hook spring を実際の force 計算へ含めるバグ修正として扱う。

## Consequences

- `triplet` と `distributed_flagellum` を比較できる。
- `distributed_flagellum` は参照論文の厳密な motor 実装ではなく、Phase 2 の診断・安定化近似である。
- Phase 2.6 の代表条件は、`time.dt_star=1.0e-4` を実行時 override とし、`motor.force_distribution=distributed_flagellum` も明示する。
- downstream の多本べん毛・束化・遊泳軌跡へ進む前に、この近似が必要十分かを Phase 2.7 以降で再評価する。
