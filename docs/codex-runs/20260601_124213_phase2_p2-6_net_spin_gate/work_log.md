# Phase 2.6 net spin gate correction

## 目的

ユーザー目視で `outputs/phase2_6_dt1e4_spin_review/2026-05-31/230124/render/swim3d.mp4` がほぼ回転していないと確認されたため、旧 gate がなぜ pass と判定したかを切り分け、Phase 2.6 のタスクを再分解する。

## 実施内容

- `step_summary.csv` を再評価し、`flag_helix_spin_phase_deg` の累積差が 0.25 s で 0.501 deg しかないことを確認した。
- `median(abs(flag_helix_spin_rate_hz))` は 24.40 Hz だが、net 回転は `0.00139` 回転であり、瞬間速度 gate がフィット jitter / 往復揺れを誤検出していたと判断した。
- `summarize_single_flagellum_helix_retention()` に `net_abs_flag_helix_spin_revolutions`, `signed_flag_helix_spin_rate_hz`, `flag_helix_spin_direction_consistency` を追加した。
- 旧 pass 代表を `motor_no_rotation` fail として pytest で固定した。
- `dt_star=1.0e-4`, `duration_s=0.05` の短時間スクリーニングを行い、net 回転が大きい条件はすべて `shape_pass_nonbody=False` になることを確認した。
- `docs/PROJECT_PLAN.md`, `docs/phase2/phase2_tasks.md`, `docs/phase2/phase2_6_helix_retention_gate.md` に原因と次タスクを記録した。

## 主要結果

- `torque=2.5e-20`, `dt_star=1.0e-4`, `local scale=(4,2,2,2)`, `duration_s=0.25`:
  - `median_abs_flag_helix_spin_rate_hz=24.400110847247923`
  - `net_abs_flag_helix_spin_revolutions=0.0013917990878730456`
  - `flag_helix_spin_direction_consistency=0.00022068721845754322`
  - 新 gate では `motor_no_rotation` fail
- 0.05 s スクリーニング:
  - net 回転が大きい条件は shape fail を伴う。
  - 形状維持と持続的な net 回転を両立する代表条件は未確立。

## 次アクション

- motor torque が root 方位の揺れや局所変形ではなく、螺旋全体の累積回転へ伝達されるかを診断する。
- 必要なら motor torque distribution または hook/root coupling の実装変更を ADR とともに行う。
- `time.dt_star=1.0e-4` を実行時 override として維持し、修正後に net 1回転以上の代表条件を再探索する。
