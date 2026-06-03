# Phase 2.6 distributed flagellar torque representative

## 目的

`duration_s` が短くて回転が見えない可能性を踏まえ、`time.dt_star=1.0e-4` を維持したまま、螺旋形状を保ちつつ単一べん毛が net 1回転以上する代表条件を探索する。

## 実施内容

- 従来の `triplet` motor では 0.25 s でも net 回転がほぼ出ないことを確認した。
- `attach-first` hook spring が topology には存在するが時間積分で未使用だったため、hook spring force を復元した。
- `motor.force_distribution=distributed_flagellum` を追加し、螺旋全体の主軸まわりへ motor torque を分布させる診断用 drive を実装した。
- `docs/adr/0001_phase2_distributed_flagellar_torque.md` に参照論文モデルとの差分を記録した。
- `torque=2.0e-20`, `dt_star=1.0e-4`, `motor.force_distribution=distributed_flagellum`, `local_spring_scale=1.2` で 0.5 s / 5000 steps の代表動画を生成した。

## 代表条件

```bash
uv run python -m scripts.01_simulate_swimming \
  --config conf/sim_swim.yaml \
  --duration-s 0.5 \
  --fps-out 100 \
  --render-flagella \
  --render-flagella-2d \
  flagella.n_flagella=1 \
  flagella.stub_mode=full_flagella \
  motor.torque_Nm=2.0e-20 \
  motor.force_distribution=distributed_flagellum \
  motor.local_hook_scale=1 \
  motor.local_spring_scale=1.2 \
  motor.local_bend_scale=1 \
  motor.local_torsion_scale=1 \
  time.dt_star=1.0e-4 \
  output_sampling.out_all_steps_3d=false \
  render.save_frames_3d=false \
  render.save_frames_2d=false \
  output.base_dir=outputs/phase2_6_distributed_spin_review
```

## 結果

- output: `outputs/phase2_6_distributed_spin_review/2026-06-01/134726/`
- `helix_retention_pass=True`
- `step_count=4999`
- `net_abs_flag_helix_spin_revolutions=1.0416129331359716`
- `signed_flag_helix_spin_rate_hz=2.0840594900681317`
- `flag_helix_spin_direction_consistency=0.9988303405443211`
- `median_abs_flag_helix_spin_rate_hz=2.142334537325277`
- `min_flag_helix_spin_fit_r2=0.9969545404781638`
- `max_hook_len_rel_err=0.4045704218609963`
- `max_flag_bond_rel_err=0.17236812212772734`
- `max_flag_bend_err_deg=0.996803008165139`
- `max_flag_torsion_err_deg=3.038765579091419`

## 残件

ユーザー目視レビューが未完了。確認対象は以下。

- `outputs/phase2_6_distributed_spin_review/2026-06-01/134726/render/swim3d.mp4`
- `outputs/phase2_6_distributed_spin_review/2026-06-01/134726/render/swim3d_final.png`
- `outputs/phase2_6_distributed_spin_review/2026-06-01/134726/render2d/projection.mp4`
- `outputs/phase2_6_distributed_spin_review/2026-06-01/134726/sim/step_summary.csv`
