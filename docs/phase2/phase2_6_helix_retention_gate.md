# Phase 2.6 single flagellum 螺旋維持 gate

## 目的

Phase 2.6 の目的は、**単一 full flagellum が motor-on 条件で螺旋形状を保ちつつ、螺旋そのものが安定的に回転して見える条件を定量 gate とユーザー目視で確認すること**である。

Phase 2.5 の高トルク break representative (`4.0e-21 N m`) では、body と hook ではなく flagellum chain の bond / bend / torsion が先に gate を超えた。Phase 2.6 では、この破綻を再現しつつ、トルクを上げたときに螺旋形状維持と螺旋スピンを両立できる条件を固定する。

本タスクでいう「長時間」は、CI で固定する 400-500 step representative と、ローカル probe で確認する 2500-step (`duration_s=0.25`, `dt_star=1.0e-4`) representative を指す。多本べん毛、後方束化、遊泳軌跡の自然さ、2D擬似顕微鏡動画としての利用可能性は Phase 2.7 以降の対象である。ただし、単一べん毛が定性的にも崩壊せず回転して見えるかは Phase 2.6 のユーザー目視レビュー対象とする。

## 仮説

単に `dt_star` を小さくすると、bond / bend / torsion の数値的な発散は抑えられる。しかし、従来の `flag_phase_rate_hz` は root azimuth 由来であり、目視の「螺旋がスピンしている」と一致しない場合がある。実際に `4.0e-21 N m`, `dt_star=2.5e-4`, `local_hook_scale=8`, `local_spring_scale=5`, `local_bend_scale=8`, `local_torsion_scale=4` は、形状は保つがユーザー目視では回転していなかった。

さらに、`flag_helix_spin_rate_hz` だけでも十分ではない。`outputs/phase2_6_dt1e4_spin_review/2026-05-31/230124/` では `median_abs_flag_helix_spin_rate_hz=24.40` だったが、`flag_helix_spin_phase_deg` の累積差は 0.25 s で 0.501 deg (`net_abs_flag_helix_spin_revolutions=0.00139`) しかなく、ユーザー目視でもほぼ回転していなかった。

そのため Phase 2.6 の回転判定は、body追従座標系で推定した螺旋位相の瞬間速度ではなく、`flag_helix_spin_phase_deg` の累積差から計算する net 回転数と方向一貫性を正本にする。現時点では、現行 motor 実装が root 方位の揺れや局所変形を生んでも、螺旋全体の持続回転へ十分に伝達していない可能性が高い。

## Gate 指標

`step_summary.csv` を入力に、`src/sim_swim/sim/helix_retention_gate.py` で判定する。

body diagnostics は現行実装では長時間条件で常に出力されないため、Phase 2.6 gate の正本は `step_summary.csv` の non-body / motor / flagellum 診断とする。

| 指標 | 役割 |
| --- | --- |
| `finite_pass` | 数値発散・NaN/Inf の検出 |
| `shape_pass_nonbody` | hook / flagellum 側の broad shape gate |
| `first_fail_category_nonbody` | non-body 側の first-fail category |
| `motor_degenerate_axis_count` | motor split の退化軸検出 |
| `motor_bond_length_clipped_count` | motor force split の basal link clipping 検出 |
| `flag_phase_rate_hz` | root azimuth 由来の旧回転 proxy。目視回転とは一致しない場合があるため参考値 |
| `flag_helix_spin_rate_hz` | body追従座標系で推定した螺旋位相の瞬間スピン速度。jitter を拾うため参考値 |
| `net_abs_flag_root_revolutions` | `flag_phase_deg` の開始・終了差から計算した root 方位の累積 net 回転数 |
| `net_abs_flag_helix_spin_revolutions` | `flag_helix_spin_phase_deg` の開始・終了差から計算した累積 net 回転数 |
| `signed_flag_helix_spin_rate_hz` | 累積 net 回転数を評価時間で割った符号付き平均回転速度 |
| `flag_helix_spin_direction_consistency` | 総位相変動に対する net 位相変動の比。往復揺れを持続回転から除外する |
| `helix_to_root_net_rotation_ratio` | root 方位の net 回転が螺旋全体の net 回転へ伝わった比。低い場合は root 近傍の運動が螺旋全体へ伝達していない |
| `flag_helix_spin_fit_r2` | 螺旋位相フィットの品質 |
| `hook_len_rel_err_max` | hook link 長の相対誤差 |
| `local_attach_first_rel_err` | attach-first hook link の相対誤差 |
| `flag_bond_rel_err_max` | flagellum bond 長の最大相対誤差 |
| `flag_bend_err_max_deg` | flagellum bend 角の最大誤差 |
| `flag_torsion_err_max_deg` | flagellum torsion 角の最大誤差 |

Phase 2.6 の hard gate では、Phase 2.5 の broad non-body gate より厳しい運用閾値を使う。

| 指標 | hard limit |
| --- | ---: |
| `hook_len_rel_err_max` | `<= 0.5` |
| `local_attach_first_rel_err` | `<= 0.5` |
| `flag_bond_rel_err_max` | `<= 0.25` |
| `flag_bend_err_max_deg` | `<= 30 deg` |
| `flag_torsion_err_max_deg` | `<= 60 deg` |
| `median(abs(flag_helix_spin_rate_hz))` | `>= 1.0` |
| `net_abs_flag_helix_spin_revolutions` | `>= 1.0` |
| `flag_helix_spin_direction_consistency` | `>= 0.5` |
| `min(flag_helix_spin_fit_r2)` | `>= 0.5` |

## 代表条件

共通条件:

- `n_flagella=1`
- `stub_mode=full_flagella`
- 決定論的条件、Brownian off

代表条件と解釈:

| 条件 | 期待結果 | 解釈 |
| --- | --- | --- |
| `torque=4.0e-21`, `dt_star=1.0e-3`, default local scales | fail (`flag`) | Phase 2.5 break representative を再現 |
| `torque=4.0e-21`, `dt_star=2.5e-4`, `local_hook_scale=8`, `local_spring_scale=5`, `local_bend_scale=8`, `local_torsion_scale=4` | fail (`motor_no_rotation`) | 形状は保つが螺旋スピンがほぼ出ない。ユーザー目視で「全く回転していない」と判定された条件 |
| `torque=8.0e-21`, `dt_star=1.25e-4`, `local_hook_scale=8`, `local_spring_scale=5`, `local_bend_scale=8`, `local_torsion_scale=4` | fail (`motor_no_rotation`) | 瞬間スピン rate は出るが、0.05 s / 400 steps で `net_abs_flag_helix_spin_revolutions=0.00090` |
| `torque=2.5e-20`, `dt_star=1.0e-4`, `local_hook_scale=4`, `local_spring_scale=2`, `local_bend_scale=2`, `local_torsion_scale=2` | fail (`motor_no_rotation`) | 0.25 s / 2500 steps で `median_abs_flag_helix_spin_rate_hz=24.40` だが `net_abs_flag_helix_spin_revolutions=0.00139`。目視でもほぼ回転しない |
| `torque=3.0e-20`, `dt_star=1.0e-4`, `local_hook_scale=4`, `local_spring_scale=2`, `local_bend_scale=2`, `local_torsion_scale=2` | fail (`flag`) | 0.10 s で `flag_bond_rel_err_max=0.314 > 0.25` |
| `torque=2.0e-20`, `dt_star=1.0e-4`, `motor.force_distribution=distributed_flagellum`, `local_spring_scale=1.2`, other local scales `1.0` | pass | hook spring force 復元後、0.5 s / 5000 steps で `net_abs_flag_helix_spin_revolutions=1.0416`, direction consistency `0.9988`, hook rel err `0.4046` |

これらの local scale は、現時点では物理最適値ではなく、Phase 2.6 の診断 representative である。参照論文モデルとの差分としては、motor-on 時の hook / spring / bend / torsion の局所補強であり、数値安定化寄りの運用パラメータとして扱う。標準 `conf/sim_swim.yaml` は local scale を `1.0` にしているため、再現する場合は CLI override で明示する。

`distributed_flagellum` は、従来の `triplet` 方式で root 近傍の揺れに消えていた motor torque を、螺旋全体の主軸回りへ分布させる診断用近似である。既定値は `triplet` のままであり、代表条件では CLI override で明示する。詳細は `docs/adr/0001_phase2_distributed_flagellar_torque.md` に記録する。

hook spring force 復元後に `triplet` を再評価した結果、root 方位と螺旋全体の net 回転はなお分離していた。例として、`torque=2.5e-20`, `dt_star=1.0e-4`, `local_scale=(4,2,2,2)` では 0.25 s で `net_abs_flag_root_revolutions=0.8796` だが、`net_abs_flag_helix_spin_revolutions=0.00139`, `helix_to_root_net_rotation_ratio=0.00158` だった。これは root 近傍の方位運動が材料ねじれとして flagellum chain へ伝搬していないことを示す。

## 標準 sweep コマンド

```bash
uv run python -m scripts.01_simulate_swimming.run_motor_scale_sweep \
  --target local_bend_scale \
  --values 2,4,8 \
  --torques 8e-21,1.2e-20,2e-20,3e-20 \
  --duration 0.05 \
  --dt-star 1.25e-4 \
  --stub-mode full_flagella \
  --output-dir outputs/phase2_6_helix_retention
```

CI では runtime を抑えるため、`duration_s=0.05`, `dt_star=1.0e-4` の 500-step 代表条件で「瞬間スピン rate が高くても net 回転がなければ fail」となることと、`distributed_flagellum` が短時間 net 回転を作ることを固定する。ローカル確認では、`torque=2.0e-20 N m`, `dt_star=1.0e-4`, `motor.force_distribution=distributed_flagellum`, `local_spring_scale=1.2` 条件で `duration_s=0.5` の 5000-step probe が pass することを確認した。

したがって、単一べん毛については、定量 gate の範囲で「5000-step までの螺旋形状維持と net 1回転以上」を確認済みである。ただし、これは単一べん毛・決定論的条件・診断用 motor 近似上の確認であり、ユーザー目視による定性検証は別途確認する。

## 定性レビュー

Phase 2.6 は、ユーザー目視レビューが完了するまで完了扱いにしない。

代表条件の再実行コマンド:

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

確認対象:

- `outputs/phase2_6_distributed_spin_review/2026-06-01/134726/render/swim3d.mp4`
- `outputs/phase2_6_distributed_spin_review/2026-06-01/134726/render/swim3d_final.png`
- `outputs/phase2_6_distributed_spin_review/2026-06-01/134726/render2d/projection.mp4`
- `outputs/phase2_6_distributed_spin_review/2026-06-01/134726/sim/step_summary.csv`

前回の目視結果:

- ユーザー確認では、`swim3d.mp4` はほとんど回転していない。
- 新 gate でも `net_abs_flag_helix_spin_revolutions=0.00139` のため fail。

新代表条件の定量結果:

- `helix_retention_pass=True`
- `net_abs_flag_helix_spin_revolutions=1.0416129331359716`
- `flag_helix_spin_direction_consistency=0.9988303405443211`
- `median_abs_flag_helix_spin_rate_hz=2.142334537325277`
- `max_hook_len_rel_err=0.4045704218609963`
- `max_flag_bond_rel_err=0.17236812212772734`
- `max_flag_bend_err_deg=0.996803008165139`
- `max_flag_torsion_err_deg=3.038765579091419`

今後の目視観点:

- 単一べん毛の螺旋形状が 0.25 s / 2500 steps の間に潰れない。
- hook 近傍で不自然な折れ曲がりや飛びが出ない。
- べん毛が停止せず、連続的に回転して見える。
- 菌体やべん毛が画面外へ飛ぶ、または急激に拡大・縮小して見えない。

## 検証

自動検証:

- `tests/test_run_state_fixed.py::test_phase26_default_break_fails_helix_retention_gate`
- `tests/test_run_state_fixed.py::test_phase26_small_dt_without_bend_scale_loses_rotation`
- `tests/test_run_state_fixed.py::test_phase26_small_dt_and_bend_scale_retain_helix`
- `tests/test_run_state_fixed.py::test_phase26_higher_torque_jitter_does_not_count_as_net_spin`
- `tests/test_helix_retention_gate.py`
- `tests/test_motor_scale_sweep.py`

Phase 2.6 は single flagellum の multi-step 定量 gate 固定であり、多本べん毛の束化や 2D 動画自然さは対象外である。ただし、単一べん毛の回転が定性的に安定して見えるかはユーザー目視レビューで確認する。
