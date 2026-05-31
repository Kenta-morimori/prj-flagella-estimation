# Phase 2.6 single flagellum 螺旋維持 gate

## 目的

Phase 2.6 の目的は、**単一 full flagellum が motor-on 条件で回転 activity を保ったまま、multi-step で bond / bend / torsion と螺旋形状を維持できる条件を定量的に固定し、ユーザー目視で定性確認すること**である。

Phase 2.5 の高トルク break representative (`4.0e-21 N m`) では、body と hook が先に破綻するのではなく、flagellum chain の bond / bend / torsion が gate を超えた。Phase 2.6 では、この破綻を再現しつつ、回転 activity を失わずに multi-step で螺旋形状を維持できる条件を固定する。

本タスクでいう「長時間」は、CI で固定する 200-step representative と、ローカル probe で確認する 1000-step (`duration_s=0.25`) representative を指す。多本べん毛、後方束化、遊泳軌跡の自然さ、2D擬似顕微鏡動画としての利用可能性は Phase 2.7 以降の対象である。ただし、単一べん毛が定性的にも崩壊せず回転して見えるかは Phase 2.6 のユーザー目視レビュー対象とする。

## 仮説

単に `dt_star` を小さくすると、bond / bend / torsion の数値的な発散は抑えられる。しかし、現行条件では回転 activity (`flag_phase_rate_hz`) も落ちるため、安定回転条件としては不十分である。

一方、`dt_star=2.5e-4` に下げたうえで、motor-on 局所補強を `local_hook_scale=8`, `local_spring_scale=5`, `local_bend_scale=8`, `local_torsion_scale=4` とすると、`4.0e-21 N m` 条件でも multi-step で bond / bend / torsion と回転 activity を両立できる。

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
| `flag_phase_rate_hz` | 回転 activity |
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
| `median(abs(flag_phase_rate_hz))` | `>= 1.0` |

## 代表条件

条件:

- `n_flagella=1`
- `stub_mode=full_flagella`
- `torque=4.0e-21 N m`
- 決定論的条件、Brownian off

| 条件 | 期待結果 | 解釈 |
| --- | --- | --- |
| `dt_star=1.0e-3`, default local scales | fail (`flag`) | Phase 2.5 break representative を再現 |
| `dt_star=2.5e-4`, `local_hook_scale=1`, `local_spring_scale=1`, `local_bend_scale=8`, `local_torsion_scale=1` | fail (`motor_no_rotation`) | 形状は保つが回転 activity を失う |
| `dt_star=2.5e-4`, `local_hook_scale=8`, `local_spring_scale=5`, `local_bend_scale=8`, `local_torsion_scale=4` | pass | bond / bend / torsion と回転 activity を両立 |

これらの local scale は、現時点では物理最適値ではなく、Phase 2.6 の診断 representative である。参照論文モデルとの差分としては、motor-on 時の hook / spring / bend / torsion の局所補強であり、数値安定化寄りの運用パラメータとして扱う。標準 `conf/sim_swim.yaml` は local scale を `1.0` にしているため、この代表条件を再現する場合は CLI override で明示する。

## 標準 sweep コマンド

```bash
uv run python -m scripts.01_simulate_swimming.run_motor_scale_sweep \
  --target local_bend_scale \
  --values 4,8 \
  --torques 4e-21 \
  --duration 0.05 \
  --dt-star 2.5e-4 \
  --stub-mode full_flagella \
  --output-dir outputs/phase2_6_helix_retention
```

CI では runtime を抑えるため、`duration_s=0.05`, `dt_star=2.5e-4` の 200-step 代表条件を固定する。ローカル確認では、同じ `dt_star` / local scale representative / `torque=4.0e-21 N m` 条件で `duration_s=0.25` の 1000-step probe も pass することを確認した。

したがって、単一べん毛については、上記の定量 gate の範囲で「1000-step までの安定的な回転」を確認済みである。ただし、これは単一べん毛・決定論的条件・現行診断 gate 上の確認であり、ユーザー目視による定性検証は未完了である。長時間の実用動画生成や多本べん毛条件での安定性を保証するものではない。

## 定性レビュー

Phase 2.6 は、ユーザー目視レビューが完了するまで完了扱いにしない。

再実行コマンド:

```bash
uv run python -m scripts.01_simulate_swimming \
  --config conf/sim_swim.yaml \
  --duration-s 0.25 \
  --fps-out 25 \
  --render-flagella \
  --render-flagella-2d \
  flagella.n_flagella=1 \
  flagella.stub_mode=full_flagella \
  motor.torque_Nm=4.0e-21 \
  motor.local_hook_scale=8 \
  motor.local_spring_scale=5 \
  motor.local_bend_scale=8 \
  motor.local_torsion_scale=4 \
  time.dt_star=2.5e-4 \
  output.base_dir=outputs/phase2_6_qualitative_review
```

今回生成した確認対象:

- `outputs/phase2_6_qualitative_review/2026-05-31/214109/render/swim3d.mp4`
- `outputs/phase2_6_qualitative_review/2026-05-31/214109/render/swim3d_final.png`
- `outputs/phase2_6_qualitative_review/2026-05-31/214109/render2d/projection.mp4`
- `outputs/phase2_6_qualitative_review/2026-05-31/214109/sim/step_summary.csv`

目視観点:

- 単一べん毛の螺旋形状が 0.25 s / 1000 steps の間に潰れない。
- hook 近傍で不自然な折れ曲がりや飛びが出ない。
- べん毛が停止せず、連続的に回転して見える。
- 菌体やべん毛が画面外へ飛ぶ、または急激に拡大・縮小して見えない。

## 検証

自動検証:

- `tests/test_run_state_fixed.py::test_phase26_default_break_fails_helix_retention_gate`
- `tests/test_run_state_fixed.py::test_phase26_small_dt_without_bend_scale_loses_rotation`
- `tests/test_run_state_fixed.py::test_phase26_small_dt_and_bend_scale_retain_helix`
- `tests/test_motor_scale_sweep.py`

Phase 2.6 は single flagellum の multi-step 定量 gate 固定であり、多本べん毛の束化や 2D 動画自然さは対象外である。ただし、単一べん毛の回転が定性的に安定して見えるかはユーザー目視レビューで確認する。
