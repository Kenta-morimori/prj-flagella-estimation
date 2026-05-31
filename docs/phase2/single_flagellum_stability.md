# Phase 2.5 single flagellum 短時間安定回転基準

## 目的

Phase 2.5 では、`body + hook + single flagellum` 条件で、motor-on の短時間安定回転を再現可能にする。

Phase 2.4 の `minimal_basal_stub` は basal hook 近傍だけを見ていた。Phase 2.5 では `stub_mode=full_flagella` に進み、hook だけでなく flagellum chain の bond / bend / torsion が motor torque 下で破綻しないかを同時に見る。

## 仮説

現行実装では、低トルクでは full flagellum が短時間の回転 activity を持ちながら形状 gate を通過する。一方、高トルクでは hook より先に flagellum 側の bond / bend / torsion が破綻する。

この仮説を以下の代表条件で固定する。

| torque | 期待結果 | first-fail category | 解釈 |
| ---: | --- | --- | --- |
| `1.2e-21 N m` | pass | `none` | 短時間の safe representative |
| `4.0e-21 N m` | fail | `flag` | flagellum 側の break representative |

## Gate 指標

`step_summary.csv` と `body_constraint_diagnostics.csv` を入力に、`src/sim_swim/sim/single_flagellum_gate.py` で短時間 gate を判定する。

主な判定軸は以下。

| 指標 | 役割 |
| --- | --- |
| `finite_pass` | 数値発散・NaN/Inf の検出 |
| `shape_pass_nonbody` | hook / flagellum 側の形状 gate |
| `first_fail_category_nonbody` | non-body 側の first-fail category |
| `body_shape_pass` | body 側 gate |
| `motor_degenerate_axis_count` | motor split の退化軸検出 |
| `motor_bond_length_clipped_count` | motor force split の basal link clipping 検出 |
| `flag_phase_rate_hz` | motor-on の回転 activity |
| `hook_len_rel_err_max` | hook link 長の相対誤差 |
| `local_attach_first_rel_err` | attach-first hook link の相対誤差 |
| `flag_bond_rel_err_max` | flagellum bond 長の最大相対誤差 |
| `flag_bend_err_max_deg` | flagellum bend 角の最大誤差 |
| `flag_torsion_err_max_deg` | flagellum torsion 角の最大誤差 |

短時間 pass には、少なくとも以下を要求する。

- `finite_pass == true`
- `shape_pass_nonbody == true`
- `body_shape_pass == true`
- `motor_degenerate_axis_count == 0`
- `motor_bond_length_clipped_count == 0`
- `median(abs(flag_phase_rate_hz)) >= 1.0`

## 代表条件

条件:

- `n_flagella=1`
- `stub_mode=full_flagella`
- `duration_s=0.02`
- `dt_s=1.0e-3`
- 決定論的条件、Brownian off

代表結果:

| torque | body gate | non-body gate | 支配的な破綻 |
| ---: | --- | --- | --- |
| `1.2e-21 N m` | pass | pass | なし |
| `4.0e-21 N m` | pass | fail (`flag`) | `flag_bond_rel_err_max`, `flag_bend_err_max_deg`, `flag_torsion_err_max_deg` |

`4.0e-21 N m` 条件では、body gate は pass のままで、hook link も先に gate を超えない。一方で flagellum 側の bond / bend / torsion が閾値を超える。そのため Phase 2.5 の break representative は、Phase 2.4 の hook-dominated failure とは別の、flagellum-chain dominated failure と解釈する。

## Phase 2.5 で実施しないこと

本タスクでは、bond / bend / torsion を維持するための物理パラメータ変更や数値安定化補助は実装しない。

理由は、Phase 2.5 の目的が「full flagellum 条件で短時間 motor-on の safe/fail representatives と first-fail 指標を固定すること」だからである。維持策を同じタスクに含めると、破綻モードの観測と対策の効果が混ざり、どの変更がどの指標を改善したか追跡しにくくなる。

Phase 2.6 では、以下を候補として分けて検証する。

- `dt_s` / `dt_star` を小さくして、bond / bend / torsion の時間発展誤差を抑えられるか。
- `local_spring_scale`, `local_bend_scale`, `local_torsion_scale` を sweep し、flagellum chain のどの拘束が支配的かを切り分ける。
- motor torque distribution が flagellum chain に過大な局所変形を与えていないかを、motor split diagnostics と併せて確認する。
- 必要であれば、物理モデル変更ではなく数値安定化として明示したうえで、最小限の補助拘束を検討する。

各候補は `flag_bond_rel_err_max`, `flag_bend_err_max_deg`, `flag_torsion_err_max_deg` の改善だけでなく、`flag_phase_rate_hz` の回転 activity を保つかも同時に見る。

## 標準 sweep コマンド

```bash
uv run python -m scripts.01_simulate_swimming.run_motor_scale_sweep \
  --target local_hook_scale \
  --values 8 \
  --torques 1.2e-21,4e-21 \
  --duration 0.02 \
  --stub-mode full_flagella \
  --output-dir outputs/phase2_5_single_flagellum
```

このコマンドは、`local_hook_scale_sweep_summary.csv` に `stub_mode`, `n_flagella`, `first_fail_category_nonbody`, `body_fail_category` を含む summary を出力する。

## 検証

自動検証:

- `tests/test_simulation.py::test_phase25_single_flagellum_torque_baseline_safe_and_first_fail`
- `tests/test_motor_scale_sweep.py::test_stub_mode_full_flagella_summary_columns`
- `tests/test_motor_forces.py`

Phase 2.5 は短時間の定量 gate 固定であり、動画の自然さや束化の評価は対象外である。そのため、この段階ではユーザー目視レビューは要求しない。
