# Phase 2.8 Issue #116 stability sweep plan

作成日: 2026-07-15

## 目的

Issue #115 の seed 固定 candidate 結果を受け，`n>=4` を dataset v1 (#119) に戻せるか判断するための少数 sweep と判定基準を固定する。

本ドキュメントは長時間 simulation の実行結果ではない。Codex は #116 の実行用 config と評価手順を追加したが，36 条件の simulation はユーザー実行に回す。

## #115 からの入力

source:

- config: `conf/phase2_multi_run/flagella_count_stability_candidates_seed00.yaml`
- output: `outputs/phase2_multi_run/flagella_count_stability_candidates_seed00/summary.csv`

主要結果:

| n | best observed condition | nonbody result | body result | interpretation |
| ---: | --- | --- | --- | --- |
| 4 | `flag_spring=2.0`, `body=1.0` | first fail なし，`max_flag_bond_rel_err=0.9075` | `body_shape_pass=true`, `body_spring_max_stretch_ratio=0.9605` | seed00 では v1 候補になり得る |
| 5 | `flag_spring=2.0`, `body=1.0` 近傍 | nonbody は改善するが body fail が残る | body fail | body 補強とのトレードオフ確認が必要 |
| 6 | `flag_spring=2.0`, `body=2.0` | `t=0.4891 s` に flag fail | `body_shape_pass=true`, `body_spring_max_stretch_ratio=0.8602` | body は改善するが flag fail が残る |

`max_flag_bond_rel_err` の最大 event は local `1-2` が多いが，`0-1`, `2-3`, `3-4`, `4-5` にも散る。既存の `motor.local_first_second_spring_scale` は #94 で主要 failure を改善しなかったため，#116 の第一 sweep では増強しない。

## 少数 sweep

追加 config:

- `conf/phase2_multi_run/flagella_count_stability_narrow_seed00.yaml`

固定条件:

- `n_flagella=[4,5,6]`
- `attach_seed=0`
- `phase_seed=0`
- `time.duration_s=1.0`
- `time.dt_star=1.0e-4`
- `motor.torque_Nm=2.0e-20`
- `motor.force_distribution=root_torque_segment_couples`
- `motor.local_first_second_spring_scale=1.0`

軸:

| axis | key | values | reason |
| --- | --- | --- | --- |
| `flag_spring` | `stiffness_scales.flag_spring` | `1.75, 2.0, 2.25` | #115 の改善中心 `2.0` を細かく確認する |
| `body_stiffness` | `stiffness_scales.body` | `1.0, 1.5, 2.0, 2.5` | `n=4` の body 1.0 と `n=6` の body 2.0 の間，および上側余裕を見る |

36 条件なので広域 heatmap ではなく seed 固定の少数 sweep として扱う。

実行コマンド:

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py config=conf/phase2_multi_run/flagella_count_stability_narrow_seed00.yaml overwrite=true

uv run python scripts/01_simulate_swimming/plot_heatmap.py config=conf/phase2_multi_run/flagella_count_stability_narrow_seed00.yaml
```

出力:

- `outputs/phase2_multi_run/flagella_count_stability_narrow_seed00/summary.csv`
- `outputs/phase2_multi_run/flagella_count_stability_narrow_seed00/plots/`

## 評価指標

優先して読む列:

- `first_fail_category_nonbody`
- `first_fail_t_s`
- `max_flag_bond_rel_err`
- `max_flag_bond_rel_err_local_0_1_per_flag`
- `max_flag_bond_rel_err_local_1_2_per_flag`
- `max_flag_bond_rel_err_local_2_3_per_flag`
- `max_flag_bond_rel_err_local_3_4_per_flag`
- `max_flag_bond_rel_err_local_4_5_per_flag`
- `body_shape_pass`
- `body_fail_category`
- `body_spring_max_stretch_ratio`
- `body_centerline_max_deviation_um`

v1-ready の最低条件:

- 対象範囲の全 `n_flagella` で `first_fail_category_nonbody` が `none`。
- 対象範囲の全 `n_flagella` で `body_shape_pass=true`。
- `max_flag_bond_rel_err` と proximal local bond 指標が #115 best より悪化していない。
- `body_spring_max_stretch_ratio` と `body_centerline_max_deviation_um` が body pass の範囲内にある。

## n=1,2,3 smoke check

候補が出た場合だけ，次の config で v0 baseline を壊していないか確認する。

- `conf/phase2_multi_run/flagella_count_stability_smoke_seed00.yaml`

default は `flag_spring=2.0`, `body=2.0` だが，少数 sweep の winner が異なる場合は runtime override で合わせる。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py config=conf/phase2_multi_run/flagella_count_stability_smoke_seed00.yaml stiffness_scales.flag_spring=2.0 stiffness_scales.body=2.0 overwrite=true

uv run python scripts/01_simulate_swimming/plot_heatmap.py config=conf/phase2_multi_run/flagella_count_stability_smoke_seed00.yaml
```

smoke pass の条件:

- `n=1,2,3` がすべて `first_fail_category_nonbody=none`。
- `body_shape_pass=true`。
- v0 と同じ seed00 条件で，明らかな body deformation または flag bond 悪化がない。

## 判定

現時点では #119 へ渡す v1-ready モデル条件は未決である。#115 だけでは `n=5,6` の failure が残っているため，#119 の dataset v1 再生成はまだ実行しない。

少数 sweep 後の分岐:

- `n=4,5,6` すべてが v1-ready なら，wide heatmap は省略し，`n=1,2,3` smoke check 後に #119 へ `n=1..6` training candidate として渡す。
- `n=4` だけが v1-ready なら，#119 へ渡す候補は `n=1..4` までに制限し，`n=5,6` は diagnostic-only とする。
- `n>=4` の v1-ready 条件が出ない場合，#119 へは進まず，proximal local bond 専用補強を別実装候補として設計する。
- failure が local `1-2` に集中し続ける場合，既存 `local_first_second_spring_scale` の再利用ではなく，flag 内 proximal bond を限定的に扱う project-specific extension を検討する。

## Visual review

この段階では自動指標を優先する。候補が v1-ready 指標を満たした場合だけ，#119 前に replay を作って body deformation，helical shape preservation，swimming-like propulsion をユーザー目視確認する。

```bash
uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py config=conf/phase2_multi_run/flagella_count_stability_narrow_seed00.yaml overwrite=true
```

確認対象:

- `outputs/phase2_multi_run/flagella_count_stability_narrow_seed00/replay/`

