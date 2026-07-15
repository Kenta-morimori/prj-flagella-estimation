# Phase 2.8 Issue #116 stability sweep plan

作成日: 2026-07-15

## 目的

Issue #115 の seed 固定 candidate 結果を受け，`n>=4` を dataset v1 (#119) に戻せるか判断するための少数 sweep と判定基準を固定する。

本ドキュメントは，#116 の実行用 config，評価手順，および user 実行済み seed 固定狭域 sweep の結果を記録する。

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

2026-07-15 の user 実行では，36 condition directory と `step_summary.csv` は生成済みだったが，root `summary.csv` が未生成だった。Codex は再シミュレーションせず，既存 condition output から generic multi-run の集約ロジックで `summary.csv` を再生成し，`plot_heatmap.py` で plots を生成した。

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

## 結果

v1-ready 最低条件を `first_fail_category_nonbody=none` かつ `body_shape_pass=true` とすると，狭域 sweep の結果は次である。

| n | v1-ready conditions | best / nearest condition | result |
| ---: | ---: | --- | --- |
| 4 | 5 | `flag_spring=2.25`, `body=2.5` | first fail なし，`max_flag_bond_rel_err=0.9145`, `body_spring_max_stretch_ratio=0.4700`, `body_centerline_max_deviation_um=0.0672` |
| 5 | 0 | `flag_spring=2.25`, `body=2.5` | `flag` fail at `t=0.9060 s`, `max_flag_bond_rel_err=1.0060`, body pass |
| 6 | 0 | `flag_spring=2.25`, `body=2.0` | `flag` fail at `t=0.4955 s`, `max_flag_bond_rel_err=1.3012`, body pass |

`n=4` は複数条件で自動指標 pass した。`n=5` は `flag_spring=2.25`, `body=2.5` で body pass かつ `max_flag_bond_rel_err=1.0060` まで近づいたが，`t=0.9060 s` に flag fail が残った。`n=6` は `body=2.0` 付近で body pass する一方，flag fail が `t=0.49 s` 付近に残った。

`n=5,6` の nearest condition でも最大 flag bond event は local `1-2` が支配的である。したがって，同じ `flag_spring/body` 軸を単純に広げる wide heatmap は優先しない。`n=5,6` を戻すには，proximal local bond `1-2` 近傍を含む別補強候補の設計が必要である。

## 判定

#119 へ渡す seed00 時点の候補は次とする。

- model condition: `stiffness_scales.flag_spring=2.25`, `stiffness_scales.body=2.5`
- candidate n range: `n=1..4`
- excluded for v1 training candidate at this stage: `n=5,6`

この条件は `n=4` の seed00 自動指標で pass し，`n=5` の nearest condition でも最も pass に近い。ただし `n=1,2,3` smoke check は未実行であるため，#119 の最初の作業として `conf/phase2_multi_run/flagella_count_stability_smoke_seed00.yaml` を `stiffness_scales.flag_spring=2.25`, `stiffness_scales.body=2.5` で実行し，v0 baseline を壊していないことを確認する。

広い heatmap は現時点では不要である。理由は，狭域 sweep で `n=5,6` の failure が同じ proximal flag bond 周辺に残り，`flag_spring/body` の単純な追加探索だけでは v1-ready 条件に届いていないためである。`n=5,6` を Phase3/4 training candidate に戻す場合は，wide heatmap より先に proximal local bond 専用補強を別 task として扱う。

## Visual review

この段階では自動指標を優先する。候補が v1-ready 指標を満たした場合だけ，#119 前に replay を作って body deformation，helical shape preservation，swimming-like propulsion をユーザー目視確認する。

```bash
uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py config=conf/phase2_multi_run/flagella_count_stability_narrow_seed00.yaml overwrite=true
```

確認対象:

- `outputs/phase2_multi_run/flagella_count_stability_narrow_seed00/replay/`
