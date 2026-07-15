# Phase 2.8 improved-model dataset v1 report

作成日: 2026-07-16

## Summary

Issue #119 では，#116 の handoff 条件を受けて改善モデル analysis dataset v1 の canonical config を追加した。

- config: `conf/phase2_multi_run/flagella_count_behavior_v1.yaml`
- dataset_id: `v1`
- raw run root: `outputs/phase2_multi_run/flagella_count_behavior_v1`
- dataset output: `outputs/phase2_analysis/flagella_count_behavior/datasets/v1`
- model_id: `flag_spring2p25_body2p5_candidate`
- dataset_scope: `nf1_2_3_4_as3_ps3_dur1p0`

v1 は `stiffness_scales.flag_spring=2.25` と `stiffness_scales.body=2.5` を使う。これは #116 の seed00 狭域 sweep で `n=4` のみ v1-ready 自動指標を満たし，`n=5,6` は残 failure があったためである。したがって v1 の candidate range は `n=1..4` とし，`n=5,6` は training candidate から外す。

## Smoke Check

#119 冒頭で必要だった `n=1,2,3` smoke check は，既存出力 `outputs/phase2_multi_run/flagella_count_stability_smoke_seed00/summary.csv` を確認した。

- command: `uv run python scripts/01_simulate_swimming/run_multi_run.py config=conf/phase2_multi_run/flagella_count_stability_smoke_seed00.yaml stiffness_scales.flag_spring=2.25 stiffness_scales.body=2.5 overwrite=true`
- run manifest commit: `b308998`
- attach_seed: `0`
- phase_seed: `0`
- result: `n=1,2,3` はすべて `final_shape_pass_nonbody=True`, `body_shape_pass=True`, `first_fail_category_nonbody=none`

| n_flagella | max_flag_bond_rel_err | body_spring_max_stretch_ratio | body_centerline_max_deviation_um |
| --- | ---: | ---: | ---: |
| 1 | 0.173480 | 0.201301 | 0.108015 |
| 2 | 0.408594 | 0.219053 | 0.158040 |
| 3 | 0.628380 | 0.308187 | 0.084332 |

同じ `attach_seed=0`, `phase_seed=0` の v0 historical dataset では `n=1,2,3` がすべて `strict_pass` であるため，seed00 smoke の範囲では v1 候補値は v0 baseline を明確には壊していない。ただし，full seed expansion は未実行である。

## Pending Full Generation

v1 の 36 sample full raw output は未生成である。長時間 run になるため，ユーザー実行に回す。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py config=conf/phase2_multi_run/flagella_count_behavior_v1.yaml overwrite=true
```

raw 生成後，Codex またはユーザーは次を実行して dataset と分布分析を作る。

```bash
uv run python scripts/02_phase2_analysis/build_dataset.py config=conf/phase2_multi_run/flagella_count_behavior_v1.yaml overwrite=true

uv run python scripts/02_phase2_analysis/plot_distributions.py --dataset-id v1 --overwrite
```

visual review が必要な場合は replay を生成する。

```bash
uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py config=conf/phase2_multi_run/flagella_count_behavior_v1.yaml overwrite=true
```

## Interpretation

現時点では v1 の canonical config と dataset ID / output path は確定したが，full v1 dataset は未生成である。そのため `n=4` を Phase3/4 training candidate に含める最終判断は未確定であり，`review_result.json` は `FAIL` のままにする。

full v1 生成後は，v0 と v1 の `n=1,2,3` について `quality_class`, `cell_mean_speed`, `cell_straightness`, `cell_angular_velocity_rms`, `hook_drift`, `first_fail_category` を比較し，`n=4` は `strict_pass` / body diagnostics / visual review を合わせて training candidate 可否を判断する。
