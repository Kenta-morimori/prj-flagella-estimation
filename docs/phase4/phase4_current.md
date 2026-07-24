# Phase 4 Current

このファイルは，Phase 4 作業の入口として読む短い現在地ドキュメントである。

## Goal

Phase 4 の目的は，Phase 3 common clip dataset を読み込み，べん毛数 `n_flagella` の学習・評価へ進むことである。

## Current Status

Phase 3 / 4 MVP の標準 clip duration は `0.5 s` に固定する。必要独立run数は，Phase 4 の loader / baseline classifier / grouped learning curve からの feedback を使って #129 で判断する。

現在の主対象:

- #146: Phase 3 common clip dataset loader smoke test。PR #147 で完了した。
- #148: common clip baseline classifier。PR #149 で完了した。
- #150: grouped learning curve evaluator。unique `group_key` を学習量として評価する。
- #128: 学習datasetへ混ぜてよい条件変更を整理する。
- #129: 1 clip の時間長と必要な独立run数を決める。`0.5 s` default は決定済みで，必要run数評価は Phase 4 feedback 待ち。
- #145: RUN-TUMBLE を含む dataset v2。Project Status は `TODO` だが，Phase 3 / 4 MVP 後に後回し。

## Input Contract

Phase 4 loader の最小入力:

- `manifest.json`
- `clip_metadata.jsonl`
- `split_summary.csv`
- `qc_summary.csv`
- `clips/<clip_id>.npy`

Loader は次の contract を検査する:

- `.npy` clip は `uint8`，shape `(T, H, W)`。
- `clip.frame_count` と `.npy` の `T` が一致する。
- `normalization.crop_size_px` と `.npy` の `(H, W)` が一致する。
- `labels.n_flagella` と `split_summary.csv` の class label が一致する。
- `track.group_key` が train / val / test をまたがない。

## Baseline Contract

#148 の baseline は最終モデル候補の採択ではなく，Phase 3 -> 4 の学習・評価導線と #129 の learning curve 入力を検証するための診断 baseline である。

- feature: clip intensity，foreground，時間差分，重心移動，radial spread の決定的な統計量
- classifier: feature standardization付き nearest centroid
- split: Phase 3 の `split_summary.csv` をそのまま使い，再分割しない
- freeze: dataset v1，`n_flagella=1,2,3`，`0.5 s`，`non_overlap`，`qc.status=pass`
- metrics: accuracy，balanced accuracy，macro F1，confusion matrix
- artifacts: `model.npz`, `metrics.json`, `predictions.csv`, `confusion_matrix.csv`, `manifest.json`, `run.log`

実pilot 18 clips / 9 independent groups の smoke result:

- output: `outputs/2026-07-24/142055/phase4_baseline_v1`
- train macro F1: `1.0`
- validation macro F1: `0.8222`
- test macro F1: `1.0`

各class 3 independent groupsのみなので，この性能値はモデル採択や必要run数の根拠としては不十分である。#129 では grouped learning curve と group-level uncertainty を評価する。

## Grouped Learning Curve

#150 evaluatorは同一 `group_key` のclip特徴を平均し，1 independent groupを1 feature vectorとして扱う。Phase 3のtest splitは最終評価用に保護し，train+val pool内でclassごとに固定数をholdoutして，残りから`k` groupsを学習する。

全27 v1 candidateの軽量pilot:

- Phase 3 output: `outputs/2026-07-24/143640/phase3_gt_passthrough_v1_full_candidates`
- Phase 4 output: `outputs/2026-07-24/144032/phase4_grouped_learning_curve_full_candidates`
- independent groups: 27 (`n_flagella=1,2,3` 各9)
- development pool: 各class 8 groups
- grouped holdout: 各class 2 groups
- protected test: 各class 1 group，learning curveでは未使用
- repeats: 各`k`で100

macro F1 mean:

| train groups / class | macro F1 mean | empirical 2.5% | empirical 97.5% |
| ---: | ---: | ---: | ---: |
| 1 | 0.846 | 0.467 | 1.000 |
| 2 | 0.964 | 0.822 | 1.000 |
| 3 | 0.987 | 0.822 | 1.000 |
| 4 | 1.000 | 1.000 | 1.000 |
| 5 | 1.000 | 1.000 | 1.000 |
| 6 | 1.000 | 1.000 | 1.000 |

現行pseudo v1と診断baseline内では`k=4`でplateauした。ただしprotected testは各class 1 groupのみであり，実動画への一般化も未評価なので，「一般に必要な独立simulation run数=4」とはまだ決めない。

## Next Actions

1. #150 grouped learning curve evaluatorをreview / mergeする。
2. #128 の freeze checklistをmachine-readable audit gateとして接続する。
3. #129で`k=4`をMVP dataset-size下限として採用するか，protected評価用runを追加してから決めるかを判断する。
4. #145 RUN-TUMBLE dataset v2はProject `TODO`のまま後回しにする。

## Key References

- Phase 4 roadmap issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/133`
- Loader smoke test issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/146`
- Baseline classifier issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/148`
- Grouped learning curve issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/150`
- Phase 3 metadata schema: `schemas/phase3_clip_metadata.schema.json`
- Phase 3 GT passthrough pipeline: `src/flagella_estimation/phase3/`
