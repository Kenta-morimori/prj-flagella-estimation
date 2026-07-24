# Phase 4 Current

このファイルは，Phase 4 作業の入口として読む短い現在地ドキュメントである。

## Goal

Phase 4 の目的は，Phase 3 common clip dataset を読み込み，べん毛数 `n_flagella` の学習・評価へ進むことである。

## Current Status

Phase 3 / 4 MVP の標準 clip duration は `0.5 s` に固定する。必要独立run数は，Phase 4 の loader / baseline classifier / grouped learning curve からの feedback を使って #129 で判断する。

現在の主対象:

- #146: Phase 3 common clip dataset loader smoke test。Phase 4 が `.npy` clip と `clip_metadata.jsonl` を再検出なしに読めることを確認する。
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

Loader は training を開始しない。最初に確認する contract は次の通り:

- `.npy` clip は `uint8`，shape `(T, H, W)`。
- `clip.frame_count` と `.npy` の `T` が一致する。
- `normalization.crop_size_px` と `.npy` の `(H, W)` が一致する。
- `labels.n_flagella` と `split_summary.csv` の class label が一致する。
- `track.group_key` が train / val / test をまたがない。

## Next Actions

1. #146 loader smoke test を merge する。
2. Phase 4 baseline classifier の小Issueを #133 配下に作成し，training loop / metrics / artifacts を最小実装する。
3. grouped learning curve を #129 に接続し，`n_flagella` ごとの必要独立run数を判断する。
4. #128 の freeze checklist を Phase 4 dataset build 前の gate として接続する。

## Key References

- Phase 4 roadmap issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/133`
- Loader smoke test issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/146`
- Phase 3 metadata schema: `schemas/phase3_clip_metadata.schema.json`
- Phase 3 GT passthrough pipeline: `src/flagella_estimation/phase3/`

