# Phase 3.4 common clip pipeline implementation plan

この文書は Issue #6 の最小実装準備メモである。#127 schema に沿った共通 clip pipeline を，擬似動画 GT passthrough から始めるための作業分解と test 境界を固定する。

## Goal

Phase 2 擬似動画と将来の実動画 detection / tracking 経路を，同じ clip / metadata schema へ流し込む。最初の実装は pseudo GT passthrough のみを対象にし，実動画 detector / tracker の採用判断は #8 / #9 の入力条件整理後に分離する。

## Minimal Pipeline

```text
Phase 2 dataset v1
  -> source video inventory
  -> GT track adapter
  -> window generator (default clip.duration_s=0.5)
  -> crop / centering / normalization
  -> clip writer
  -> #127 metadata writer
  -> grouped split summary / QC summary
```

## Work Breakdown

| step | scope | output | test boundary |
| ---: | --- | --- | --- |
| 1 | input inventory | source video path, fps, frame count, provenance | small JSON fixture |
| 2 | GT track adapter | `track_id`, bbox/center frame sequence, `group_key` | synthetic track fixture |
| 3 | window generator | `0.5 s` non-overlap default; `0.25 s` / `1.0 s` and overlap comparison windows | pure function pytest |
| 4 | crop planner | `crop_xywh_px`, boundary padding policy | image-free geometry pytest |
| 5 | clip writer | `.npy` or frame sequence path and shape | tmpdir smoke test |
| 6 | metadata writer | #127 JSON schema valid records | schema contract pytest |
| 7 | grouped split summary | no `group_key` crosses split | leakage pytest |
| 8 | CLI wrapper | manifest, run.log, output paths | lightweight subprocess smoke only after library tests |

## Inputs

Pseudo GT passthrough minimum input:

- Phase 2 analysis dataset root: `outputs/phase2_analysis/flagella_count_behavior/datasets/v1`
- raw run directory or replay directory for each `run_id`
- render video path or frame sequence path
- frame rate and source frame count
- GT label `n_flagella`
- provenance fields: `model_id`, `dataset_version`, `run_id`, `render_id`
- default clip config: `clip.duration_s=0.5`, `clip.window_policy=non_overlap`

The first implementation may use a small checked-in fixture instead of real generated videos. It should not require heavy simulation, sweep, render, or real AVI processing.

## Outputs

Minimum output tree:

```text
outputs/YYYY-MM-DD/HHMMSS/phase3_common_clip/
├─ run.log
├─ manifest.json
├─ clip_metadata.jsonl
├─ split_summary.csv
├─ qc_summary.csv
└─ clips/
   └─ <clip_id>.npy
```

`clip_metadata.jsonl` contains one #127-compatible metadata object per clip. `manifest.json` records config path, overrides, input paths, output paths, git commit, environment, and schema version when available.

## Library Interfaces

Initial module candidates under `src/flagella_estimation/phase3/`:

- `inventory.py`: read source descriptors and provenance.
- `tracks.py`: normalize GT / detection track records into a common track object.
- `windows.py`: generate frame windows from fps, source frame count, duration, and policy.
- `crops.py`: compute crop windows and padding.
- `metadata.py`: build #127 metadata objects.
- `splits.py`: grouped split helper by `track.group_key`.
- `cli.py`: orchestration boundary for scripts.

Implementation should keep real-video detection behind an adapter interface so that pseudo GT passthrough can land first.

## Real Video Boundary

Do not finalize these in the first #6 implementation PR:

- detector model or segmentation method
- tracker association policy
- scale normalization based on real body length
- actual use of `data/20250716data1.avi` / `data/20250716data2.avi`
- visual acceptance of real video crops

Real video work should wait for #8 / #9. The current real AVI files remain git-ignored and must not be added or pushed.

## Initial Tests

Recommended first tests:

- `test_phase3_window_generation_defaults_to_0p5s_non_overlap`
- `test_phase3_window_generation_supports_0p25s_1p0s_and_overlap`
- `test_phase3_grouped_split_rejects_cross_split_group_key`
- `test_phase3_gt_passthrough_metadata_matches_schema`
- `test_phase3_crop_planner_clips_or_pads_boundaries`
- `test_phase3_manifest_records_schema_and_inputs`
- `test_phase3_training_freeze_rejects_torque_variation_for_mvp`

These are library-level and should run without generated videos. CLI subprocess tests can be added after the library boundary is stable.

## User-Run Required Later

Once the lightweight pipeline exists, a pilot run over Phase 2 dataset v1 can be requested. Codex should not run this automatically if it requires heavy video IO or long render.

Exact command draft:

```bash
uv run python scripts/03_phase3/build_clip_dataset.py \
  config=conf/phase3/gt_passthrough_v1.yaml \
  input_dataset=outputs/phase2_analysis/flagella_count_behavior/datasets/v1 \
  output_dir=outputs/2026-07-23/phase3_gt_passthrough_v1 \
  clip.duration_s=0.5 \
  clip.window_policy=non_overlap
```

Expected outputs:

```text
outputs/2026-07-23/phase3_gt_passthrough_v1/run.log
outputs/2026-07-23/phase3_gt_passthrough_v1/manifest.json
outputs/2026-07-23/phase3_gt_passthrough_v1/clip_metadata.jsonl
outputs/2026-07-23/phase3_gt_passthrough_v1/split_summary.csv
outputs/2026-07-23/phase3_gt_passthrough_v1/qc_summary.csv
outputs/2026-07-23/phase3_gt_passthrough_v1/clips/
```

Decision points:

- `clip_metadata.jsonl` が #127 schema に通るか。
- all split summary で `group_key` leakage がないか。
- crop が body center を含み，source boundary で破綻しないか。
- `0.5 s` non-overlap の clip 数が #129 評価に十分か。
- pseudo GT passthrough で Phase 4 が再検出なしに読める形になっているか。
