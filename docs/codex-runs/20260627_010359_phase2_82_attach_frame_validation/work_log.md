# 2026-06-27 Phase 2 Issue #82 Attach-Frame Validation

## Context

Issue #82 / PR #90 の hook過伸長対策について，実装後の Stage A，dt sweep，Stage B，後方条件0.5 s定性評価を整理した。

## Findings

- Stage A: `local_attach_frame_position_scale=3`, `local_attach_frame_tangent_scale=1.5` は baseline の hook過伸長を大きく抑えた。
- dt sweep: `dt_star=1.0e-4, 5.0e-5, 2.5e-5` で破綻種別は解消せず，`dt_star` は本質原因ではない。
- Stage B: `local_first_second_spring_scale=1..3` を追加しても `flag_bond_rel_err_max` は改善しなかった。
- 後方条件0.5 s定性評価: ユーザー確認では hook根元挙動は問題なし。ただし長時間安定性は未確認。
- body-flagella 貫通可能性は Issue #93 へ分離した。

## Output Policy

定量・定性評価後も，`step_summary.csv`，`trajectory.csv`，`state_archive.npz`，`manifest.json` は再現に必要なため削除対象にしない。削除対象は途中停止run，重複run，報告に使わない再生成可能な動画・frame出力に限定する。

## Next

merge前に長時間後方条件 sweep を実行し，summary と根拠figを Issue / PR に残す。
