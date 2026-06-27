# 2026-06-27 Phase 2 Issue #82 Attach-Frame Validation

## Context

Issue #82 / PR #90 の hook過伸長対策について，実装後の Stage A，dt sweep，Stage B，後方条件0.5 s定性評価を整理した。

## Findings

- Stage A: `local_attach_frame_position_scale=3`, `local_attach_frame_tangent_scale=1.5` は baseline の hook過伸長を大きく抑えた。
- dt sweep: `dt_star=1.0e-4, 5.0e-5, 2.5e-5` で破綻種別は解消せず，`dt_star` は本質原因ではない。
- Stage B: `local_first_second_spring_scale=1..3` を追加しても `flag_bond_rel_err_max` は改善しなかった。
- 後方条件0.5 s定性評価: ユーザー確認では hook根元挙動は問題なし。ただし長時間安定性は未確認。
- 長時間3D定性評価: `outputs/phase2_82/qualitative_long_flag_bond_review/frame_fp3_ft1p5_fs1p5/2026-06-27/134540` では first fail `t=0.4363 s`, category `flag`, `hook_len_rel_err_max=0.0157`, `flag_bond_rel_err_max=1.0006`。final `t=1.9999 s` では `hook_len_rel_err_max=0.0229`, `flag_bond_rel_err_max=2.0500`, `flag_bond_len_max_over_b=1.7690` となり，破綻は hook ではなく flagellum bond 過伸長である。
- `root_torque_segment_couples` は反作用トルクを消していない。body 側へ反対向き torque を入れているが，attach-frame補強が body-root 間の相対運動を抑えるため，一体回転に近く見える。
- 全 local scale default は `1.0` に据え置く。`fp=3, ft=1.5` は hook過伸長を抑える診断候補だが，長時間 flag bond 過伸長が残るため標準defaultへ昇格しない。
- body-flagella 貫通可能性は Issue #93 へ分離した。

## Output Policy

定量・定性評価後も，`step_summary.csv`，`trajectory.csv`，`state_archive.npz`，`manifest.json` は再現に必要なため削除対象にしない。削除対象は途中停止run，重複run，報告に使わない再生成可能な動画・frame出力に限定する。

## Next

merge前に長時間後方条件 sweep の summary と根拠fig，および `frame_fp3_ft1p5_fs1p5` の3D定性評価結果を Issue / PR に残す。

Issue #82 の sub-issue 候補として `[Phase2] attach-frame補強後のflag bond過伸長を診断・安定化する` を切り出す。最初の目的は，`flag_bond_rel_err_max` が発生する `flag_id` と bead pair を per-step で特定し，root近傍か下流helixかを判定すること。
