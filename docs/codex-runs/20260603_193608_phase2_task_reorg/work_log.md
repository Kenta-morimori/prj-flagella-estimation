# Phase 2 task reorganization

## Summary

2026-06-03 の方針に基づき、P2-7-006 と旧 P2-8-007 を統合した。Phase 2.7 は、複数べん毛で「崩壊しないこと」と「後方束化すること」を同じ条件探索で扱う。

遊泳挙動の運動指標検証は P2-8-008 として整理した。`dt_star=1.0e-4` の長時間実行で3D動画出力枚数が過大になる問題は、P2-9-009 として `output_sampling.fps_out_3d` 導入タスクに分離した。

## Verification

- `git diff --check`
- commit hook:
  - `uv run ruff format --check .`
  - `uv run ruff check .`
  - `uv run pytest -q` (`149 passed`)

## Notes

- P2-7-006 の主条件は `motor.force_distribution=material_twist_local_couple` とする。
- `distributed_flagellum` は複数べん毛でも診断用比較modeとして残す。
- 後方束化は完全な二値判定にせず、partial bundle / independent flagellum を分類できるようにする。
- body姿勢変化、姿勢揺らぎ、bundle axis と body axis のずれは、束化条件が得られた後の P2-8-008 で評価する。
