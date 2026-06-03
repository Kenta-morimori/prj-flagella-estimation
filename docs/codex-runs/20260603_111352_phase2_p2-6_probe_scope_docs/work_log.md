# P2-6 probe scope docs work log

## 目的

P2-6-007のmerge前整理として、本PRをprobe成功までのscopeに限定し、完全物理実装を次タスクへ分離する。

## 変更内容

- `axial_torque_flux_probe` と `local_twist_transmission_probe` を ADR 0002 (`docs/adr/0002_phase2_torque_transmission_probes.md`) に統合した。
- 既存 ADR 0003 は履歴ADRとして `superseded` にした。
- 完全物理実装の方針を ADR 0004 として新規作成した。
- `docs/phase2/phase2_6_triplet_twist_dof_design.md` を、probe結果と完全物理実装方針が分かれるように修正した。
- `docs/phase2/phase2_tasks.md` に P2-6-008 を proposed として追加した。
- `docs/PROJECT_PLAN.md` で、P2-6-007はprobe成功、完全物理実装は次タスクであることを明記した。
- 既存の P2-6-007 レビュー記録を、Phase 2.6完全完了ではなくprobe成功として読めるように修正した。

## 判断

本PRでは `local_twist_transmission_probe` により、root側orientation activityが先端側へ伝われば単一べん毛でnet 1回転以上とshape gate PASSを両立できることを示した。ただし bead force 変換はまだ近似であり、完全な material frame / segment twist 物理実装ではない。

そのため、merge条件はprobe成功までとし、material frame / segment twist / 局所force coupleによる安定回転はP2-6-008の条件にする。

## 確認

- `uv run python -m json.tool docs/codex-runs/20260602_221638_phase2_p2-6-007_local_twist_probe/review_result.json`
- `uv run python -m json.tool docs/codex-runs/20260603_111352_phase2_p2-6_probe_scope_docs/review_result.json`
- `uv run ruff check .`
- commit hook: `uv run ruff format --check .`: PASS
- commit hook: `uv run ruff check .`: PASS
- commit hook: `uv run pytest -q`: PASS, 141 passed
