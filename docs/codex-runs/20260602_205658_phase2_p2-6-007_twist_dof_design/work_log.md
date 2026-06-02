# Work log

## 対象

- Phase 2.6 P2-6-007
- triplet motor でのねじれ・回転自由度不足の明確化

## 実施内容

- `docs/phase2/phase2_6_triplet_twist_dof_design.md` を追加し、現行 bead-position-only model に不足している状態量を整理した。
- material frame / segment twist / axial torque flux の意味と、現行 torsion force との差分を文書化した。
- root torque を flagellum chain へ伝える候補として、material frame + segment twist、segment torsional torque flux、quasi-rigid helical body approximation、distributed_flagellum を比較した。
- `distributed_flagellum` を diagnostic upper-bound として残し、最終解にはしない方針を明記した。
- helix retention gate で root-to-helix 伝達比を評価するため、`flag_phase_deg` と `flag_phase_rate_hz` を必須列に追加した。
- root phase 欠損時に `nonfinite` fail する pytest を追加した。
- `docs/phase2/phase2_tasks.md` で P2-6-007 を完了扱いにし、次タスク P2-6-008 を proposed として追加した。

## 判断

- P2-6-007では物理モデル変更を実装しないため、新規ADRは作成しない。
- 次に物理モデルを変更する場合は、material frame / segment twist 導入のADRを先に作成する。
- `triplet + hook spring fix` だけでは、root 方位の net 回転が螺旋全体の net 回転へ十分に伝わっていない。

## 検証

- `uv run pytest tests/test_helix_retention_gate.py -q`
- `uv run ruff format --check src/sim_swim/sim/helix_retention_gate.py tests/test_helix_retention_gate.py`
- `uv run ruff check src/sim_swim/sim/helix_retention_gate.py tests/test_helix_retention_gate.py`
- `uv run pytest -q`
