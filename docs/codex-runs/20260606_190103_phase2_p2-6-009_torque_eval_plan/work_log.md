# P2-6-009 torque transmission model evaluation planning

## 実施内容

- 最新 `main` を取得した。
- `feature/phase2-54-torque-transmission-eval` を作成した。
- Issue #54 を、後方束化前に行う単一べん毛モデル詳細評価として整理した。
- `docs/phase2/phase2_6_torque_transmission_model_evaluation.md` を追加した。
- `docs/phase2/phase2_tasks.md` に `P2-6-009` を追加した。

## 方針

最初に `local_*_scale=1.0` で torque sweep を行い、`material_twist_local_couple` の安定境界を確認する。

その後、`local_spring_scale`, `local_bend_scale`, `local_torsion_scale`, `local_hook_scale` の one-factor sweep を行い、支配的な局所項を見つける。必要な軸だけ heatmap 化する。

多べん毛・後方束化・遊泳検証は Issue #58 に分離し、Issue #54 の結果で代表条件を決めてから進める。

## 検証

- `git diff --check`
- commit hook:
  - `uv run ruff format --check .`
  - `uv run ruff check .`
  - `uv run pytest -q` -> 149 passed
