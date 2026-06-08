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

2026-06-06の追加方針:

- torque sweep は初期上限を `4.0e-20 N m` とする。
- 全条件を `duration_s=0.5` で評価する。
- `distributed_flagellum` と probe 系 mode は必須比較に含めず、`material_twist_local_couple` を単独評価する。
- torsion force OFF 診断は必須にしない。
- 代表PASS条件では、遊泳に必要な駆動力の前段評価として菌体重心変位・平均速度・body axis角度変化も確認する。

## 検証

- `git diff --check`
- commit hook:
  - `uv run ruff format --check .`
  - `uv run ruff check .`
  - `uv run pytest -q` -> 149 passed
