# Issue #255: 初期hookとmotor反作用比較の実行契約

## 今回の到達点

`2010_hex_project`と`2010 project`の初期形状を確認し、1τ short screenを設計する。simulation、cs10接続、予約、job開始は今回行わない。両モデルともbody–flagella spring-segment排除をOFFにし、flagella–flagella反発は維持する。既存profile既定値、canonical model、dataset、strict QC閾値は変更しない。

## 初期形状

campaignで`flagella.initial_hook_force_neutral=true`と`flagella.initial_hook_body_axis_perpendicular=true`を指定する。付着ビーズとべん毛第1ビーズを固定したまま、各べん毛鎖全体を菌体長軸まわりに剛体回転する。根元接線をhookに垂直にする二候補から、非付着bodyビーズおよび他べん毛との最小距離が最大の組合せを決定的に選び、付着ビーズとの非重複も確認する。初期位相seedは回転前形状を決め、実際の軸回り位相はこの幾何補正で変わる。

preflightではhook–長軸角度誤差、hook角度誤差、hook長、hook初期力、非付着ビーズ間距離を記録する。非重複配置が得られなければ実行前に停止する。hex 13配置とproject n=1..6の固定camera画像・各値は、simulationを行わない`preview_initial_geometry.py`のmanifestで確認する。画像は初期形状だけを示し、時間発展の安定性は示さない。

2026-10-05の実装後preflightは38/38 armで通過した。hook–長軸角度誤差は0°、hook角度誤差の最大は`2.84e-14°`、最小非付着ビーズ間距離はhex `0.440 µm`・project `0.362 µm`、各付着ビーズとの最小距離は`0.250 µm`であり、ビーズ直径`0.200 µm`を上回る。固定cameraの19形状画像とSHA-256をローカル`outputs/2026-10-05/120858/issue255_initial_geometry_preview/manifest.json`に記録した。

画像の再生成は、MacのPython環境で次を実行する。`--output-dir`には新しいJST時刻のdirectoryを指定する。

```bash
.venv/bin/python scripts/01_simulate_swimming/preview_initial_geometry.py \
  --config conf/phase2_multi_run/2010_hex_project_motor_reaction_1tau_issue255.yaml \
  --config conf/phase2_multi_run/2010_project_motor_reaction_1tau_issue255.yaml \
  --output-dir outputs/YYYY-MM-DD/HHMMSS/issue255_initial_geometry_preview
```

## 1τ short screen

| モデル | 形状 | motor反作用 | 独立condition |
| --- | ---: | --- | ---: |
| 2010 hex | #245の13 attachment配置 | 現行軸方向／全vector | 26 |
| 2010 project | `seeded_surface`、n=1..6、attach seed 0 | 現行軸方向／全vector | 12 |

各対は初期位置、seed、torque、`dt_star`、相互作用、出力間隔を一致させる。`T=2.5e-20 N m/flagellum`、reference torque同値、`dt_star=1e-4`、1τ=10,000 internal steps、`phase_seed=0`、Brownian/switching OFF、compact archiveを固定する。既存の全vector反作用は**全body beadへ分配する実装**であり、局所hook反作用の物理的採択を意味しない。

共通`development_evaluation`を使い、finite・body・hook length・flag・motorを別々にPASS/FAIL判定する。hook angleは診断に残す。同一stepのmotor force/torque分解、`run_summary.json`、manifest、condition別archiveを保存する。反作用方式間のトルク残差と形状を対で比較し、合算PASSだけで差を隠さない。

cs10用jobは`conf/phase2_parallel/issue255_motor_reaction/`のhex/project各YAMLを使う。`max_workers: auto`、`worker_policy: cs10_qualified`、実効8 workers、数値ライブラリ各1 thread、condition別output分離と全条件geometry preflightを必須とする。#245のhex 13条件・1τの暫定約42分を単純外挿するとhex 26条件だけで約84分相当であり、project 12条件を含む見積りは実測で更新する。

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py --config conf/phase2_parallel/issue255_motor_reaction/hex_1tau_job.yaml --dry-run
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py --config conf/phase2_parallel/issue255_motor_reaction/project_1tau_job.yaml --dry-run
```

dry-run確認は予約・simulation開始を意味しない。cs10接続、enqueue、dispatcher/tmux、job開始・停止は操作ごとにユーザー明示承認を要する。

## 2sへの移行

1τ完了後にQCと固定camera replayをレビューする。finite・body・hook長・flag形状が安定した**対のみ**を候補にし、2s config/job、condition数、wall timeをその時点で確定する。現行反作用のmotor FAILが残る対を実行する場合は、strict FAILを保持した診断専用とし、特徴量評価やモデル採択に用いない。新たな反作用方式、閾値変更、局所hook物理解釈は別判断とする。
