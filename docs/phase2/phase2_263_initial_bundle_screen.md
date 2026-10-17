# Issue #263：初期束配置と短期screen

## 比較契約

`conf/phase2_multi_run/2010_hex_nf3_bundle_initial_screen_issue263.yaml`の3条件を比較する。`nf03__slots024__C`は#255四方式比較のCと同じ初期座標・物理設定。`bundle`はopt-inの後方コンパクト配置、`bundle_phase5`は同配置のflag 1だけを後方軸まわりに5°回す。3条件はnf3、slots024、T=2.5e-20 N m/本、dt_star=1e-4、Brownian/switching OFF、菌体–べん毛排除OFF、べん毛同士反発ONを共有する。

初期配置の手順・制限はADR 0027。preflightでは束半径・後方軸角・べん毛間軸角、hook長・角度・力、最小ビーズ間距離を記録する。非重複・hook長が成立しなければsimulationを始めない。t=0の幾何成立を遊泳の成立と扱わない。

## 実行前の確認と実行境界

Issue targetは`cs10_user_run`、独立3条件、parallel jobは`conf/phase2_parallel/issue263_initial_bundle/screen_job.yaml`。`cs10_qualified`、`max_workers: auto`、条件別出力を使用する。Macの長時間simulationは30分超と見積もり、実行しない。短期screenの実測後に長時間条件の所要時間を更新する。

Macではtargeted test、初期配置preview、`run_parallel.py --config conf/phase2_parallel/issue263_initial_bundle/screen_job.yaml --dry-run`まで確認する。cs10接続、queue予約、dispatcher/job開始は操作ごとに固定commit・範囲を示したユーザー明示承認が必要。過去の#255予約・許可は使用しない。

ローカルdry-run planは`outputs/2026-10-17/193744/dry_run_plan.json`に保存した。実効3 workers、数値ライブラリthread各1、3条件の出力分離を確認済み。これはMac上の計画確認であり、cs10の実効資源・実行許可を示すものではない。

screen成果物は共通`development_evaluation`のfinite/body/hook length/flag/motor QCを先に確認し、固定cameraの3D／2D replayと`bundle_rotation`の時系列・周期診断を確認する。Cとの差分が初期幾何だけであることを保持する。QC未達なら幾何起因を切り分け、QC通過でも遊泳未達なら#255のA/B/C/D保存済み証拠と照合してCの物理前提を一項目ずつ再検討する。長時間stageの条件・必要回転周期数・挙動判定はscreenレビュー後に確定する。

## ローカル初期幾何確認

初期の探索では表面内向きhookを持つ配置が現れたため、候補をADR 0027の外向き・接線方向に限定した。制限後の画像とpreflightは`outputs/2026-10-17/193744/issue263_initial_geometry_preview/manifest.json`に保存した。初期束半径はC 1.041 µm、bundle 0.813 µm、phase5 0.807 µm。最小非付着ビーズ間距離は順に0.440、0.250、0.250 µmで、直径0.200 µmを超えた。hook方向は菌体外向き径方向から90°、hook長誤差は最大5.3e-23 m、hook初期力normは最大4.0e-28 N。これは静的確認であり、短期screenは未実行。
