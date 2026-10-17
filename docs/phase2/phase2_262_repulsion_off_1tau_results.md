# Issue #262：排除体積全OFFの1τ short screen

2026-10-17にcs10 queue #27で`2010_hex_project`の13配置を実行した。固定commitは`a444a1bd45422d9b5b911df44ee61b397a30f809`、jobは`conf/phase2_parallel/issue262_repulsion_off/hex_1tau_job.yaml`。#255局所one-ring全ベクトル反作用モデルを維持し、body–flagellaとflagella–flagellaのsegment反発をともにOFFにした。#263の候補やC条件は含めない。

## 完走、来歴、QC

jobは19:57:16～20:22:44 JST、**25分28秒**で`succeeded`。13/13配置が各10,000 step、1τ = 0.04 sを完走した。各`run_summary.json`はstep連続性とonline finite/shapeゲートがPASS。共通`development_evaluation`のstrict finite/body/hook length/flag/motor判定も**13/13 PASS**で、first failureはない。全条件のmotor force残差比最大は`1.16e-15`、torque残差比最大は`2.80e-15`。既定上限`1e-8`、`0.02`を変更していない。hook長相対誤差最大`0.0394`、flag bond相対誤差最大`0.0367`、body spring stretch最大`0.0618`。hook角ずれ最大`3.81°`は診断値で、合否閾値に含めない。

各子campaignの`run_manifest.json`を照合した。全13配置で`body_flagella_enabled=false`、`flagella_flagella_enabled=false`、`body_flagella_active=flagella_flagella_active=active_total=0`。candidate数は非ゼロであり、OFFはpair選択の結果として記録された。13 manifestのSHA-256は全件一致。aggregate portable artifact 122件もリモートとローカルのファイル集合・SHA-256が全件一致した。operational logとcredentialは同期していない。

## geometryとreplay

13配置のgeometry preflightは実行前にPASS。共通評価器が初期配置図をn=1～6別に生成し、固定cameraの3D/2D replayを各6本、計12本生成した。全動画41フレームをdecode確認した。3Dのn=3とn=6の終端frameを目視確認し、この短時間内に明白な菌体崩壊は見られない。2Dは菌体像が主に見える描画であり、べん毛の空間形状は3Dで確認する。**1τのPASSは2秒の安定性や菌体軸揺れの判定ではない。**

- 元campaign：`outputs/2026-10-17/202338/issue262_hex_repulsion_off_1tau/`
- 共通QC・初期配置図・replay：`outputs/2026-10-17/202338/issue262_1tau_evaluation/`
- SHA-256照合：`outputs/2026-10-17/202338/sync_verification.json`と`condition_manifest_verification.json`
- 子manifest：`outputs/2026-10-17/202338/issue262_condition_manifests/`
- 13配置のQC表と固定camera動画をまとめたレビュー入口：`outputs/2026-10-17/202338/issue262_1tau_review.html`

## 次の判断

ユーザーがstrict QC、初期geometry、固定camera replayをレビューするまで2秒jobは予約・開始しない。13配置すべてが原則上の2秒候補だが、短時間結果から「定常／揺れあり」の診断はしない。1τ実測を単純に50倍すると2秒の目安は約21時間13分であり、長時間の性能や完走を保証しない。#255局所ON対照の2秒は約41時間だった。2秒のOFF結果が揃ってから同時刻の軸方向時系列、0.1秒窓、wobble、replayを対照と比較する。

queueはenqueue時に既存dispatcherが稼働しており、#27が自動開始した。開始の個別承認より先に動いたことをユーザーへ直ちに報告し、その後に実行継続の明示承認を得た。今後は**enqueue前にdispatcher稼働状態を確認**し、稼働中ならenqueueと自動job開始の両操作について先に明示承認を得る。
