# Phase 2 Issue #184: 2015 project nf1–6 実行契約

PR #242では10τ完走を追わず、同一topologyの短時間並列probeから計算費用を測定する。cs10手順は`docs/phase2/phase2_184_2015_nf1_6_cs10_runbook.md`を正本とする。

| 項目 | PR #242 probe |
| --- | --- |
| model | 2015 project、`seeded_surface`、nf1–6、attach/phase seed 0 |
| motor/time | `2.5e-20 N m`、motor=reference=force torque、tracking-reference、`dt_star=1e-5` |
| motion | RUN固定、switchingなし、Brownian OFF |
| duration | 各0.01実秒=0.25τ=25,000 steps。configは丸め誤差を避けるため`0.25 tau`で指定 |
| execution | `cs10_qualified`、単一parallel-job予約、最大3 worker、compact checkpoint 2,500 steps |
| purpose | 0.5実秒=12.5τ=1,250,000 stepsの計算費用外挿。物理的PASS・dataset採択・profile昇格・canonical選定は対象外 |

#61の再解析済み3 torque・1τ decisionをpreflightでauditし、既知の3/3 strict FAILを記録する。decision欠落・不整合時は起動を拒否する。各conditionのgeometry構成はoutput root作成前に確認する。

旧reservation 5の`seeded_center_layer` nf1–3完走結果は長時間performanceとの比較資料として保持する。旧nf6とdirect nf4/nf5は未完走のまま停止し、対象確認後にchild出力だけを削除する。旧・新rootを混ぜた暫定nf1–6物理集計は行わない。将来10τ screenが必要なら、全6条件を同一`seeded_surface` topology、同一commit、clean campaignで実行する。

2026-09-23、probeの6条件実測と不確かさ、nf1–3の既存10τ実測との差を確認したユーザー判断により、2015 project modelは今回のdataset候補に不採用となった。速度probeの完走は物理的PASSではなく、#61 strict FAILとnf5/nf6のonline hook違反も解消していない。今後のclean 10τ campaignは本判断の前提ではなく、再検討する場合だけ新しい契約と承認を要する。2015 paper profileやStage A実装の維持には影響しない。
