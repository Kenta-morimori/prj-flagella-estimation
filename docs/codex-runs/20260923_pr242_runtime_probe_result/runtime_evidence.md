# PR #242 / Issue #184: 2015 project計算費用の実測

## 入力と検証

- cs10予約8、固定commit `8278c82`、`cs10_qualified`、最大3 worker。`seeded_surface`・seed 0・`2.5e-20 N m`・tracking-reference・`dt_star=1e-5`・RUN固定・Brownian OFF。
- 6条件は各0.01実秒=0.25τ=25,000 internal stepsを完了。`job_manifest.status=succeeded`、`aggregation.status=completed`、`campaign_completion.status=completed`、child exit codeは全て0。6本の相対symlinkとsummary/performance/compact archiveを確認した。
- NASの予約8 rootから83ファイル、旧reservation 5のnf1–3比較rootから34ファイルをローカルへ同期した。運用logは同期せず、両rootの全同期ファイルのSHA-256がcs10と一致した。reservation 7の途中artifactは入力に含めない。
- 元root: `outputs/2026-09-22/203850/parallel/queue-00008-d1ca2234__480d09541b15/`。旧比較root: `outputs/2026-09-08/044446/parallel/queue-00005-14c7fb53__ccff7345439f/`。費用表: 元rootの`analysis/runtime_projection/runtime_projection.csv/json`。

| 主なartifact | SHA-256 |
| --- | --- |
| 予約8 `job_manifest.json` | `2cba8660eb04e3fb3458163c6140e034207a1a3db2ff872dc75ea1299ea857f1` |
| campaign `run_manifest.json` | `644136b1a9e728da8976baa3963fcb42ed79ce191bdd3921dd131cf81dd821be` |
| campaign `summary.csv` | `12e9dff453972ecbc74a68011a2848d8ed270de0107fbd4f6e566e8313071729` |
| `runtime_projection.csv` | `f72cdccfd109fb8fe66a1a5d496e3d5526487da7fb08057688d0331611b1fc3a` |
| `runtime_projection.json` | `b68176e58754bca2af3221245771bda2e43341ba14c180d0c3d1c01b7069596a` |

## 費用表

0.5実秒=12.5τ=1,250,000 stepへの外挿係数は50。日数はcs10のcondition wall timeに50を掛けた値で、Macの実測ではない。

| 本数 | 0.01秒実測（時間） | steps/s | 0.5秒外挿（日） | 旧10τ実測からの0.5秒外挿（日） | online nonbody shape |
| --- | ---: | ---: | ---: | ---: | --- |
| nf1 | 0.84 | 8.278 | 1.75 | 2.47 | 違反なし |
| nf2 | 1.92 | 3.618 | 4.00 | 4.67 | 違反なし |
| nf3 | 3.33 | 2.087 | 6.93 | 7.41 | 違反なし |
| nf4 | 4.95 | 1.402 | 10.32 | — | 違反なし |
| nf5 | 6.90 | 1.007 | 14.37 | — | hook、step 0 |
| nf6 | 9.35 | 0.743 | 19.48 | — | hook、step 0 |

probe jobの実測wall timeは12.68時間。6条件を同じ順序で3 workerへ割り当てた0.5秒の外挿makespanは26.42日（condition計算時間の総和は56.85 worker日）。旧nf1–3の10τ実測との外挿比はそれぞれ0.71/0.86/0.94であり、この短時間probeからの単純外挿は旧実測より約6〜29%短い。旧runはattachment topologyが異なるため、比は較正の参考にとどめる。nf4–6には長時間実測による較正がない。

有限性とbody shapeのオンラインgateは全6条件で違反なし。nf5/nf6のhook系nonbody shape違反は最初の内部step（`t_s=4e-7 s`）に記録された。`runtime_projection`のstrict QC欄は`not_evaluated`であり、速度測定の成功を物理的PASSと扱わない。長時間のthroughput、I/O、並列時の資源競合は変化し得る。2015を採用候補から外すかはユーザー判断待ちである。
