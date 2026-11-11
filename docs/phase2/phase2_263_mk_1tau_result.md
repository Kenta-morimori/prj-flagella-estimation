# Issue #263：M/Kの短期QCと1秒比較への移行

2026-11-11。Refs #263 / PR #264。[承認済みロードマップ](phase2_263_autonomous_roadmap.md)に沿った第1回のshort。**4/4共通strict QC PASS、固定camera replayの代理レビューPASS。1秒比較#44を同じsimulation commitで開始した。束形成や定常遊泳の採択ではない。**

## 固定実行と成果物

- short予約#43、simulation/解析source `e32c4ad033411c4bbdfb3fcc6583eab10b93e08e`。nf3/slots024のDJ・M・K・MK、同じJ座標、各0.04秒・10,000連続step。実終端0.03999999999999624秒。
- short job SHA256 `cabf77fef5a73c45e5731f3d665be2267e1ac7dd167b88f4e9a0014414844314`、campaign SHA256 `587470b09539c69280c941b1659afd7a9df964c5ba4ab8a4601251873322f636`。4ID・固定commit・入力と完走を照合。
- 11/11 23:02 JST開始、最遅条件は23:38 JSTに完走。実効4workers・数値thread各1。既存dispatcher/FIFOを利用し、他予約を変更していない。
- portable **37 files / 1,016,482 bytes**。明示したportable集合だけを同期し、remoteの転送前後とlocalの集合・件数・全SHA-256が一致。operational log・credential・partialは同期していない。
- 各`run_summary.json`を先に読み、共通`development_evaluation`のstrict QC・heatmap・固定camera 3D/2D replay、QC-PASS archiveの共通`bundle_rotation`を実施した。独自のQC gateは追加していない。

ローカル正本：`outputs/2026-11-11/234018/issue263_mk_short43_review/`。`portable_index.json`、`run_summary_first.json`、`evaluation/`、`bundle_rotation/`、`saved_contact_body.json`、`short_quantitative_summary.json`、`replay_verification.json`、`review_result.json`を保持する。[実行手順](phase2_263_mk_runbook.md)。

## 短期の定量結果

| 条件 | strict QC | 各鎖の相対材料点spin [回] | 菌体のlabel付き材料点自転 [回] | 最大曲げ偏差 [°] | 最大torsion偏差 [°] |
| --- | --- | --- | ---: | ---: | ---: |
| DJ | PASS | .1850 / .1836 / .1849 | .05478 | 2.566 | .598 |
| M | PASS | .1813 / .1792 / .1807 | .05109 | 3.131 | 4.303 |
| K | PASS | .1848 / .1833 / .1850 | .05467 | 1.242 | .714 |
| MK | PASS | .1811 / .1787 / .1807 | .05084 | 1.516 | .845 |

Kを加えると、この初期窓では曲げ偏差が小さい。M単独はtorsion偏差がDJより大きく、Mを加えただけで全てが改善したとはいえない。これは全鎖・全内部stepの最大偏差であり、近位だけのRMS効果とは区別する。M/Kの整列・前進への効果や交互作用は1秒比較で確認する。

hook相対伸びの最大はDJ/Kで約.035、M/MKで約.042。motorの合力・合トルク残差は丸め誤差の範囲。Mの解析力は実行中も指定vectorを再現し、旧D専用metricのnullは未適用として扱う。QC閾値は変更していない。

保存41状態で、runtimeと同じactive線分pair/production距離を用いた最短距離はBFが全条件.213832 µm、非body pairは.575255〜.575967 µm。2a=.20 µmを下回る保存状態はない。ただし保存点の確認を未保存の内部stepの非貫通保証とはしない。

各鎖はまだ約.18回転、菌体も約.05自転であり、完全周期・持続回転・定常前進を判定する窓ではない。online PCA-frame rollをlabel付き材料点自転として使っていない。

## replayの代理レビュー

共通視野・固定cameraの3D/2D MP4は各41frame・25fps、物理coverage 0〜.04秒。ffprobeのframe/fpsとSHAを照合し、3Dの開始・中間・終端、2Dの中間を視認した。同じ初期座標、全geometryの共通視野内表示、明白な崩壊・交差・実装/描画不一致がないことを確認した。2Dは既存の菌体投影で、べん毛形状は3Dと保存geometryで確認する。広いJ初期配列を完成した密束とは扱わない。

代理レビューはユーザーの明示委任範囲であり、追加承認待ちは設けない。[短期review_result](../codex-runs/20261111_issue263_mk_short_result/review_result.json)。

## 承認済み1秒比較#44

shortの入口確認後、同じ固定commit `e32c4ad033411c4bbdfb3fcc6583eab10b93e08e` の `nf3_mk_1s_job.yaml` をMac/cs10でsimulationなしdry-runし、1回だけenqueueした。予約#44は11/11 23:46 JSTに4条件running/failed0で開始。入力・4ID・固定commit・実効4workers/thread1を照合した。

- 1秒job SHA256 `f369b6fba4b0c8a6451cfa4cd86628c3d555478a258f6d4149a69e4731a8713a`、campaign SHA256 `0834a3f481eee142891e0ff39a2c80d9934a02c694ca64889c599e0e085a7d2f`。
- 各1秒・公称250,000step、既存ceil計画は250,001step。公称終端1.000004秒、主比較は1.0秒まで。差を隠すためのdt/計画実装変更はしない。
- T/reference T=2.5e-20 Nm/本、dt4µs、粘性1mPa s、FF/BF ON、自然形状と反発曲線・既存QCはshortと同じ。変更は実行時間だけ。
- 短期実測からwall約15時間、初期終了見込みは11/12 14:40 JSTごろ。1秒自身のcompact checkpointから更新する。共有負荷と形状変化で速度は変わり得る。
- 台帳は**2/5予約使用**、short1/2、1秒1/2・4/8条件、公称1,040,000step・runtime1,040,004step。失敗も枠に算入する。

比較窓`.04–.20 / .20–.65 / .65–1.00秒`で、同じ遠位bead3以降の軸・幅、根元内部変形、材料点spin、符号付き前進を比較する。入口を満たさない候補を時間延長で救わない。3.5秒は合意済みの候補条件を満たす場合だけで、現時点では未予約。追加位置復元code/scale削除はmerge前残件、ADR0029はProposed、科学的採択・mergeは行っていない。

## 待機中の静的補助診断

`outputs/2026-11-11/232848/issue263_conical_axis_reference/`には、現在のJ根元を固定した対称な直線軸族の計算を保持する。根元半径.7462 µm、遠位区間1.6677〜5.0115 µmの族では、相互角20°以内の最小の全区間最大広がりが.4118 µmで、通常らせん半径.2865 µmより大きい。**この限定族では**全体を内側へ傾けるだけで全長を束にできず、近位曲がりや軸の並進を区別する必要がある。

一般の付着点まわり回転・非対称運動・柔軟な近位部・動的形状・実際のPCA評価の不可能証明ではない。接触・弾性成立も判定しておらず、実軌道の入口判定には用いない。初期frozen RPYではKの初期力の差は丸め誤差、Mのmotor仕事率はDJの3.696e-18 Wから4.262e-18 Wへ変化する。これは承認済みの仕事率差の確認であり、瞬間rateから動的な改善を断定しない。
