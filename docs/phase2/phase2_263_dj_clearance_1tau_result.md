# Issue #263：非重複D/J 1τ screenと3.5秒予約

後続記録：main#42は11/7に4/4完走。[3.5秒結果](phase2_263_dj_3p5s_result.md)（11/11評価）を現在状況の正本とし、以下の予約/ETA/継続中表記は11/4時点の記録として保持する。

2026-11-04。修正後の予約#41は17:31 JSTごろに4/4条件が完走した。**共通strict QC 4/4 PASS、固定camera 3D/2D replayと保存点幾何を確認し、承認済み3.5秒4条件を予約#42として開始した。** モデル採択・定常遊泳の判定は未実施。

## 固定実行と照合

short/mainとも固定実行commitは`59890c05460ba3032d50813afaa488ed3440f271`。nf3/slots024、R_clearance_control・D_clearance_vector・J_transition・DJ_combinedの4条件。T=2.5e−20 N m/本、dt4 µs、粘性1 mPa s、FF/BF ON、Brownian/切替OFF、数値thread各1/4 workersを維持した。

#41は全条件10000連続step、実終端0.03999999999999624秒。portable37 files・946942 bytesの集合/件数/全SHA一致を確認した。run_summaryを詳細解析前に読み、全内部stepでfinite/body/nonbodyにFAIL履歴なし。共通development_evaluation strict QC 4/4 PASS、nf3 heatmapを生成。旧#40や距離バグを含む過去archiveを新しいQCへ混在させていない。

short job SHA256は`c0794b9f57e13679439e66ac8efb840288cae035af88bc344c777a3714748f8b`、campaignは`49dd38d724ac7cb31fa5f792ddc251405f6a260a6b10d60bf6f102a9186c2992`。

固定camera共通視野3D/2D replayは各41frame/25fps、全4条件0〜.04秒を収録。ffprobe/manifest/frame/fps/SHAを照合し、3D中間/終端と2D中間を視認した。目立つ崩壊・脱離・交差・視野切れはない。2Dは従来通り菌体投影。

## 初期過渡の定量

| 条件 | べん毛0の菌体相対回転 [回] | 菌体材料点自転 [回] | 前方変位 [µm] | 終端束半径指標 [µm] | 保存点BF最小距離 [µm] |
| --- | ---: | ---: | ---: | ---: | ---: |
| R 非重複compact | .116562 | .044640 | −.020500 | .958631 | .199968 |
| D 非重複compact・有効motor | .137440 | .050161 | −.016794 | .976971 | .204788 |
| J 根元移行部 | .120831 | .059128 | −.002761 | .748119 | .200133 |
| D+J | .184970 | .054779 | +.001198 | .758047 | .213832 |

flag位相は内部step積分したlabel付き鎖の剛体fitから菌体運動を除いた値。body材料点自転は保存1msごとのlabel付き菌体fitを長軸へ投影した値で、PCA方位角のjumpではない。両者の軸方向の定義は異なり、符号を単純に加算しない。全条件1回転未満で、初期2周期を除いた同位相比較はできない。小さいD+J前進だけで持続推進を宣言しない。完成した密束の維持を証明した結果でもない。

最大motor横トルク/TはR=1.355、D≈5.25e−15、J=88.775、DJ≈6.16e−15。最大flag torsion偏差はR/D/J/DJ=3.078/1.030/59.122/.598°。Jの大きい不要横トルクと過渡変形を保持し、長期での停止・変形との対応を観測する。D/DJの力/トルク仕様成立と、目標遊泳の成立は別である。

## geometry/replay gate

初期のR/D非隣接BF距離は.204788 µmで、旧root-on-body-edgeの距離0は解消した。保存41状態でもBF距離はすべて正で、中心線交差は検出されなかった。

Rの保存点最小距離.1999682068 µmは2a=.20 µmより約.0000318 µm（32 pm）小さい。これはsoft反発が働く動的接触圧縮として記録し、以前の初期中心線一致とは区別した。**初期preflightの非重複条件は厳格に保持し、動的なhard非貫通を証明したとはしない。** 1ms保存点の監査で、未保存内部stepの非貫通までは保証しない。共通QCを緩めず、明白な崩壊・貫通/実装異常が見られないことを固定camera replayと合わせてレビューし、long gateをPASSとした。

この判断はユーザーが承認したshort QC/replay→3.5秒の範囲内で、条件縮小・物性/閾値調整・自動再実行を行っていない。

## main予約#42

- 2026-11-04 17:40:36 JSTにpriority0でenqueue、17:40:50 JSTに正常開始。事前のrunning/queued予約は0、既存dispatcherを利用。4条件すべてrunning・failed0、fixed commit/入力SHA/4workers/thread1を照合した。
- job：`conf/phase2_parallel/issue263_initial_bundle/nf3_dj_3p5s_job.yaml`。各3.5秒/875000step、同じ4条件。job SHA256 `12f8176e93ddade147c45555e79d1c7ca7f87bf6bf663f15cee52022b5f25b58`、campaign SHA256 `a71a51c4aec35716b31aa80189e6ee78197cecd8bc525a396b38555e4372056d`。Mac/cs10で直前dry-runを再確認した。
- main出力：`/net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs/2026-11-04/174050/parallel/queue-00042-62fbb2d2__e2f77772d674`。
- 短期実測のwallはR2062.95/D2008.76/J2081.73/DJ2005.25秒。875000/10000倍の外挿で各48.74〜50.60時間、最遅の終了見込みは**11/6（金）20:17 JSTごろ（報告目安20:20）**。共有負荷・長期形状で変わる暫定値で、長期自条件の進捗速度はまだ未取得。
- fixed commitを直接指定した予約であり、後続PR headの記録commitに追従/置換しない。既存予約のcancel/置換/優先度変更はしていない。

今回のフォロー目的は3.5秒予約の正常性確認まで。自動確認`issue263-d-j-qc-3-5`は目的達成後に削除した。main計算は継続中で、完走QC・定性/定量評価・必要周期/摂動・数値感度の判断は残る。時間長だけで定常性・canonical freezeを宣言しない。追加位置復元処理/scale削除はmerge前残件。

ローカルshort正本は`outputs/2026-11-04/173618/issue263_clearance_short_review/`。portable_index、run_summary_first、evaluation/replay、bundle_rotation、short_diagnostics、replay_verification、short_screen_review、long dry-run、再現script、manifest/review_result/run.logを保全する。main予約証拠は`outputs/2026-11-04/174036/issue263_clearance_main_reservation/`。

[初期配置修正と実行手順](phase2_263_dj_clearance_runbook.md)／[受入案・機序窓](phase2_263_model_repair_direction.md)／[ADR0029 Proposed](../adr/0029_phase2_effective_torque_and_root_transition_direction.md)。
