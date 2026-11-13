# Issue #263：実線分接触対照の短期結果と接続経路の静的監査

2026-11-13。Refs #263 / PR #264。[実行手順](phase2_263_connector_execution_runbook.md)の現J対照shortを代理レビューし、同じ固定sourceの1秒比較へ進めた。科学的採択・束形成の達成報告ではない。

## #45：2条件とも短期QC・replay PASS

固定simulation/analysis commit `823f81873aa601e7a12ad026c80234db16f7fc2d`。nf3/slots024、現Jの座標・自然参照、D drive/Kなしを共通とし、菌体–べん毛反発BFだけをON/OFFにした。べん毛間反発FF ON、T/refT2.5e−20Nm/本、dt4µs、η1mPa sを固定。各10000連続step、実終端.03999999999999624秒、cs10_qualified実効2workers/thread各1。

| 項目 | BF ON | BF OFF |
| --- | ---: | ---: |
| 共通strict QC | PASS | PASS |
| max bend誤差 [°] | 2.56614 | 2.56614 |
| max torsion誤差 [°] | .598079 | .598079 |
| max hook相対伸長 | .0351205 | .0351205 |
| max body相対stretch | .0295756 | .0295756 |
| 実BF最短距離 [µm] | .213831995 | 対象pairなし |
| 実BF cutoff内時間 [s] | 0 | 対象pairなし |
| 実interFF最短距離 [µm] | .727881129 | .727881129 |
| child wall [s] | 2134.420 | 564.096 |

新しい実線分診断は全10000 force evaluationで記録された。BF ONは3342候補pairで距離が2a=.20µmを上回り、初接触時刻のsentinel−1・cutoff内時間0。OFFの距離−1は対象外のsentinelで物理的な距離ではない。全41保存点のON/OFF座標・時刻arrayはbit単位で一致した。従って、この短期窓ではまだBFの作用後を比較できない。

portable **21files/537390bytes**の集合・件数・全SHA-256をremote前後とlocalで照合し、両`run_summary.json`を詳細分析より先に読んだ。共通`development_evaluation` strict QC・本数別heatmap・標準固定camera 3D/2D replay、共通`bundle_rotation`を実施。replayは41frame/25fps、物理時刻0–.04秒、動画長1.64秒。ffprobe/coverage/SHAを検証し、3Dの0/20/40frameと2D終端を視認した。明白な崩壊・交差は見えなかったが、保存点から内部全stepの非貫通は保証しない。初期配置は広いJであり、密束ではない。

ローカル正本：`outputs/2026-11-13/190235/issue263_contact_short45_review/`。`review_result.json: PASS`は代理短期ゲートのPASS、`scientific_acceptance:false`。replayは同directory内`evaluation/replay/dt1e-04/nf03/grid_swim3d.mp4`と`grid_projection2d.mp4`。

## #46：同じ固定commitの1秒比較を開始

#45レビュー完了後、既存queue/dispatcherと入力を再照合し、11/13 19:07:00 JSTに予約#46、19:07:22に正常開始した。他のactive予約はなく、既存dispatcherを利用。他予約や実行条件を変更していない。

- job：`conf/phase2_parallel/issue263_initial_bundle/nf3_connector_contact_1s_job.yaml`
- job SHA-256：`c4b1a9b841d528909938e5013b1e65041e4b7cf0753a74c54c39b30bcef0c863`
- campaign SHA-256：`f7d8beea326db711f28972273125d33ec0417adb033c37b04e89013666869b6f`
- IDs：`nf03__slots024__J_bfon_control` / `nf03__slots024__J_bfoff_control`
- 各250001連続step、公称1秒・予定終端1.000004秒。比較窓は1.0秒まで。既存ceilによる1step差を隠さない。

short自身の平均速度を用いる初期予測はOFFが11/13 23:02 JST、ONが11/14 09:57 JST前後。長時間形状・共有負荷によって変わるため、実測checkpointと台帳を更新する。予定時刻は終了保証ではない。

両条件の実接触開始・根元の内部変形/全体傾斜・束幅/相互角・材料spin・菌体自転・符号付き前進を同時刻で比較する。BF接触が初動原因か、開いた後の増幅/結果かを分ける。完走後のQC/replay・機序評価は未実施。今回2/7予約を使用した。

## 静的監査：初期配置だけでは回転成立を示せない

ここは時間積分を伴わない有限範囲の幾何検算である。自然参照・弾性/散逸・連続回転の合格証拠ではなく、特定の探索失敗を全配置の不可能性とも扱わない。証拠は準備outputの`G1_*` / `G2_*` JSONとscript・NPZに残した。

| 参照案 | 接続link/contour | 有限範囲の結果 | 未確定の点 |
| --- | --- | --- | --- |
| G1 | 現行3×.58=1.74µm、11bead/本 | 全3鎖を満たす初期案なし。個別の菌体非接触経路を組み合わせても鎖0/2で接触 | 有限family外の可能性、連続spin/仕事/受動応答 |
| G2A | 4×.435=1.74µm、12bead/本 | 選択した鎖0/2の168組が全て未達。最良でも最短.139612µm<.20µm | bead増加だけではこのfamilyの障害を解消しなかった |
| G2B（初期幾何の参照のみ） | 4×.58=2.32µm、12bead/本 | 初期全active線分clearanceと近束幾何を満たす一例を得た | 全長/配置変更、連続spin、固定自然参照、材料/流体/drive対応は未検証 |

現Jの遠位8点・局所phaseを保って軸を対称に寄せる限定familyでは、軸半径.10〜.286543µmの試験点に遠位らせん同士の重複があった。一方、軸半径0/.05µmは初期遠位区間だけなら非重複だった。後者へ接続する個別経路を調べたところ、G1の例ではBF最短.204788µmを確保してもFF最短.022960µmとなった。菌体を避けるだけでは全鎖経路を満たさない。

G2Aの衝突witnessは全168組で鎖0の最初の遠位segmentと鎖2の最後の接続segmentだった。単に接続link同士だけを検査しても見落とす。各鎖の経路選択を独立に最適化するのではなく、接続–遠位を含む全組み合わせのclearanceが必要である。

G2Bはhook.25µmと遠位8点を保ち、共通横移動.30µm、後方shift0の一例。初期BF最短.204788µm、FF最短.240456µm、遠位幅.084616µm、軸間最大角4.093°、後方軸角.429°。ただし全filament長が5.80→6.38µm、bead数11→12に変わり、G2Aからは遠位の後方配置も変えた。成立をlink追加だけの効果と帰属しない。

**G2Bを動的候補や予約可能なprofileとして実装していない。** 長さによる抗力/RPY・受動材料応答・全鎖PCAによるD drive軸と仕事の対応、固定自然参照の連続spin経路を静的に検証してから、必要な構造変更と比較条件を具体化する。G1へ黙って置き換えず、現Jの1秒対照を進めながら設計根拠を積む。

## 残作業

#46の完走検証・共通QC/replay・同時刻機序評価、G1/G2の回転/受動応答設計が残る。1秒の後半で既定の束入口を満たさない候補を3.5秒へ延長しない。11/29までの自律範囲と予約上限は[実行手順](phase2_263_connector_execution_runbook.md)を維持する。ADR0029 Proposed。追加位置復元code/scale削除はmerge前残件。
