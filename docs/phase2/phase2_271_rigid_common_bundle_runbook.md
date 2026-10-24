# Issue #271：全べん毛後方束RUNの形状保存候補・予約前runbook

## 目的と境界

対象は「全べん毛が後方束を作るRUN」という一条件であり、全遊泳状態や本数差の一般化ではない。2010 hex（30 body beads、各flag 11 beads）の`paper_table1`全鎖を剛体変換し、同じ真の離散らせん軸へ向ける。有限長PCA軸は形状生成には使わない。後方束を**初期に与える**条件付き検証であり、自発的束形成の証拠ではない。

このrunbookの実装基点は未mergeのPR #264 head `82962a3cd1e8a467e968c7b8b84fb645021a58a4`。#271 PRは同branch宛てに積み、#264側の変更やmergeがあれば、予約前に実効設定・座標hash・dry-runを再照合する。既定profileは変更せず、`rigid_common_axis_rear`をopt-inにする。

## 形状、接続近似と静的成立範囲

`conf/phase2_initial/2010_hex_rigid_common_bundle_v1.json`を座標の正本とし、SHA-256は`8844f37f5495ab66d19ab0a1c174e51c4bd42b8f1e6c837a85b3952981ca991d`。生成CLIは`design_rigid_common_bundle.py`。normal離散鎖の真のらせん半径は0.286543 µm。bead 0を中央ring上の各slotから方位+19.797585°、共通後方軸上の同じ長手位置へ置くとhook長0.25 µmとなる。n=3 slots024の摂動はflag 0全鎖をbead 0を通るhook軸の周りに1°だけ回す。内部形状とhook長は維持され、比較対象の他2本は動かない。

菌体との貫通許容はhookと最初の3 shaft segment（bead 0→3）だけ。そこでは最短body-edge距離が0.148160 µmまで下がるため、当該候補は菌体–べん毛反発OFFと対にする。bead 3以降のshaft segmentはbody-edgeから0.448376 µm以上離れる。べん毛同士はhookを含む有限segment間で直径/反発cutoff 0.2 µm以上を必須とする。内向きhookは外向き半径方向から約157.16°であり、接続の物理近似として結果の解釈に残す。

13代表配置（n1=0、n2=01/02/03、n3=012/013/014/024、n4=0123/0124/0134、n5=01234、n6=012345）と摂動1件の静的auditを`outputs/2026-10-23/155900/issue271_static_geometry/`へ保存した。14件すべてで内部pairwise距離誤差は最大1.8e-15 b、hook長誤差は最大1e-16 b、自己非隣接segment距離は0.58 µm以上。別べん毛segment間の最小値は0.20151677 µm（隣接slotsを含む配置）、n2 slots02は0.38425648 µm。初期normal曲げ・torsion誤差は浮動小数点丸め内で、C力学の初期弾性力最大は1.3e-26 N未満、初期弾性エネルギーは1.8e-47 N m未満。各配置のbead/segment距離、hook、力・エネルギーは同rootの`static_geometry_audit.json`と`manifest.json`を参照する。余裕の小さい隣接slotsは動的QCを保証しない。

静的な固定camera 3D・軸方向投影の5条件previewは`outputs/2026-10-23/160222/issue271_initial_preview/`に保存した。これはt=0の可視化であり、動的束維持の証拠ではない。

## 旧FAILの再現と原因

旧予約#33は固定commit`50445ec4b6c8385fc7bc08c354cbc24ac0e3ea03`、nf02 slots02 ON/OFFとも8 µs（2 step）にstrict flag QC FAIL。Macの現行#264基点でも同じbond 3.027227、bend 113.195619°、torsion 163.454691°を再現した。旧座標のflag 1の根元側曲げ角は初期179.224953°、1 step後179.787060°で、torsion planeの最小sineが約0.003717まで減る。解析torsion力のbead最大は次stepで8.394e-10 Nに増え、最大変位1.5465 µmとなる。motor力は約7.46e-14 N、反発力はON/OFFとも0。したがって旧接続部のほぼ共線な四点と残留torsion応力が、解析torsion力の鋭い増幅を起こす直接の数値的原因である。

切分けでは、同じ4 µs内部刻み・同じ弾性スケールでmotorだけOFFにしても2 step目にFAIL。dtを半分の2 µsにした旧配置は同じ8 µsまでstrict形状QCを通り、torsion剛性だけ0にした診断対照も2 stepを通る。いずれも採用する物理変更ではない。対策は接続部を変形させない全鎖剛体配置で、同じnf02 slots02・C力学・4 µs刻みに戻して10 stepを確認した。新配置の初期最大曲げ角は142°、10 step時の最大bond相対誤差は0.000996で、5短期候補すべて最初の10 stepはstrict形状QCを通過した。これは予約後の1τ全時間QCの代用ではない。一次記録は`outputs/2026-10-23/160000/issue271_initial_steps/manifest.json`と各条件の`initial_step_diagnostic.json`。

## 予約対象と確認点

短期campaignは`conf/phase2_multi_run/2010_hex_rigid_bundle_1tau_issue271.yaml`、parallel jobは`conf/phase2_parallel/issue271_rigid_bundle/short_1tau_job.yaml`。各0.04秒/10,000 stepで、独立5条件は`nf01__slots0__base`、`nf02__slots02__base`、`nf03__slots024__base`、`nf06__slots012345__base`、`nf03__slots024__tilt1deg`。T=2.5e-20 N m/本、`dt_star=1e-4`、Brownian/switching OFF、Cの局所one-ring全vector反作用・解析torsion・追加frame位置復元なし・腕重み/横トルク補正なしを使う。菌体–べん毛反発OFF、べん毛間反発ON。共通`development_evaluation`のfinite/body/hook length/flag/motor strict QCを維持し、hook角は診断表示のみとする。

Macのsimulationなしdry-runは`outputs/2026-10-23/160001/issue271_dry_run/`に保存。5条件すべてgeometry preflight PASS、独立output、実効5 workers、数値ライブラリthread各1。job SHA-256 `931359e7c2ec2cd7da2533795a74577d0c2d4855614ecf54002c74289a193f1b`、campaign SHA-256 `896fee4f861ad8d5b39560d91a7b5236cc55ec504dba0687ca4be6d032efff77`。#263のnf6・2秒実測約21.7時間を単純に50分の1へ換算すると、cs10短期nf6は約26分が下敷きとなる。5 worker並列でも共有負荷・新幾何の影響を含むため、初期wall計画枠は約30〜60分、Macでの全campaignは30分超と見込む。短期実測で更新する。

予約前にはPRの最終commitと上記3入力SHA、Mac/cs10のdry-run条件・output対応、cs10のNAS mount/容量・queue状態を照合する。cs10接続、queue enqueue、dispatcher/tmux、job開始はそれぞれユーザーの当該操作への明示承認が必要。予約のfixed-commitはbranch更新で差し替えない。

## 2026-10-24 Mac 10,000-step確認と予約状態

ユーザーの追加依頼により、上記5条件をMacで各10,000 internal step（各1τ、0.04秒）まで実行した。固定commitは`5ac720ce9de2f86b023e0c03d25541f8d04fbfc5`で、PR #264 headは`82962a3cd1e8a467e968c7b8b84fb645021a58a4`のまま。座標・campaign・jobのSHA-256は上記から変わっていない。一次記録は`outputs/2026-10-24/141900/issue271_mac_10k/manifest.json`、条件別`run_summary.json`、集約campaign、`model_development_evaluation/evaluation_manifest.json`。合計50,000 stepの全時間集計でfinite/body/flagに観測FAILはなく、共通`development_evaluation`も**5/5 PASS**（hook角のみ診断）だった。

| condition ID | Mac wall time | 最大hook長相対誤差 | 最大shaft結合相対誤差 | strict QC |
| --- | ---: | ---: | ---: | --- |
| `nf01__slots0__base` | 154秒 | 0.0523 | 0.0342 | PASS |
| `nf02__slots02__base` | 243秒 | 0.0551 | 0.0395 | PASS |
| `nf03__slots024__base` | 357秒 | 0.0576 | 0.0383 | PASS |
| `nf06__slots012345__base` | 880秒 | 0.0632 | 0.0508 | PASS |
| `nf03__slots024__tilt1deg` | 358秒 | 0.0576 | 0.0382 | PASS |

5条件並列のMac wall timeは約15分。固定cameraの3D replay終端ではnf2/nf6の後方らせん束が視認できた。2D replayは既存renderer仕様により菌体投影のみを描くため、べん毛束の2D可視性を立証しない。8本のreplay動画（本数別3D/2D）と入力・成果物SHA-256を同manifestに固定した。短期QC PASSは複数周期の束維持・菌体自転・符号付き前進の判定を意味しない。

ユーザーはMac確認後のcs10予約を許可した。最初のSSH試行はMac既定のユーザー名`kentatakemori`を使ったため拒否されたが、ユーザー指摘を受けて正しい`Ktakemori@cs10`で接続できた。2026-10-24の予約前照合ではqueueはpauseなし、予約#35がrunning。NAS mountとproject出力先の書込み権限は正常で、空き容量約47 TB、cs10 Python 3.11.15とGitHub CLI認証も確認した。PR #264 headは`82962a3cd1e8a467e968c7b8b84fb645021a58a4`で不変、予約直前のPR #272 headは`fc320bd0b2c7b84b64dfa35c8477675ce59a549e`。3入力SHAも上記の記録と一致した。

cs10のsimulationなしdry-runは`outputs/2026-10-24/143900/issue271_cs10_preflight/`に保存した。全5条件のgeometry preflightとoutput分離、実効5 workers、数値ライブラリthread各1を確認した。dry-run planのSHA-256は`6ac2e23cb43939739788319ee3e26c6a12334cdef95cc29a7bf17165fffdc016`。`queue.py enqueue`で**予約#36**（priority 0）を作成し、固定commitは`fc320bd0b2c7b84b64dfa35c8477675ce59a549e`、configは`conf/phase2_parallel/issue271_rigid_bundle/short_1tau_job.yaml`。直後のqueue照合では#35がrunning、#36がqueued、pauseなし。dispatcher/tmux起動とjob開始はCodexから実施していない。以後の文書commitでPR branchが進んでも予約#36のcommitは保持し、差し替え・取消・pause/resumeは別途明示許可が必要。

完走後はportable artifactをローカルへ同期して件数・SHA-256・QCを確認し、各`run_summary.json`を先に読む。共通evaluatorの全時間strict QCと本数別heatmap、固定camera 3D/2D replayをレビューしてから長時間stageを設計する。複数完全周期では初期緩和を挙動評価から分けてもQCの全時間FAIL履歴を保持し、全べん毛の束幅・軸角、位相差、べん毛回転、菌体材料点自転、周期平均の符号付き前進を同時刻で評価する。短期だけで遊泳再現を採択しない。
