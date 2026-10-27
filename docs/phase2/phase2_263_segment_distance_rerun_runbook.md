# Issue #263：線分間距離修正の実行計画と結果

## 現在の状態・承認範囲

**予約#37の1τ screenは2条件とも完走し、2026-10-27にportable照合・共通strict QC 2/2 PASS・固定camera replay/定量診断を完了した。**結果の正本は[1τ screenの結果](phase2_263_segment_distance_1tau_result.md)。1秒stage・26条件2秒は未予約で、QC/replayのユーザーレビュー待ち。

2026-10-25に距離ソルバー修正と実行準備、10月26日に予約/正常開始を確認した。同日13:40 JST予定の1回確認は22:40 JSTに起動したが接続不可で、automationは削除して追加確認しなかった。10月27日のVPN復旧連絡・cs10操作許可に基づき、今回完了確認と評価を実施した。

実行targetはIssue #263の`cs10_user_run` / `execution:cs10`。予約#37の固定commitは **`0017fbbd922997f32dbad4b40aa5af03da974c92`**、対象は`nf6_segment_distance_1tau_job.yaml`。予約前にremoteでもsimulationなしdry-runと入力SHA一致を確認した。記録用のPR headが進んでも、この予約commitを置換しない。

## 予約#37と実行範囲

- ユーザーの2026-10-26の予約指示を根拠に、12:36:29 JST、priority 0でenqueueした。直前のqueueはunpaused、他のrunning/queuedは0だった。
- 既存dispatcherを利用し、既存予約の取消・置換・並べ替え、dispatcherの起動/停止/pause/resumeは行っていない。launch時刻は12:36:53 JST。
- fixed commit、2条件のID、実効2 workers・各数値thread 1、2出力への分離をremoteで確認した。開始確認時のjobはrunning、failed configs 0。まだcheckpoint進捗がないため、全stepのQCは未評価。
- 予約時の完了目安は事前費用見積りから**同日13:40 JST前後**（実測未取得）。動的形状・共有負荷で変わる。同日13:40 JSTに1回だけ確認し、完走なら照合・QC・replayまで進む。未完走/接続不可/FAILでも状況を報告して確認を終了し、自動再スケジュールはしない。
- remote root：`/net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs/2026-10-26/123653/parallel/queue-00037-09668d5a__03cd2c703b5b`。
- ローカル証拠：`outputs/2026-10-26/123600/issue263_distance_screen_reservation/`。予約・launch・remote dry-runを保存した。後続の記録commitは実行commitと区別する。

## 目的と変更

[原因監査](phase2_263_off_spin_cause_audit.md)で、nf6の菌体–べん毛反発OFFにも全鎖の回転鈍化を認めた。近位segment接触・弾性抵抗に加え、距離計算の単位依存を再現した。今回の目的は、正しい距離計算で接触・減速がどう変わるかと、時間刻みへの依存を分けて調べること。

`src/sim_swim/dynamics/forces.py::_closest_points_on_segments`を変更した。

- `min ||p1+s(q1-p1)-p2-t(q2-p2)||, 0<=s,t<=1`を解く。
- segment長で座標差を正規化し、mとµmで分岐が変わる絶対閾値を除く。
- 外積による内点解で、ほぼ平行な場合の`a*c-b*b`の桁落ちを避ける。
- 内点解が範囲外の場合は4境界それぞれを最小化する。従来のs/t独立clampによる端点側の誤りも修正した。
- 点に縮退したsegment、平行・共線segmentも扱う。総力・総トルク保存と、回転/平行移動/端点順序の整合をtestした。

このhelperは全profileのsegment反発で共有されるため、修正は菌体–べん毛とべん毛間の両方に適用される。旧solverを切り替えるruntime optionは追加していない。過去の再現には固定commitを使用する。係数・cutoffの式、hook、motor、torsion、初期配置は変えていない。cutoffの不連続とsegment頂点付近の重複加算は、残る別の検証事項である。

## 比較条件

比較は3条件、**新規実行は修正後の2条件**。全てnf6 / slots012345、従来`compact_rear`、菌体–べん毛反発OFF、べん毛間反発ON。

| 比較条件 | 距離計算 | dt | 実行方法 |
|---|---|---:|---|
| 旧対照 | 修正前 | 4 µs | #35の保存済み0〜1秒を再利用 |
| `nf06__slots012345__distance_fixed_dt4us` | 修正後 | 4 µs | 新規 |
| `nf06__slots012345__distance_fixed_dt2us` | 修正後 | 2 µs | 新規 |

旧対照は固定`71ed7529d362474235b931494af6cafcf3fe2edf`、`2010_hex_compact_count13_body_contact_2s_issue263.yaml`の`nf06__slots012345__bfoff`。26/26 strict QC PASSと全portable照合済み。過去2秒の完走結果を、同じ時刻窓に切って比較する。新campaignのQCへ別durationの旧manifestを混ぜたり、IDを書き換えて共通evaluatorへ適合させたりしない。

旧対照とのconfig比較では観測時間の入力表現とdt以外が同一、同じbuilderでの初期座標は完全一致。cs10保存座標とMac生成座標の最大差は約1.69e-21 mで、浮動小数点の丸め程度だった。Brownian/motor切替OFF、各motor torque `2.5e-20 N m`、diffusive重み、全ベクトルone-ring反作用、解析torsion、腕重み/横トルク補正/追加位置復元なし、local phase offset各0を維持する。

## Stageとparallel config

| Stage | 共通contract stage | 観測時間 | 4 µs / 2 µsのstep数 | job YAML |
|---|---|---|---|---|
| 最初のscreen | `short_screen` | 1τ ≈0.04秒 | 10,000 / 20,000 | `conf/phase2_parallel/issue263_initial_bundle/nf6_segment_distance_1tau_job.yaml` |
| 機序比較 | `long_duration` | 25τ ≈1秒 | 250,000 / 500,000 | `conf/phase2_parallel/issue263_initial_bundle/nf6_segment_distance_25tau_job.yaml` |

campaignはそれぞれ`conf/phase2_multi_run/2010_hex_nf6_segment_distance_{1tau,25tau}_issue263.yaml`。τは約0.04秒、25τのfloat上の終端は0.9999999999999999秒。新規profileではτ指定で上記step数を確認し、既存runtimeのceil処理は変更しない。

両jobは`execution.max_workers: auto`・`worker_policy: cs10_qualified`。2条件なので実効2 workers、OMP/OPENBLAS/MKL各1 thread。`geometry_all_conditions`のpreflightと、条件別の子run出力をdry-runで確認した。baselineは新規jobに入らない。

**1τの共通QC・固定camera replayをレビューしてから、1秒stageを予約・開始する。** 1τでは0.4〜0.8秒の減速機序は観測できないため、その結果だけで減速解消とは判定しない。後続の2秒×26条件の再実行は、このnf6対照のレビュー後に設計する。

## 準備時のローカル確認結果（動的screen前）

- 修正前の新規幾何test：9 FAIL / 18 PASS。SIの内点接触見逃し、端点clamp、ほぼ平行、座標変換を検出した。
- 修正後：新規幾何27ケース＋既存potential25ケース＋campaign3ケース、計**55 PASS**。
- 1τ・25τのparallel-job dry-run：両方`planned`、2 workers、条件IDと出力分離を確認。
- 共通geometry preflight/preview：2条件の初期座標が同じ。hook長誤差0、hook初期力0、非付着bead間の最小中心距離0.25 µm。
- production engineの`step(0)`：初期反発0、初期弾性力ノルム2.94e-26 N、べん毛間segment最小距離0.4221 µm（cutoff 0.2 µm）。位置更新なし。
- hook角の旧診断偏差51.77°は記録する。基本hookの片側力がゼロになる配置であり、共通contractどおりhook角をQCへ戻さない。
- 共通previewの3D像を確認。この準備時点では動的実行・QC/replayは未実施だった。完走後の評価は上記の結果文書を参照。

## 費用見積り（予約前の外挿）

旧nf6 OFFのcs10実測平均は5.5548 step/s。Macで4保存状態×3回、位置を進めずproduction evaluatorを呼び、旧helper/修正helperの費用比を比較した。median 0.0752 / 0.0772秒、修正/旧の比約1.027。これを旧cs10速度へ掛け戻した外挿は約5.41 step/s。

| Stage | 4 µs | 2 µs | 2条件並列のwall目安 |
|---|---:|---:|---:|
| 1τ | 約0.51時間 | 約1.03時間 | **約1時間** |
| 25τ | 約12.8時間 | 約25.7時間 | **約26時間** |

cs10実測benchmarkではない。共有負荷・修正後の形状/接触で変わり、予算幅は外挿の0.5〜2倍（1秒stage約13〜52時間）とする。待ち行列待ちとレビュー時間は含めない。Macでの長時間simulationは行わず、開始後の実測checkpointから更新する。

1τ完走後の実測は4 µsで6.237、2 µsで6.073 step/s。1秒への単純外挿は並列wall約23時間だが、接触発生後の速度は未測定。[結果・費用の限界](phase2_263_segment_distance_1tau_result.md)を参照。

## 評価手順

1. 固定commit・2条件・config SHA・予定step数を照合。portable同期のファイル集合・全SHAを確認し、各`run_summary.json`を詳細分析より先に読む。
2. 共通`development_evaluation`のfinite/body/hook length/flag/motor strict QC、heatmap、固定camera 3D/2D replay。既存閾値を維持し、FAIL時に閾値を緩めたり自動再実行したりしない。
3. QC-PASSのみ、各鎖の材料点spin/phase、菌体材料点自転、符号付き前進、束幅/軸角、根元横腕、近位segment接触/距離/反発力を同じ物理時刻で比較。
4. 旧4 µsと修正4 µsで距離計算の影響、修正4 µsと2 µsでdt依存を分ける。主要窓は0.30〜0.40秒、0.55〜0.65秒、0.90〜1.00秒。0.2〜1秒内の最小100 ms平均spin、1秒までの積算回転、全鎖のばらつきも報告する。
5. 1 ms保存状態の接触/力再計算を、4/2 µs内部stepの時間平均とは扱わない。回転は全internal step積算の共通recorderを使う。旧力の再現は旧commit、新solverで旧座標を評価した場合は静的な対照として区別する。
6. 固定時間窓の回転数・時系列・形状差を報告し、2つのdtだけで収束次数や定常遊泳成立を宣言しない。1秒は接触機序の比較であり、回転周期数が十分かは結果から別途判断する。

新たなQC evaluator、torsion方式比較、モデル採択、canonical freeze、mergeは含めない。マージ前の追加位置復元処理/scale削除は既存の残作業として保持する。

## 残タスクと判断の順序

| 順序 | 残タスク | 次へ進む条件・判断点 |
|---|---|---|
| 1 | **評価完了・ユーザーレビュー待ち**：#37完走照合、portable同期、共通strict QC 2/2 PASS、固定camera 3D/2D replay | 2条件ともQCを確認し、初期過渡に破綻がないことをレビューする。1τだけで減速解消とは判定しない |
| 2 | 同じ固定commitの25τ・新規2条件を予約・実行 | 1τのQC/replayレビュー後、対象jobと許可範囲を確認する。暫定wall約26時間 |
| 3 | 旧4 µs／修正4 µs／修正2 µsの機序比較 | 同じ0.55〜0.65秒のspin・接触・弾性抵抗、回復時刻を比較する。距離修正の効果とdt依存を分ける |
| 4 | 26条件×2秒の再実行条件を確定 | 採用するdt、修正後の共通モデル、本数/付着13配置×菌体–べん毛反発ON/OFFを確定し、費用を再見積りする。旧実行の承認から新モデルの予約を推測しない |
| 5 | 束維持・持続回転・推進安定性の成立を評価 | 各べん毛のspin、菌体材料点自転、符号付き前進を同時刻で確認する。初期緩和を除く完全周期数を数え、観測期間が不足すれば延長条件を相談する |
| 6 | マージ前整理・検証・レビュー | 合意済みの追加位置復元処理/scaleを削除し、参照config/test/docsを整える。過去証拠は固定commitで再現できるよう保持。required checksとreview gate、ユーザーの物理採択判断を経る |

**今回の1τ予約に追加の要件確認は不要**。上記の条件と許可範囲で確定している。事前に詰める必要があるのは、1秒比較を受けた後の採用dt・26条件への展開・必要観測周期数である。

距離修正後も減速が残る場合は、cutoffの不連続、頂点近傍のsegment対加算、hook/motor/流体相互作用を次の独立仮説として整理し、追加変更前に根拠と対照条件を提示する。一度に複数を変えない。torsion方式の追加比較は不要という合意を維持し、重み3方式比較は低優先度で保留。密集・同位相配置は別Issue #271の管理範囲とする。

Issueには重要な判断・実行状況・結論・阻害要因だけを短く残し、同じ話題は既存コメントを更新する。詳細なチェック・SHA・操作記録は本runbookとローカルmanifestを参照する。

## 証拠と再現コマンド

ローカル正本：`outputs/2026-10-25/223700/issue263_distance_rerun_preparation/`。
`short_dry_run.json`、`mechanism_dry_run.json`、`preview/manifest.json`、`preparation_evidence.json`、`prepare_evidence.py`、`manifest.json`、`run.log`を保存する。source/設定/testのSHAと準備commitはmanifestへ記録する。

以下は**ローカルdry-runだけ**を行う。

```bash
python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/issue263_initial_bundle/nf6_segment_distance_1tau_job.yaml \
  --dry-run --output-root outputs/<JST-date>/<time>/issue263_distance_short_plan

python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/issue263_initial_bundle/nf6_segment_distance_25tau_job.yaml \
  --dry-run --output-root outputs/<JST-date>/<time>/issue263_distance_mechanism_plan
```

予約承認後に固定commitでremote側の同じdry-runとqueue状態を確認する。予約・開始の操作はこの文書のコマンドから自動実行しない。
