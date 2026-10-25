# Issue #263：線分間距離修正の再実行準備

## 現在の状態・承認範囲

2026-10-25のユーザー依頼「再実行が必要。準備。予約は許可してから」に基づく。距離ソルバーの修正、targeted test、simulationなしの初期preflight、preview、parallel-job dry-runまで準備した。**cs10への接続・予約・開始はこの準備では行っていない。全jobは未予約**。監視automationもPAUSEDのままとする。

実行targetはIssue #263の`cs10_user_run` / `execution:cs10`。予約前にこの準備を含む**固定commitと対象stage**を提示し、ユーザーの明示許可を受ける。許可前にenqueueやdispatcher操作をしない。現在のPR head候補はPRコメント・ローカルmanifestの`prepared_commit`を正本とする。将来branchが進んでも予約commitを推測で置換しない。

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

## ローカル確認結果

- 修正前の新規幾何test：9 FAIL / 18 PASS。SIの内点接触見逃し、端点clamp、ほぼ平行、座標変換を検出した。
- 修正後：新規幾何27ケース＋既存potential25ケース＋campaign3ケース、計**55 PASS**。
- 1τ・25τのparallel-job dry-run：両方`planned`、2 workers、条件IDと出力分離を確認。
- 共通geometry preflight/preview：2条件の初期座標が同じ。hook長誤差0、hook初期力0、非付着bead間の最小中心距離0.25 µm。
- production engineの`step(0)`：初期反発0、初期弾性力ノルム2.94e-26 N、べん毛間segment最小距離0.4221 µm（cutoff 0.2 µm）。位置更新なし。
- hook角の旧診断偏差51.77°は記録する。基本hookの片側力がゼロになる配置であり、共通contractどおりhook角をQCへ戻さない。
- 共通previewの3D像を確認。科学条件の時間発展、cs10でのbenchmark、QC/replayの実行結果はまだない。

## 費用見積り

旧nf6 OFFのcs10実測平均は5.5548 step/s。Macで4保存状態×3回、位置を進めずproduction evaluatorを呼び、旧helper/修正helperの費用比を比較した。median 0.0752 / 0.0772秒、修正/旧の比約1.027。これを旧cs10速度へ掛け戻した外挿は約5.41 step/s。

| Stage | 4 µs | 2 µs | 2条件並列のwall目安 |
|---|---:|---:|---:|
| 1τ | 約0.51時間 | 約1.03時間 | **約1時間** |
| 25τ | 約12.8時間 | 約25.7時間 | **約26時間** |

cs10実測benchmarkではない。共有負荷・修正後の形状/接触で変わり、予算幅は外挿の0.5〜2倍（1秒stage約13〜52時間）とする。待ち行列待ちとレビュー時間は含めない。Macでの長時間simulationは行わず、開始後の実測checkpointから更新する。

## 完走後の評価

1. 固定commit・2条件・config SHA・予定step数を照合。portable同期のファイル集合・全SHAを確認し、各`run_summary.json`を詳細分析より先に読む。
2. 共通`development_evaluation`のfinite/body/hook length/flag/motor strict QC、heatmap、固定camera 3D/2D replay。既存閾値を維持し、FAIL時に閾値を緩めたり自動再実行したりしない。
3. QC-PASSのみ、各鎖の材料点spin/phase、菌体材料点自転、符号付き前進、束幅/軸角、根元横腕、近位segment接触/距離/反発力を同じ物理時刻で比較。
4. 旧4 µsと修正4 µsで距離計算の影響、修正4 µsと2 µsでdt依存を分ける。主要窓は0.30〜0.40秒、0.55〜0.65秒、0.90〜1.00秒。0.2〜1秒内の最小100 ms平均spin、1秒までの積算回転、全鎖のばらつきも報告する。
5. 1 ms保存状態の接触/力再計算を、4/2 µs内部stepの時間平均とは扱わない。回転は全internal step積算の共通recorderを使う。旧力の再現は旧commit、新solverで旧座標を評価した場合は静的な対照として区別する。
6. 固定時間窓の回転数・時系列・形状差を報告し、2つのdtだけで収束次数や定常遊泳成立を宣言しない。1秒は接触機序の比較であり、回転周期数が十分かは結果から別途判断する。

新たなQC evaluator、torsion方式比較、モデル採択、canonical freeze、mergeは含めない。マージ前の追加位置復元処理/scale削除は既存の残作業として保持する。

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
