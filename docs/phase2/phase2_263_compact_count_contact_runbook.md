# Issue #263：従来compact・13付着配置×菌体–べん毛反発

## 合意した範囲

2026-10-20、ユーザーは密集・同位相初期配置の検証を子Issue #271へ切り分け、従来`compact_rear`による13付着配置×菌体–べん毛反発ON/OFFを本PR #264で予約まで進めると指示した。#271の対象モデル・profile・実装条件は保留であり、#263の完了依存にしない。#271は`execution:triage`で実装・実行しない。旧試作と予約#33の証拠は保存し、今回の予約を差し替えない。

今回予約するのは**26条件・各1τのshort screen**。2秒×26条件の設定も用意するが、短期QC・固定camera replayのユーザーレビューまで予約・開始しない。本スレッドのcs10接続・予約・開始の承認を用い、固定commitと範囲を開始前に明示する。Macではsimulationなしのpreflight/dry-runと短いtargeted testだけを行う。

## 条件と対照

| nf | 中央layerの付着slots |
|---|---|
| 1 | 0 |
| 2 | 01、02、03 |
| 3 | 012、013、014、024 |
| 4 | 0123、0124、0134 |
| 5 | 01234 |
| 6 | 012345 |

回転のみを同一視し、反射を区別した13代表。各ON/OFF対で初期座標・位相・その他物理設定を一致させ、`potentials.spring_spring_repulsion.body_flagella_enabled`だけを変える。べん毛間反発はON固定。本数が増えると総駆動トルクも増えるため、本数間の差を付着位置や総トルクから独立な効果とは解釈しない。nf1の束・べん毛間指標は該当なし。

共通力学は既存の3.5秒機序比較の`nf03__slots024__bundle`（compact単独）を引き継ぐ。`compact_rear`、hook中立化・初期直交化あり、各べん毛の局所位相offset=0。局所offset=0は空間的な同位相化を意味しない。密集座標ファイルは読み込まない。

motorは拡散segment重み（`fixed_weights`）と局所one-ringへの全ベクトル反作用。腕重み・横トルク補正・追加位置復元なし。解析torsion、各べん毛T=2.5e-20 N m、基準Tも同値、dt*=1e-4（4 µs）、τ=0.04秒。Brownian・motor切替OFF。追加位置復元scale=1は現行互換設定で、その処理とscale自体の削除は既に合意したマージ前作業として残す。今回その削除によるモデル差を混ぜない。torsion追加比較はしない。

## 静的確認

共通geometry preflightは26/26成立。13対の物理設定・初期座標が一致し、nf3 slots024 OFFは過去compact単独と時間以外の全解決済み設定・座標が一致することを確認した。初期曲げ・二面角誤差は丸め誤差、初期弾性力のbead最大は7.72e-27 N以下。ONでは初期反発力の全配列ノルムが0.597〜1.400 pN、OFFでは0であり、hook中立化を総力ゼロとは解釈しない。この反発の初期過渡も短期screenの比較対象に含める。

simulationなし固定camera 3D・軸投影で13形状を確認した。証拠は`outputs/2026-10-20/213900/issue263_compact_count_contact/`の`preflight_audit.json`、`short_dry_run.json`、`long_dry_run.json`、`preview/manifest.json`。新しい初期状態を作らず、既存compactの条件展開として扱う。

## stageと実行契約

| stage | campaign | parallel job | 条件数・step |
|---|---|---|---|
| short screen | `2010_hex_compact_count13_body_contact_1tau_issue263.yaml` | `compact_count13_contact_1tau_job.yaml` | 26×10,000 step（0.04秒） |
| long checkpoint | `2010_hex_compact_count13_body_contact_2s_issue263.yaml` | `compact_count13_contact_2s_job.yaml` | 26×500,000 step（2秒） |

campaignは`conf/phase2_multi_run/`、jobは`conf/phase2_parallel/issue263_initial_bundle/`。両stageの物理設定・初期配置は同一で、時間と評価stageだけを変える。全26条件geometry preflight、独立output、`cs10_qualified`の`max_workers: auto`（最大8 workers）、数値library thread各1を使う。

短期wallは待ち時間を除いて約1〜3時間を粗い計画枠とし、各条件の実測で更新する。2秒は総13,000,000 stepの数日規模を見込み、短期の本数・ON/OFF別速度から再見積りする。既存#33がrunningの場合、新予約はその後続とし、既存予約・dispatcherを変更・停止しない。

## QC・replay・解釈

固定commit、26 condition ID、入力SHA-256を照合し、portable成果物をローカル同期してファイル集合・件数・SHA-256全件一致を確認する。各`run_summary.json`を先に読み、共通`development_evaluation`のfinite/body/hook length/flag/motorのstrict QCを行う。hook角の診断除外と既存閾値を維持し、終端回復を全時間PASSとは扱わない。FAILなら原因を記録して2秒を保留する。

本数別heatmap、固定camera・ON/OFF共通視野の3D/2D replayを生成する。QC-PASSでは束幅・べん毛間軸角・菌体軸角・接近/接触、べん毛回転・菌体材料点自転、符号付き前進、後退割合・持続時間、首振りと横変位を記録する。位相PCAの見かけの回転と材料点自転を区別する。shortは初期過渡・数値安定性のscreenであり、定常性を判定しない。

2秒checkpointでは全区間・共通の後半1〜2秒・旧事象を含む0.64〜0.80秒を比較し、初期2回転を除いた周期平均も併記する。観測周期数不足は未判定として記録する。nf3の既存compact周期約0.37秒から約5回転を見込めるが、他本数へ外挿して成立としない。完全対称条件の軸一定や回転一致を自然な揺れ・自発同期の証拠にしない。科学的採択、canonical freeze、dataset採択、mergeは含めない。

## 短期予約 #34

2026-10-20 21:43 JST、固定commit `71ed7529d362474235b931494af6cafcf3fe2edf`、`compact_count13_contact_1tau_job.yaml`で予約#34を作成し、queuedを確認した。Mac/cs10 dry-runは26 ID・実効8 workers・thread各1・job/campaign SHA-256が一致。幾何値はplatform差の丸め誤差（rtol/atol=1e-12）内で一致した。job SHA-256は`864d8773425520eb0d8999a1b9fd707c2d4b933efd25ffa06241a2805b45d9f9`、campaign SHA-256は`745360e95d0e4d9894605ceaad5e9e2d764336bb616774bcc92b7bfddc8e4362`。

同時刻の先行#33は23/26完走、最遅nf6 ONが2,500/10,000 step・2.414 step/sで、平均速度外挿の完了見込みは22:21 JSTごろ。#34は先行終了後の開始待ちであり、22:21開始を保証しない。既存dispatcherを用い、#33・queue・実行条件は変更しない。この後のbranch更新は予約記録だけで、#34の固定commitを差し替えない。次回22:40 JSTごろ開始状態と自条件の速度を確認し、完了見込みに応じて同じ監視を再設定する。2秒stageは未予約。

## 短期screen完了

予約#34は2026-10-21 00:10:43 JSTに26/26完走。portable291ファイル・8,974,401 bytesのSHA-256全件一致、共通strict QC26/26PASS、14本の固定camera3D/2D replayとarchive診断を完了した。結果は[`phase2_263_compact_count_contact_short_result.md`](phase2_263_compact_count_contact_short_result.md)。完全周期0、全条件で初期の前進変位は負。2秒stageは未予約でユーザーレビュー待ち。
