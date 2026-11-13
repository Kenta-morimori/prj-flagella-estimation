# Issue #263：接続経路・菌体接触の自律検証手順

2026-11-13。Refs #263 / PR #264。ユーザーは「実験構成まで進め、予約まで可、不在中の検証と再来週までの予約を自律的に」と明示した。前段の[計画](phase2_263_connector_repair_plan.md)の予約保留は今回の範囲で解除された。前フェーズM/Kは未達報告で完了済み。今回の予約台帳は別に作り、旧#43/#44を新規使用枠へ混在させない。

## 委任の運用範囲

- 再来週は11/23〜11/29 JSTと解釈。11/29までに、成立し必要な比較を順に構成・予約・代理QC/replayレビュー・分析する。日付まで無条件にsimulationを埋める運用ではない。
- cs10接続・固定commit enqueue・今回必要なdispatcher起動/job開始・short PASS後の1秒移行は許可済み。同じ範囲の再承認待ちを作らない。他予約のcancel/置換/pause/resume/優先度変更/停止はしない。
- 接続targetは`Ktakemori@cs10`を明示。VPNが必要ならmacOS標準の既存接続設定を使い、credentialを保存/変更しない。11/13はnative service表示はDisconnectedだったが、SSHと`id -un`は成功したため再接続/設定変更を行わなかった。表示だけでVPN障害と断定しない。
- 自分の実行上限を**最大7予約**とする：現J対照short+1秒の2、G1比較short+1秒の2、根拠ある同仮説内修理short+1秒の2、成立候補＋対照3.5秒の1。失敗も算入。1予約最大4条件、3.5秒は最大2条件。公称最大4,350,000step＋既存ceil差。不要な枠は使わない。
- G2の有限hook/離散化は具体的設計と静的検算まで。抗力/材料/PCA drive軸に及ぶ新しい構造をG1へ黙って混ぜない。科学的採択、merge、canonical freeze、T/dt/粘性sweep、3.5秒超の予約はしない。
- 新job追跡automationを増やさず、継続Goal内でfollow-upする。既存のread-onlyメール監視設定は変更しない。Mac/Codex停止や接続不可を越えた無条件の常駐は保証しない。

## 最初に進める現Jの2条件

G1は全3本の静的成立性がまだ確認できていない。その待ち時間にも、既存の正しいJで菌体接触の寄与を切り分ける。先にGを作れたと見なして4条件を予約しない。

| 条件ID | 幾何/drive | BF | FF |
| --- | --- | --- | --- |
| `nf03__slots024__J_bfon_control` | 現J / D / Kなし | ON | ON |
| `nf03__slots024__J_bfoff_control` | 上と座標・自然参照・全設定一致 | OFF | ON |

共通：nf3/中央slots024、11bead/本、hook自然長.25µm、局所phase0、T/refT2.5e−20Nm/本、dt4µs、η1mPa s、Brownian/切替OFF、解析torsion、遠位材料/ばね/既存指数反発固定。BF以外は一致し、OFFは診断で採択候補ではない。新しい力は加えない。

| Stage | campaign | parallel-job | 各条件 |
| --- | --- | --- | --- |
| short | `conf/phase2_multi_run/2010_hex_nf3_connector_contact_1tau_issue263.yaml` | `conf/phase2_parallel/issue263_initial_bundle/nf3_connector_contact_1tau_job.yaml` | .04秒・10000step |
| contact mechanism | `conf/phase2_multi_run/2010_hex_nf3_connector_contact_1s_issue263.yaml` | `conf/phase2_parallel/issue263_initial_bundle/nf3_connector_contact_1s_job.yaml` | 1秒・runtime250001step（公称250000との差1step） |

両stageとも`cs10_qualified`/max_workers:auto、実効2workers・数値thread各1、条件別output。Macではshortと1秒のsimulationなしdry-run/全条件初期geometry検査を確認する。Macのsimulationは32step/条件の回帰checkのみ（各.000128秒）。旧実行e32の7保存状態×BF ON/OFFと32step終端、計30arrayについて、診断ONの新実装の物理力・座標がbit単位で一致した。benchはON約3.95秒/32step、OFF約1.02秒/32step（pytest同時負荷下）で、Macの長時間計算はしない。

## 診断を物理式から分離する

`output.segment_contact_diagnostics: true`で実線分のBF/interFF/selfFFを同じproduction force passから記録する。[共通summary contract](phase2_run_summary_contract.md#opt-in実線分接触診断)に未適用/時刻/単位/first contact/cutoff内時間を定義した。既定false。ビーズ距離proxyを実反発と扱わず、QC閾値を変えない。

最初のshort完走では、新列のavailable、記録step数、pre-step時刻、候補pair数を照合する。明白な描画崩壊・貫通/実装異常がなくcommon strict QC PASSなら、同じ固定source/input契約で1秒を予約する。1τの周期不足だけで移行を止めない。

## fixed commitでの予約と完走レビュー

1. テスト・source回帰・共通preflight/preview・Mac planを確認し、物理source/configをcommit/push。具体的commit、job/campaign SHA、2ID、step数をlocal予約manifestとIssue/PR要約へ残す。
2. SSH identity、queueと既存dispatcher、固定commitを確認。隔離sourceでcs10のshort/1秒dry-runを再確認する。
3. 既存FIFOを使い`queue.py enqueue --branch <fixed SHA> --config <short job>`。不要なdispatcher重複起動をしない。必要な今回dispatcher/job開始は許可済み。
4. 予約ID/固定commit/正常開始/入力SHA/実測ETAを照合。応答が不明なら同じ予約を確認して重複enqueueを避ける。branch head移動を予約の置換と解釈しない。
5. 完走後はportableのみlocal JST outputsへ同期し、集合・件数・SHA全件一致、`run_summary.json`先読、共通development_evaluation strict QC/heatmap/固定3D/2D replay、frame/fps/coverage/SHA/視認を確認。QC-PASS archiveだけで機序を分析。
6. short PASSなら1秒を予約。全条件の実segment first contact/力/占有と、近位角、全体tilt、遠位幅/軸角、材料spin、菌体材料点自転、符号付き前進を`.04–.20/.20–.65/.65–1秒`で同時刻比較する。非保存内部stepの診断と保存形状の区別を保つ。

旧#44の実測ON4条件は約15〜16h/1秒だった。今回OFFはpair数が少なく費用が異なるため、これを終了保証とせず短期/自条件の実測で更新する。正常な中間ポーリングを固定周期で繰り返さず、次の完了見込みから確認時刻を選ぶ。通知は開始・評価完了・重要変更・介入に限る。

## G1と後続の分岐

G1は[静的設計条件](phase2_263_connector_repair_plan.md#具体案g後方をかわす経路と回転できる接続を同時に設計する)に従い、既存探索未成立を認めたbounded監査を進める。初期形状だけでなく固定自然形状で連続spin、hook長/非重複、受動応答、端力moment/drive仕事を確認する。G1成立時だけJ/G1×BF ON/OFFの4条件short→1秒を構成する。既存Jデータを再利用する場合はsourceと物理契約の同等性を検証する。

1秒の最後.65–1秒で既定の束幅/角度3条件の同時達成90%以上、各鎖正味材料spinとblock前進正、全期間strict QC PASS・replay正常なら、候補＋対照各3.5秒を予約できる。未達を時間延長だけで救わない。初期近束の維持と広い配置からの自発形成、定常性を区別する。

G1不成立なら、経路長/自由度/材質/drive軸のどこが不足かを示し、G2の具体設計と次条件を提出する。M/K係数sweep、QC緩和、body固定、初期位置復元力で救済しない。追加位置復元code/scale削除はmerge前残件。ADR0029 Proposed。

## 証拠と台帳

準備正本は`outputs/2026-11-13/180742/issue263_connector_execution_preparation/`。旧M/K台帳は保持し、今回の`execution_ledger.json`へ新規予約/condition/step・接続/レビュー/方針変更を記録する。Issue既存要約は重要3〜5項目だけ、詳細は本書/PRとlocal review_resultへ置く。
