# Issue #262：segment反発全OFF診断の実行契約

## 固定モデルと解釈範囲

#255の`2010_hex_project`局所one-ring全ベクトル反作用を13 attachment配置の固定対照とする。対照は`conf/phase2_multi_run/2010_hex_project_local_reaction_2s_issue255.yaml`で実行済みで、body–flagella反発OFF、べん毛同士・同一べん毛内の反発ON。#262はこの設定の`potentials.spring_spring_repulsion.flagella_flagella_enabled`だけをfalseにし、hookを含むflagella segment候補のactive数もゼロにする。body–bodyは従来の拘束ベース挙動を維持する。

両armとも旧有限差分torsionとlegacy frameを使用する。解析torsionの#255四方式C、#263の開発候補、#245・#257のarchiveとは物理条件が異なる。今回の結果はこの固定モデルの診断であり、canonical modelやdatasetの採択根拠にしない。#255対照は13/13完走・共通strict QC PASSだが、束化・軸整列には懸念が残り、旧torsionの有限差分誤差も解釈上の限界である。

## 条件とstage gate

`conf/phase2_multi_run/2010_hex_project_repulsion_off_1tau_issue262.yaml`と`2010_hex_project_repulsion_off_2s_issue262.yaml`を使う。1～6本のC6回転のみを同一視した13配置、`phase_seed=attach_seed=0`、各べん毛と基準トルク`2.5e-20 N m`、`dt_star=1e-4`、Brownian/switching OFF、初期hook中立・菌体軸直交、archive間隔1 msは#255と同一。screenは1τ・10,000 step、長時間stageは50τ・2秒・500,000 step。旧#255のON対照は新たにsimulationしない。

parallel-jobは`conf/phase2_parallel/issue262_repulsion_off/hex_1tau_job.yaml`と`hex_2s_job.yaml`。13条件を別outputに保存し、全条件geometry preflight、`max_workers: auto`、`worker_policy: cs10_qualified`を使う。dry-runで実効8 workers、各数値ライブラリthread=1、13の一意な出力を確認する。1τの共通`development_evaluation`によるfinite/body/hook length/flag/motor strict QC、初期geometry、固定camera 3D/2D replayをユーザーがレビューするまで2秒は開始しない。strict FAIL条件は原則進めず、診断専用の例外には条件と理由を明記する。hook角は診断値であり、新しい合否閾値は追加しない。

実行先は`cs10_user_run`／`execution:cs10`。Macでheavy simulationを行わない。#255局所方式の実測参考値は13条件の1τが42分41秒、2秒が約41時間。OFFの所要時間はscreen実測から再見積りする。cs10接続、enqueue、dispatcher/tmux、job開始・停止、予約の変更は操作ごとにユーザーの明示承認を得る。予約後のcommit移動を予約変更の許可とみなさない。

2026-10-17の#262 OFF screenは25分28秒で完了した。結果と2秒の暫定所要時間は`phase2_262_repulsion_off_1tau_results.md`を参照する。既存dispatcherが稼働中ならenqueue直後に自動開始するため、次のreservationではenqueue前に稼働状態を確認し、自動開始の承認も先に得る。

## archiveの検証と診断

対照は`outputs/2026-10-08/171556/issue255_hex_local_2s/`のcompleted archiveを使う。13 condition ID、profile、config、実行commit `e94c7cd4f0e31f99169b28ff5025e58e3b02b7cc`、各500,000 step、`outputs/2026-10-08/171556/sync_verification.json`の122ファイルSHA-256を照合する。対照と新規OFFは各々の元campaign configで共通評価し、設定一致を検証するまで評価器入力を結合しない。

completed archiveの菌体長軸3D単位方向、初期方向からの角度、0.1秒ごとの方向RMS・最大偏差・窓端間ドリフト、`body_axis_wobble_rms_deg`、固定camera replayを配置別に示す。`sim_swim.analysis.body_axis_diagnostic`はarchiveを読むだけで再simulationや新たなQC gateは行わない。同時刻ON/OFF窓比較は13配置と反発flag以外の設定が一致した場合のみ作る。各配置の「定常／揺れあり／判定不能」は時間帯、定量時系列、replay、QCを根拠に判定し、全体値だけで一時的な乱れを隠さない。未完走・archive欠損・判定に影響するQC FAILは定常としない。OFF単独結果から原因を断定しない。

cs10での各stage後は`run_summary.json`を先に確認し、ユーザー確認用artifact・manifest・summary・必要なreanalysis archiveをローカル同期して件数、SHA-256、QCを照合する。operational logとcredentialは同期しない。

## 実装時のローカル確認（2026-10-17）

旧#255対照の`run_manifest.json`・`job_manifest.json`・`run_summary.json`を確認し、13/13 completed、各500,000 step、固定commit`e94c7cd4f0e31f99169b28ff5025e58e3b02b7cc`、同期記録の122/122ファイルSHA-256一致を再検証した。`outputs/2026-10-17/195156/issue262_reference_axis_audit/`に対照のみの軸方向時系列・各配置20個の0.1秒窓・図・manifest・run.logを保存した。これは新規OFF結果やON/OFF差の判断ではない。
