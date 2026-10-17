# Issue #255：nf3の束形成・軸整列・持続回転探索

> 現在地（2026-10-17）：本書はqueue #24/#25の実行契約の履歴。1τ・1sは完了し、追加2sは見送り。束維持・定常回転・推進のモデル開発は#263へ引き継ぐ。本書は新規実行の承認ではない。

## 目的と境界

現在のnf03_slots024の後方平行・hook中立初期配置から、べん毛が接近し、束軸を菌体軸へ揃えて持続回転する条件、またはそれを妨げる要因を特定する。反作用の釣り合いと望ましい挙動は分けて評価する。目標状態になることは事前に保証しない。

排除体積ON/OFFは#262へ分離する。今回、菌体–べん毛排除OFF、べん毛間・自己接触の既存反発ONを維持する。nf4/5、hook曲げ剛性・長さ、初期位置・位相、T、dt、RPY、Brownian/switchingの変更は含めない。canonical採択・dataset変更・merge・2s延長は自動で行わない。

## torsionの事前検証と実装

#23のrun_summary先読後、nf3/4/5各5状態（0/.2/.6/.82/1s）の同じエネルギーを再評価した。有限差分幅1e-9/1e-10/1e-11mで最大総トルクは1.081e-20/1.076e-22/1.076e-24 N m、最大総力は1.126e-14/1.123e-16/1.123e-18 N。刻み依存が二次収束している。解析勾配では総トルク最大3.685e-34 N m、総力最大2.779e-28 N。

`potentials.torsion.force_method: analytic`を追加する。既定`finite_difference`と既存profileの結果を維持し、同じwrapped quadraticポテンシャルを正規化・射影・atan2まで解析微分する。縮退した二面角とwrapped-energyの微分不能点では明示エラーとし、射影による残差相殺やfallbackはしない。今回の4条件で共通にanalyticを選ぶ。#20/#23とはtorsionが異なるため、過去結果は参考であり新規4条件間を主な因果比較にする。

証拠：`outputs/2026-10-15/160836/issue255_nf3_bundle_preparation/torsion_archive_audit.json`。energy gradient、FD収束、並進・回転共変性、総力・総トルク、縮退・cusp、既定不変をtargeted testで確認する。

## 4条件

| ID末尾 | frame位置復元scale | motor横方向補正 |
|---|---:|---|
| A | 1.25（現行） | none |
| B | 1.25 | minimum_norm |
| C | 1.0（追加分なし） | none |
| D | 1.0 | minimum_norm |

全件`motor.attach_frame_reaction: energy_gradient`。scale=1.0は追加位置復元力だけを外し、hook長・接続・基本曲げ抵抗は維持する。菌体側は実際のべん毛motor全トルクへのattach_one_ring反作用、fallbackなし。横方向補正OFFは保存則を外す意味ではなく、べん毛側の駆動分布の選択である。

全件nf=3、slots024、seed=0、T=2.5e-20 N m/本、dt_star=1e-4（4µs）、reference τ=.04s。初期後方平行・hook中立・菌体軸直交、runtime90°拘束なし。phase_seed=0、Brownian/switching OFF。4条件の初期位置が完全一致し、2つの比較設定以外が一致することを検証する。

campaignは`2010_hex_project_nf3_bundle_1tau_issue255.yaml`と`2010_hex_project_nf3_bundle_1s_issue255.yaml`。parallel-jobは`conf/phase2_parallel/issue255_motor_reaction/nf3_bundle_1tau_job.yaml`と`nf3_bundle_1s_job.yaml`。condition IDは`nf03__slots024__A/B/C/D`。cs10_qualified、auto=独立4条件の実効4 workers、各thread1、出力分離、全条件geometry preflight。

## 回転計測と挙動診断

`output.rotation_diagnostics: true`（既定false）で`rotation_diagnostics.csv`と定義manifestを追加する。各内部stepで菌体の並進・最良剛体回転を除き、ラベル付きべん毛鎖（hook第1ビーズ除外）の最良剛体回転を前時刻の向き付きPCA軸へ射影して積算する。これは変形体のbest-fit axial spinであり、らせん軸の首振りや厳密な弾性rodのmaterial twistとは区別する。

積算位相はmoduloを取らず、1ms出力の間に何回転しても積算値を保持する。内部stepで0.25rad以上の回転・直線退化は判別不能として以降の位相を無効化し、リセットして見かけの同期を作らない。非剛体fit残差の区間最大も保存する。既知の剛体回転・並進・菌体共回転・出力間20rad回転を検証する。motor/全力のF(prestep)・実変位から全ビーズの累積仕事も保存する（lab frame、Brownian OFFの診断）。

共通development_evaluationでfinite/body/hook length/flag/motorを先に判定し、QC通過completed archiveを再simulationなしで別の再利用可能な`sim_swim.analysis.bundle_rotation`へ渡す。Issue専用QC evaluatorは追加しない。

診断は束の広がり、べん毛間軸角、正規化輪郭長25/50/75/100%の横断的広がり、束軸対後方方向、遊泳方向対菌体軸、速度変動、各べん毛の菌体相対回転数・角速度、累積位相差の変動を扱う。位相差は初期物理位相ではなく積算回転の差であり、bounded driftから回転協調を検討する。初期位相が同じであることだけを同期の証拠にしない。

初期緩和を含む全時間と0～.2/.2～.5/.5～.75/.75～1sの窓を示す。束の平均軸が揃っても個々の軸が開いていれば束形成とはしない。少なくとも終盤の複数周期で形状・軸角・回転速度・位相差が維持されるかを判断する。周期が不足・fitが不確実なら定常性は判定不能とする。新しい挙動FAIL閾値を機械的に採用せず、比較時系列とreplayを根拠にユーザーレビューする。

## 段階実行と見積り

1τ=10,000 stepを4条件で実行し、共通QC・初期形状・局所solver・力の極値・固定camera3D/2D replayをレビューする。束化未達だけでscreenから除外しない。QC未達や判別不能な回転計測が残れば、原因を報告し自動で長時間stageへ進めない。

レビュー通過後に同じ4方式の1s=250,000 stepを行い、根拠付き比較を行う。新しいtorsion実装は計算費用が変わるため、#22 nf3の1τ約12分/#23 nf3の1s約5.5時間は参考値であり、所要時間と完了予定時刻は新screen実測で更新する。Macでheavy runは行わず、Mac wall timeは30分超として扱う。長時間共有機監視はユーザー指定がなければ追加しない。

固定commitと操作範囲を示してcs10接続・enqueue・dispatcher/job開始の承認を確認する。既存予約を変更しない。remote成果物をローカル同期し、件数・SHA-256、QC、video frame数を照合する。Issue #255・PR #261・local reviewへ証拠を記録する。目標を満たさない場合も成功したとせず、どの比較がどの仮説を支持/反証したかと未解決点を報告する。
