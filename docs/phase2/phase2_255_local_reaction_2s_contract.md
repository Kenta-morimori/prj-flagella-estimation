# Issue #255：局所・全ベクトル反作用2s実行契約

2010 hexの13 attachment配置と2010 projectのn=1..3、計16条件をcs10で2s実行する。1τで全16条件が共通short-screen QCを通過した結果は`phase2_255_local_reaction_1tau_results.md`に記録済み。project n=4..6、軸方向反作用、全body反作用の新規2s実行は対象外。canonical modelは採択しない。

## 条件と実行前FAIL

1τから変更する物理条件は実行時間だけ。`T=2.5e-20 N m/flagellum`、`dt_star=1e-4`、phase/attach seed 0、Brownian・switching OFF、body–flagella排除OFF、flagella同士の反発ON、局所one-ring全ベクトル反作用、compact archive 1 ms間隔を維持する。実行時間は`50τ`で指定する。reference `τ≈0.04 s`なので2.0 s相当、浮動小数点のceilによる余分な1 stepを避けて**500,000 step**となる。

両parallel jobは`geometry_all_conditions`で次を全条件に要求する。一件でもFAILなら予約・開始へ進まない。

- 第2ビーズ以降で推定した各らせん軸と菌体後方との角度が`1e-6°`以下。縮退軸はFAIL。
- hook曲げ角と平衡角90°の誤差、hookと菌体長軸の直交誤差が各`1e-6°`以下、hook長誤差が`1e-15 m`以下、hook初期曲げ力のnormが`1e-18 N`以下。ポテンシャルは有効のまま平衡位置から始める。
- hookが外向きで、非付着body bead・他flagellumとの中心距離と付着body beadからの最小距離がビーズ直径以上。

1τと同じ`development_evaluation`のfinite/body/hook length/flag/motorを全stepのPASS/FAILに使う。motor force residualの上限は`1e-8`、torque residualの上限は`0.02`。時間発展中のhook角とattach→first対菌体長軸角は**診断専用**とし、新しい90°閾値を設けない。局所solver失敗はcondition FAILとし、全bodyへのfallbackはしない。

|job|campaign|条件数|実効workers|
|---|---|---:|---:|
|`hex_local_2s_job.yaml`|`2010_hex_project_local_reaction_2s_issue255.yaml`|13|8|
|`project_local_2s_job.yaml`|`2010_project_local_reaction_2s_issue255.yaml`|3|3|

両jobは`conf/phase2_parallel/issue255_motor_reaction/`、campaignは`conf/phase2_multi_run/`にある。condition別outputを分離する。1τのcs10実測からの単純比例見積りはhex約35.6時間、project約8.1時間、計約43.6時間。長時間での実効速度は未確認。Mac wall timeは30分超、実行先はIssue #255の`cs10_user_run`・`execution:cs10`。

## 予約・実行・結果確認

targeted test、両job dry-run、local review PASS、commit/push後、cs10でremote commit・queue・NAS容量を確認する。確定commitからhex→projectの順にqueueへ予約し、許可済みdispatcherを起動して両jobを順次開始する。enqueue後のbranch更新を理由に予約を差し替えない。失敗・pause時は状況を報告し、cancel・置換・resumeは自動で行わない。

完了後はrun summaryを先に読み、必要artifactをローカルへ同期して件数・SHA-256を照合する。共通long-duration evaluatorで全step QC、first failure、motor、局所support/fallback、hook角時系列を確認する。固定cameraの3D/2D replayと保存座標の初期・終端軸方向2D配置図をレビューする。結果はIssue #255、PR #261、local reviewに記録する。特徴量評価・canonical採択は別判断とする。
