# Issue #255：局所反作用1τ実行契約

Issueのexecution targetは`cs10_user_run`、labelは`execution:cs10`。
対象はhex 13配置とproject n=1..3の3形状、計16条件。
既存の全body・全ベクトル16条件（固定commit `6592e14`）と軸方向16条件を対照とする。
べん毛側driveと初期形状を対内で一致させ、菌体側supportだけを変更する。

|job|campaign|条件|実効workers|
|---|---|---:|---:|
|hex|`conf/phase2_multi_run/2010_hex_project_local_reaction_1tau_issue255.yaml`|13|8|
|project|`conf/phase2_multi_run/2010_project_local_reaction_1tau_issue255.yaml`|3|3|

job YAMLは`conf/phase2_parallel/issue255_motor_reaction/`の
`hex_local_1tau_job.yaml`と`project_local_1tau_job.yaml`。
両方`cs10_qualified`、数値ライブラリ各1 thread、geometry全条件preflight、
condition別出力を要求する。Macでのwall time見積りは30分超。
旧1τのcs10実測はhex26条件63.06分、project6条件9.97分。
新jobの単純見積りはhex約32分、project約10分、合計約42分だが、
局所solverの費用と新軌道による変動は未測定である。

共通条件は`T=2.5e-20 N m/flagellum`、`dt_star=1e-4`、
1τ=10,000 step、phase seed 0、Brownian/switching OFF、
body–flagella排除OFF、flagella同士の反発ON、1 ms間隔のcompact archive。
共通`development_evaluation`のfinite/body/hook length/flag/motorを個別に判定し、
motor force residual上限`1e-8`、torque residual上限`0.02`を維持する。
hook対菌体長軸角度は診断値であり、新しい閾値を設けない。

固定commitとlocal review PASS後にcs10のqueueへ予約する。
接続、各enqueue、dispatcher/job開始には操作ごとのユーザー明示承認を確認する。
完了後はrun summaryを先に読み、必要成果物をローカル同期して
件数とSHA-256を照合し、QC・角度・固定camera replayをレビューする。
2sは開始しない。

実行結果は`docs/phase2/phase2_255_local_reaction_1tau_results.md`に記録した。
実測wall timeはhex 42分41秒、project 9分41秒、合計約52分22秒。
