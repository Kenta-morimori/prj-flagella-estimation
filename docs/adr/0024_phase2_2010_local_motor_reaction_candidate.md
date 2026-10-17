# ADR 0024: 2010系の局所・全ベクトルmotor反作用を評価候補として追加する

- Status: Proposed
- Date: 2026-10-05
- Issue: #255

## Context

Issue #255の1τ比較では、軸方向反作用のmotor QCは5/16条件のみPASSし、
全body beadへの全ベクトル反作用は16/16条件PASSした。両方式のべん毛側driveは同一である。
全bodyへの反作用分配がhook近傍の局所負荷と同じ運動を生むとは限らない。
2015 paper profileの`hook_coupled_body_reaction`は局所反作用を持つが、
べん毛側driveも基部3ビーズ方式へ変えるため、今回の対照には用いない。

## Decision

`root_torque_segment_couples`の菌体側に、opt-in設定
`motor.body_reaction_support: attach_one_ring`を追加する。
`motor.body_reaction_full_vector: true`との組合せでのみ有効とし、
各flagellumの付着body beadとring/vertical edgeで直接接続されたbody beadsへ、
べん毛側で実際に発生したtorque全ベクトルの逆符号を合力ゼロのforce coupleとして分配する。
べん毛側のdrive、segment重み、初期geometry、時間積分とQC閾値は変えない。
既定の`all_body`は従来通りとする。

局所supportが不足またはsolverが縮退した場合、conditionを失敗として記録する。
全bodyへの自動fallbackはしない。manifestへsupport種別、観測ビーズ数、
fallbackなしを記録する。これはpending候補の評価であり、canonical採択ではない。

## Evaluation

hexの13配置とproject n=1..3の3形状で1τ short screenを行う。
既存の軸方向および全body・全ベクトル1τ成果物を同一条件の対照に用いる。
共通`development_evaluation`のfinite/body/hook length/flag/motor判定、
hook対菌体長軸角度、固定camera replayを比較する。
2sの対象と採択は1τの結果をレビュー後に別途判断する。
