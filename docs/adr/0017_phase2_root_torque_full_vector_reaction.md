# ADR 0017: root torque segment couplesの全ベクトル反作用

- Status: Superseded in part by ADR 0019
- Date: 2026-08-13
- Issue: #61

## Context

2010 projectの`root_torque_segment_couples`は、べん毛の各segmentへ局所force coupleを分配し、菌体へ逆向きの反作用を分配するproject実装である。従来は、べん毛側で実際に発生した回転力のうち、主軸方向の成分だけを菌体側で相殺していた。

Issue #61のtorque-linked `1 tau` screenでは、合力は機械精度で釣り合った一方、主軸に直交する回転力が残り、motor torque balance residualが約3--5%になった。これは`dt_star`比較の前提であるaction-reaction gate（`<=1e-8`）を満たさない。

## Decision

Issue #61のexperiment conditionでは、`motor.body_reaction_full_vector=true`を明示する。これにより`root_torque_segment_couples`は、各flagellumへ実際に加えた回転力の**全ベクトル**を計算し、全body beadsへ最小ノルムのzero-net-force force coupleとして逆符号で分配する。

これにより、軸方向だけでなく横方向の回転力もswimmer全体で相殺する。flagellum側のsegment weight、目標torque、時間scale、material coefficientは変更しない。このopt-inを持たない既存2010 baseline、および`triplet`・`root_torque_axis_projection`・2015 project profileの既存挙動は変更しない。

既存runの数値結果は再解釈しない。反作用修正前に作成した#61 initial screenは比較・採択に使わず、新しいcommitで再実行する。

## Consequences

- #61のmotor action-reaction gateを、閾値緩和ではなく力学実装の整合で満たせる。
- 固定τcontrol、dataset採択、2015 profileの採否は変更しない。新規runのτ policyはADR 0019に従う。
- run manifestには`root_torque_segment_couples`の全body・全ベクトル反作用を記録する。
- full-vector reactionは力の分配を変えるため、既存baseline screenの数値値と直接同一視しない。#61の1 tau再実行後に形状・遊泳・dt一致度を評価する。
