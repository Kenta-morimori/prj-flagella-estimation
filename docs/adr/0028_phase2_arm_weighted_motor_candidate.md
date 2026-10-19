# ADR 0028: 根元腕長に頑健なsegment torque分配候補

Status: Proposed

## Context

#262の旧局所反作用モデルでは、根元segmentの主軸垂直腕が約0.007 µmへ縮み、固定重みのsegment force coupleで大きな力・横方向トルクが生じた。反発ON/OFFの両方で菌体自転の鈍化と同時期だった。ただし#263の開始条件である四方式Cの1秒archiveでは、同じ16τ付近に腕の崩壊は見られない。今後のcompact配置で再発するかは未検証である。

## Decision

既定の`motor.segment_torque_distribution: fixed_weights`を維持し、`arm_weighted_minimum_norm`をopt-inの研究候補として追加する。正規化済みの既存segment重みを`w_i`、segmentのべん毛主軸への垂直腕を`r⊥,i`とし、軸トルクの配分を

```text
τ_i = τ w_i |r⊥,i|² / Σ_j w_j |r⊥,j|²
```

にする。従来の局所force coupleをこの`τ_i`で計算する。個別の腕がゼロへ近づいても、他の腕が有限なら力は発散しない。全腕が事実上ゼロの場合は失敗を明示する。局所one-ringの全ベクトル反作用は維持し、既存の`minimum_norm`横トルク補正との同時指定は認めない。

## Consequences

- 正味の軸トルク、局所force coupleのゼロ合力、菌体との全力・全トルク釣り合いをtargeted testで確認する。既定経路の数値列は変えない。
- segment間のtorque分配が姿勢に依存するため、力学的に別の候補である。横方向トルクは残り得る。遊泳改善や定常性をこの式だけから主張しない。
- 共通`development_evaluation` QC・固定camera replayを通した後、compact対照と同じ初期座標・トルク・刻みで、少なくとも16τを含む時間窓を比較する。候補の採択とcanonical変更は別判断とする。
