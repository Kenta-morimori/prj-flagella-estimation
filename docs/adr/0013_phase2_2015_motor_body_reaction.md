# ADR 0013: 2015 motor dynamicsを局所body反作用付き近似として導入する

- Status: Accepted
- Date: 2026-08-02
- Issue: #167

## Context

Kong et al. (2015)のrefined profileでは、flagellar motorの駆動torqueとbodyへの
反作用を対として扱う必要がある。既存の`triplet`はattach/first/secondの局所
force coupleで駆動torqueを作るが、独立したbody反作用を持たない。
`root_torque_segment_couples`と`root_torque_axis_projection`は全body beadsへ反作用を
分配するproject実装であり、hook近傍の局所反作用とは区別が必要である。

#166の120-bead refined geometryは未実装であるため、#167ではgeometryに依存しない
motor dynamicsとprovenanceだけを導入する。

## Decision

`motor.force_distribution=hook_coupled_body_reaction`を追加する。
`sim_swim_2015_paper.yaml`だけがこのmodelを使い、project比較profileの
`sim_swim_2015.yaml`は既存`root_torque_segment_couples`を維持する。

各flagellumについて、first-second方向をmotor軸とし、flagellum基部の最初の3 beadsへ
合力ゼロかつtorque全ベクトルが`motor.torque_Nm * motor_axis`となる最小ノルム力を
分配する。body側にはflagellum側で実際に生じたtorque全ベクトルの逆符号を与える。
これにより、軸方向成分だけでなく軸直交成分もswimmer全体で相殺する。body supportは
attach beadと、`body_ring_edges`または
`body_vertical_edges`で直接接続されたone-ring beadsとする。このsupportで軸torqueを
生成できない場合だけ、全body beadsを使うzero-net-force torque coupleへfallbackする。

このモデルは論文中の離散force配置を完全再現するものではないため、provenanceを
`paper_inspired_approximation`とする。manifestには次を保存する。

- `dynamics.implementation_status: implemented`
- `dynamics.force_distribution`
- `dynamics.provenance`
- body reaction modelとfallback model
- run中に観測したreaction support bead数
- run中に一度でもfallbackしたか
- `geometry.implementation_status`
- `simulation.implementation_status`とblocker
- `paper_reference.parameters`の項目別source classification

2015 project / paper profilesは、#166と#168が完了するまで
`model_profile.implementation_status: pending`、`geometry: pending`、
`simulation: blocked`を維持する。

## Consequences

flagellum側とbody側は個別に合力ゼロとなり、motor軸まわりtorqueが相殺される。局所
support縮退時も全body fallbackで反作用を維持でき、その使用はmanifestから監査できる。

一方、120-bead geometry上の安定性、長時間の形状保持、論文再現性は本ADRでは主張しない。
motor-off `0.1 tau`、motor-on `1 tau`とそれ以降の評価は#168で扱う。
parameter対応表とsupported昇格条件は
`docs/phase2/phase2_167_2015_paper_conditions.md`を正本とする。
