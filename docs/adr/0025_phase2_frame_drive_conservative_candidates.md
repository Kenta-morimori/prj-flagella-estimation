# ADR 0025：付着frameとmotor駆動の釣り合い修正を評価候補として追加する

- Status: Proposed
- Date: 2026-10-08
- Issue: #255

## 背景

局所全vector反作用の#20はmotor QCを通過したが、ユーザーの望む後方束化・菌体軸整列・滑らかな遊泳を満たしていない。motor反作用の保存則だけでは、継承したframe復元力とべん毛駆動の妥当性を保証できない。

## 評価候補

1. `motor.attach_frame_reaction: energy_gradient`でframe位置依存を含む全ポテンシャル勾配を用いる。既定は`legacy`。
2. `motor.segment_torque_correction: minimum_norm`で既存segment駆動への最小ノルム補正を加え、合力ゼロ・駆動軸方向Tの全トルクを満たす。既定は`none`。

両者を独立に選択できる。縮退frame/solverは失敗とし、代替モデルへfallbackしない。局所motor反作用・初期中立配置・QC閾値を維持する。時間発展中の90°拘束は加えない。

## 検証と限界

数値勾配、保存則、回転/並進整合性、縮退と既存方式への非影響を検証する。実行候補はユーザー指定のnf3/4/5各1配置・両修正の3条件で1τ→1sとする。計画変更により単独修正の動的寄与は識別できない。最小ノルム補正は横方向トルクを抑えても大きな局所力を必ず解消しない。

詳細な条件・実行境界は`docs/phase2/phase2_255_frame_drive_contract.md`。このADRは候補の実装記録であり、物理モデルの採択、canonical変更、cs10開始、2sへの展開を承認しない。
