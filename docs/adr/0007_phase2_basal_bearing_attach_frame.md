# ADR 0007: Phase 2 basal bearing attach-frame tangent mode

- status: accepted
- date: 2026-07-08
- scope: Phase 2 / Issue #103

## Context

Issue #94 / PR #95 で導入した attach-frame 補強は，後方条件の hook 過伸長を抑えた。一方で，`local_attach_frame_position_scale=3` と `local_attach_frame_tangent_scale=1.5` を併用した条件では，べん毛が螺旋軸中心に回るというより，body を含めた剛体回転に近く見える。

既存 `local_attach_frame_tangent_scale` は，body表面局所frameで見た `first -> second` 根元接線ベクトルを初期ベクトルへ戻す。これは hook/root 近傍の姿勢を安定化するが，root の軸まわり spin 自由度も同時に抑える可能性がある。

## Decision

`motor.local_attach_frame_tangent_mode` を追加する。

- `vector`: 既存挙動。`first -> second` ベクトル全体を body attach frame の target へ戻す。
- `basal_bearing`: `attach -> first` 軸まわりの方位角を拘束せず，`first -> second` の軸方向成分と垂直半径だけを保つ。

default は `vector` とし，既存条件の挙動を変えない。`basal_bearing` は Issue #103 の診断・候補比較用 extension として扱う。

## Consequences

- `basal_bearing` は参照論文モデルそのものではなく，現行 bead-spring 実装で hook安定化と root axial spin 自由度を分離するための project-specific numerical stabilizing approximation である。
- `local_attach_frame_position_scale` は attach-root 位置安定化として維持し，tangent mode だけを切り替えることで `fp` と `ft` の寄与を分離できる。
- `body_roll_phase_deg` と `axis_center_body_relative_phase_deg` を出力し，body roll と螺旋軸中心 spin を同じ summary で比較する。

## Verification

今回のPRでは長時間 sweep / render は実行しない。短時間 smoke と単体テストで以下を確認する。

- `basal_bearing` force は `attach -> first` 軸まわりの純回転に反応しない。
- 軸方向成分または垂直半径が target からずれた場合は復元力を出す。
- `basal_freedom_diagnostic.yaml` は #103 用の5条件を展開する。
- `summary.csv` と `flag_helix_axis_diagnostics.csv` に body roll / body-relative spin 指標が出る。
