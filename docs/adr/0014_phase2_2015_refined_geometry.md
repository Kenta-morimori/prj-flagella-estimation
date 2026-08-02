# ADR 0014: 2015 refined geometryを六角形5層の120-bead modelとして実装する

- Status: Accepted
- Date: 2026-08-02
- Issue: #166

## Context

Kong et al. (2015)はrefined modelをbody 30 beadsと30 beads/flagellum x 3本の
total 120 beadsとする一方、本文ではbody断面をpentagonと記載する。Issue #166では
Fig. 1、30 body beads、5 layersとの整合を優先し、正六角形x5層をproject inferenceとして
採用することが決定済みである。

## Decision

- bodyは軸方向5 layers、各layer 6 beads、layer spacing `0.5 b`、length `2.0 b`とする。
- body radiusはdefault `0.5 b`（width `1.0 b`）。paper比較のwidth `0.7 b`は
  `body.prism.radius_over_b=0.35` overrideで表す。
- ring、vertical、diagonal edgesを区別し、2015 profilesでは
  `body.prism.diagonal_braces_enabled=false`とする。2010互換defaultはtrueとする。
- 1本のflagellumは30 beadsと29本の`0.29 b` bondを持つ。body attachからflagellum bead 0
  までをHookとし、`hook.length_over_b=0.25`をimplementation assumptionとして使う。
- `placement_mode=seeded_center_layer`は中央層の120度間隔slotsを選び、seed 0/1は
  `[0,2,4]` / `[1,3,5]`となる。`attach_seed`と`phase_seed`は独立とする。
- manifestはnominal構成に加えてactual bead数、body寸法、edge数、Hook/bond長、
  attachment layer/slotを記録する。

## Evaluation state

2015 profilesは採否未確定のため`model_profile.implementation_status=pending`を維持する。
geometryとdynamicsが実装済みになったためsimulationは`evaluation_ready`とし、#168用の
CLI/campaign/replay実行を許可する。manifestの残blockerは#168のみとする。

## Consequences

paper本文のpentagonを完全再現したとは主張しない。六角形採用は
`figure_and_bead_count_inference`、Hook長は`implementation_assumption`、brace OFFは
`paper_condition`として区別する。motor-off `0.1 tau`、motor-on `1 tau`、長時間安定性、
brace採否、supported昇格は#168へ残す。
