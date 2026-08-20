# ADR 0018: 2010 paper profileのtorqueと積分刻みを原著値へ訂正する

## Context

`conf/sim_swim_2010_paper.yaml`はWatari & Larson (2010, DOI: 10.1016/j.bpj.2009.09.044)を参照するprofileである。しかし、motor torqueには過去互換の`-1` sentinel、積分刻みには`dt_star=1e-4`が残っており、原著の数値条件を明示していなかった。

## Decision

2010 paper profileの各flagellum motor torqueとreference torqueを`1.2e-18 N m`、積分刻みを`Δt/τ = dt_star = 1e-3`に設定する。`τ`は`η b^3/T`から導出する。

`-1` sentinelの読取り互換は`SimulationConfig`に残すが、2010 paper profileの標準configには使わない。新規runのτ policyはADR 0019に従い、過去の固定τcontrol、#61のcampaign、既存runは変更・再解釈しない。

## Consequences

- paper profileのmanifestは明示的なphysical torque、`τ`、`Δt`を記録する。
- 2010 paper profileの長時間runは、project profileの11-bead/分散torque/後方軸条件とは別の再現検証として扱う。
- #61の既存8条件は2010 projectの無次元相似性診断であり、paper profileの検証結果には用いない。
