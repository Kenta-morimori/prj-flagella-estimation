# ADR 0011: 2010 project profileの初期べん毛軸を後方整列defaultにする

- status: accepted
- date: 2026-08-02
- scope: Phase 2 / Issue #171

## Context

`conf/sim_swim_2010.yaml` は従来、べん毛付着点の局所法線に基づくside-attach初期化を互換defaultとしていた。一方、RUN状態のproject profileでは3本の螺旋中心軸を菌体後方へ揃えた初期条件が、定性評価と本数別比較の基準として明確である。

paper profileは論文参照性を保つ必要があり、pendingの2015 project profileはrefined geometryとdynamicsが未実装である。そのため、全4 profileを一律に変更するのではなく、実行・評価済みの2010 projectだけを対象とする。

## Decision

`conf/sim_swim_2010.yaml` の `flagella.initial_helix_axis_from_rear_deg` を `0` とし、第2ビーズ以降の螺旋軸を菌体後方へ揃える初期条件をproject defaultにする。

`conf/sim_swim_2010_paper.yaml`、`conf/sim_swim_2015_paper.yaml`、`conf/sim_swim_2015.yaml` は `null` を維持する。`null` は今後もside-attach参照条件を明示するruntime overrideとして利用できる。

## Evidence

user実行run `outputs/phase2_project_posterior_default/2010_project/2026-08-02/181138` を評価した。条件は `attach_seed=0`、`phase_seed=0`、`duration_s=0.5`、`dt_star=1.0e-4`、後方軸 `0 deg` である。

- 5,000 / 5,000 stepで `finite_pass`、`shape_pass_nonbody`、`shape_pass_nonbody_strict` がtrue
- max hook relative error: `0.08801`
- max flag bond relative error: `0.08626`
- max bend / torsion error: `15.61 deg` / `32.41 deg`
- 個別helix-axis rearward projection min: `0.78186`
- bundle rearward projection min: `0.96679`
- final bundle axis vs rear: `14.62 deg`
- axis fit R2 min: `0.94857`、degenerate count: `0`
- axis-center/body-roll ratio mean: `2.354`
- user qualitative review: 問題なし

個別軸は全期間最大 `38.57 deg` まで開くが、bundle軸は後方を維持する。standard 0.5 s runはbody diagnostics CSVを出力しないため、body変形の判定はuser visual reviewを根拠とする。

## Boundaries

- Issue #165の `tau` / `s` duration schemaと実時間換算は変更しない。
- Issue #166の2015 refined geometry、#167の2015 paper dynamics、#168の安定性評価には踏み込まない。
- 2015 projectの後方default化は、pending解除後の別の実行・採用判断とする。

## Consequences

- defaultの2010 project runは、明示overrideなしで後方整列初期条件を使う。
- 従来のside-attach参照は `flagella.initial_helix_axis_from_rear_deg=null` で再現できる。
- paper profileの初期化とpendingの2015 project候補の意味は変更しない。
