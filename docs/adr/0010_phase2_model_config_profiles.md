# ADR 0010: Phase 2 model config profile と provenance

- status: accepted
- date: 2026-08-02
- scope: Phase 2 / Issue #164

## Context

従来の `conf/sim_swim.yaml` は Watari & Larson (2010) を基礎とする一方、11 beads/flagellum、body対角brace、分散トルク、attach-frame補強などのproject固有変更を含む。論文条件とproject defaultを同じ曖昧なconfig名で扱うと、runの由来と実装可能範囲をmanifestから判別できない。

また、2010年paper profileは15 body beadsと15 beads/flagellum、2015年refined profileは30 body beadsと30 beads/flagellumを使う。`11 beads` は現行projectの1本あたりのべん毛bead数、`120 beads` は3本条件での全体bead数であり、同じ `resolution` 名に混在させない。

## Decision

simulation configを次の4 profileへ分け、曖昧な旧入口 `conf/sim_swim.yaml` は残さない。

| Config | year | variant | resolution | implementation_status | nominal composition |
| --- | ---: | --- | --- | --- | --- |
| `conf/sim_swim_2010.yaml` | 2010 | `project` | `legacy_project` | `supported` | body 15 + flagellum 11 x 3 = 48 |
| `conf/sim_swim_2010_paper.yaml` | 2010 | `paper` | `coarse` | `supported` | body 15 + flagellum 15 x 3 = 60 |
| `conf/sim_swim_2015.yaml` | 2015 | `project` | `refined` | `pending` | body 30 + flagellum 30 x 3 = 120 |
| `conf/sim_swim_2015_paper.yaml` | 2015 | `paper` | `refined` | `pending` | body 30 + flagellum 30 x 3 = 120 |

`model_profile` は上表の `year`、`variant`、`resolution`、`implementation_status` に加え、`body_beads`、`flagellum_beads_per_filament`、`nominal_flagella_count`、`nominal_total_beads` を持つtyped provenanceとする。bead数は正の整数とし、`nominal_total_beads = body_beads + flagellum_beads_per_filament * nominal_flagella_count` を満たす。`nominal_*` はprofileの代表構成を表し、実行時に許可された `n_flagella` overrideを行ったrunの実bead数そのものではない。

canonical profileの `model_profile` はsource configの識別情報であり、CLIやcampaignからoverrideしない。metadataを持たない既存fixtureや過去configはparser互換のため読み込みを許可するが、4つのcanonical profileでは必須とする。

single-run manifestとgeneric multi-runのcampaign-level `run_manifest.json` には、正規化したtop-level `model_profile` と、読み込んだconfigを示すtop-level `source_config_path` を保存する。既存の `input.config` や `base_config` は互換性のため維持する。

CLI、sweep、multi-run、replay fallbackのdefaultは `conf/sim_swim_2010.yaml` とする。このprofileは旧 `conf/sim_swim.yaml` の有効設定値を、`model_profile` 追加と名称変更を除いて維持する。

## Paper と project の区別

`sim_swim_2010.yaml` は既存runとの互換を優先するproject profileである。`sim_swim_2010_paper.yaml` は15 beads/flagellum、bond長 `0.58 b`、14 bondsに対応するcontour length `8.12 b` と、ADR 0009で採用した `fene_fraenkel` を用い、project固有のattach-frame補強を無効相当へ戻す。

ただし `implementation_status: supported` は現行実装でsimulationを開始できることを示し、論文の完全再現を保証しない。現行builderが常時生成するbody diagonal braceや、paperと完全一致しないtorque/counter-torque実装などの制限はconfigコメントと本ADRで追跡する。

2015年の2 profileは #154で整理したrefined model値を機械可読な候補として保持する。`sim_swim_2015_paper.yaml` はpaper値、`sim_swim_2015.yaml` は今後のproject採用候補であり、Issue #164時点では両者に同じ初期候補値を置いても採用判断とはみなさない。

## Pending profile と後続Issueの境界

`implementation_status: pending` のprofileはYAML parse、profile表示、manifestを作らないdry-runまで許可するが、実simulationはoutput directory作成前に明示的に拒否する。

- Issue #165: `tau` / `s` duration schema、参照トルクからの `tau_s` 導出、正式な実時間換算を扱う。Issue #164では現行time schemaを変更しない。
- Issue #166: 正六角形断面 x 5層の30-bead body、30-bead flagellum、中央層attachment、brace切替を実装する。Issue #164ではgeometry builderを変更しない。
- Issue #167: 2015 paperのpotential、motor、counter-torque、局所補強OFFなどを実装・検証し、paper/project差分を確定する。Issue #164では2015 profileを実行可能とはしない。

この境界により、2015 configに値が存在することと、その値を現行solverが論文どおり解釈・実行できることを区別する。

## Consequences

- config名とmanifestだけで、year、paper/project、resolution、実装状態、source pathを判別できる。
- 2010 project defaultは名称変更後も既存挙動を維持する。
- 2010 paper profileを現行solverの実装制限付き参照条件として実行できる。
- 2015 profileを後続実装前に誤って本計算へ使うことを防ぐ。
- 過去manifestは新fieldなしでも読み込めるが、profile provenanceは復元できない場合がある。

## Verification status

4 profileのparse・metadata整合、2010 projectの旧config同値性、default path、manifest provenance、pending実行拒否、metadata override拒否を自動testで確認する。長時間simulationと動画目視はIssue #164の完了条件に含めない。
