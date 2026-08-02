# Phase 2 Issue #167: Kong et al. (2015) paper条件対応表

## Scope

Issue #167は2015 profileのparameter schema、motor dynamics、反作用、provenanceを
実装した。Issue #166で120-bead geometryも実装済みであり、2015 project / paper profileは
#168の評価実行が可能である。ただしsupported昇格と安定性採否は#168で扱う。

## Profile status

| profile | dynamics | geometry | simulation | motor model |
| --- | --- | --- | --- | --- |
| `sim_swim_2015.yaml` | implemented | implemented | evaluation ready | `root_torque_segment_couples` |
| `sim_swim_2015_paper.yaml` | implemented | implemented | evaluation ready | `hook_coupled_body_reaction` |

paper motorはKong et al. (2015)のmotor/body counter-torqueを参考にした
`paper_inspired_approximation`であり、論文の離散force配置の完全再現ではない。
project motorは既存のmaterial twist伝達比較用実装である。

## Parameter correspondence

| config key / condition | value | provenance | implementation status |
| --- | ---: | --- | --- |
| `scale.bead_radius_a_over_b` | `0.06` | paper Table 1 | implemented |
| body width | `1.0 b` | paper default | profile default |
| body width comparison | `0.7 b` | paper comparison | `body.prism.radius_over_b=0.35`でparse可能 |
| `body.prism.n_prism` | `6` | figure、30 body beads、5 layersからのproject推論 | implemented |
| `flagella.bond_L_over_b` | `0.29` | paper Table 1 | implemented |
| `potentials.spring.formulation` | `fene_fraenkel` | #163採用結果 | implemented |
| `potentials.spring.H_over_T_over_b` | `1000` | paper Table 1 | implemented |
| `potentials.spring.s` | `0.1` | FENE-Fraenkel implementation assumption | implemented |
| `potentials.bend.kb_over_T` | `30` | paper Table 1 | implemented |
| normal bend / torsion | `161 deg` / `-30 deg` | paper Table 1 | implemented |
| semicoiled bend / torsion | `138 deg` / `32.5 deg` | paper Table 1 | implemented |
| curly bend / torsion | `142 deg` / `60 deg` | paper Table 1 | implemented |
| `potentials.torsion.kt_over_T` | `5` | paper Table 1 | implemented |
| `motor.torque_Nm` | `1.2e-18 N m` | 1200 pN nm reference torque | implemented |
| `motor.reference_torque_Nm` | `1.2e-18 N m` | #165 time-scale contract | implemented |
| `time.integration.dt_star` | `1e-5` | paper integration step | implemented |
| `brownian.enabled` | `false` | paper condition | implemented |
| Hook length | `0.25 b` | implementation assumption | implemented |
| body diagonal brace | OFF | paper condition | implemented |
| center-layer attachment slots | `[0,2,4]` / `[1,3,5]` | figure-based geometry interpretation | implemented |
| project local stiffness extensions | OFF (`1.0`) | paper/project separation policy | implemented |

`manifest.json` / campaign manifestの`paper_reference.parameters`は、主要parameterを
`paper_table_1`、`paper_condition`、`figure_and_bead_count_inference`、
`implementation_assumption`、`paper_inspired_approximation`、
`project_comparison_model`に分類して保存する。

## Motor balance

`hook_coupled_body_reaction`はflagellum基部3 beadsへ合力ゼロかつ指定torque全ベクトルを
与える。body側はattach beadとring/vertical one-ringへその逆torque全ベクトルを与え、
局所supportが縮退した場合だけ全body beadsへfallbackする。refined hexagonal geometryでも
中央層attach bead、ring隣接2 beads、vertical隣接2 beadsの計5 beadsを使用し、fallbackなしで
swimmer全体のnet forceとnet torqueが数値許容範囲内であることを自動検証する。

## Promotion to supported

2015 profilesを`implementation_status: supported`へ変更するには、少なくとも次が必要である。

1. #166で30-body + 3x30-flagella geometry、brace OFF、Hook長、attachment slotsを実装する。
2. refined geometry上でparameter/provenanceとnominal 120-bead構成を検証する。
3. #168でmotor-off `0.1 tau`とmotor-on `1 tau`を実行し、有限性、body/flagellum形状、
   net force/torque、fallback、collapse/fly-awayを評価する。
4. project/paper profileごとに安定性結果と採否を記録する。

Issue #166完了後は次のcommandで#168評価を開始できる。ただし長時間条件は#168でのみ実行する。

```bash
uv run python -m scripts.01_simulate_swimming config=conf/sim_swim_2015_paper.yaml
```
