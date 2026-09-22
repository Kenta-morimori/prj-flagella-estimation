# Phase 2 Issue #61: 2015 project 1τ torque stability 契約

## 目的と境界

2015 project profileのtracking-reference条件を、3つの物理torqueで1τまで実行し、strict physical QCとthroughputを記録する。これは先行safety evidenceであり、10τ安定性、`dt_star`採択、dataset採択、canonical model freeze、supported profile昇格を行わない。

## 固定条件

| item | value |
| --- | --- |
| profile | `conf/sim_swim_2015.yaml` のprojectのみ |
| motion / Brownian | RUN固定、switching/reversalなし、OFF |
| geometry / seeds | `n_flagella=3`、attach/phase seed `0` |
| integration | `dt_star=1e-5`、`duration_tau=1` |
| torque | per flagellum `1e-21`, `2.5e-20`, `1e-19 N m` |
| scale policy | tracking-reference: motor = reference = force torque、`reference_torque` |

各conditionは `tau_s`、`dt_internal_s`、`total_steps`、material coefficients、wall time、steps/sをcampaign manifestへ保存する。Stage A evidenceはコピーしない。実行時にsource manifestのSHA-256を計算し、`reference_evidence`へ保存する。Task Dは実行済みoutputが存在する場合だけ同様に参照し、未実行のTask Dをevidenceとして偽装しない。evidence fileがrun開始時に指定されなかった場合は`not_recorded_at_run_start`と記録し、sourceが存在しなかったとは推論しない。

## 判定

locked Stage A thresholdによるbody/non-body shape、hook/bond、bend/torsion、helix pitch/radius、motor action-reaction、finite/completionを1τの全記録stepで確認する。trajectoryとstate archiveがbody/flagella motionの記録として存在することも集計する。1条件でもFAILなら最初のcriterionを保存し、2015 profileの昇格とIssue #184へのhandoffを禁止する。3条件すべてPASSなら、後続評価に渡せるだけである。

複数のgate / thresholdが違反した場合、`first failure`はsummary列の順序ではなく、raw diagnosticから得られる最初のobserved crossingを時刻・step順に選ぶ。時刻・stepを観測できない違反は、観測済みの早期違反より優先しない。

## 実測済み1τ結果（2026-09-23訂正）

不変campaign `outputs/2026-09-06/220250/parallel/issue61-2015-1tau__3654804140d9/campaign` をstreaming再解析した結果、`1e-21`、`2.5e-20`、`1e-19 N m`はいずれもstrict FAILである。全3条件の最初の観測済み違反は`max_motor_torque_balance_residual_ratio`で、step 0（それぞれ`1e-5`、`4e-7`、`1e-7 s`）だった。helix pitchも後続stepで閾値を超えるが、first failureではない。訂正結果は`outputs/2026-09-23/023008/analysis/issue61_2015_1tau_corrected/`に保存する。locked thresholdは変更しない。

fixed-reference・同一実時間の効率比較は本契約の対象外であり、tracking条件間のwall timeを同一物理系の直接比較として解釈しない。

## 1.2e-18 N m supplemental evidence

既存3 shardを再実行せず、2015 project・`n_flagella=3`・attach/phase seed `0`・
`dt_star=1e-5`・1τ・tracking-referenceの単独conditionとして`1.2e-18 N m`を追加する。
これは論文対応torque値のproject-model診断であり、既存3条件のcampaign manifestへ混在させない。
同じlocked QCとthroughputを比較表に記録するが、profile昇格、#184 handoff、dataset採択の根拠にはしない。
