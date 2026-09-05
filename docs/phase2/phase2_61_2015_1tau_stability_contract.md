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

各conditionは `tau_s`、`dt_internal_s`、`total_steps`、material coefficients、wall time、steps/sをcampaign manifestへ保存する。Stage A evidenceはコピーしない。実行時にsource manifestのSHA-256を計算し、`reference_evidence`へ保存する。Task Dは実行済みoutputが存在する場合だけ同様に参照し、未実行のTask Dをevidenceとして偽装しない。

## 判定

locked Stage A thresholdによるbody/non-body shape、hook/bond、bend/torsion、helix pitch/radius、motor action-reaction、finite/completionを1τの全記録stepで確認する。trajectoryとstate archiveがbody/flagella motionの記録として存在することも集計する。1条件でもFAILなら最初のcriterionを保存し、2015 profileの昇格とIssue #184へのhandoffを禁止する。3条件すべてPASSなら、後続評価に渡せるだけである。

fixed-reference・同一実時間の効率比較は本契約の対象外であり、tracking条件間のwall timeを同一物理系の直接比較として解釈しない。
