# Issue #204 Feature Study Reference Run

Issue #204の3D / 2D運動特徴量・時系列plot整理では，以下を固定診断参照runとする．

```text
outputs/phase2_multi_run/flagella_count_behavior_v1_r2/reference/
2010_project_tau_linked_2s_nf1_4_as3_ps3_2026-08-18/
```

## Provenance

- 実行開始: 2026-08-18 17:48:28 JST
- 入力config: `conf/phase2_multi_run/flagella_count_behavior_v1_r2_tau_linked_2s.yaml`
- model: 2010 project，tau-linked，`T=2.5e-20 N m`，`dt_star=1e-3`，2.0 s RUN固定
- 条件: `n_flagella=1,2,3,4` × attach seed `0,1,2` × phase seed `0,1,2`（36 independent runs）
- root artifacts: `manifest.json`，`run_manifest.json`，`summary.csv`，`run.log`
- Git provenance: `fedb810dcb436400971985f7673f3aec6bce8a45`

旧path `outputs/2026-08-18/174828` は上記への互換symlinkである．過去のmanifestに記録されたpathと既存の解析導線はこのsymlinkで解決する．

## Positioning

これはIssue #204のfeature解析入力として固定する診断参照である．datasetのfreeze，物理modelの採択，ML training candidateの採択を自動的に意味しない．`n_flagella=4` はdiagnostic-onlyであり，現行training scopeには含めない．

raw artifactは約24 GBのためGit管理しない。`data/`は入力data用かつGit ignoreであり，本simulation outputの保管先には使用しない．

## Analysis contract

`scripts/03_dataset_building/analyze_motion_features.py --config conf/phase2_analysis/issue204_motion_feature_study.yaml` はgeneric multi-runのraw runを直接読む共通解析である．3D / 2Dのbody mean speed，body-axis angular velocity，mean flagella-axis angular velocity，body--flagella axis angleをframe時系列と`0.25 / 0.5 / 1.0 s` non-overlap windowへ保存する．

2D body axisはcanonical body-only silhouetteから得るpixel-observable量である。一方、tracking-center前のprojected centroid速度、mean flagella axis、body--flagella relationはsimulation GT由来のprojected latent featureとしてmanifestへ明記する。投影面はYAMLの直交2基底で指定し、既定はXY面である。
