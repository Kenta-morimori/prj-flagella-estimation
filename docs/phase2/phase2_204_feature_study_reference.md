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

3D GTは`state_archive.npz`の全simulation stepを使う。2Dはcampaign manifestの`output_sampling.fps_out_2d`で選んだ出力映像frame時刻を使い、pixel-observable body axisとprojected latent featureを同じ観測時刻で比較する。`frame_rate_hz`を3D GTの間引きに使うことは禁止し、manifestから2D frame rateを一意に得られないcampaignは明示設定なしでは解析を失敗させる。windowは`[start_s, end_s)`の実時間非重複区間であり、2.0 s runでは`0.25 / 0.5 / 1.0 s`がそれぞれ`8 / 4 / 2`窓となる。CSVにはsampling sourceと各windowのsample countを残す。

解析結果は参照run配下の`analysis/motion_features/`へ置く。`time_series/`には元時系列の`3D_*.png` / `2D_*.png`、`windows/`には3つのwindow幅を横並びにした同名形式のplotを置く。個別runは半透明の破線、`n`ごとの平均は太い実線で描く。高周波のflagella-axis angular velocityだけは可読性のためplot上で20 ms bin平均とし、raw CSVは間引かない。各CSVは同rootに保持し、strict QCに失敗した`n=4`条件はCSVとmanifestには残すが、両plot群から除外する。

2D body axisはcanonical body-only silhouetteから得るpixel-observable量である。一方、tracking-center前のprojected centroid速度、mean flagella axis、body--flagella relationはsimulation GT由来のprojected latent featureとしてmanifestへ明記する。投影面はYAMLの直交2基底で指定し、既定はXY面である。
