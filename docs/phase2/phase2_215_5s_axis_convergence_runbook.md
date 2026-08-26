# Issue #215 5.0 s axis-angle convergence runbook

## Scope

Issue #204 の固定診断参照runと同じ2010 project、tau-linked、RUN fixed、Brownian OFF、
`T=2.5e-20 N m` per flagellum、`dt_star=1e-3`を、2.0 sから5.0 sへ延長する。
`n_flagella=1,2,3,4`とattach / phase seed各3条件の36 independent runsを扱う。これは
body--flagella axis angleの長時間推移を確認する診断であり、physical model、dataset、ML policyを
採択・変更しない。`n=4`はdiagnostic-onlyである。

Heavy/runtime execution targetは`cs10_user_run`（label: `execution:cs10`）である。Codexはcs10への
接続、tmux起動、job開始・停止を行わない。以下はユーザーがcs10上で実行する。

## Qualification

PRのcommitを取得した対話cs10 sessionで、まず全36 shardを`0.001 s`へ短縮してqualificationする。

```bash
cd ~/src/prj-flagella-estimation
git pull --ff-only origin codex/issue-215-5s-axis-convergence
git rev-parse --short HEAD

.venv-cs10/bin/python scripts/cs10/parallel_tmux.py start \
  --config conf/phase2_parallel/issue215_5s_axis_convergence/qualification_job.yaml \
  --session issue215-qualification --label issue215_5s_qualification
```

```bash
.venv-cs10/bin/python scripts/cs10/parallel_tmux.py status \
  --control-dir outputs/YYYY-MM-DD/HHMMSS/cs10_parallel/issue215_5s_qualification
```

`job=succeeded`、`failed_configs=[]`、aggregate=`completed`、
`campaign.run_summary_count=36`を満たす場合だけ本実行へ進む。

## 5.0 s campaign and analysis

同じworker policy（最大8 worker、各BLAS thread=1）で実行する。

```bash
.venv-cs10/bin/python scripts/cs10/parallel_tmux.py start \
  --config conf/phase2_parallel/issue215_5s_axis_convergence/job.yaml \
  --session issue215-5s --label issue215_5s_axis_convergence
```

完了時の`status`で上記4条件を再確認する。表示された`campaign_root`を`$CAMPAIGN`に設定し、
#204と同じmotion-feature解析をcs10で実行する。

```bash
CAMPAIGN=outputs/YYYY-MM-DD/HHMMSS/parallel/issue215_tau_linked_5s_axis_convergence__UUID/campaign

.venv-cs10/bin/python scripts/03_dataset_building/analyze_motion_features.py \
  --config conf/phase2_analysis/issue204_motion_feature_study.yaml \
  run_dir="$CAMPAIGN" output_dir="$CAMPAIGN/analysis/motion_features" overwrite=true
```

## Review points and transfer

`run_summary.json`を各conditionのfirst-read artifactとし、strict非PASS、first-fail時刻・category、
欠損artifactを一覧化する。`analysis/motion_features/`の3D / 2D時系列と`0.25 / 0.5 / 1.0 s` window
plotから、body--flagella axis angleの2.0--5.0 s推移をnとseed別に確認する。未完了・失敗conditionを
PASSへ補完したり、plotの平均へ混入させたりしない。

Macへ移送するのはreviewに必要なmanifest、`run_summary.json`、motion-feature CSV / plotなどの
compact artifactに限定する。raw archiveの一括移送、dataset freeze、model採択は本Issueのscope外である。
