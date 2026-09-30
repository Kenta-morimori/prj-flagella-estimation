# Issue #257: cs10 parallel runbook

Issue #257は`execution:cs10`のUser-run診断である。Codexはjobをenqueue、dispatcher/tmuxを起動、停止しない。各操作はユーザーの明示承認を要する。

## Jobs

| stage | job | conditions |
| --- | --- | ---: |
| screen | `hex_screen_job.yaml` | 26 |
| screen | `project_screen_job.yaml` | 6 |
| 2 s diagnostic | `hex_2s_job.yaml` | 26 |
| 2 s diagnostic | `project_2s_job.yaml` | 6 |

全jobは`geometry_all_conditions`、`max_workers: auto`、`worker_policy: cs10_qualified`を固定する。hexとprojectはbase configが異なるため別jobである。

## Dry-run

```bash
uv run python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/issue257_body_flagella_contact/hex_screen_job.yaml dry_run=true
```

他の3 jobにも同じ操作を行い、planned status、全condition、childごとのoutput_dir、8 workerを確認する。

## Execution and review

screen完了後、manifest、summary、run_summary、performance、state archiveをローカル同期し、SHA-256とQCを照合する。motor residualはstrict FAILとして残す。finite/body/hook length/flagを確認した後、ユーザーが2秒jobを明示承認した場合だけ実行する。

2秒完了後は共通evaluatorでQC heatmap、window QC、stall summary/time-series、fixed-camera replayを生成する。partial archiveとcs10 operational logは採択入力・同期対象にしない。
