# Issue #257: cs10 parallel runbook

Issue #257は`execution:cs10`のUser-run診断である。jobのenqueue、dispatcher/tmuxの起動・停止は、各操作についてユーザーの明示承認を要する。screenの予約・起動承認は2秒診断へ継承しない。

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
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  --config conf/phase2_parallel/issue257_body_flagella_contact/hex_screen_job.yaml --dry-run
```

他の3 jobにも同じ操作を行い、planned status、全condition、childごとのoutput_dir、8 workerを確認する。

## cs10 runtime の固定

cs10では通常の`.venv`、`uv run`、PATH上の`python`を実行に使わない。必ずcheckout rootの`.venv-cs10/bin/python`を明示する。これは`requirements/cs10.txt`から構築されたCentOS 7向けの隔離runtimeである。

非対話SSHでは`~/.local/bin`がPATHへ継承されず、`uv`が見つからない場合がある。`uv`を使う必要があるruntime再構築は、対話sessionで`bash scripts/cs10/setup_environment.sh`を実行する。既存`.venv-cs10`の再構築は環境を変更するため、実行予約とは別にユーザー承認を要する。

実行前には次を確認する。失敗時は`.venv`へfallbackせず、runtime setupの問題として止める。

```bash
.venv-cs10/bin/python --version
.venv-cs10/bin/python -c 'import yaml, numpy; print("cs10 runtime imports: ok")'
```

今回のfixed commit reservationでは、queue dispatcherにも`CS10_RUNTIME_PYTHON=/home/people/Ktakemori/src/prj-flagella-estimation/.venv-cs10/bin/python`を渡し、workerが同じruntimeを使うことを固定する。

## Execution and review

screen完了後、manifest、summary、run_summary、performance、state archiveをローカル同期し、SHA-256とQCを照合する。motor residualはstrict FAILとして残す。finite/body/hook length/flagを確認した後、ユーザーが2秒jobを明示承認した場合だけ実行する。

2秒完了後は共通evaluatorでQC heatmap、window QC、stall summary/time-series、fixed-camera replayを生成する。partial archiveとcs10 operational logは採択入力・同期対象にしない。
