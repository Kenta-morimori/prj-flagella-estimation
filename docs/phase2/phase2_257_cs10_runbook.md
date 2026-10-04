# Issue #257: cs10 parallel runbook

Issue #257は`execution:cs10`のUser-run診断である。jobのenqueue、dispatcher/tmuxの起動・停止は、各操作についてユーザーの明示承認を要する。screenの予約・起動承認は2秒診断へ継承しない。

## Jobs

| stage | job | conditions |
| --- | --- | ---: |
| screen | `hex_screen_job.yaml` | 26 |
| screen | `project_screen_job.yaml` | 6 |
| 2 s diagnostic | `hex_2s_job.yaml` | 13（OFFのみ） |

残す3 jobは`geometry_all_conditions`、`max_workers: auto`、`worker_policy: cs10_qualified`を固定する。hexとprojectはbase configが異なるため別jobである。project 2秒jobは実行せず、PR #259から削除する。

実行済みhexの2秒比較では、Issue #245のcompleted 2秒archive（commit `f41a693`、反発ON、13 attachment配置）を固定参照として再利用し、反発OFFの13条件だけを新規実行した。当時の共通evaluatorは旧sourceのselector省略を既定値`true`として正規化した。この#257専用の結合特例は現行評価器から除く。過去診断の再解析には固定commit `de05985` と当時のcompleted archiveを使う。hook角度のraw first-failはdiagnostic-only、motor residualのstrict FAILは維持し、#257結果もdiagnostic-onlyとする。

screen reservation #13/#14はfixed-commit contractであり、本変更に伴って停止、取消、差替えしない。projectの2秒比較は実施しない。

PR #259の追加候補`hex_hook_neutral_screen_job.yaml`は別のpending 26条件screenで、既存job・archive・reservationを変更しない。`flagella.initial_hook_force_neutral=true`をこの新設定だけに適用し、geometry preflightで初期hook力ゼロとbead clearanceを確認する。dry-runまではローカルで実施可能だが、cs10接続、enqueue、dispatcher起動は各操作への明示承認を受けるまで行わない。short screenの結果をユーザーがreviewするまで、2秒実行へ進めない。

## Dry-run

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  --config conf/phase2_parallel/issue257_body_flagella_contact/hex_screen_job.yaml --dry-run
```

残るjobも同じ操作を行い、planned status、全condition、childごとのoutput_dir、8 workerを確認する。

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

実行当時は共通evaluatorでQC heatmap、window QC、stall summary/time-series、fixed-camera replayを生成した。現行HEADでは#257専用stall・pair処理を除去し、completed `attachment_pattern` archiveの初期形状図を標準出力する。partial archiveとcs10 operational logは採択入力・同期対象にしない。

## 実行済みhex 2秒診断の判断

#245のON 13条件と#257のOFF 13条件のcompleted archiveをローカルでSHA-256照合し、共有workspaceの`outputs/2026-10-04/131754/issue257_hex_2s_on_off_visualization/`へ比較bundleを同期した。13ペアすべてで40 ms stall候補はON/OFFとも0窓のため、この定義では改善の有無を判定できない。hookの生記録は初期配置時点の角度誤差によるdiagnostic-onlyであり、26 armのstrict FAILはmotor torque residual超過による。詳細な結果と未解決の初期hook設計整合性は比較契約と`phase2_tasks.md`のP2-D26を参照する。
