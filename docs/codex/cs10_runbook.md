# cs10 Runbook

cs10 は CentOS 7 / glibc 2.17 の CPU 実行環境である。Mac の `pyproject.toml`、`uv.lock`、`.venv` を canonical environment とし、cs10 はこの runbook と `requirements/cs10.txt` に従う `.venv-cs10` を使う。

## Setup

GitHub SSH 鍵を agent へ登録してから、Issue 用 branch を checkout する。

```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519_github_cs10
cd ~/src/prj-flagella-estimation
git pull --ff-only
scripts/cs10/setup_environment.sh
```

setup は uv 管理の `$HOME/.local/bin/python3.11` を使い、source build を禁止して `.venv-cs10` に runtime + pytest を install する。`pyproject.toml` の Mac canonical dependency は install しない。失敗時は package を手動で source build せず、エラー全文と `outputs/.../cs10_setup_probe/` を共有する。

## Qualification

環境を作成後、最初に runtime probe を記録する。

```bash
.venv-cs10/bin/python scripts/cs10/probe_runtime.py \
  --output-dir outputs/2026-08-23/HHMMSS/cs10_runtime_probe
```

続けて fixed workload benchmark を実行する。これは任意 config の並列 launcher ではなく、Issue #209 の worker default を決めるための qualification 専用 runner である。最初に約 10-step の smoke を実行する。

```bash
.venv-cs10/bin/python scripts/cs10/benchmark.py \
  --worker-counts 1 --repetitions 1 --duration-s 0.00004
```

worker default を確認する screen は、次の 251-step 条件を使う。cs10 での実測所要時間は全 12 条件で約 6 分 32 秒（各並列 batch は約 30--36 秒）だった。

```bash
.venv-cs10/bin/python scripts/cs10/benchmark.py \
  --worker-counts 1,2,4,6,8,10 --repetitions 1 --duration-s 0.001
```

通常の短縮 benchmark は `shape_stability_grid.yaml` の 1 条件を、`n_flagella=3`、`duration_s=0.05`（実効 12,501 steps）、`dt_star=1e-4`、固定 seed、state archive 有効で実行する。workers `1,2,4,6,8,10` と、thread 数未設定 / `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` を各 5 回比較する。`--duration-s` で測定時間を明示的に変更できる。manifest には config から算出した `tau_s`、内部刻み `dt_s`、`total_steps` も記録される。

`duration_s=0.05` の full benchmark は cs10 では未実測である。251-step screen の約50倍の積分 step を含むため、短時間の qualification として直接開始せず、長時間 campaign 前の再qualificationとして実行する。

`summary.csv` と `manifest.json` を確認し、全 5 試行成功かつ最大 throughput の setting を候補とする。共有利用のため 2 physical core を留保し、通常 default は最大 8 workers とする。10 workers は上限確認の記録のみである。短縮 benchmark の結果は暫定 default とし、0.5秒以上の長時間 campaign を移す直前に同じ worker 数を再確認する。

## Qualified default (2026-08-23)

cs10 実機では、251-step（`duration_s=0.001`）screen を workers `1,2,4,6,8,10`、各 1 回で実行した。全 job が成功し、`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` での aggregate throughput は、8 workers で `0.2493 jobs/s`、10 workers で `0.3033 jobs/s` だった。peak RSS は約 55 MiB/job である。

通常 default は共有利用のため **8 workers + thread 数各 1** とする。10 workers は2 physical coreの留保方針に反するため default にしない。screen artifact は cs10 上の `outputs/2026-08-23/174000/cs10_worker_screen/` にある。

## 並列 sweep job（Issue #209）

複数の既存 sweep profile はコピーせず、`conf/phase2_parallel/<job_name>/job.yaml` から参照する。cs10 qualification を使う job には、`execution.max_workers: auto` と `execution.worker_policy: cs10_qualified` を指定する。この policy は最大8 workers と `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` を各子 process に設定する。

まず simulation を起動しない dry-run で、worker 数、command、出力 namespace を確認する。

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/example_stage_a_validation/job.yaml \
  dry_run=true
```

実行時は job ごとに `outputs/YYYY-MM-DD/HHMMSS/parallel/<job_name>__<uuid>/` が作られる。`job_manifest.json` の `status`、`failed_configs`、各 config の `exit_code` と `wall_time_s` を確認する。失敗時も他 config は回収されるため、該当 directory の `stderr.log` を確認する。

```bash
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  config=conf/phase2_parallel/example_stage_a_validation/job.yaml
```

所要時間は選択する profile の条件数・積分時間に依存する。`duration_s=0.5` 以上の長時間 campaign を開始する前には、同じ worker policy で representative workload を再qualificationする。

## 長時間 parallel job の標準運用

Issueの`execution:cs10`は、独立conditionが8以上またはMac見積りwall timeが30分超の
heavy/runtime execution targetとして選ぶ。target未記載・`execution:triage`のIssueでは、
cs10 runtimeを開始しない。

`duration_s >= 0.5`、または複数の独立conditionを実行する場合は、まず
`conf/phase2_parallel/<job>/job.yaml` の利用を検討する。serial実行は、condition間に
依存がある場合など、Issue runbookに明示した理由がある場合だけ許可する。cs10では
`execution.worker_policy: cs10_qualified`（最大8 worker、BLAS thread各1）を使う。

Git更新はGitHub SSH agentを持つ**対話**sessionで完了させる。非対話SSHではagentや
`$HOME/.local/bin`が継承されないことがあるため、`git pull`、tmux起動、長時間jobの
開始を一つのSSH one-linerへ混在させない。`tmux`はまず`command -v tmux`、次に
`$HOME/.local/bin/tmux`を確認する。user-local tmuxをビルドする必要がある場合は、
アクセス不能なinclude/lib pathを持つ`C_INCLUDE_PATH`と`LIBRARY_PATH`をunsetしてから行う。

起動・marker・statusは`parallel_tmux.py`を唯一の入口にする。helperはclean worktree、
`.venv-cs10`のruntime import、tmux、`cs10_qualified` policyを検査し、実行前に確定した
output rootとcontrol directoryを記録する。手書きの`RUN_ID`、`JOB_ROOT=$(...)`、`exit`を
含むtmux commandは使わない。

重い simulation artifact の既定出力先は NAS の
`/net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs/YYYY-MM-DD/HHMMSS/parallel/...`
とする。`Ktakemori`直下はプロジェクト単位で分け、他プロジェクトの出力と混在させない。
`launch.json`、launcher の stdout/stderr、exit marker はリポジトリ内の
`outputs/YYYY-MM-DD/HHMMSS/cs10_parallel/...` に残る。開始前にhelperがNASのプロジェクト
directoryを作成して書込みを検査し、利用できなければtmuxを起動せずに失敗する。NAS容量は
各実行前に次で確認する。

```bash
df -h /net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs
```

```bash
.venv-cs10/bin/python scripts/cs10/parallel_tmux.py start \
  --config conf/phase2_parallel/<job>/job.yaml \
  --session <unique-session> --label <job-label>

.venv-cs10/bin/python scripts/cs10/parallel_tmux.py status \
  --control-dir outputs/YYYY-MM-DD/HHMMSS/cs10_parallel/<job-label>

.venv-cs10/bin/python scripts/cs10/parallel_tmux.py attach \
  --session <unique-session>
```

`status`で job `succeeded`、`failed_configs=[]`、aggregate `completed`、campaignの
`run_summary_count`が期待condition数と一致することを確認する。canonical `conditions/`は
child artifactへのsymlinkであるため、通常の`find`の件数を合否判定に使わない。helperは
Pythonの`Path.is_file()`でsymlinkを辿って数える。`grep "status"`もGit provenanceなどを
混在させるため、合否判定に使わない。

Codexは既定でcs10へ接続・操作しない。接続、tmux起動、long job開始・停止を行うのは、
その操作についてUserが明示許可した場合だけとする。User実行のjobは、command、想定出力、
確認点、最小転送artifactをIssue runbookに記録する。

## Scope and requalification

RTX 3090 は hardware record のみであり、CUDA、PyTorch、GPU benchmark はこの scope に含まれない。OS/glibc、CPython、cs10 requirements、主要 simulation 実装、hardware が変わった場合は setup、probe、benchmark を再実行する。

## Issue #210: serial / parallel qualification

parallel launcher の実 simulation qualification は、同じ clean Git commit と同じ既存 sweep config を使い、次の順で実施する。Mac canonical environment、既存 serial CLI、既存 sweep config、物理解釈は変更しない。比較は bitwise identity を要求しない。

1. Mac で `tests/test_phase2_parallel_job.py`、`tests/test_phase2_sweep_profiles.py`、`tests/test_cs10_qualification.py`、`tests/test_phase2_qualification.py` を実行する。`shape_stability_grid.yaml` と `torque_distribution_grid.yaml` の dry-run で、2010 project profile を確認する。
2. Mac serial と cs10 serial で、各 profile の先頭1 conditionを`duration_s=0.001`へ短縮して実行する。`torque_distribution_grid`には`torque_nm=2.5e-20`も指定し、両 config を251 stepsへ揃える。
3. cs10 で `conf/phase2_parallel/issue210_2010_project/job.yaml` を `run_parallel.py` から実行する。`job_manifest.json` の status が `succeeded`、`failed_configs` が空、全 child の exit code が 0 であることを確認する。
4. 下記比較器を Mac vs cs10 serial の config ごと、および cs10 serial vs parallel job に対して実行する。

```bash
uv run python scripts/01_simulate_swimming/compare_qualification.py \
  --left outputs/.../mac_serial_campaign \
  --right outputs/.../cs10_serial_campaign \
  --output-dir outputs/YYYY-MM-DD/HHMMSS/qualification/mac_vs_cs10_shape_stability

uv run python scripts/01_simulate_swimming/compare_qualification.py \
  --parallel-job-manifest outputs/.../parallel/.../job_manifest.json \
  --serial conf/phase2_sweeps/shape_stability_grid.yaml=outputs/.../cs10_serial_shape_stability \
  --serial conf/phase2_sweeps/torque_distribution_grid.yaml=outputs/.../cs10_serial_torque_distribution \
  --output-dir outputs/YYYY-MM-DD/HHMMSS/qualification/cs10_serial_vs_parallel
```

251 stepsは、#208 の実機 worker screenと同じ最低 execution qualification 条件である。物理的な0.5秒 stabilityやdataset採択の根拠には使わない。比較器は compact な `summary.csv` と`run_manifest.json`、各 condition の`run_summary.json`を読み、Git commit、clean worktree、effective override、expected/completed steps、contiguous step index、finite / shape gate、failure category、主要 stability/residual metric を確認する。数値は `abs(a-b) <= max(1e-9, 1e-6 * max(abs(a), abs(b)))` とする。provenance の不一致または partial failure は数値が近くても FAIL とする。

共有する最小 artifact は、runtime probe の `manifest.json`、各 serial / parallel child の `summary.csv`・`run_manifest.json`・`manifest.json`・各 condition の `run_summary.json`、parallel の `job_manifest.json` である。失敗時だけ対応する `stdout.log`、`stderr.log`、`failure.json`（存在すれば）を追加する。raw `step_summary.csv` と state archive は共有不要である。

serial が失敗した場合は environment / simulation core を先に切り分け、parallel へ進まない。serial が PASS して parallel のみ失敗した場合は child command、output namespace、thread override、stderr を確認する。数値差のみの場合は failure gate と最大差分を report で確認し、根拠なく tolerance を変更しない。
