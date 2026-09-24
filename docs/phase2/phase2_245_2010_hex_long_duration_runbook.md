# Issue #245: 2010 hex 13 attachment pattern runbook

## Scope

対象は**`2010_hex_project`のみ**である。2010 project modelのnf5/nf6は含めない。初期らせん位相は`phase_seed=0`に固定し、六角中心環のattachment slotだけを変更する。回転（C6）のみを同一視し、反射は区別する。

| n | canonical slot集合 |
| ---: | --- |
| 1 | `0` |
| 2 | `01`, `02`, `03` |
| 3 | `012`, `013`, `014`, `024` |
| 4 | `0123`, `0124`, `0134` |
| 5 | `01234` |
| 6 | `012345` |

condition IDは`nf03__slots013`形式で全13条件である。n=6の全slot占有は`full_ring_rotation_equivalent: true`としてmanifest・summaryへ保存する。`T=2.5e-20 N m / flagellum`、`dt_star=1e-4`、RUN、Brownian/switching OFF、archive interval `0.001 s`、archive保存ONを固定する。

screenは1τ（10,000 steps）、mainは2.0 s（reference torqueで50τ、500,000 steps）である。cs10 targetは`cs10_user_run`（`execution:cs10`）。新規jobは`max_workers: auto`と`cs10_qualified`を使い、実効8 workers、`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1`である。1τ実測に基づく暫定見積りはscreen約42分、main約35時間、合計約36時間であり、screen実測後に更新する。

## 予約前確認と連続実行

予約・dispatcher・tmux・job開始・停止はすべて操作ごとのユーザー明示許可が必要である。予約・起動の前に、ユーザーはreview済みPR commit、clean worktree、13 shard dry-run、NAS容量、
queue/GitHub認証、既存reservationとの衝突なしを確認する。以下は予約を作成しないpre-reservation
checklistである。

```bash
cd ~/src/prj-flagella-estimation
git fetch origin
git show --quiet --format='%H %s' origin/codex/issue-245-2010-hex-long-duration
git status --porcelain
.venv-cs10/bin/python scripts/01_simulate_swimming/run_parallel.py \
  --config conf/phase2_parallel/issue245_2010_hex_long_duration/attachment_screen_job.yaml \
  --dry-run
df -h /net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs
gh auth status --hostname github.com
.venv-cs10/bin/python scripts/cs10/queue.py status
```

ユーザーはscreenとmainを**連続してreservation / dispatcher実行してよい**と指定している。screen FAILの場合もmainを停止せず実行できるが、mainのarchive、heatmap、replayはstrict failure診断専用であり、特徴量評価・採択・canonical化の根拠にしてはならない。キャンセル等の変更はこの例外からは許可されず、別操作として明示許可を要する。

```bash
cd ~/src/prj-flagella-estimation
git fetch origin
.venv-cs10/bin/python scripts/cs10/queue.py enqueue \
  --branch origin/codex/issue-245-2010-hex-long-duration \
  --config conf/phase2_parallel/issue245_2010_hex_long_duration/attachment_screen_job.yaml \
  --priority 0
```

```bash
.venv-cs10/bin/python scripts/cs10/queue.py status
```

screenとmainのenqueueおよびdispatcher起動は、実行時にそれぞれユーザーが明示許可した場合だけ行う。予約はfixed commit contractであり、branch更新だけで差し替えない。`progress.json`、`diagnostic_samples.csv`、`state_archive.partial.npz`はstrict failure時のdiagnostics/replay専用であり、completed archiveや採択根拠として使わない。

## 2 s campaign, sync, and aggregation

screen FAILでもユーザー指定の連続mainは実行可能である。ただしstrict failureが一件でもあれば、main成果物はdiagnostic-onlyであり、特徴量評価・採択・canonical化へ進めない。

```bash
.venv-cs10/bin/python scripts/cs10/queue.py enqueue \
  --branch origin/codex/issue-245-2010-hex-long-duration \
  --config conf/phase2_parallel/issue245_2010_hex_long_duration/job.yaml \
  --priority 0
```

完了したcampaignをMacへ同期する際は、再解析に必要な各conditionの`state_archive.npz`、`run_summary.json`、`performance.json`、root manifestとsummaryを保持し、SHA-256を検証する。

```bash
.venv/bin/python scripts/cs10/sync_reference_from_cs10.py \
  --host cs10 --remote-dir "$CAMPAIGN" \
  --local-dir outputs/YYYY-MM-DD/HHMMSS/issue245_2010_hex_2s

.venv/bin/python scripts/03_dataset_building/analyze_dataset.py \
  --analysis-kind model-development-evaluation \
  --config conf/phase2_multi_run/2010_hex_project_long_duration_2s_issue245.yaml \
  --run-dir outputs/YYYY-MM-DD/HHMMSS/issue245_2010_hex_2s \
  --output-dir outputs/YYYY-MM-DD/HHMMSS/model_development_evaluation \
  --render-replay
```

集約器は`manifest.json`、`run.log`、統合summary、`window_qc.csv`、first failure、wall time、steps/s、artifact SHA-256を作る。n行・canonical slot pattern列のsparse heatmapでは該当しないセルをmaskする。replayはn/pattern別にページングした固定camera 3D/2Dで、screen/main、PASS/FAIL、diagnostic-onlyをmanifestと画面へ記録する。hook angleは値とfirst failureを保存するdiagnostic-onlyであり、finite/body/hook length/flag/motorだけがstrict PASS/FAILである。
