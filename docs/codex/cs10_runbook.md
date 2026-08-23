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

続けて fixed workload benchmark を実行する。これは任意 config の並列 launcher ではなく、Issue #209 の worker default を決めるための qualification 専用 runner である。

```bash
.venv-cs10/bin/python scripts/cs10/benchmark.py
```

benchmark は `shape_stability_grid.yaml` の 1 条件を、`n_flagella=3`、`duration_s=0.5`、`dt_star=1e-4`、固定 seed、state archive 有効で実行する。workers `1,2,4,6,8,10` と、thread 数未設定 / `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` を各 5 回比較する。

`summary.csv` と `manifest.json` を確認し、全 5 試行成功かつ最大 throughput の setting を候補とする。共有利用のため 2 physical core を留保し、通常 default は最大 8 workers とする。10 workers は上限確認の記録のみである。

## Scope and requalification

RTX 3090 は hardware record のみであり、CUDA、PyTorch、GPU benchmark はこの scope に含まれない。OS/glibc、CPython、cs10 requirements、主要 simulation 実装、hardware が変わった場合は setup、probe、benchmark を再実行する。
