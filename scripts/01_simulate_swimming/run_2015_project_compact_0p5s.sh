#!/usr/bin/env bash
set -euo pipefail
# Issue #186.  Run on a separate machine/worktree; this intentionally starts
# no render or replay job.  Archive cadence is 1 ms (1000 fps maximum replay).
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  conf/phase2_multi_run/2015_project_compact_0p5s.yaml
