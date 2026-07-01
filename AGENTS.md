# AGENTS.md

This file defines the minimum repository-level instructions for Codex.

## Project Overview

This repository develops a pipeline to estimate bacterial flagella counts from swimming microscopy videos.

Phases:

1. Phase 1: repository, CLI, config, logging, reproducibility foundations.
2. Phase 2: 3D physical simulation and 2D pseudo-microscopy video generation.
3. Phase 3: cell detection and per-cell clip generation.
4. Phase 4: flagella-count model training and evaluation.
5. Phase 5+: prediction visualization and real-data analysis support.

Current focus: **Phase 2**.

## Source Of Truth

Read only what is needed for the task.

Priority:

1. Repository-wide rules: `AGENTS.md`
2. Current Phase 2 entry point: `docs/phase2/phase2_current.md`
3. Accepted task status: `docs/phase2/phase2_tasks.md`
4. Overall project map: `docs/PROJECT_PLAN.md`
5. Codex workflow details: `docs/codex/codex_workflow.md`
6. ADRs and task-specific docs under `docs/`
7. Historical context: Git history, issue/PR history, and `docs/codex-runs/*/review_result.json`

There is no `prompts/` source of truth in the current repository. Do not recreate prompt-based workflow files unless a new task explicitly reintroduces them.

## Context Reading Policy

* Always start from this file, the user's latest request, and the target issue or PR when provided.
* For Phase 2 work, read `docs/phase2/phase2_current.md` before longer project documents.
* Use `rg -n` before opening long Markdown files, logs, CSVs, or generated outputs.
* Read only relevant sections by default.
* For prior Codex runs, read `review_result.json` before `work_log.md`.
* Do not read large generated files under `outputs/` unless summaries and docs are insufficient.

## Language Policy

* Communicate with the user in Japanese unless explicitly requested otherwise.
* Write user-facing project documents in Japanese by default, especially under `docs/phase*/` and `docs/codex-runs/`.
* Keep technical identifiers, config keys, filenames, column names, commands, and accepted domain terms in their original form when translation would reduce precision.
* Prefer concise explanations with clear evidence and verification results.

## Repository Rules

* Do not work directly on `main` or `master`.
* Check current branch and `git status` before making changes.
* Keep changes within the requested task scope.
* Do not make broad refactors unless required by the task.
* Do not add dependencies without explaining why they are necessary.
* Do not commit secrets, API keys, tokens, private data, or generated credentials.
* Do not run remote scripts directly, such as `curl ... | sh`, unless explicitly approved by the user.
* Do not create GitHub issues unless explicitly requested.
* Do not merge pull requests.
* When creating a PR for an issue, link the original source issue in the PR body using an issue-closing keyword when the PR is intended to complete it.
* Target the branch specified by the task or issue. If no target branch is specified, target the repository default branch.
* Do not mark a task complete without a local `review_result.json` with `"status": "PASS"`.

## Directory Responsibilities

* `scripts/`: user-facing CLI entrypoints; keep orchestration, config loading, output setup, logging, and manifest handling here.
* `src/`: reusable implementation and core algorithms.
* `src/sim_swim/`: Phase 2 simulation implementation.
* `src/flagella_estimation/`: Phase 3-4 detection and ML implementation.
* `tools/codex/`: Codex-specific workflow helpers.
* `docs/phase*/`: user-facing phase docs; new files should start with a phase prefix such as `phase2_6_...md`.
* `docs/codex/`: Codex workflow details that are too long for this file.
* `docs/codex-runs/`: Codex run logs and review results.

## Output And Reproducibility

* Store project outputs under `outputs/YYYY-MM-DD/HHMMSS/` using JST.
* Each applicable run should produce `run.log` and `manifest.json`.
* `manifest.json` should record config paths, overrides, seeds, input/output paths, Git commit information, and environment details when available.
* Phase 2 simulation logs use `step_summary.csv` as the unified per-step summary.
* Do not reintroduce `step_summary_full.csv` unless a task explicitly changes this decision.

## Phase 2 CLI Command Convention

For Phase 2 simulation, sweep, and heatmap commands, prefer `KEY=VALUE` arguments in new user-facing examples:

`uv run python -m scripts.01_simulate_swimming config=conf/sim_swim.yaml time.duration_s=0.5 time.dt_star=1.0e-4`

`uv run python scripts/01_simulate_swimming/run_sweep.py config=conf/phase2_sweeps/hook_overstretch.yaml dry_run=true sample_limit=3`

Option-style arguments such as `--config`, `--duration-s`, and `--fps-out` remain only for legacy compatibility.

## Testing And Review

* Add or update pytest tests when changing physics, geometry, output format, or pipeline behavior.
* Prefer library-level tests over slow subprocess-based CLI tests.
* Stability-sensitive Phase 2 behavior may require multi-step tests; one-step tests are often insufficient.
* Do not rely only on manual video inspection.
* If tests or simulations cannot run, record what was not run, why, and what alternative checks were performed.
* Completion, review_result schema, commit, push, PR, ADR, and Cloud review details live in `docs/codex/codex_workflow.md`.

## Phase 2 Visual Review

Phase 2 may require user visual review when automated metrics cannot judge behavior fully, including:

* posterior or side bundling,
* collapse or fly-away behavior,
* body deformation,
* helical flagellar shape preservation,
* swimming-like propulsion,
* naturalness of 2D pseudo-microscopy video,
* usability for ML training.

When visual review is required, report the command, expected output directory, files to inspect, evaluation points, automated checks already passed, and what cannot be determined automatically.

Standard command pattern:

`uv run python -m scripts.01_simulate_swimming time.duration_s=0.5 time.dt_star=1.0e-4 ...`

If visual review is required but not completed, the review result remains `FAIL`.

## Phase 2 Modeling Policy

When working on Phase 2, distinguish:

* reference paper model,
* current repository implementation,
* project-specific extensions,
* numerical stabilizing approximations.

Important modeling choices must be documented in task logs or ADRs when they affect interpretation, reproducibility, or downstream ML data generation.

Failures such as collapse, fly-away, unstable hook motion, body deformation, or loss of helix geometry can be meaningful diagnostics. Preserve reproducible failure conditions when useful.

Phase 2 baseline configs should omit `stiffness_scales` when all values are parser defaults (`1.0`). If non-default stiffness scales are kept, document the reason.

For Phase 2 simulation runs, do not change the default config just to set the integration step. Use explicit runtime overrides such as `time.dt_star=1.0e-4` unless the task requires another value.

## Final Response

Final response should be sent only after the task state has been committed and pushed when remote access is available.

After each task, report:

* summary,
* changed files,
* tests/checks run,
* user visual review required or not,
* review result,
* ADR created or not created,
* commit hash,
* push status,
* remaining issues.
