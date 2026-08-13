# AGENTS.md

This file defines the minimum repository-level instructions for Codex.

## Project Overview

This repository develops a pipeline to estimate bacterial flagella counts from swimming microscopy videos.

Phases:

1. Phase 1: repository, CLI, config, logging, and reproducibility foundations.
2. Phase 2: 3D physical simulation and 2D pseudo-microscopy video generation.
3. Phase 3: cell detection and per-cell clip generation.
4. Phase 4: flagella-count model training and evaluation.
5. Phase 5+: prediction visualization and real-data analysis support.

## Context Routing

Read only what is needed for the task, in this order:

1. `AGENTS.md`
2. The user's latest request and the target Issue or PR
3. `docs/phaseX/phaseX_current.md`
4. `docs/phaseX/phaseX_guide.md` when it exists and phase-specific rules are relevant
5. The relevant section of `docs/phaseX/phaseX_tasks.md` when past decisions are needed
6. Live schemas, contracts, configs, tests, and active validation documents
7. ADRs
8. Issue / PR history, Git history, and `docs/codex-runs/*/review_result.json`

Use `rg -n` before opening long Markdown files, logs, CSVs, or generated outputs.
For prior Codex runs, read `review_result.json` before `work_log.md`.
Do not read large files under `outputs/` unless compact summaries and manifests are insufficient. For Phase 2 diagnostics, read `run_summary.json` before using the bounded `inspect_step_summary.py` CLI; never load `step_summary.csv` in full for routine analysis.

## Language

* Communicate with the user in Japanese unless explicitly requested otherwise.
* Write user-facing project documents in Japanese by default.
* Keep technical identifiers in their original form when translation reduces precision.

## Repository Rules

* Do not work directly on `main` or `master`.
* Check the current branch and `git status` before making changes.
* Keep changes within scope; do not make broad refactors unless required.
* Do not add dependencies without explaining why they are necessary.
* Do not commit secrets, tokens, credentials, private data, or generated authentication files.
* Do not run remote scripts such as `curl ... | sh` without explicit user approval.
* Target the branch specified by the task or Issue; otherwise target the default branch.
* Track related work with GitHub-native relationships: make a bounded child task a sub-issue of its parent, and add `blocking` / `blockedBy` only for a real completion dependency. Do not create a dependency edge merely because work is related.
* Link the source Issue from a PR. Use `Closes #<issue>` only when merge completes that Issue; otherwise use a non-closing reference and state what remains (for example, a user-run experiment or result review).
* Do not mark a task complete without a local `review_result.json` whose status is `PASS`.
* Do not merge unless required checks and `codex-review-gate` pass.
* Do not merge changes to physical interpretation, dataset adoption, phase boundaries, output contracts, or ML policy without explicit user approval.
* After merge, sync the default branch and delete the completed task branch unless intentionally retained.

## Directory Responsibilities

* `scripts/`: user-facing CLI entrypoints and orchestration.
* `src/`: reusable implementation and core algorithms.
* `conf/`: reproducible runtime configuration.
* `schemas/`: machine-readable contracts.
* `docs/phase*/`: phase state, decisions, contracts, and active validation documents.
* `docs/adr/`: important design decisions.
* `docs/codex/`: Codex workflow and documentation policy.
* `docs/codex-runs/`: run logs and review results.
* `tools/codex/`: Codex workflow helpers and skills.

## Phase Documentation

* Follow `docs/codex/phase_document_policy.md`.
* Treat `phaseX_current.md` as the current-state entry point.
* Treat `phaseX_tasks.md` as the compact decision record.
* Keep task-specific documents only when they remain live contracts, active validations, or reusable reports.
* Update current, tasks, and ADRs when a semantic decision changes.
* Before deleting a document, preserve its decision-bearing information and update all references.
* Use the `phase-document-maintenance` skill for consolidation, migration, or deletion.

## Output And Reproducibility

* Store project outputs under `outputs/YYYY-MM-DD/HHMMSS/` using JST.
* Applicable runs should produce `run.log` and `manifest.json`.
* Record configs, overrides, seeds, paths, Git information, and environment details when available.
* Preserve full diagnostics when needed, but prefer compact summaries for routine analysis.
* Phase 2 uses `step_summary.csv`; do not reintroduce `step_summary_full.csv` without an explicit decision.

## Testing And Review

* Keep default pre-commit checks lightweight.
* Start with targeted tests and expand only when required.
* Do not require full pytest for docs-only, planning-only, or workflow-only changes.
* Add or update tests when changing physics, geometry, schemas, output formats, or pipeline behavior.
* Prefer library-level tests over slow subprocess tests.
* Do not rely only on visual inspection when automated checks are possible.
* Record tests or simulations that were not run and the reason.
* Do not run long simulations, sweeps, training jobs, or renders unless the user explicitly asks.
* For user-executed long runs, provide the command, expected outputs, evaluation points, and checks already passed.
* Documentation deletion or consolidation must include stale-reference checks.
* Request Cloud review only for the merge-ready final candidate by default.
* Resolve actionable review threads before merge.

## Completion Report

Report the summary, changed files, tests and checks, unexecuted long runs, user-review status, review result, documentation updates, ADR status, commit hash, push status, and remaining issues.
