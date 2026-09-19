---
name: model-development-evaluation
description: Use when adding or materially changing a pending Phase 2 physical-model candidate, including geometry, discretization, attachment, motor, interactions, or time integration. Do not use for feature-only analysis or an already accepted model.
---

# Model Development Evaluation

For a pending Phase 2 candidate, use the reusable `development_evaluation`
contract instead of creating Issue-specific analysis scripts.

1. Read the target Issue execution label and the candidate profile.  Confirm
   that the source profile has a `development_evaluation` declaration.
2. Define a `short_screen` with a bounded matrix that exposes the modified
   parameter and relevant numerical controls.  Preserve archives and produce
   a manifest, per-count PNG heatmaps, and fixed-camera replay MP4s from the
   generic evaluator.
3. Treat finite, body, hook length, flag, and motor QC as PASS/FAIL.  Keep
   diagnostics that are explicitly excluded by the contract visible but out
   of the gate.
4. Do not start the `long_duration` stage until the short screen is reviewed.
   Respect the Issue execution target; `execution:cs10` runs are user-run.
5. Keep swimming-feature analysis out of this workflow.  A later task may
   consume PASS long-duration archives through the feature registry.

Existing results may be reused only after validating condition coverage,
profile identity, source manifests, and Git provenance.  The generic
`model-development-evaluation` analysis consumes completed run directories;
it must not re-simulate them.
