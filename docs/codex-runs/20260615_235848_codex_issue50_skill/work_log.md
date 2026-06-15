# Codex Run Log

## Task

Issue #50「SKILL.mdの追加」に対応し、issue 単位作業で Codex の context / tool / simulation resource を削減するための repository-local skill を追加する。

## Scope

- Add a repository-reviewed skill under `tools/codex/skills/flagella-issue-workflow/`.
- Keep `SKILL.md` concise and move detailed routing / completion / resource policy into references.
- Record whether SKILL.md can reduce resources and list additional reduction ideas for future issues.

Out of scope:

- Install the skill into `~/.codex/skills`.
- Change Phase 2 physical simulation code.
- Run Phase 2 simulations.

## Implementation

- Created `tools/codex/skills/flagella-issue-workflow/SKILL.md`.
- Added `references/context-routing.md` to route required project documents by task type.
- Added `references/completion-policy.md` to summarize review_result, commit, push, and PR completion handling.
- Added `references/resource-reduction.md` to document expected resource savings and future issue template ideas.
- Fixed generated `agents/openai.yaml` default prompt so it references `$flagella-issue-workflow`.

## Verification

- `python3 /Users/kentatakemori/.codex/skills/.system/skill-creator/scripts/quick_validate.py tools/codex/skills/flagella-issue-workflow`
- `rg -n "TODO|\\[TODO|Use -issue" tools/codex/skills/flagella-issue-workflow`

## Result

The skill validates successfully. This is a workflow/docs-only change, so Phase 2 simulations and full pytest were not run.
