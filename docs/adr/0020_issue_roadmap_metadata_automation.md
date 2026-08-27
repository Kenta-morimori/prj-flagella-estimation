# ADR 0020: Issue Roadmap metadata の自動同期

## Status

Adopted

## Context

GitHub Project #8 の Roadmap は Milestone、`Start date`、`Target date` を使うが、
手動入力だけでは新規Issueの分類漏れや終了日の欠損が起こる。GitHub Issue Formは
MilestoneやProject fieldを直接設定できず、user-owned Projectは`GITHUB_TOKEN`から
更新できない。

## Decision

- 新規Issueは必須の`Roadmap category (Milestone)`を選択する。選択肢は現行の
  8 Milestoneと一対一に対応させる。
- `issue-roadmap-sync` workflowがIssue lifecycleをProject #8へ同期する。
  Start dateは作成日（JST）、Planned target dateは任意の予定日、未設定の
  Target dateはclose時の終了日（JST）とする。
- 既存Target dateはclose時に上書きしない。
- Form外作成・不正なmetadataは`roadmap:triage`、reopen後は
  `roadmap:needs-review`で明示し、修正前の実装着手を禁止する。
- Project更新にはrepository secretの`PROJECT_AUTOMATION_TOKEN`を使う。これは
  user-owned Project #8を更新できるclassic PAT（`repo` + `project` scope）である。

## Consequences

- Issue Formのカテゴリ変更とMilestone構成変更は同じPRで同期する必要がある。
- PAT未登録時はworkflowが明示的に失敗し、Project metadataを部分更新しない。
- 予定日と実績終了日を同一`Target date`に保存するため、予定が存在したclose
  Issueの実績終了時刻はGitHub Issueの`closed_at`を正本とする。
