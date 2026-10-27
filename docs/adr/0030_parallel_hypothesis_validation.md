# ADR 0030: 複数仮説の検証計画と cs10 の実行順

- Status: accepted
- Date: 2026-10-27
- Source: Issue #273

## 背景

原因が未解明で有力な仮説が複数あるとき、仮説を一件ずつ検討すると、共有 cs10 の実行待ちの間に進められる調査と実験設計も止まる。cs10 では job 内の独立 condition は並列化できるが、異なる固定 commit の job は逐次 queue で扱う。

## 判断

仮説ごとの予測・反証条件・識別指標・短時間 screen・判定基準を並行して準備する。複数仮説を区別できる共通実験を優先し、screen の結果から長時間 job の必要性と順序を見直す。結果は支持・棄却・未判定に分け、未判定を採択根拠にしない。cs10 の parallel job 契約、fixed-commit queue、操作ごとの明示承認は維持する。

具体的な手順の正本は `docs/codex/codex_workflow.md` とする。
