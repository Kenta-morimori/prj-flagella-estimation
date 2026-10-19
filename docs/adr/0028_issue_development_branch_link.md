# ADR 0028: Issue 用 branch と GitHub Development の連携

- Status: accepted
- Date: 2026-10-19
- Source: Issue #267（親 Issue #266）

## 背景

従来の workflow は task-specific branch と PR 本文での source Issue 参照を求めていたが、branch 作成時の Issue Development 欄との連携確認を定めていなかった。#262 / PR #265 では事後対応が必要になった。#263 / PR #264 は作業中の事例として扱い、既存関係を変更しない。

## 判断

Issue 用 branch の作成前に linked branch を確認し、未連携なら `gh issue develop` で task / Issue の target branch（指定がなければ default branch）から branch を作る。既存・削除済み branch は現存する Development 欄と PR の関係を先に確認し、重複した連携や履歴の復元を避ける。PR 本文にも source Issue を記し、完了を意図しない PR に `Closes` / `Fixes` を付けない。linked PR の merge による Issue 自動 close の可能性があるため、merge 後は残受入条件と Issue 状態を照合する。

具体的な手順の正本は `docs/codex/codex_workflow.md` とする。
