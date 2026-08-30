# Issue専用スクリプト lifecycle

Issue検証のために追加した script は、merge前に `promoted`、`deleted`、
`retained-temporarily` のいずれかへ分類する。`promoted` は共通CLIまたは`src`の
再利用可能なAPIへ移す。`deleted` は active docs/config/tests/scripts/src の参照をすべて
移行してから削除する。`retained-temporarily` は削除Issue・期限・利用者をPR本文と
`review_result.json`へ記録する。

通常の検証出力は `outputs/YYYY-MM-DD/HHMMSS/` をrun rootとし、派生物は
`dataset/<id>/`、`replay/<kind>/<timestamp>/`、`analysis/<kind>/<timestamp>/`へ出す。
固定output pathはcanonical datasetなど明示的に再利用するものだけに使う。

`tools/codex/check_issue_script_lifecycle.py --base <base>` は、削除された`scripts/`配下の
pathを active consumer が参照していないことを確認する。`docs/codex-runs/`は歴史的記録のため除外する。

`scripts/cs10/queue.py` はIssue #219で追加するcs10共通運用CLIであり、`promoted`として維持する。
予約の固定commit、worktree隔離、逐次dispatch、通知はcs10 runbookを正本とする。
