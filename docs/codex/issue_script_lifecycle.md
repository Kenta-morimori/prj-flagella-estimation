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
