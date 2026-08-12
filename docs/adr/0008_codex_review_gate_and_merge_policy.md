# ADR 0008: Codex review gate and scoped merge policy

## Status

Accepted

Amended on 2026-08-06 to add GitHub Copilot review as an explicit fallback.

## Context

ADR 0006により，Codex review自体はChatGPT/Codex Cloud connectorへ一本化しており，repository-managed `openai/codex-action` workflowは再導入しない方針である。

`codex-review-gate`は，allowlisted Codex connectorによるreview実施とreview thread対応をrequired statusとして検査するために導入された。`main` rulesetでは`test`と`codex-review-gate`をrequired status checksとしている。

一方，Codexの利用制限によってCloud connector reviewを実行できない場合がある。Codex reviewだけを必須にすると，GitHub Copilot reviewで指摘を受け，すべて対応済みであってもmergeできない。

review resourceを抑えるため，review後の修正commitごとに再reviewを要求しない。PR履歴中のcommitに対して，CodexまたはGitHub Copilotによる有効なreviewを最低1回完了していればよい。

## Decision

Codex reviewを既定のreview手段として維持し，利用できない場合はGitHub Copilot reviewを正式なfallbackとして認める。

repository側には，PR metadata，comments，reviews，reactions，review threadsだけをGitHub APIから検査するlightweight gateを置く。status context名はrulesetとの互換性のため`codex-review-gate`のまま維持する。

### Common policy

- `.github/workflows/codex-review-gate.yml`はPR branchをcheckoutしない。
- PR branchのcode，script，action，生成物を実行しない。
- `OPENAI_API_KEY`や`openai/codex-action`は使わない。
- workflowはPR head SHAにcommit status `codex-review-gate`を付与する。
- PR comments，PR reviews，PR commits，review threadsは必要な範囲でページネーションする。
- reviewer loginは完全一致allowlistだけを許可する。
- review対象commitが現在のPR commit履歴に含まれることを要求する。
- review後に新しいcommitが追加されても，再reviewは要求しない。
- force-pushやrebaseによってreview対象commitが現在のPR履歴から消えた場合，そのreviewは無効とする。
- CodexまたはCopilotが作成した，未解決かつoutdatedでないthreadが1件でもあればfailureとする。
- CodexとCopilotの両方が存在する場合も，どちらかの有効な未解決threadを無視しない。
- 有効なCodex reviewまたは有効なCopilot reviewのどちらか一方が完了し，trusted threadがすべてresolvedまたはoutdatedならsuccessとする。
- review本文，PR本文，comment本文，commit message，変更ファイルはuntrusted inputとして扱い，codeとして実行しない。

### Codex review

- Codex reviewは原則としてmerge-readyなfinal candidateに対して1回だけ要求する。
- request commentは`@codex review <PR-commit-sha>`形式とし，指定SHAが現在のPR commit履歴に含まれることを要求する。
- Cloud connector loginは`chatgpt-codex-connector`または`chatgpt-codex-connector[bot]`の完全一致だけを許可する。
- request commentへのallowlisted connectorのthumbs-up reactionは，request更新時刻以後のreactionだけを有効とする。
- PR review responseはallowlisted login，submitted state，review対象commit，request後の時刻をAPIから検証する。
- Cloud connectorがGitHubの正式reviewではなくPR commentで応答する場合は、`@codex review <SHA>`要求（編集された要求は`updated_at`、未編集時は`created_at`以後）に、allowlisted loginが現在のPR履歴にある同じ一意のcommit SHAを`Reviewed commit`として明記し、`Didn't find any major issues`を返した場合だけ成功シグナルとして扱う。
- `APPROVED` review，`Final verdict PASS`，またはrequestへのthumbs-upも、同じreview対象commit・時刻検証を満たす場合の成功シグナルとして扱う。
- feedbackを伴うreviewは，current Codex-authored threadがすべてresolvedまたはoutdatedになった時点で完了とする。
- feedback修正後にPR headが変わっても，再度`@codex review`を要求しない。

### GitHub Copilot review fallback

- GitHub Copilot reviewはCodexを利用できない場合の明示的なfallbackとする。
- reviewer loginは`copilot-pull-request-reviewer`または`copilot-pull-request-reviewer[bot]`の完全一致だけを許可する。
- reviewがsubmitted済みであり，`DISMISSED`でないことを要求する。
- REST APIの`review.commit_id`が現在のPR commit履歴に含まれることを要求する。
- Copilot review本文の文言はPASS判定に使用しない。
- Copilot reviewの存在，reviewer login，review commit，threadの作成者，`isResolved`，`isOutdated`をGitHub APIから検証する。
- Copilot review後に追加された修正commitだけを理由に再reviewを要求しない。
- current Copilot-authored threadがすべてresolvedまたはoutdatedならreview完了とする。threadが作成されなかったreviewも，有効なreview submissionが確認できれば完了とする。

### Trigger and security boundary

- `pull_request_target`，PR上の`issue_comment`，`workflow_dispatch`，`schedule`で再評価する。
- `pull_request_review`はPR merge commit側のworkflow定義を実行し得るため，write permissionを持つgateのtriggerには使用しない。
- `pull_request_target`のsecurity boundaryを維持し，default branch上のworkflow定義だけで評価する。
- review submissionやthread resolve後の即時再評価が必要な場合は`workflow_dispatch`を使う。
- scheduled scanはopen PRをページネーションして再評価する。

### Merge policy

Codexは，小タスクに限り以下を満たす場合だけmergeしてよい。

- `docs/codex-runs/<run-id>/review_result.json`が`PASS`
- CIがpass
- `codex-review-gate`がpass
- ユーザー定性評価や重大な設計判断が不要

Cloud reviewをcommitごとのlintとして使わない。指摘が出た場合はactionable threadを一括取得し，まとめて修正し，merge-final self-checkを再実行してから全current actionable threadsをresolveする。修正不要と判断してresolveする場合は理由commentを残す。

通常commitのhookはADR 0007のlightweight checksのままとし，merge-final self-checkでdocs，schema，Issue状態，必要testを確認する。

重大判断やユーザー定性評価が必要な場合は，`review_result.json`を`FAIL`とし，必要な実行command，出力先，評価観点，未判断点を提示して止める。

## Consequences

Codex reviewを既定として維持しながら，Codex利用制限時もGitHub Copilot reviewでmerge gateを満たせる。

review後の修正commitごとの再reviewを不要とするため，review resourceと待ち時間を抑えられる。

review対象commitがPR履歴から消えた場合はreviewを再利用できないため，review済みPRへmainを取り込む場合は原則としてrebaseではなくmergeを使う。

`pull_request_review`による即時再評価は行わないため，review submissionまたはthread resolveからstatus更新まで，manual dispatchまたはscheduled scanを待つ場合がある。

物理モデル，dataset採用条件，Phase境界，ML学習条件，動画の自然さなどの判断は自動mergeしない。
