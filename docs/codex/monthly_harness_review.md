# Codex 月次ハーネス棚卸し

## 目的と境界

この手順は，Codex のハーネス，`AGENTS.md`，repository skill，
`.codex/config.toml`，および関連運用文書を月次で見直すためのものである．
ユーザーが相談スレッドを開始したときだけ実施し，定期自動実行や自動Issue作成は行わない．

通常の実装taskはこの文書を読む必要がない．月次棚卸し，モデル更新の提案，または
AGENTS / skill の大きな再編時だけ参照する．

## 入力と読み方

対象期間は前回の月次親Issueの作成日翌日から棚卸し日までとする．前回がなければ，
初回ベースラインとして現行状態を記録する．

次を根拠として使う．

1. `AGENTS.md`，`.agents/skills/*/SKILL.md`，`.codex/config.toml`，
   `docs/codex/` の正本運用文書．
2. 対象期間のGitHub Issue，PR，Git履歴，CI / review gate と
   `docs/codex-runs/*/review_result.json`．run recordは必ず
   `review_result.json` を先に読み，必要な場合だけ `work_log.md` を読む．
3. 当該repositoryのCodex task一覧と要約．詳細本文は，棚卸しの論点を裏付ける
   必要があるtaskだけに限定する．他projectのtask，私的会話，無関係な会話は読まない．
4. モデル，skills，Codex の運用提案には，その時点の公式OpenAI documentationを使う．

大きなログや生成物は読まない．既存のcompact summary，manifest，Issue / PR記録を優先する．

## 評価表

各観点につき，少なくとも一つの根拠を添えて次の欄を埋める．問題がなければ
「維持」と理由を記録する．

| 観点 | 良かった点 | 悪かった点・兆候 | 根拠 | 推奨action | 優先度 | 想定効果 | 実装対象 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 開発効率 | 着手からPRまでの反復，調査重複，テスト粒度，待機・手戻り |  |  |  | P0--P3 | 時間・反復の削減 | doc / config / skill / code |
| トークン・コンテキスト効率 | 常時指示量，文書読込，skill選択，ツール往復 |  |  |  | P0--P3 | context / 往復の削減 | AGENTS / doc / skill |
| 品質・再現性 | review result，テスト，manifest，PR review，実行証跡の整合 |  |  |  | P0--P3 | 失敗・再解析の削減 | contract / test / workflow |
| 安全性・統制 | 権限境界，cs10操作，破壊的操作，秘密情報，外部変更 |  |  |  | P0--P3 | 安全な自律性 | AGENTS / runbook / workflow |
| 保守性・知識衛生 | 重複，矛盾，古さ，正本の位置 |  |  |  | P0--P3 | 読み違いの削減 | AGENTS / doc / skill |
| 利用者体験・協働 | 進捗，確認依頼，ユーザー実行手順 |  |  |  | P0--P3 | 判断・操作の明確化 | prompt / workflow / runbook |
| モデル適合性 | 過剰制約，完遂性，理由付け，テスト過剰 |  |  |  | P0--P3 | 品質・cost・latency | config / AGENTS / skill |

実トークン数または料金を取得できる場合は記録する．取得不能な場合は，`AGENTS.md` と
読み込んだskillの行数，読んだ文書数，ツール往復数，重複した調査・checkを代理指標とし，
代理指標であることを明記する．推測した数値を記録してはならない．

## AGENTS と skill の判定

`AGENTS.md` には，全taskに必要な不変の安全境界，repository規約，最小のcontext routingだけを
置く．task固有の手順，Phase固有の詳細，長い例，完了済みtaskの履歴は，対象文書または
用途限定skillの正本へ置く．同じ規則は一つの正本だけに置き，他方からリンクする．

各skillは次を確認する．

* descriptionが短く，適用される作業を狭く明示していること．
* root `SKILL.md` が必要最小限のルーターであり，詳細はreferencesやscriptに遅延読込できること．
* 重複・矛盾する指示，または対象期間に再利用実績がないskillがないこと．
* model固有の回避策を，現在のモデルでも必要な根拠なしに保持していないこと．

削除・統合・移動は，実際のtask evidenceと参照検索を確認してから子Issueで実施する．

## モデル比較と移行判定

既定モデルは月次棚卸しだけでは変更しない．候補（原則として公式OpenAI documentationが
推奨する現行general-purpose model）を更新する場合は，採択後の独立子Issueで現行設定と比較する．

比較corpusは，対象期間の代表taskから次の4種を選ぶ．各taskは秘密情報を除いた再実行可能な
promptと期待結果を記録する．

1. 小規模な文書・workflow変更．
2. 限定的なコード変更と対象test．
3. Issue / PR 作成・更新を含む運用task．
4. Phase 2のcs10統制を含む，runtimeを開始しない計画task．

各候補を同じcorpusと同じ成功条件で評価し，成功率，手戻り，無用な確認停止，必要な検証の
欠落，ユーザー修正，取得できるtoken / cost / latency，および前節の代理指標を比較する．
重大な安全・再現性の退行がなく，現行より総合的に劣らないことを移行の最低条件とする．
移行子Issueには比較表，更新するconfig・AGENTS・skill，移行後の対象回帰checkを含める．

GPT-6系を候補にする場合は，skillsの短い適用条件とprogressive disclosure，必要な文書だけを
読むrouting，task相応のtesting，明示された完遂条件を特に再評価する．理由なく旧モデル向けの
強い逐次手順や一律test要求を移植しない．

## 月次Issue化

棚卸しの最後に，`execution:none` の月次親Issueを作る．Issue Formのruntime欄は
`no_runtime`，独立condition数とMac wall timeは`N/A`，Roadmap categoryは
`Project & Operations` とする．親Issueには次を記録する．

* 監査日，対象期間，現行モデルと設定，参照した根拠．
* 評価表，採択・保留・却下の判断，保留の再評価条件．
* 採択した子Issueへのリンクと，未採択案を実装しない理由．

採択した改善だけをGitHub sub-issueとして親Issueに追加する．子Issueには目的，変更対象，
受入条件，リスク，必要なexecution target，独立condition数，Mac wall time見積りを記載する．
単なる関連性ではblocking関係を作らない．

この手順そのものを導入する初回親Issueでは，現行モデル，`AGENTS.md`，詳細workflow，
project skill数，直近の代表taskをベースラインとして記録する．実装候補が採択されるまで
子Issueを作らない．
