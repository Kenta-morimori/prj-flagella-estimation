# Issue #263：M/K short予約#43

2026-11-11。固定e32c4ad033411c4bbdfb3fcc6583eab10b93e08eの実装をpush後、ユーザーが明示委任した接続・fixed-commit enqueue・既存dispatcher利用/必要時開始の範囲で実施した。1回目はcs10名前解決不可。FortiClient画面が空白で状態を断定せず、復旧依頼を一度提示。再照合でSSH復旧し、ユーザーからもVPN再接続の返答を受領した。再承認待ちは作っていない。

queue未pause・active予約なしと既存dispatcherを確認。中央checkoutのbranch/codeを変更せずfetchし、隔離sourceでdry-run。4条件の初期clearance PASS、入力SHA、実効4workers/thread1とNFS namespaceを確認。単一M条件のdt=0静的力でruntime互換性も確認した。新しいdispatcherは起動せず既存を利用した。

台帳へ先に1枠を割り当て、enqueueを1回だけ実行。#43を23:02:22 JSTにqueuedとして受領。23:02:43 JSTに既存dispatcherが開始。queueレコードのoutput_root/control_dirは完了後に埋まる実装なので、running中は当該worktreeのlaunch metadataからNFS job manifestを読み、固定commit・4ID・job/campaign SHA・running4/failed0を照合した。初回2500step checkpointは23:07時点では未取得。

対象はnf3/slots024・同一J配置のDJ/M/K/MK、各1τ=.04秒/10000step。job SHA256 cabf77fef5a73c45e5731f3d665be2267e1ac7dd167b88f4e9a0014414844314、campaign SHA256 587470b09539c69280c941b1659afd7a9df964c5ba4ab8a4601251873322f636。実行中にsimulation/queue/code/条件を変更していない。初期ETAは旧実績の約45分という予算であり、自条件checkpointで更新する。

ローカル証拠はoutputs/2026-11-11/230216/issue263_mk_short_reservation/。NFS rootは /net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs/2026-11-11/230243/parallel/queue-00043-d8090358__bf8ed5038049。operational log/credentialは同期せず、必要なmetadataだけを記録した。

完走時はfixed source /private/tmp/issue263-fixed-e32c4adで解析し、portable集合/件数/全SHA、run_summary先読、common QC/heatmap/固定replay・回転/束/推進/根元を代理レビューする。合意済み条件を満たせば1秒stageへ進む。予約予算1/5使用、1秒はまだ未予約。Goalはactiveのままで、今回の手続PASSを科学的成功としない。
