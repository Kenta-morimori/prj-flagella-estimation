# Issue #263：承認済みM/K解析力・比較準備

2026-11-11。前Goal turnのロードマップ・明示委任・Goal開始は進捗に分類。実行台帳は使用0で、新規cs10操作をまだしていないことを確認した。専用PR264 worktreeを使いprimary checkoutは変更していない。

Mは相対frame位相の仮想仕事から解析pullbackを導出。proper fitのpolar微分とPCAの固有vector微分を含め、変形時に最小ノルムwrenchへ代用しない。独立のSVD/rigid-fit/quaternion位相oracleで非剛体仮想仕事・全座標差分を確認し、objectivity/units/gauge/1回転仕事・逆軸/縮退拒否をtest。差分力をproductionに入れていない。

Kは固定近位自然角を保ち、9 tripletの角度接線応答だけを正規化する係数配列として実装。保存力のエネルギー差分を確認。遠位、ばね、hook、反発、torsion、T/dt/粘性を変えていない。新modeはopt-in、既定profileは維持。

共通診断にdrive contract version、指定vector/軸投影、入力/出力角/gain/conditionとχ仕事を追加。旧D専用全vector誤差をMでは未適用とし、NaN/nullが物理有限性gateを変えないことを共通QCで確認。step仕事の時刻はpost-step、力/torqueはpre-stepとしてsummary極値時刻も区別した。

基準59890c0のsrc/confを隔離抽出し、D+Jの初期/64step終端と13保存状態のmotor/total力を新sourceと比較。各配列がbit一致。比較のため旧prototypeを読み込んだ際に更新された以前のstatic JSONは、元のwall_sを戻してmanifestの元SHAと完全一致を確認した。以降のreproducerは旧outputを書き換えない。

既存run_summary先読後にportable37件の集合/件数/全SHAを再確認。新sourceの13保存状態で独立方向仕事の差は最大~1.06e−8(T/µm)、couple誤差は丸め誤差。Macの64step費用は.108〜.122s/step（各物理.000256秒のみ）。短期・1秒profile/jobを作り共通preflightと実効4worker/thread1・分離outputのdry-run、同一初期形状のpreviewを確認した。

targeted142 PASS、full1137 PASS(293秒、既存layout warning1)。その後のpost-step時刻/manifest整理は関連83と最終phase/DJ28で確認。最初のtest harnessはbody gateへ正規body recorderを渡していなかったため修正し、invalid override testもnested APIに合わせた。物理QC閾値は変更していない。

1秒は既存のfloat/ceilで250001step。contractへ明示し、1.000004sのうち比較は1.0秒までを基本とする。shortは10000step。コード/設定/文書をcommit/pushしてfixed commitを確定後、許可済み接続・queue確認・remote dry-run・short予約へ進む。userの代理レビューは許可済みだが、科学的採択/merge・新構造への拡張は含まない。
