# Issue263 D/J実装・予約準備

ユーザーのD/J承認と実行予約指示に従い、既存PR264 branchでopt-in実装。Dは同じSVD軸の純vector forceをゼロから作りone-ringへ返す。Jは固定prefix自然角をモデルに保存しengine/共通step-summary/preflightが共有。legacy/default profileは保持。

全条件のpreflight、短期/長期の4 worker dry-run、targeted38/full1108、production受動Hessian2幅、初期previewを確認。Jは幅0.762、currentは0.813、元Cは1.041 µm。初期配置差・currentの1.039pN反発とJの0を明示し、密束採択としない。

cs10読み取り接続でqueue unpaused・active0を確認。予約はfixed commit push/remote dry-run後に行い、実ID・開始状況はローカルreservation manifestとIssue/PRへ反映する。
