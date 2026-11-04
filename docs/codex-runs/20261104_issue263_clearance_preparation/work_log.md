# 非重複R/Dの修正・準備

新mode compact_clearance_rearは旧shapeを剛体並進し、hook35°でBF最短.204788um、body-hull下界.204498umを満たす。自然角/物性/力学は維持。共通pair selectorをruntimeとpreflightで共有し、新4条件contractで全active線分clearanceを必須化。旧root-on-edgeは予約前に拒否する回帰testを追加。

Targeted29/full1117 PASS。Macは限定32stepと検査だけ。geometry preview/preflight、2差分幅受動安定、同一geometry座標一致/旧J-DJ初期不変を確認。short/main dry-runは4workers/各thread1。固定commitと実予約IDはpush後のlocal reservation manifestへ記録する。
