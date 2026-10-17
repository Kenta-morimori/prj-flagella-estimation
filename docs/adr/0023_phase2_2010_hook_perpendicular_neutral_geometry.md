# ADR 0023: 2010系候補のhook直交・中立初期配置

## Context

#257で導入した`initial_hook_force_neutral`は、hookとべん毛根元接線を直交させ、初期hook力とビーズ重なりを避ける。しかし菌体長軸との直交は保証せず、従来の側面付着配置とは異なる。#255では、両2010系モデルでhookを菌体長軸にも垂直にしたうえでmotor反作用を比較する。

## Decision

既存の`initial_hook_force_neutral`挙動は残し、`initial_hook_body_axis_perpendicular`を既定OFFのopt-in設定として追加する。両方をONにした場合、付着ビーズと第1べん毛ビーズを固定し、べん毛全体を菌体長軸まわりに剛体回転する。根元接線とhookを直交させる二候補から、他ビーズとの最小距離が最大の組合せを選び、非重複配置が得られなければ構築を失敗させる。

#255の候補campaignではbody–flagella spring-segment排除をOFFにし、flagella–flagella反発を維持する。既存の2010 projectと2010 hex profileの既定値は変更しない。

## Consequences

- 初期hookは長軸・根元接線に直交し、hook長と内部べん毛形状を保つ。seedで生成した初期位相は回転前形状を規定し、回転後の物理位相と同一ではない。
- body–flagella貫通を許すのは#255候補runに限る。初期ビーズ重なりは許さない。
- t=0の幾何学的成立は動的安定性やmotor反作用の物理的妥当性を証明しない。1τのQC/replay後に2s条件を判断する。
