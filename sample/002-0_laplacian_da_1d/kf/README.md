# kf/ — Kalman Filter によるデータ同化

## 概要

このフォルダは `oi/` の次のステップとして、**Kalman Filter (KF)** を実装する予定です。

## OI との違い

| 項目 | OI（`oi/`） | KF（このフォルダ） |
|---|---|---|
| 背景誤差共分散 | 固定 B | 時変 P（毎ステップ伝播） |
| 状態遷移行列 A_d | 不要 | 必要（OF の線形化近似） |
| ゲイン K | 毎回同じ | P から毎ステップ更新 |
| 実装難度 | 低 | やや高め |

## 実装予定

1. `A_d`（離散状態遷移行列）を Python 側で Euler 近似
2. 予測ステップ: `P = A_d P A_d^T + Q`
3. 更新ステップ: `K = P H^T (H P H^T + R)^{-1}`、`x_kf = x_pred + K(y - Hx_pred)`、`P = (I - KH) P`

## フォルダ構成（実装後）

```
kf/
├── case/           # KF 用 OF 作業ディレクトリ（case_base/ からコピー）
├── da_main.py
├── run.sh
└── results/
```

## 参考

- `../docs/theory.md` — KF の理論式と OpenFOAM との接続
- `../docs/da_taxonomy_and_refs.md` — KF の位置づけと参考文献
