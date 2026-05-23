# IoT 温度測定プロジェクト

K型熱電対による固体表面温度測定。最終的にOpenFOAM計算結果とのデータ同化を目指す。

## フォルダ構成

```
IoT/
├── docs/
│   ├── memo.md        # 機材・配線・トラブルシュートメモ
│   └── plan.md        # プロジェクト計画（フェーズ別タスクリスト）
├── arduino/
│   └── sketches/
│       └── thermocouple_basic/   # Phase 1: Arduino UNO 基本スケッチ
├── python/                        # Phase 2以降: CSV保存・グラフ・データ同化
├── data/
│   ├── raw/                       # 生CSVデータ
│   └── processed/                 # 処理済みデータ
└── README.md
```

## フェーズ概要

| Phase | 内容 | 状態 |
|---|---|---|
| 1 | Arduino UNO + MAX31856 + K型熱電対 → シリアルモニタ表示 | 準備中 |
| 2 | Python で CSV保存・リアルタイムグラフ | 未着手 |
| 3 | ESP32 で Wi-Fi Web 表示 | 未着手 |
| 4 | OpenFOAM 計算結果とのデータ同化 | 未着手 |

詳細は `docs/plan.md` を参照。
