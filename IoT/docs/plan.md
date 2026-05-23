# IoT 温度測定 × OpenFOAM データ同化 プロジェクト計画

更新日: 2026-05-08

---

## プロジェクト概要

K型熱電対で固体表面温度を実測し、最終的にOpenFOAMの熱流体計算結果とデータ同化を行う。

---

## フェーズ構成

```
Phase 1  Arduino UNO でシリアルモニタ温度表示    ← 現在
Phase 2  Python で CSV 保存・リアルタイムグラフ
Phase 3  ESP32 で Wi-Fi Web 表示
Phase 4  OpenFOAM 計算結果とのデータ同化
```

---

## Phase 1: Arduino UNO + シリアルモニタ表示

### ゴール
Arduino UNO でK型熱電対の温度を取得し、シリアルモニタに1秒ごとに表示する。

### タスクリスト

- [ ] MAX31856搭載モジュールを購入する
- [ ] ジャンパーワイヤ、カプトンテープ等を購入する
- [ ] Arduino UNO と MAX31856 を SPI で配線する
- [ ] K型熱電対を MAX31856 に接続する
- [ ] Arduino IDE に `Adafruit MAX31856` ライブラリをインストールする
- [ ] `arduino/sketches/thermocouple_basic/` にスケッチを保存する
- [ ] Arduino UNO に書き込み、シリアルモニタで室温が表示されることを確認する
- [ ] 熱電対の先端を指でつまみ、温度が上がることを確認する（極性チェック）
- [ ] 固体に熱電対を固定し、表面温度を測定する

### 成功基準
- シリアルモニタに `熱電対温度: XX.XX ℃` が1秒ごとに表示される
- 熱電対先端を温めると値が上昇する
- Fault が出ない

---

## Phase 2: Python CSV 保存・リアルタイムグラフ

### ゴール
ArduinoのシリアルデータをPythonで受信し、CSVに保存しながらリアルタイムグラフ表示する。

### タスクリスト

- [ ] Python 環境を準備する（pyserial, matplotlib, pandas）
- [ ] `python/serial_logger.py` を作成する（シリアル受信 + CSV保存）
- [ ] `python/realtime_plot.py` を作成する（リアルタイムグラフ）
- [ ] `data/raw/` に日付付きCSVが保存されることを確認する
- [ ] 測定データのサンプルを `data/raw/` に保存する

### 成功基準
- `data/raw/YYYYMMDD_HHMMSS.csv` にタイムスタンプ付き温度データが保存される
- グラフで温度変化がリアルタイムに見える

---

## Phase 3: ESP32 Wi-Fi Web 表示

### ゴール
ESP32 を使って、ブラウザから温度をリアルタイム確認できるようにする。

### タスクリスト

- [ ] ESP32 開発ボードを購入・入手する
- [ ] Arduino IDE に ESP32 ボードサポートを追加する
- [ ] ESP32 + MAX31856 の配線を行う（SPI）
- [ ] `arduino/sketches/esp32_webserver/` にスケッチを保存する
- [ ] Wi-Fi に接続し、簡易 Web サーバで温度を JSON 配信するスケッチを作成する
- [ ] ブラウザからアクセスして温度が見えることを確認する

### 検討事項
- Arduino UNO R4 WiFi も候補（Wi-Fi付きArduino）
- データをクラウド（MQTT / InfluxDB / Grafana）に流すことも将来検討

---

## Phase 4: OpenFOAM 計算結果とのデータ同化

### ゴール
IoT で計測した固体表面温度と、OpenFOAM の熱流体計算結果を組み合わせてデータ同化を行い、精度の高い熱流体状態を推定する。

### 参照 OpenFOAM ケース
`D:\work\002_CAE\openfoam\20260505_datadoka\`（本プロジェクトの兄弟ディレクトリ）

### タスクリスト

- [ ] OpenFOAM 側の計算ケースを整理する（境界条件、メッシュ）
- [ ] IoT 計測点とOpenFOAMメッシュの対応を確認する
- [ ] データ同化手法を選定する（EnKF / 4DVar / 最小二乗法等）
- [ ] `python/data_assimilation/` にデータ同化スクリプトを作成する
- [ ] 同化結果を可視化する

### 検討事項
- 計測点数が少ないため、表面温度1点 or 数点からの推定が現実的
- まずはシンプルな最小二乗フィッティングから始める
- OpenFOAMの境界条件を計測値で更新するアプローチが現実的

---

## フォルダ構成

```
IoT/
├── docs/
│   ├── memo.md        # 機材・配線・トラブルシュートメモ
│   └── plan.md        # 本ファイル（プロジェクト計画）
├── arduino/
│   └── sketches/
│       ├── thermocouple_basic/   # Phase 1: 基本スケッチ
│       └── esp32_webserver/      # Phase 3: ESP32 Webサーバ
├── python/
│   ├── serial_logger.py          # Phase 2: シリアル受信 + CSV保存
│   ├── realtime_plot.py          # Phase 2: リアルタイムグラフ
│   └── data_assimilation/        # Phase 4: データ同化
├── data/
│   ├── raw/                      # 生CSVデータ
│   └── processed/                # 処理済みデータ
└── README.md
```

---

## 現在の状態（2026-05-08）

- Phase 1 着手前（機材購入待ち）
- MAX31856モジュール・ジャンパーワイヤ等の購入が次のアクション
