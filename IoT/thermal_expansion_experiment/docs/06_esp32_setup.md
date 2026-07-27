# 06. ESP32 開発環境セットアップ

更新日: 2026-07-27

## 6.1 なぜArduino UNOではなくESP32か

- GPIO数に余裕があり、MAX31856×3枚＋Digimatic入力を同時に配線しやすい。
- Digimatic信号処理と温度記録を同時に行うのに十分な処理性能・メモリがある。
- Arduino IDEで開発でき、USBシリアルでPCへデータを送信できる。
- ロジックは3.3V系で統一でき、MAX31856（3.3V対応品）と電圧を合わせやすい。

なお、ESP32はArduino UNOの完全な上位互換ではない（ピン配置・ADC特性・一部ライブラリの挙動が異なる）。既存の`IoT/arduino/sketches/thermocouple_basic`（Arduino UNO向けMAX31856サンプル）はロジックの参考にはなるが、そのまま流用はできない。

## 6.2 開発方式: PlatformIOとArduino IDEの両対応

本プロジェクトのファームウェアは `firmware/` 以下にPlatformIO形式（`platformio.ini` + `src/*.cpp,*.h`）で格納している。

- **PlatformIO（推奨）:** VSCode拡張のPlatformIO IDEを使い、`firmware/`フォルダを開いてビルド・書き込みするだけでよい。`lib_deps`に必要ライブラリが記載されているため、初回ビルド時に自動でダウンロードされる。
- **Arduino IDEでビルドする場合:** Arduino IDEは「スケッチフォルダ名と同名の.inoファイル」を必要とする。以下の手順で読み替える。
  1. `firmware/src/`の内容を、新しいスケッチフォルダ（例: `thermal_expansion_logger/`）にコピーする。
  2. `main.cpp` を `thermal_expansion_logger.ino` にリネームする（中身はそのままでよい。Arduino IDEは`.ino`をC++として扱う）。
  3. 他の`.h`/`.cpp`ファイルはそのままスケッチフォルダに置く（Arduino IDEはタブとして複数の.cpp/.hファイルを扱える）。
  4. ボードマネージャで「ESP32」を追加し、ライブラリをライブラリマネージャからインストールする（6.4節）。

## 6.3 ボード設定（Arduino IDE / PlatformIO共通の考え方）

| 項目 | 設定値 |
|---|---|
| ボード | ESP32 Dev Module（ESP32-DevKitC-32E） |
| Flash Size | 4MB（ボードに応じる） |
| Upload Speed | 921600（書き込み不安定なら460800や115200に下げる） |
| シリアルモニタ速度 | 115200 bps |

Arduino IDEの場合、事前に「ファイル > 環境設定 > 追加のボードマネージャのURL」にESP32用URL（`https://raw.githubusercontent.com/espressif/arduino-esp32/gh-pages/package_esp32_index.json` 等、Espressif公式ドキュメントの最新情報を確認）を登録し、ボードマネージャから「esp32 by Espressif Systems」をインストールする。

## 6.4 必要ライブラリ

| ライブラリ名 | 用途 | インストール方法 |
|---|---|---|
| Adafruit MAX31856 library | MAX31856読み取り（K型熱電対、CJ補償、フォルト取得） | Arduino IDE: ライブラリマネージャで検索・インストール / PlatformIO: `platformio.ini`の`lib_deps`に記載済み |
| Adafruit BusIO | 上記の依存ライブラリ | 同上（Adafruit系ライブラリインストール時に自動提案される） |
| SPI | ESP32標準搭載（追加インストール不要） | - |

Digimatic/SPC信号の読み取りは、汎用ライブラリに頼らず`firmware/src/digimatic_reader.cpp`で直接実装している（`05_digimatic_protocol.md`参照）。既存のEspDRO等のコードを参考にした場合は、参考にした部分とライセンス表記を`references/references.md`に記載すること。

## 6.5 初回動作確認（ステップ1〜3、`README.md`の段階的構築ガイドと対応）

**ステップ1: ESP32単体動作確認**
1. Arduino IDE/PlatformIOをセットアップする。
2. 空のスケッチ（`Serial.println("hello");`をloopで1秒ごとに出す程度）を書き込む。
3. シリアルモニタ(115200bps)に文字列が表示されることを確認する。

成功条件: 書き込みが通り、シリアルモニタに周期的に文字列が出る。
よくある失敗: ドライバ未インストールでCOMポートが認識されない → CP2102/CH340等のUSBシリアルドライバを確認する。書き込み時に「Connecting....」で止まる → ESP32のBOOTボタンを押しながら書き込み開始する機種がある。

**ステップ2: MAX31856を1枚接続**
1. `04_electrical_wiring.md`のピン割当でMAX31856 #1のみ配線する。
2. `firmware/src/max31856_manager.cpp`のロジックを流用し、1ch分だけ読み取る簡易スケッチで動作確認する。
3. 熱電対1本を接続し、室温付近の値が表示されることを確認する。

成功条件: readCelsius相当の値がNaNでなく、室温相当（例: 20〜30℃）の値になる。フォルトフラグが立たない。
よくある失敗: 配線ミス（特にMOSIの配線忘れ。MAX31856はMAX31855と異なりMOSIも必要）、CS配線忘れ、熱電対の極性逆（値が異常に低い/高い、または変化が逆）。

**ステップ3: MAX31856を3枚へ増設**
1. CS2, CS3を配線し、`NUM_THERMOCOUPLE_SENSORS`を3に設定する。
2. 3ch同時に値が表示されることを確認する。

成功条件: 3ch分の温度がすべて妥当な値で、CSVの列数が3温度分正しく出る。
よくある失敗: 全チャンネルのSCK/MISO/MOSIの結線忘れによる特定chのみ動作しない状態。CSピンの取り違え。

以降のステップ（インジケータ単体確認、Digimatic波形確認、統合、PC保存、治具評価、加熱予備試験、本実験）は、`README.md`の「段階的構築ガイド」および`05_digimatic_protocol.md`、`07_experiment_procedure.md`を参照。

## 6.6 シリアル出力の運用ルール

- 起動時ログ（初期化成功/失敗メッセージ等）と、CSVデータ行は明確に区別する。
- データ行の途中にデバッグ文字列を混在させない（CSV解析側でのパース事故を防ぐため）。デバッグ・診断メッセージは`[INFO]`/`[ERROR]`等のプレフィックスを付け、CSVデータ行とは別の行として出力する。
- 通信エラー（MAX31856読み取り失敗等）が起きてもプログラム全体を停止させない設計とする（`firmware/src/main.cpp`参照）。
