# SUS304丸棒 熱膨張実験システム

更新日: 2026-07-27
状態: **設計段階**。多くの電気的・機械的詳細が未確定であり、各所に「要確認」「要実機確認」「要信号波形確認」「要メーカー仕様確認」と明記している。実機組立・信号確認を進めながら本ドキュメント一式を更新していくことを前提とする。

## 目的

SUS304丸棒をシリコンラバーヒーターで局所加熱し、以下を自動計測する。

- 丸棒の軸方向熱膨張量（自由端変位、ミツトヨ デジマチックインジケータ）
- 丸棒温度3点（K型熱電対 + MAX31856 ×3）
- 経過時間
- 各センサの正常/異常状態

得られた実験結果を、理論式 `ΔL = αLΔT` および `ΔL = ∫α(T)[T(x,t)-T0(x)]dx`、FEM/熱構造解析結果と比較し、OpenCAE関連の発表・技術資料に使える実験システムを構築する。

本プロジェクトは既存の `IoT/` フォルダ配下に独立サブプロジェクトとして構築している。既存の `IoT/arduino/sketches/esp32_thermal_logger`（住居内温湿度ロガー、BME280+microSD）とはセンサ構成・目的が異なるため独立させている。

## 全体構成（概要）

```
K型熱電対 T1 → MAX31856 #1 ─┐
K型熱電対 T2 → MAX31856 #2 ─┼→ ESP32 → USBシリアル → PC → CSV保存
K型熱電対 T3 → MAX31856 #3 ─┘
ミツトヨ デジマチックインジケータ → Digimatic/SPC → ESP32(同上)
```

ESP32が温度3点・変位1点・経過時間を1つのCSV行に集約してPCへ送る。センサが個別にCSVを送るわけではない。詳細は `docs/01_system_overview.md`、`diagrams/system_block_diagram.md`。

## フォルダ構成

```
thermal_expansion_experiment/
├─ README.md                      (本ファイル)
├─ docs/                          技術資料一式(01-11)
├─ firmware/                      ESP32用ファームウェア(PlatformIO形式)
│  ├─ platformio.ini
│  └─ src/                        main.cpp, config.h, 各managerクラス
├─ python/                        PC側: 記録・解析・比較スクリプト
├─ config/                        実験条件・材料物性のJSON設定
├─ data/
│  ├─ raw/                        生CSV(serial_logger.py出力先)
│  ├─ processed/                  解析結果(png, summary.json等)
│  └─ sample/                     サンプルCSV(合成データ、動作確認用)
├─ diagrams/                      ブロック図・配線図・治具図(Mermaid+ASCII)
└─ references/                    参考文献・参考実装の出典
```

## ドキュメント一覧

| # | ファイル | 内容 |
|---|---|---|
| 01 | `docs/01_system_overview.md` | システム概要、測定項目、比較対象 |
| 02 | `docs/02_parts_list.md` | 部品表・価格目安・購入時注意 |
| 03 | `docs/03_mechanical_fixture.md` | 治具設計(支持方式、たわみ計算、摩擦確認) |
| 04 | `docs/04_electrical_wiring.md` | ピン割当、AC100V部の隔離 |
| 05 | `docs/05_digimatic_protocol.md` | Digimatic/SPC信号の扱い方と安全策(要実機確認多数) |
| 06 | `docs/06_esp32_setup.md` | 開発環境構築、ライブラリ、書き込み手順 |
| 07 | `docs/07_experiment_procedure.md` | 準備〜無加熱試験〜加熱試験〜繰返し試験 |
| 08 | `docs/08_calibration.md` | インジケータ・熱電対・治具の校正 |
| 09 | `docs/09_error_budget.md` | 誤差要因の定量化 |
| 10 | `docs/10_safety.md` | 安全上の注意 |
| 11 | `docs/11_fem_comparison.md` | 理論式・FEM比較の考え方 |

## CSV出力形式

```
elapsed_ms,T1_C,T2_C,T3_C,indicator_mm,disp_mm,spc_valid,fault1,fault2,fault3,cold_junction_1_C,cold_junction_2_C,cold_junction_3_C
0,23.812,23.750,23.781,5.0002,0.0000,1,0,0,0,23.9,24.0,23.8
1000,24.125,23.906,23.844,5.0012,0.0010,1,0,0,0,23.9,24.0,23.8
```

- `#`で始まる行はメタデータ・コメント（起動時のゼロ点校正結果等）で、CSVパーサ側は無視できる（`pandas.read_csv(..., comment="#")`）。
- 冷接点温度列(`cold_junction_*_C`)は`firmware/src/config.h`の`ENABLE_COLD_JUNCTION_LOG`で有効/無効を切替可能。

サンプルデータ: `data/sample/sample_log.csv`（40分・2400行の**合成データ**。実測値ではなく動作確認用）。

## クイックスタート（PC側だけを試す場合）

```bash
cd python
pip install -r requirements.txt

# サンプルCSVで解析スクリプトを試す
python3 plot_results.py ../data/sample/sample_log.csv
python3 analyze_experiment.py ../data/sample/sample_log.csv

# 実機接続時
python3 serial_logger.py --port COM5   # Windowsの例
python3 serial_logger.py --port /dev/ttyUSB0   # Linuxの例
```

## 段階的構築ガイド（初心者向け、25節対応）

完全初心者が一度に全部作るのは難しいため、以下の順で進めることを推奨する。各ステップの詳細（成功条件・よくある失敗と対処）は該当ドキュメントを参照。

| ステップ | 内容 | 詳細ドキュメント |
|---|---|---|
| 1 | ESP32単体動作確認(Arduino IDE/PlatformIOセットアップ、シリアル出力) | `docs/06_esp32_setup.md` 6.5節 |
| 2 | MAX31856を1枚接続、K型熱電対1本を読む | `docs/06_esp32_setup.md` 6.5節 |
| 3 | MAX31856を3枚へ増やし、温度3点をCSV出力 | `docs/06_esp32_setup.md` 6.5節 |
| 4 | インジケータ単体動作確認、表示値の繰返し確認 | `docs/08_calibration.md` 8.1節 |
| 5 | Digimaticケーブル信号確認(ロジックアナライザで波形確認) | `docs/05_digimatic_protocol.md` 5.4節 |
| 6 | ESP32で変位のみ取得(RAWキャプチャモードで確認) | `docs/05_digimatic_protocol.md` 5.6節 |
| 7 | 温度3点と変位を統合 | `firmware/src/main.cpp` |
| 8 | PCでCSV保存 | `python/serial_logger.py` |
| 9 | 治具を無加熱で評価 | `docs/07_experiment_procedure.md` 7.2節、`docs/03_mechanical_fixture.md` 3.6節 |
| 10 | ヒーターを付けて低温予備試験 | `docs/10_safety.md`、`docs/07_experiment_procedure.md` |
| 11 | 本実験 | `docs/07_experiment_procedure.md` 7.3-7.4節 |

各ステップの成功条件・よくある失敗は`docs/06_esp32_setup.md`（ステップ1-3）、`docs/05_digimatic_protocol.md`（ステップ4-6）、`docs/07_experiment_procedure.md`（ステップ9-11）に記載している。

## 必要なArduino/PlatformIOライブラリ

- Adafruit MAX31856 library
- Adafruit BusIO（上記の依存ライブラリ）
- SPI（ESP32標準搭載）

インストール手順は`docs/06_esp32_setup.md`6.4節を参照。Digimatic/SPC受信は外部ライブラリを使わず自前実装している。

## 未確定事項（サマリ）

- インジケータ型番(543-790B-10)・ケーブル型番(905338)の最終適合
- Digimatic/SPC信号の電気仕様、保護回路の要否、正確なビットフレーム構成
- ヒーター容量・目標最高温度・温調器/SSR型番
- 熱電対固定方法・支持部材質の最終選定
- FEM側の境界条件

詳細は各`docs/*.md`内の該当箇所、および`docs/01_system_overview.md`1.8節を参照。これらは「仮案」であり、実機確認後に更新すること。

## OpenFOAM/OpenCAE活用との接続

本システムで得られる実測データ(温度3点・変位1点の時系列)は、以下のような形でOpenCAE関連の取り組みに活用できる見込み。

- 単純なFEM/理論式との比較による、実験誤差要因の定量的な理解(`docs/09_error_budget.md`)
- 本リポジトリ内の他のデータ同化の取り組み(`study/practice/`)と同様に、実測温度から変位を推定する、あるいは逆に変位から熱源・境界条件を推定するデータ同化の題材としての活用（将来課題）
