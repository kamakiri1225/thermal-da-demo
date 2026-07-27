/*
 * config.h
 * ============================================================
 * SUS304丸棒 熱膨張実験ロガー 設定ファイル
 *
 * ここに書かれたピン割当・定数は docs/04_electrical_wiring.md,
 * docs/05_digimatic_protocol.md の内容と対応している。
 * 実機確認・信号波形確認の結果に応じて随時更新すること。
 * ============================================================
 */
#pragma once

#include <Arduino.h>

// ============================================================
// 測定周期
// ============================================================
// 1秒周期を基本とする（要件7節）。変更可能。
static const unsigned long MEASURE_INTERVAL_MS = 1000;

// ============================================================
// MAX31856 (K型熱電対アンプ) 設定
// ============================================================
// 実装しているセンサ数（丸棒T1-T3で3。基台温度等を追加する場合は4まで拡張可）
static const int NUM_THERMOCOUPLE_SENSORS = 3;

// 共通SPIバス (docs/04_electrical_wiring.md)
static const int PIN_SPI_SCK  = 18;
static const int PIN_SPI_MISO = 19;
static const int PIN_SPI_MOSI = 23;

// CSピン（個体ごと）。4番目は基台温度等の将来拡張用（未使用ならそのままでよい）。
static const int MAX31856_CS_PINS[4] = {25, 26, 27, 33};

// センサ名（CSVヘッダ・診断メッセージ用）
static const char* THERMOCOUPLE_LABELS[4] = {
  "T1_heated", "T2_center", "T3_free_end", "T4_frame_opt"
};

// ============================================================
// Digimatic/SPC (変位計) 設定
// ============================================================
// 入力専用ピン。内部プルアップ/プルダウンなし（docs/05_digimatic_protocol.md 要確認）。
static const int PIN_SPC_DATA  = 34;
static const int PIN_SPC_CLOCK = 35;

// 一部機種で必要とされる要求(REQ)信号線。未使用ならUSE_SPC_REQ_LINEを0のままにする。
// 要信号波形確認: 実機でREQが必要と判明した場合のみ1にしてピンを配線する。
#define USE_SPC_REQ_LINE 0
#if USE_SPC_REQ_LINE
static const int PIN_SPC_REQ = 32;  // 仮案。REQの極性・要否は要確認。
#endif

// 1フレームとみなす最大ビット数（要実機確認、5.5節の一般的な構成を仮定した初期値）
static const int SPC_MAX_FRAME_BITS = 64;

// このミリ秒間クロックの変化がなければ1フレーム受信完了とみなす（要調整）
static const unsigned long SPC_FRAME_TIMEOUT_MS = 50;

// このミリ秒間、有効なフレームを受信できなければ spc_valid=0 とする（通信断とみなす）
static const unsigned long SPC_STALE_TIMEOUT_MS = 5000;

// RAWビット列をデバッグ出力するか（第一段階のブリングアップで1にする。docs/05_digimatic_protocol.md 5.6節）
#define SPC_RAW_DEBUG_OUTPUT 1

// インジケータの表示単位。trueならinch表示なのでmmへ変換する（要実機確認）。
static const bool INDICATOR_IS_INCH_MODE = false;
static const double INCH_TO_MM = 25.4;

// ---- BCDデコードのビット割当（仮の暫定値。実測波形で必ず調整すること） ----
// docs/05_digimatic_protocol.md 5.5-5.6節の通り、公開情報でも細部の記述が
// 一致しないため、ここでは「よく説明される構成」を暫定値として置いている。
// SPC_RAW_DEBUG_OUTPUT の生ビット列とインジケータの表示値を突き合わせて
// 実測で確定させること（要実機確認）。
static const int SPC_EXPECTED_FRAME_BITS = 52;  // これと一致した場合のみデコードを信頼する
static const bool SPC_BIT_ORDER_MSB_FIRST = false; // false=LSBファースト(仮), true=MSBファースト
static const int SPC_SIGN_BIT_INDEX = 0;         // 符号ビットの位置（仮）
static const int SPC_DATA_START_BIT = 4;         // BCDデータ開始ビット（仮）
static const int SPC_DATA_NUM_DIGITS = 6;         // BCD桁数（仮）
static const int SPC_DECIMAL_POINT_DIGITS_FROM_RIGHT = 3; // 小数点位置（仮、下3桁が小数部）

// ============================================================
// ゼロ点校正 (17節)
// ============================================================
static const unsigned long ZERO_CAL_DURATION_S = 30;   // 測定開始前に何秒間ゼロ点を取るか
static const unsigned long ZERO_CAL_INTERVAL_MS = 1000; // ゼロ点取得中のサンプリング周期
// ゼロ点取得中の標準偏差がこれを超えたら「不安定」として警告する（要確認、仮の目安値）
static const double ZERO_CAL_STDDEV_WARN_MM = 0.003;

// ============================================================
// CSV / シリアル出力
// ============================================================
// 冷接点温度(CJ)を追加列として出力するか
static const bool ENABLE_COLD_JUNCTION_LOG = true;

static const unsigned long BAUD_RATE = 115200;
