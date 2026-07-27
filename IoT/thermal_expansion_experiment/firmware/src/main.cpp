/*
 * main.cpp
 * ============================================================
 * SUS304丸棒 熱膨張実験ロガー (ESP32)
 *
 * 測定項目:
 *   - 丸棒温度3点 (K型熱電対 + MAX31856 x3)
 *   - 軸方向変位1点 (ミツトヨ デジマチックインジケータ, Digimatic/SPC)
 *   - 経過時間
 *   - 各センサの正常/異常
 *
 * 出力: USBシリアル(115200bps)へCSV1行/測定周期。
 *       PC側は python/serial_logger.py で受信・保存する。
 *
 * 詳細は docs/ 以下、特に
 *   docs/04_electrical_wiring.md (配線)
 *   docs/05_digimatic_protocol.md (Digimatic信号、要実機確認多数)
 *   docs/06_esp32_setup.md (ビルド・書き込み手順)
 * を参照。
 * ============================================================
 */
#include <Arduino.h>
#include <SPI.h>

#include "config.h"
#include "max31856_manager.h"
#include "digimatic_reader.h"
#include "zero_calibration.h"
#include "csv_output.h"

static unsigned long s_startMillis = 0;
static unsigned long s_lastMeasureMs = 0;
static double s_zeroMeanMm = NAN;

static void measureAndLog() {
  unsigned long elapsedMs = millis() - s_startMillis;  // unsigned演算によりmillis()オーバーフローも安全に扱える

  ThermocoupleReading readings[3];
  for (int i = 0; i < NUM_THERMOCOUPLE_SENSORS && i < 3; i++) {
    readings[i] = g_thermocouples.read(i);
  }
  for (int i = NUM_THERMOCOUPLE_SENSORS; i < 3; i++) {
    readings[i].temperatureC = NAN;
    readings[i].coldJunctionC = NAN;
    readings[i].faulted = true;  // 未使用チャンネルは「無効」として扱う
    readings[i].faultRaw = 0;
  }

  g_digimatic.update();
  double indicatorMm = g_digimatic.getDisplacementMm();
  bool spcValid = g_digimatic.isValid();

  double dispMm = NAN;
  if (spcValid && !isnan(indicatorMm) && !isnan(s_zeroMeanMm)) {
    dispMm = indicatorMm - s_zeroMeanMm;
  }

  csvPrintRow(elapsedMs, readings, indicatorMm, dispMm, spcValid);
}

void setup() {
  Serial.begin(BAUD_RATE);
  delay(1500);
  Serial.println();
  Serial.println("=== SUS304丸棒 熱膨張実験ロガー 起動 ===");
  Serial.printf("測定周期: %.1f 秒 / 熱電対点数: %d\n",
                MEASURE_INTERVAL_MS / 1000.0, NUM_THERMOCOUPLE_SENSORS);

  SPI.begin(PIN_SPI_SCK, PIN_SPI_MISO, PIN_SPI_MOSI, -1);

  g_thermocouples.begin(NUM_THERMOCOUPLE_SENSORS);
  g_digimatic.begin();

  // ゼロ点取得（17節: 1点でなく一定期間の平均を使う）
  ZeroCalibrationResult zc = g_zeroCalibration.run();
  s_zeroMeanMm = zc.meanMm;  // 失敗時はNaNのままdisp_mmはNaN表示になる

  Serial.printf("# zero_cal_ok=%d zero_cal_mean_mm=%.4f zero_cal_stddev_mm=%.4f zero_cal_n=%d\n",
                zc.ok ? 1 : 0, zc.meanMm, zc.stdDevMm, zc.sampleCount);
  if (!zc.ok) {
    Serial.println("[WARN] ゼロ点が不安定なまま計測を開始します。disp_mmの信頼性に注意してください。");
  }

  csvPrintHeader();

  s_startMillis = millis();
  s_lastMeasureMs = 0;
  Serial.println("---- 測定開始 ----");
}

void loop() {
  unsigned long now = millis();
  if (now - s_lastMeasureMs >= MEASURE_INTERVAL_MS || s_lastMeasureMs == 0) {
    s_lastMeasureMs = now;
    measureAndLog();
  } else {
    // 測定周期の合間もSPCフレーム処理を継続する(取りこぼし低減)
    g_digimatic.update();
  }
}
