#include "zero_calibration.h"
#include "config.h"
#include "digimatic_reader.h"

#include <math.h>

ZeroCalibration g_zeroCalibration;

ZeroCalibrationResult ZeroCalibration::run() {
  ZeroCalibrationResult result;
  result.ok = false;
  result.meanMm = NAN;
  result.stdDevMm = NAN;
  result.sampleCount = 0;

  const unsigned long totalMs = ZERO_CAL_DURATION_S * 1000UL;
  const unsigned long stepMs = ZERO_CAL_INTERVAL_MS;
  const int maxSamples = (int)(totalMs / stepMs) + 1;

  double sum = 0.0;
  double sumSq = 0.0;
  int n = 0;

  Serial.printf("[INFO] ゼロ点取得を開始します（%lu秒間、%lu msごと）。丸棒・治具を動かさないでください。\n",
                (unsigned long)ZERO_CAL_DURATION_S, stepMs);

  unsigned long startMs = millis();
  unsigned long nextSampleMs = startMs;

  while (millis() - startMs < totalMs) {
    // 待機中もSPCフレーム処理を継続する
    g_digimatic.update();

    if (millis() >= nextSampleMs) {
      double v = g_digimatic.getDisplacementMm();
      if (!isnan(v) && g_digimatic.isValid()) {
        sum += v;
        sumSq += v * v;
        n++;
      }
      nextSampleMs += stepMs;
    }
    delay(10);
  }

  result.sampleCount = n;

  if (n < 3) {
    Serial.println("[ERROR] ゼロ点取得: 有効サンプルが不足しています(Digimatic未接続/未デコードの可能性)。"
                   "disp_mmはNaNのまま記録されます。");
    return result;
  }

  double mean = sum / n;
  double variance = (sumSq / n) - (mean * mean);
  if (variance < 0) variance = 0;  // 数値誤差対策
  double stddev = sqrt(variance);

  result.meanMm = mean;
  result.stdDevMm = stddev;
  result.ok = (stddev <= ZERO_CAL_STDDEV_WARN_MM);

  Serial.printf("[INFO] ゼロ点取得完了: mean=%.4f mm, stddev=%.4f mm, n=%d\n", mean, stddev, n);
  if (!result.ok) {
    Serial.printf("[WARN] ゼロ点の標準偏差(%.4fmm)が閾値(%.4fmm)を超えています。"
                  "治具の振動・スティックスリップ・接触不良の可能性があります。\n",
                  stddev, ZERO_CAL_STDDEV_WARN_MM);
  }

  return result;
}
