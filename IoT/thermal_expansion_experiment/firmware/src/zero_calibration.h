/*
 * zero_calibration.h
 * ============================================================
 * 変位の初期値を「測定開始前の一定期間の平均値」として取得する(17節)。
 *   u(t) = x(t) - mean(x0..xN)
 * 単純に1点をゼロとしない。
 * ============================================================
 */
#pragma once

#include <Arduino.h>

struct ZeroCalibrationResult {
  bool   ok;             // サンプルが十分取れ、安定していると判断できたか
  double meanMm;          // 初期値平均 [mm]
  double stdDevMm;        // 初期値の標準偏差 [mm]
  int    sampleCount;      // 有効サンプル数
};

class ZeroCalibration {
public:
  // ZERO_CAL_DURATION_S 秒間サンプリングし、平均・標準偏差を求める。
  // 呼び出し中も g_digimatic.update() を回してSPCフレーム処理を継続する。
  ZeroCalibrationResult run();
};

extern ZeroCalibration g_zeroCalibration;
