/*
 * max31856_manager.h
 * ============================================================
 * MAX31856 (K型熱電対アンプ) 複数個の初期化・読み取りをまとめる。
 * 共通SPIバス + 個別CSピン構成 (docs/04_electrical_wiring.md)。
 * ============================================================
 */
#pragma once

#include <Arduino.h>

struct ThermocoupleReading {
  double temperatureC;   // 熱電対温度。読み取り失敗時はNAN
  double coldJunctionC;  // 冷接点温度。読み取り失敗時はNAN
  bool   faulted;        // 何らかのフォルトが立っているか
  uint8_t faultRaw;      // フォルトビットの生値（診断用）
};

class Max31856Manager {
public:
  // numSensors: 使用するセンサ数 (1-4)
  void begin(int numSensors);

  // idx番目のセンサを読み取る（0始まり）
  ThermocoupleReading read(int idx);

  // idx番目のセンサが初期化に成功しているか
  bool isReady(int idx) const;

private:
  int _numSensors = 0;
};

extern Max31856Manager g_thermocouples;
