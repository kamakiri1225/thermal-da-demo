#include "max31856_manager.h"
#include "config.h"

#include <SPI.h>
#include <Adafruit_MAX31856.h>

Max31856Manager g_thermocouples;

// ハードウェアSPI(共有バス)使用。CSのみ個別。
static Adafruit_MAX31856* s_sensors[4] = {nullptr, nullptr, nullptr, nullptr};
static bool s_ready[4] = {false, false, false, false};

void Max31856Manager::begin(int numSensors) {
  _numSensors = numSensors;
  if (_numSensors > 4) _numSensors = 4;

  for (int i = 0; i < _numSensors; i++) {
    // CS指定のみのコンストラクタ = ハードウェアSPI(共有&SPIバス)を使用する。
    // mosi/miso/sckを明示するコンストラクタはこのライブラリではソフトウェアSPIに
    // 切り替わってしまうため使わない。ハードウェアSPIのピン自体は
    // main.cpp の SPI.begin(SCK, MISO, MOSI, -1) でGPIO18/19/23へ remap 済み。
    s_sensors[i] = new Adafruit_MAX31856(MAX31856_CS_PINS[i]);
    if (s_sensors[i]->begin()) {
      s_sensors[i]->setThermocoupleType(MAX31856_TCTYPE_K);
      s_ready[i] = true;
      Serial.printf("[OK] MAX31856 #%d (%s, CS=GPIO%d) 初期化完了\n",
                    i + 1, THERMOCOUPLE_LABELS[i], MAX31856_CS_PINS[i]);
    } else {
      s_ready[i] = false;
      Serial.printf("[ERROR] MAX31856 #%d (%s, CS=GPIO%d) 初期化失敗。配線を確認してください。\n",
                    i + 1, THERMOCOUPLE_LABELS[i], MAX31856_CS_PINS[i]);
    }
  }
}

bool Max31856Manager::isReady(int idx) const {
  if (idx < 0 || idx >= 4) return false;
  return s_ready[idx];
}

ThermocoupleReading Max31856Manager::read(int idx) {
  ThermocoupleReading r;
  r.temperatureC = NAN;
  r.coldJunctionC = NAN;
  r.faulted = false;
  r.faultRaw = 0;

  if (idx < 0 || idx >= 4 || !s_ready[idx] || s_sensors[idx] == nullptr) {
    r.faulted = true;
    return r;
  }

  double tc = s_sensors[idx]->readThermocoupleTemperature();
  double cj = s_sensors[idx]->readCJTemperature();
  uint8_t fault = s_sensors[idx]->readFault();

  r.temperatureC = tc;
  r.coldJunctionC = cj;
  r.faultRaw = fault;

  if (fault != 0 || isnan(tc)) {
    r.faulted = true;
    Serial.printf("[ERROR] 熱電対 #%d (%s) フォルト 0x%02X: ", idx + 1, THERMOCOUPLE_LABELS[idx], fault);
    if (fault & MAX31856_FAULT_CJRANGE) Serial.print("CJ範囲外 ");
    if (fault & MAX31856_FAULT_TCRANGE) Serial.print("TC範囲外 ");
    if (fault & MAX31856_FAULT_CJHIGH)  Serial.print("CJ高温 ");
    if (fault & MAX31856_FAULT_CJLOW)   Serial.print("CJ低温 ");
    if (fault & MAX31856_FAULT_TCHIGH)  Serial.print("TC高温 ");
    if (fault & MAX31856_FAULT_TCLOW)   Serial.print("TC低温 ");
    if (fault & MAX31856_FAULT_OVUV)    Serial.print("過電圧/低電圧 ");
    if (fault & MAX31856_FAULT_OPEN)    Serial.print("断線(OPEN) ");
    if (fault == 0 && isnan(tc))        Serial.print("原因不明(配線・ノイズ?)");
    Serial.println();
  }

  return r;
}
