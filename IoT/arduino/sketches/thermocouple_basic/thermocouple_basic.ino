#include <SPI.h>
#include "Adafruit_MAX31856.h"

#define MAX31856_CS 10

Adafruit_MAX31856 maxthermo = Adafruit_MAX31856(MAX31856_CS);

void setup() {
  Serial.begin(9600);

  Serial.println("MAX31856 K-type Thermocouple Test");

  if (!maxthermo.begin()) {
    Serial.println("MAX31856が見つかりません。配線を確認してください。");
    while (1);
  }

  maxthermo.setThermocoupleType(MAX31856_TCTYPE_K);

  Serial.println("測定開始");
}

void loop() {
  double tempC = maxthermo.readThermocoupleTemperature();
  double cjTempC = maxthermo.readCJTemperature();

  Serial.print("熱電対温度: ");
  Serial.print(tempC);
  Serial.println(" ℃");

  Serial.print("基板温度: ");
  Serial.print(cjTempC);
  Serial.println(" ℃");

  uint8_t fault = maxthermo.readFault();
  if (fault) {
    Serial.print("Fault: 0x");
    Serial.println(fault, HEX);

    if (fault & MAX31856_FAULT_OPEN) {
      Serial.println("熱電対が接続されていない、または断線の可能性があります。");
    }
    if (fault & MAX31856_FAULT_OVUV) {
      Serial.println("電圧異常の可能性があります。");
    }
    if (fault & MAX31856_FAULT_TCRANGE) {
      Serial.println("熱電対の測定範囲外の可能性があります。");
    }
  }

  Serial.println("--------------------");
  delay(1000);
}
