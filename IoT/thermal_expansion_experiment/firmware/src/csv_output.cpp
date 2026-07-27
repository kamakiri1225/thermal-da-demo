#include "csv_output.h"
#include "config.h"

#include <math.h>

static String formatValue(double v, int decimals) {
  if (isnan(v)) return "NaN";
  return String(v, decimals);
}

void csvPrintHeader() {
  String header = "elapsed_ms,T1_C,T2_C,T3_C,indicator_mm,disp_mm,spc_valid,fault1,fault2,fault3";
  if (ENABLE_COLD_JUNCTION_LOG) {
    header += ",cold_junction_1_C,cold_junction_2_C,cold_junction_3_C";
  }
  Serial.println(header);
}

void csvPrintRow(
    unsigned long elapsedMs,
    const ThermocoupleReading readings[3],
    double indicatorMm,
    double dispMm,
    bool spcValid
) {
  String line = String(elapsedMs);

  for (int i = 0; i < 3; i++) {
    line += "," + formatValue(readings[i].temperatureC, 3);
  }

  line += "," + formatValue(indicatorMm, 4);
  line += "," + formatValue(dispMm, 4);
  line += "," + String(spcValid ? 1 : 0);

  for (int i = 0; i < 3; i++) {
    line += "," + String(readings[i].faulted ? 1 : 0);
  }

  if (ENABLE_COLD_JUNCTION_LOG) {
    for (int i = 0; i < 3; i++) {
      line += "," + formatValue(readings[i].coldJunctionC, 3);
    }
  }

  // データ行はここでのみ出力する。デバッグ文字列を混ぜない([INFO]/[ERROR]等は
  // max31856_manager.cpp/digimatic_reader.cpp側で別行として出している)。
  Serial.println(line);
}
