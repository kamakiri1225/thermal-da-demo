/*
 * csv_output.h
 * ============================================================
 * CSVヘッダ生成・1行分の組み立て・Serial出力を担当。
 * 列構成は docs/01_system_overview.md / README.md のCSV仕様と対応。
 *   elapsed_ms,T1_C,T2_C,T3_C,indicator_mm,disp_mm,spc_valid,fault1,fault2,fault3
 *   (+ ENABLE_COLD_JUNCTION_LOG時: cold_junction_1_C,2_C,3_C)
 * ============================================================
 */
#pragma once

#include <Arduino.h>
#include "max31856_manager.h"

// 起動時に1回だけヘッダを出力する
void csvPrintHeader();

// 1回分の測定データをCSV1行としてSerialへ出力する
void csvPrintRow(
    unsigned long elapsedMs,
    const ThermocoupleReading readings[3],
    double indicatorMm,
    double dispMm,
    bool spcValid
);
