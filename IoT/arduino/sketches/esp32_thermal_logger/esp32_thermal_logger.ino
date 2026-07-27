/*
 * esp32_thermal_logger.ino
 * ============================================================
 * ESP32 温度・湿度・表面温度ロガー（第1段階: microSD保存）
 *
 * 測定項目:
 *   - 固体表面温度 : K型熱電対 + MAX31855（1〜4点、拡張可能）
 *   - 雰囲気温度   : BME280
 *   - 相対湿度     : BME280
 *   - 気圧         : BME280
 *
 * 保存先:
 *   - microSDカード（CSV形式、ヘッダ自動書き込み）
 *   - シリアルモニタにも同時出力（115200 bps）
 *
 * 配線の詳細は docs/wiring.md を参照してください。
 *
 * ── 使い方 ──────────────────────────────────────────────
 * 1. 下の「設定」セクションで測定間隔・センサ数を確認する
 * 2. Arduino IDE でボード「ESP32 Dev Module」を選び書き込む
 * 3. シリアルモニタ(115200bps)を開いて動作を確認する
 * 4. microSDの LOG.CSV にデータが溜まる
 *
 * ── 時刻について ────────────────────────────────────────
 * USE_WIFI_NTP を 1 にして Wi-Fi 情報を書くと、NTPから実時刻
 * を取得して "2026-07-17 10:00:00" 形式で記録します。
 * 0 のままなら「起動からの経過秒」を記録します（Wi-Fi不要）。
 * まず動かすだけなら 0 のままでOKです。
 * ============================================================
 */

#include <SPI.h>
#include <SD.h>
#include <Wire.h>
#include <Adafruit_MAX31855.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_BME280.h>

// ============================================================
// 設定（ここだけ変更すれば運用できます）
// ============================================================

// 測定間隔 [ミリ秒]。10000 = 10秒ごと。
const unsigned long MEASURE_INTERVAL_MS = 10000;

// 表面温度（MAX31855）の点数。1〜4で変更。
// 2点にしたいときは 2 にして、2個目のモジュールのCSを
// MAX31855_CS_PINS の2番目のピンへ配線するだけです。
const int NUM_SURFACE_SENSORS = 1;

// MAX31855 の CS(チップセレクト) ピン割り当て（最大4点分）
const int MAX31855_CS_PINS[4] = {25, 26, 27, 33};

// microSD モジュールの CS ピン
const int SD_CS_PIN = 5;

// ログファイル名（8.3形式が無難）
const char* LOG_FILENAME = "/LOG.CSV";

// Wi-Fi + NTP で実時刻を使うか（0=使わない, 1=使う）
#define USE_WIFI_NTP 0

#if USE_WIFI_NTP
#include <WiFi.h>
#include <time.h>
const char* WIFI_SSID     = "あなたのSSID";
const char* WIFI_PASSWORD = "あなたのパスワード";
// 日本標準時 (UTC+9)
const long  GMT_OFFSET_SEC      = 9 * 3600;
const int   DAYLIGHT_OFFSET_SEC = 0;
const char* NTP_SERVER          = "ntp.nict.jp";
#endif

// ============================================================
// グローバル変数
// ============================================================

// MAX31855はハードウェアSPI(共有バス)を使用。CSだけ個別。
Adafruit_MAX31855* thermocouples[4] = {nullptr, nullptr, nullptr, nullptr};

Adafruit_BME280 bme;          // I2C接続 (SDA=21, SCL=22)
bool bmeFound = false;        // BME280が見つかったか
bool sdReady  = false;        // SDカードが使えるか

unsigned long lastMeasureMs = 0;

// ============================================================
// 初期化
// ============================================================
void setup() {
  Serial.begin(115200);
  delay(1500);  // シリアルモニタを開く時間の余裕
  Serial.println();
  Serial.println("=== ESP32 温度・湿度ロガー 起動 ===");

  // --- MAX31855（表面温度） ---
  for (int i = 0; i < NUM_SURFACE_SENSORS; i++) {
    thermocouples[i] = new Adafruit_MAX31855(MAX31855_CS_PINS[i]);
    if (!thermocouples[i]->begin()) {
      Serial.printf("[エラー] MAX31855 #%d (CS=GPIO%d) の初期化に失敗。配線を確認してください。\n",
                    i + 1, MAX31855_CS_PINS[i]);
    } else {
      Serial.printf("[OK] MAX31855 #%d (CS=GPIO%d) 初期化完了\n", i + 1, MAX31855_CS_PINS[i]);
    }
  }
  delay(500);  // MAX31855の安定待ち

  // --- BME280（雰囲気温度・湿度・気圧） ---
  // アドレスは 0x76 が多いが、モジュールにより 0x77 の場合もある
  if (bme.begin(0x76)) {
    bmeFound = true;
  } else if (bme.begin(0x77)) {
    bmeFound = true;
  }
  if (bmeFound) {
    Serial.println("[OK] BME280 初期化完了");
  } else {
    Serial.println("[エラー] BME280 が見つかりません。I2C配線(SDA=21, SCL=22)と");
    Serial.println("        アドレス(0x76/0x77)を確認してください。測定は続行します。");
  }

  // --- microSD ---
  if (SD.begin(SD_CS_PIN)) {
    sdReady = true;
    Serial.println("[OK] microSD 初期化完了");
    writeCsvHeaderIfNeeded();
  } else {
    Serial.println("[エラー] microSD の初期化に失敗。カードの挿入と配線(CS=GPIO5)を");
    Serial.println("        確認してください。シリアル出力のみで続行します。");
  }

#if USE_WIFI_NTP
  // --- Wi-Fi + NTP（実時刻） ---
  Serial.printf("Wi-Fi '%s' に接続中", WIFI_SSID);
  WiFi.begin(WIFI_SSID, WIFI_PASSWORD);
  for (int i = 0; i < 20 && WiFi.status() != WL_CONNECTED; i++) {
    delay(500);
    Serial.print(".");
  }
  Serial.println();
  if (WiFi.status() == WL_CONNECTED) {
    Serial.println("[OK] Wi-Fi 接続完了。NTPで時刻同期します。");
    configTime(GMT_OFFSET_SEC, DAYLIGHT_OFFSET_SEC, NTP_SERVER);
  } else {
    Serial.println("[警告] Wi-Fi に接続できません。経過秒で記録します。");
  }
#endif

  Serial.printf("測定間隔: %.1f 秒 / 表面温度: %d 点\n",
                MEASURE_INTERVAL_MS / 1000.0, NUM_SURFACE_SENSORS);
  Serial.println("---- 測定開始 ----");
}

// ============================================================
// メインループ
// ============================================================
void loop() {
  unsigned long now = millis();
  if (now - lastMeasureMs >= MEASURE_INTERVAL_MS || lastMeasureMs == 0) {
    lastMeasureMs = now;
    measureAndLog();
  }
}

// ============================================================
// 1回分の測定 → シリアル出力 + SD追記
// ============================================================
void measureAndLog() {
  // --- タイムスタンプ ---
  String timestamp = getTimestamp();

  // --- 表面温度（MAX31855, 複数点） ---
  double surfaceTemp[4];
  for (int i = 0; i < NUM_SURFACE_SENSORS; i++) {
    surfaceTemp[i] = readThermocouple(i);
  }

  // --- 雰囲気温度・湿度・気圧（BME280） ---
  float ambientTemp = NAN, humidity = NAN, pressure = NAN;
  if (bmeFound) {
    ambientTemp = bme.readTemperature();   // [°C]
    humidity    = bme.readHumidity();      // [%]
    pressure    = bme.readPressure() / 100.0F;  // [hPa]
  }

  // --- CSVの1行を組み立てる ---
  // 例: 2026-07-17 10:00:00,28.5,27.1,61.2,1004.5
  String line = timestamp;
  for (int i = 0; i < NUM_SURFACE_SENSORS; i++) {
    line += "," + formatValue(surfaceTemp[i]);
  }
  line += "," + formatValue(ambientTemp);
  line += "," + formatValue(humidity);
  line += "," + formatValue(pressure);

  // --- シリアルモニタへ ---
  Serial.println(line);

  // --- microSDへ追記 ---
  if (sdReady) {
    File f = SD.open(LOG_FILENAME, FILE_APPEND);
    if (f) {
      f.println(line);
      f.close();  // 毎回閉じる: 電源断でもデータが残りやすい
    } else {
      Serial.println("[エラー] ログファイルを開けません。SDカードを確認してください。");
    }
  }
}

// ============================================================
// 熱電対1点の読み取り（エラー内容も表示）
// ============================================================
double readThermocouple(int idx) {
  if (thermocouples[idx] == nullptr) return NAN;

  double tc = thermocouples[idx]->readCelsius();
  if (isnan(tc)) {
    uint8_t err = thermocouples[idx]->readError();
    Serial.printf("[エラー] 熱電対 #%d 読み取り失敗: ", idx + 1);
    if (err & MAX31855_FAULT_OPEN)      Serial.print("断線(OPEN) ");
    if (err & MAX31855_FAULT_SHORT_GND) Serial.print("GNDショート ");
    if (err & MAX31855_FAULT_SHORT_VCC) Serial.print("VCCショート ");
    if (err == 0)                       Serial.print("原因不明(配線・ノイズ?)");
    Serial.println();
  }
  return tc;
}

// ============================================================
// CSVヘッダをファイルが無いときだけ書く
// ============================================================
void writeCsvHeaderIfNeeded() {
  if (SD.exists(LOG_FILENAME)) {
    Serial.println("既存のログファイルに追記します。");
    return;
  }
  File f = SD.open(LOG_FILENAME, FILE_WRITE);
  if (!f) {
    Serial.println("[エラー] ヘッダの書き込みに失敗しました。");
    return;
  }
  // ヘッダをセンサ数に合わせて動的に生成
  // 例: timestamp,surface_temp_1_C,ambient_temp_C,humidity_percent,pressure_hPa
  String header = "timestamp";
  for (int i = 0; i < NUM_SURFACE_SENSORS; i++) {
    header += ",surface_temp_" + String(i + 1) + "_C";
  }
  header += ",ambient_temp_C,humidity_percent,pressure_hPa";
  f.println(header);
  f.close();
  Serial.println("CSVヘッダを書き込みました: " + header);
}

// ============================================================
// タイムスタンプ文字列
//   USE_WIFI_NTP=1: "2026-07-17 10:00:00"
//   USE_WIFI_NTP=0: 起動からの経過秒 "123.4"
// ============================================================
String getTimestamp() {
#if USE_WIFI_NTP
  struct tm t;
  if (getLocalTime(&t, 100)) {
    char buf[24];
    strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", &t);
    return String(buf);
  }
  // NTP未同期のときは経過秒にフォールバック
#endif
  return String(millis() / 1000.0, 1);
}

// ============================================================
// 数値→文字列（NaNは "NaN" と書く。後処理のPython側で除外できる）
// ============================================================
String formatValue(double v) {
  if (isnan(v)) return "NaN";
  return String(v, 2);
}
