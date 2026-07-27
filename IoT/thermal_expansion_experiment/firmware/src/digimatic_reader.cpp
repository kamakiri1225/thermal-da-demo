#include "digimatic_reader.h"
#include "config.h"

DigimaticReader g_digimatic;

// ============================================================
// ISR共有状態（volatile必須）
// ============================================================
static volatile uint64_t s_bitBuffer = 0;
static volatile int s_bitCount = 0;
static volatile unsigned long s_lastEdgeMs = 0;

static double s_lastDisplacementMm = NAN;
static bool s_lastValid = false;
static unsigned long s_lastValidMs = 0;

// ------------------------------------------------------------
// クロックエッジ割込み。
// 要信号波形確認: 現状は「立ち上がりで1ビットサンプル」としているが、
// 実測波形次第でCHANGE/FALLINGへの変更が必要になる可能性がある。
// ------------------------------------------------------------
static void IRAM_ATTR onSpcClockEdge() {
  bool bit = digitalRead(PIN_SPC_DATA);
  if (s_bitCount < SPC_MAX_FRAME_BITS) {
    s_bitBuffer = (s_bitBuffer << 1) | (bit ? 1ULL : 0ULL);
    s_bitCount++;
  }
  s_lastEdgeMs = millis();
}

void DigimaticReader::begin() {
  pinMode(PIN_SPC_DATA, INPUT);   // 入力のみ。OUTPUT HIGHにしない(docs/05参照)
  pinMode(PIN_SPC_CLOCK, INPUT);  // 入力のみ
#if USE_SPC_REQ_LINE
  // REQが必要と判明した場合のみ有効化する。要信号波形確認。
  pinMode(PIN_SPC_REQ, OUTPUT);
  digitalWrite(PIN_SPC_REQ, LOW);
#endif
  attachInterrupt(digitalPinToInterrupt(PIN_SPC_CLOCK), onSpcClockEdge, RISING);
  Serial.printf("[OK] Digimatic/SPC 入力初期化 (DATA=GPIO%d, CLOCK=GPIO%d)\n",
                PIN_SPC_DATA, PIN_SPC_CLOCK);
  Serial.println("[INFO] Digimaticの詳細ビット割当は未確定です。docs/05_digimatic_protocol.md を参照してください。");
}

// ------------------------------------------------------------
// bits(chronological, index0=最初に受信したビット) を数値へデコードする。
// 暫定実装: config.h の SPC_* 定数は仮値。実測結果で調整すること。
// ------------------------------------------------------------
static bool decodeBcdFrame(const uint8_t* bitsChrono, int bitCount, double& outValueMm) {
  if (bitCount != SPC_EXPECTED_FRAME_BITS) {
    return false;  // フレーム長が想定と違う間はデコードしない(誤検出防止)
  }

  // ビット順の解釈: SPC_BIT_ORDER_MSB_FIRSTがtrueなら受信順そのまま、
  // falseなら値としては逆順(LSBファースト受信)とみなして並べ替える。
  uint8_t bits[SPC_MAX_FRAME_BITS];
  if (SPC_BIT_ORDER_MSB_FIRST) {
    for (int i = 0; i < bitCount; i++) bits[i] = bitsChrono[i];
  } else {
    for (int i = 0; i < bitCount; i++) bits[i] = bitsChrono[bitCount - 1 - i];
  }

  int dataEndBit = SPC_DATA_START_BIT + 4 * SPC_DATA_NUM_DIGITS;
  if (SPC_SIGN_BIT_INDEX >= bitCount || dataEndBit > bitCount) {
    return false;  // ビット位置設定が現在のフレーム長と矛盾している
  }

  bool negative = bits[SPC_SIGN_BIT_INDEX] != 0;

  long rawInt = 0;
  for (int d = 0; d < SPC_DATA_NUM_DIGITS; d++) {
    int nibbleStart = SPC_DATA_START_BIT + d * 4;
    uint8_t nibble = 0;
    for (int b = 0; b < 4; b++) {
      nibble = (nibble << 1) | (bits[nibbleStart + b] ? 1 : 0);
    }
    if (nibble > 9) {
      // BCDとして不正な値。ビット割当がまだ合っていない可能性が高い。
      return false;
    }
    rawInt = rawInt * 10 + nibble;
  }

  double value = (double)rawInt;
  for (int i = 0; i < SPC_DECIMAL_POINT_DIGITS_FROM_RIGHT; i++) value /= 10.0;
  if (negative) value = -value;
  if (INDICATOR_IS_INCH_MODE) value *= INCH_TO_MM;

  outValueMm = value;
  return true;
}

void DigimaticReader::processFrame(uint64_t bits, int bitCount) {
  if (bitCount == 0) return;

  uint8_t chrono[SPC_MAX_FRAME_BITS];
  for (int i = 0; i < bitCount; i++) {
    chrono[i] = (uint8_t)((bits >> (bitCount - 1 - i)) & 1ULL);
  }

#if SPC_RAW_DEBUG_OUTPUT
  Serial.printf("[SPC_RAW] bits=%d value=0x%llX  raw=", bitCount, (unsigned long long)bits);
  for (int i = 0; i < bitCount; i++) Serial.print(chrono[i]);
  Serial.println();
#endif

  double value;
  if (decodeBcdFrame(chrono, bitCount, value)) {
    s_lastDisplacementMm = value;
    s_lastValid = true;
    s_lastValidMs = millis();
  } else {
    // フレーム長不一致 or BCD不正 → 無効値として扱う(誤ったもっともらしい値を出さない)
    s_lastValid = false;
    Serial.printf("[INFO] SPCフレームをデコードできませんでした (bits=%d, 想定=%d)。"
                  "config.h の SPC_* 定数を実測波形に合わせて調整してください。\n",
                  bitCount, SPC_EXPECTED_FRAME_BITS);
  }
}

void DigimaticReader::update() {
  unsigned long now = millis();
  bool timedOut;
  uint64_t bufSnapshot;
  int countSnapshot;

  noInterrupts();
  timedOut = (s_bitCount > 0) && (now - s_lastEdgeMs >= SPC_FRAME_TIMEOUT_MS);
  if (timedOut) {
    bufSnapshot = s_bitBuffer;
    countSnapshot = s_bitCount;
    s_bitBuffer = 0;
    s_bitCount = 0;
  }
  interrupts();

  if (timedOut) {
    processFrame(bufSnapshot, countSnapshot);
  }

  // 長時間フレームを受信できていない場合は無効扱いにする
  if (now - s_lastValidMs > SPC_STALE_TIMEOUT_MS) {
    s_lastValid = false;
  }
}

double DigimaticReader::getDisplacementMm() const {
  return s_lastValid ? s_lastDisplacementMm : NAN;
}

bool DigimaticReader::isValid() const {
  return s_lastValid;
}

unsigned long DigimaticReader::getLastValidMs() const {
  return s_lastValidMs;
}
