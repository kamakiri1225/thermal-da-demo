/*
 * digimatic_reader.h
 * ============================================================
 * ミツトヨ Digimatic/SPC 信号の読み取り。
 *
 * 重要: この実装は「生ビット列を確認しながら仕様を確定させる」ことを
 * 前提にした暫定実装。docs/05_digimatic_protocol.md を必ず参照。
 * デコード結果(getDisplacementMm)は SPC_EXPECTED_FRAME_BITS と一致した
 * フレームのみ有効とみなす設計にして、誤ったもっともらしい値を
 * 記録しないようにしている。
 * ============================================================
 */
#pragma once

#include <Arduino.h>

class DigimaticReader {
public:
  void begin();

  // メインループから毎周期呼ぶ。フレームタイムアウトの判定・デコードを行う。
  void update();

  // 最後にデコードできた変位値[mm]。無効なら NAN。
  double getDisplacementMm() const;

  // 直近のフレームが有効か（SPC_EXPECTED_FRAME_BITSと一致し、かつ
  // 一定時間内に受信できているか）
  bool isValid() const;

  // 最後に有効なフレームを受信した時刻(millis)
  unsigned long getLastValidMs() const;

private:
  void processFrame(uint64_t bits, int bitCount);
};

extern DigimaticReader g_digimatic;
