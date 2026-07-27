#!/usr/bin/env python3
"""
serial_logger.py
============================================================
ESP32(SUS304丸棒 熱膨張実験ロガー)のUSBシリアル出力を受信し、
CSVファイルへ保存する。Windows/Linux/macOS で動作する想定。

使い方:
    python3 serial_logger.py --port COM5
    python3 serial_logger.py --port /dev/ttyUSB0 --baud 115200
    python3 serial_logger.py --port COM5 --outdir ../data/raw --reconnect

主な機能 (docs/06_esp32_setup.md, README.md の要件19節に対応):
    - COMポート/ボーレート指定 (--port, --baud)
    - 自動ファイル名作成 (thermal_expansion_YYYYMMDD_HHMMSS.csv)
    - UTF-8で保存
    - 受信データのリアルタイム表示
    - 1行ごとにflush（電源断でもデータが残りやすい）
    - Ctrl+Cで安全終了
    - 不正な行(列数不一致・数値変換失敗)は別ログへ保存
    - 測定開始日時などのメタデータをJSONへ保存
    - 接続断時のエラー処理、任意で自動再接続(--reconnect)
============================================================
"""
import argparse
import json
import os
import sys
import time
from datetime import datetime

try:
    import serial
except ImportError:
    print("[エラー] pyserial がインストールされていません。 pip install -r requirements.txt を実行してください。")
    sys.exit(1)

# ESP32側 main.cpp / csv_output.cpp が出力する基本ヘッダ(先頭列名)。
# 列数チェックの基準として使う。ENABLE_COLD_JUNCTION_LOG=trueなら13列、falseなら10列。
EXPECTED_MIN_COLUMNS = 10


def build_filenames(outdir: str):
    os.makedirs(outdir, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base = f"thermal_expansion_{stamp}"
    csv_path = os.path.join(outdir, base + ".csv")
    invalid_path = os.path.join(outdir, base + ".invalid.log")
    meta_path = os.path.join(outdir, base + ".meta.json")
    return csv_path, invalid_path, meta_path, stamp


def is_comment_or_info_line(line: str) -> bool:
    """'#'で始まるメタデータ行、'[OK]'/'[ERROR]'/'[INFO]'/'[WARN]'等の診断行、
    '===...==='形式の起動バナー、'----...----'形式の区切りはCSVデータ行ではない。"""
    stripped = line.strip()
    if not stripped:
        return True
    if stripped.startswith("#"):
        return True
    if stripped.startswith("[") :
        return True
    if stripped.startswith("=") or stripped.startswith("-"):
        return True
    return False


def validate_data_line(line: str, expected_columns: int):
    """列数と数値変換を確認する。問題なければ(True, None)、問題があれば(False, 理由)。"""
    parts = line.strip().split(",")
    if len(parts) < expected_columns:
        return False, f"列数不足 ({len(parts)} < {expected_columns})"

    # 先頭(elapsed_ms)と、NaN許容の数値列を軽くチェックする。
    # 完全な列名チェックは行わず、致命的な壊れ方(非数値の混入等)だけ弾く。
    for i, p in enumerate(parts):
        p = p.strip()
        if p == "" or p.upper() == "NAN":
            continue
        try:
            float(p)
        except ValueError:
            return False, f"列{i}が数値/NaNとして解釈できません: '{p}'"
    return True, None


def open_serial(port: str, baud: float, timeout: float = 2.0):
    return serial.Serial(port, baudrate=int(baud), timeout=timeout)


def main() -> int:
    parser = argparse.ArgumentParser(description="ESP32熱膨張実験ロガー受信・CSV保存スクリプト")
    parser.add_argument("--port", required=True, help="シリアルポート (例: COM5, /dev/ttyUSB0)")
    parser.add_argument("--baud", type=float, default=115200, help="ボーレート (既定: 115200)")
    parser.add_argument("--outdir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "raw"),
                         help="CSV保存先ディレクトリ (既定: ../data/raw)")
    parser.add_argument("--expected-columns", type=int, default=EXPECTED_MIN_COLUMNS,
                         help="1行あたりの最小列数(既定10。冷接点温度列を含む場合は13を指定)")
    parser.add_argument("--reconnect", action="store_true", help="接続断時に自動再接続を試みる")
    parser.add_argument("--reconnect-interval", type=float, default=3.0, help="再接続の再試行間隔[秒]")
    args = parser.parse_args()

    csv_path, invalid_path, meta_path, stamp = build_filenames(args.outdir)

    meta = {
        "start_datetime": datetime.now().isoformat(),
        "port": args.port,
        "baud": args.baud,
        "csv_file": os.path.basename(csv_path),
        "invalid_log_file": os.path.basename(invalid_path),
    }
    with open(meta_path, "w", encoding="utf-8") as mf:
        json.dump(meta, mf, ensure_ascii=False, indent=2)
    print(f"メタデータを保存しました: {meta_path}")

    csv_f = open(csv_path, "w", encoding="utf-8", newline="")
    invalid_f = open(invalid_path, "w", encoding="utf-8", newline="")

    header_written = False
    n_rows = 0
    n_invalid = 0

    ser = None
    try:
        ser = open_serial(args.port, args.baud)
        print(f"接続しました: {args.port} @ {args.baud}bps")
        print(f"保存先: {csv_path}")
        print("Ctrl+C で終了します。")

        while True:
            try:
                raw = ser.readline()
            except serial.SerialException as e:
                print(f"[エラー] シリアル通信エラー: {e}")
                if not args.reconnect:
                    break
                ser.close()
                print(f"{args.reconnect_interval}秒後に再接続を試みます...")
                time.sleep(args.reconnect_interval)
                try:
                    ser = open_serial(args.port, args.baud)
                    print("再接続に成功しました。")
                    continue
                except serial.SerialException as e2:
                    print(f"[エラー] 再接続失敗: {e2}")
                    continue

            if not raw:
                continue  # タイムアウト(データなし)

            line = raw.decode("utf-8", errors="replace").rstrip("\r\n")
            if not line:
                continue

            print(line)  # リアルタイム表示

            if is_comment_or_info_line(line):
                # メタデータ行(# zero_cal_...)・起動ログ・診断メッセージはCSV本体には書かず、
                # そのまま保存ファイルへコメント行として残す(pandas等はcomment='#'で無視可能)。
                if line.strip().startswith("#"):
                    csv_f.write(line + "\n")
                    csv_f.flush()
                continue

            if not header_written:
                # 最初の非コメント行をヘッダとみなす
                csv_f.write(line + "\n")
                csv_f.flush()
                header_written = True
                continue

            ok, reason = validate_data_line(line, args.expected_columns)
            if ok:
                csv_f.write(line + "\n")
                csv_f.flush()
                n_rows += 1
            else:
                invalid_f.write(f"{datetime.now().isoformat()}\t{reason}\t{line}\n")
                invalid_f.flush()
                n_invalid += 1
                print(f"[警告] 不正な行をスキップしました: {reason}")

    except serial.SerialException as e:
        print(f"[エラー] シリアルポートを開けません: {e}")
        return 1
    except KeyboardInterrupt:
        print("\nCtrl+Cを検知しました。安全に終了します。")
    finally:
        if ser is not None and ser.is_open:
            ser.close()
        csv_f.close()
        invalid_f.close()
        print(f"正常行: {n_rows}  不正行: {n_invalid}")
        print(f"CSV: {csv_path}")
        print(f"不正行ログ: {invalid_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
