# 参考文献・参考実装

更新日: 2026-07-27

## 材料物性値

- 日本機械学会・日本金属学会等の材料便覧に見られるSUS304(AISI304相当)の代表的物性値
- NIMS(物質・材料研究機構) 材料データベース(MatNavi)の公開代表値
- ASM Handbook Vol.1 (Properties and Selection: Irons, Steels, and High-Performance Alloys) のSUS304/AISI304相当データ
- 各社ステンレス鋼材カタログの物性値表

物性値は代表値であり温度依存性がある。詳細は `config/material_properties_sus304.json` のコメントを参照。実際に使用する試験片のミルシートがあれば、それを優先すること。

## センサ・部品

- MAX31856データシート（Maxim Integrated / Analog Devices）: フォルトフラグ・CJ補償の仕様確認用（**要メーカー仕様確認**）
- Adafruit MAX31856 library ドキュメント（GitHub: adafruit/Adafruit_MAX31856）
- K型熱電対 JIS C1602 規格（許容差クラスの確認用）
- ミツトヨ デジマチックインジケータ（543-790B-10候補）の製品仕様（**購入前に最新メーカー資料で要確認**）
- ミツトヨ Digimatic/SPC接続ケーブル（905338候補）の製品仕様（**購入前に最新メーカー資料で要確認**）

## Digimatic/SPCプロトコルに関する公開情報（参考、未検証）

以下は公開されている実装例・解説記事として一般に参照されるものの例。**これらの存在は、本プロジェクトで使用する具体的な型番の組合せでの動作を保証するものではない**（`docs/05_digimatic_protocol.md`参照）。

- MGX3D / EspDRO プロジェクト（ミツトヨDigimatic/SPC信号読み取りの公開実装例）
- Roger-random / mitutoyo リポジトリ（GitHub上の公開実装例）
- Arduino Digimatic関連の公開記事・フォーラム投稿（EEVblog等）多数
- Mitutoyo Digimatic SPC protocol の非公式な解説記事

本プロジェクトの `firmware/src/digimatic_reader.cpp` は、上記のような一般的な説明を参考にしつつ、特定のソースコードを直接移植してはいない（生ビットキャプチャ+要調整デコードという独自の暫定実装）。もし今後、上記や他の公開実装から具体的なコードを移植する場合は、当該プロジェクトのライセンス条件を確認し、ここに出典とライセンス名を追記すること。

## OpenCAE・データ同化関連（本プロジェクトの応用先）

- 本リポジトリ内 `study/practice/slides_opencae_intro.html`（データ同化の解説スライド、関連する先行研究の紹介を含む）
- Kalman Filter-Driven State Observer 関連の先行研究（`docs/Kalman Filter-Driven State Observer for.pdf`、Manufacturing Letters誌掲載、ETH Zurich）

## 未確定事項に関する注記

本ドキュメント一覧に挙げた資料の多くは「一般的な代表値・一般に説明される情報」であり、実機・実データでの検証が済んでいない。各`docs/*.md`内の「要確認」「要実機確認」「要信号波形確認」「要メーカー仕様確認」の表記を優先し、実験を進める中で本ファイルおよび関連ドキュメントを更新すること。
