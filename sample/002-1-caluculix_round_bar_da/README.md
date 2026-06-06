# 002-1 CalculiX round bar data assimilation

`002-1_laplacian_da_round_bar` と同じ丸棒データ同化を、CalculiX の熱伝導解析で行うための準備フォルダです。

フォルダ名は依頼に合わせて `002-1-caluculix_...` としています。ソルバ名としては CalculiX、実行コマンドは通常 `ccx` です。

現時点のこの環境では `ccx` が PATH に無いため、CalculiX 実行ケースは未実行です。まずはインストール方法と実装方針を整理しています。

## 目的

OpenFOAM 版で行っている次の流れを CalculiX に置き換えます。

```text
丸棒温度場を CalculiX で1ステップ予測
  ↓
同化用センサ温度で OI 補正
  ↓
補正後温度場を次ステップ初期温度として書き戻す
  ↓
DAなし / OIあり / truth を比較してグラフ化
```

## 丸棒条件

OpenFOAM 版 `sample/002-1_laplacian_da_round_bar` と同じ条件を使う予定です。

| 項目 | 値 |
|---|---:|
| 丸棒長さ | 0.30 m |
| 丸棒直径 | 0.01 m |
| 軸方向分割 | 30 |
| 時間刻み | 5 s |
| 計算時間 | 900 s |
| 熱拡散率 | 6.6e-5 m2/s |
| 同化センサ | N3, N12 |
| 検証センサ | N21, N27 |

## 次に作るもの

| ファイル | 役割 |
|---|---|
| `case/round_bar.inp` | CalculiX の熱伝導入力ファイル |
| `ccx_interface.py` | 初期温度を書き、`ccx` を1ステップ実行し、温度結果を読む |
| `oi/da_main.py` | OpenFOAM 版と同じ OI ループ |
| `oi/results/` | グラフと CSV の出力先 |

## インストール

CalculiX の導入方法は [INSTALL.md](INSTALL.md) を参照してください。

## 作業ログ

今回行ったことは [WORK_LOG.md](WORK_LOG.md) に記録しています。
