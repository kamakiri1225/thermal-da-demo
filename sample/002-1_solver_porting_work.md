# 002-1 solver porting work

## 作成したフォルダ

| フォルダ | 内容 |
|---|---|
| `sample/002-1-FrontISTR_round_bar_da` | FrontISTR 版の準備フォルダ |
| `sample/002-1-caluculix_round_bar_da` | CalculiX 版の準備フォルダ |

## 確認したこと

この環境で次のコマンドを確認した。

```bash
command -v fistr1
command -v ccx
```

どちらも PATH 上に見つからなかったため、FrontISTR / CalculiX の実解析とグラフ生成はまだ実行していない。

## 記載したこと

- FrontISTR のインストール候補を `sample/002-1-FrontISTR_round_bar_da/INSTALL.md` に記載した。
- CalculiX のインストール候補を `sample/002-1-caluculix_round_bar_da/INSTALL.md` に記載した。
- 丸棒データ同化として使う解析条件を各 `README.md` に記載した。
- 今回行ったことと未実施項目を各 `WORK_LOG.md` に記載した。

## 次の作業

1. `fistr1` と `ccx` が使える環境を準備する。
2. FrontISTR / CalculiX の熱伝導1ステップ実行ケースを作る。
3. Python インターフェースで初期温度の書き込み、実行、温度読み取りを実装する。
4. OpenFOAM 版と同じ OI ループで `results_of_da.png` と `summary_rmse.csv` を出す。
