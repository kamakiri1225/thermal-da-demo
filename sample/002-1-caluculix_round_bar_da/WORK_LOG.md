# Work log

## 2026-06-06

- `sample/002-1-caluculix_round_bar_da/` を作成した。
- この環境で `ccx` を確認したが、PATH 上に見つからなかった。
- CalculiX 版は未実行。まずインストール手順を `INSTALL.md` に整理した。
- OpenFOAM 版 `sample/002-1_laplacian_da_round_bar` と同じ丸棒条件、同化センサ、検証センサを使う方針を `README.md` に記載した。

## 未実施

- CalculiX 用 `.inp` 熱伝導モデルの作成。
- Python から `ccx` を1ステップずつ実行する `ccx_interface.py` の実装。
- OI ループの実行と `results_of_da.png` の生成。
