# Work log

## 2026-06-06

- `sample/002-1-FrontISTR_round_bar_da/` を作成した。
- この環境で `fistr1` を確認したが、PATH 上に見つからなかった。
- Codex から `sudo apt update` を実行したが、sudo パスワード入力が必要で非対話実行できなかった。
- ユーザーが WSL ターミナルで実行できるよう、`INSTALL.md` に apt 確認、依存パッケージ、ソースビルド候補を追記した。
- ユーザーが `~/src/FrontISTR` に公式リポジトリを clone した。
- Codex 側で CMake configure と build を実行し、`/home/kamakiri/local/frontistr/bin/fistr1` を生成した。
- `/home/kamakiri/local/frontistr/bin/fistr1 -v` で FrontISTR 5.9 の起動を確認した。
- 公式リポジトリ内の既存熱伝導サンプル `tests/analysis/heat/exP/P361.*` を `/tmp/frontistr_heat_smoke` で実行し、FrontISTR の熱伝導解析が動くことを確認した。
- OpenFOAM 版 `sample/002-1_laplacian_da_round_bar` と同じ丸棒条件、同化センサ、検証センサを使う方針を `README.md` に記載した。
- FrontISTR 用の `fistr_interface.py` を追加した。
  - φ10 mm x 300 mm の丸棒を、中央ブロック + 外周4ブロックの六面体メッシュとして生成する。
  - `!INITIAL_CONDITION,TYPE=TEMPERATURE` で節点温度を毎ステップ書き戻す。
  - `!FIXTEMP` で左右端温度境界を与える。
- `oi/da_main.py` と `oi/run.sh` を追加した。
  - truth / FrontISTR-only / FrontISTR+OI の3ケースを実行する。
  - OpenFOAM 版と同様に OI でセンサ節点の観測を同化する。
  - OpenFOAM 版と同じ 180 step = 900 s の条件で実行する。
  - 外周節点を円弧上に配置し、FrontISTR 用の六面体近似丸棒メッシュとした。
- 画像出力先を `oi/results/img/` に統一した。
  - `results_frontistr_da.png`、`measurement_points_frontistr_da.png`、`all_axial_nodes_timeseries_frontistr_da.png`、`axial_nodes_heatmap_frontistr_da.png`、`axial_nodes_error_heatmap_frontistr_da.png` をここへ保存する。
  - README のリンク先も `oi/results/img/` に更新した。

## 未実施

- `oi/run.sh` の完走確認。
- `oi/results/img/results_frontistr_da.png` の生成確認。
