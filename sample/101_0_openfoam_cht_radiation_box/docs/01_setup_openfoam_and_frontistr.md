# 01. セットアップ手順(OpenFOAM + FrontISTR、GitHubからclone後)

更新日: 2026-07-29

このリポジトリをGitHubから`git clone`した状態から、`101_0_openfoam_cht_radiation_box`(OpenFOAM)と`101_1_frontistr_cht_box_thermal_expansion`(FrontISTR)を実際に動かすまでの手順をまとめる。**環境はLinux(Ubuntu系)を想定**している(WSL2でも同様)。

## 0. 全体の流れ

```
git clone <このリポジトリ>
  ├─ 1. OpenFOAM(ESI/OpenCFD版)をインストール
  ├─ 2. FrontISTRをソースからビルド
  ├─ 3. 101_0を実行(OpenFOAMのCHT+輻射計算)
  └─ 4. 101_1を実行(101_0の温度分布→FrontISTRで熱膨張計算)
```

## 1. OpenFOAMのインストール(ESI/OpenCFD版)

**重要**: OpenFOAMには系統が2つある(`docs/00_tutorial_derivation.md`参照)。本ケースは**ESI/OpenCFD版(openfoam.com系)** 用に書かれている。OpenFOAM Foundation版(openfoam.org系、`openfoam13`等)では動作しない可能性がある。

### 1.1 公式aptリポジトリを登録してインストール(推奨)

```bash
# 公式リポジトリの鍵とソースリストを追加
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash

# パッケージ一覧を更新
sudo apt-get update

# 使いたいバージョンをインストール(例: v2512。他バージョンも同様の書式)
sudo apt-get install -y openfoam2512-default
```

インストール可能なバージョン名は`apt-cache search openfoam`で確認できる(`openfoam2406`, `openfoam2506`, `openfoam2512`等)。**本ケースは`openfoam2512`で動作確認済み**(`openfoam2406`でも構築・検証済み)。

### 1.2 環境の有効化

OpenFOAMはターミナルを開くたびに環境変数を読み込む必要がある。

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
```

**毎回打つのが面倒な場合**は`~/.bashrc`の末尾に追記しておくとよい(ただし、複数バージョンを使い分けたい場合は都度sourceする方が安全)。

### 1.3 インストール確認

```bash
which blockMesh chtMultiRegionFoam splitMeshRegions topoSet changeDictionary postProcess
blockMesh -help | head -3
```

コマンドのパスが表示され、`-help`がエラーなく実行できればOK。

## 2. FrontISTRのインストール(ソースビルド)

FrontISTRはaptパッケージが無いため、**ソースからビルド**する。ここでは並列計算なし(MPI無効)の最小構成でビルドする手順を示す(101_1のケースは小規模なので並列化は不要)。

### 2.1 必要パッケージ

```bash
sudo apt-get install -y build-essential gfortran cmake git libblas-dev liblapack-dev
```

- `build-essential`: C/C++コンパイラ
- `gfortran`: Fortranコンパイラ(FrontISTRはFortran90で書かれた部分が多い)
- `cmake`: ビルドシステム
- `libblas-dev` / `liblapack-dev`: 一部の解析機能で使用(`WITH_LAPACK`)

**MPI並列やMETISによる領域分割を使いたい場合**(大規模計算向け、本ケースでは不要)は、追加で以下も入れる。

```bash
sudo apt-get install -y libopenmpi-dev libmetis-dev
```

### 2.2 ソース取得とビルド

```bash
git clone https://github.com/FrontISTR/FrontISTR.git
cd FrontISTR

mkdir build && cd build
cmake ..
make -j$(nproc)
```

ビルドが成功すると`build/fistr1/fistr1`(実行ファイル)ができる。

**要確認**: MPI/METIS/Trilinos-ML等を有効にしたい場合は、`cmake ..`の代わりに`cmake -DWITH_MPI=ON -DWITH_METIS=ON ..`のようにオプションを付ける(詳細は`FrontISTR/README.md`の「System Requirements & Dependencies」節、および公式マニュアル https://manual.frontistr.com/ja/ を参照)。本リポジトリの動作確認は`WITH_MPI=OFF`, `WITH_METIS=OFF`, `WITH_LAPACK=ON`(デフォルトに近い最小構成)で行った。

### 2.3 実行ファイルへのパスを通す

101_1側のPythonスクリプト(`fistr_case.py`)は、既定で`/home/kamakiri/local/frontistr/bin/fistr1`のようなパスを見るのではなく、**環境変数`FISTR1`で指定したパス、無指定なら`fistr1`をそのままコマンドとして呼ぶ**ようになっている(`fistr_case.py`の`FISTR_BIN`参照)。以下のいずれかの方法でパスを通す。

```bash
# 方法A: 環境変数で明示的に指定する(推奨、都度指定)
export FISTR1=/path/to/FrontISTR/build/fistr1/fistr1

# 方法B: PATHに追加する(ログインシェル起動時に毎回必要ならbashrcへ)
export PATH="/path/to/FrontISTR/build/fistr1:$PATH"
```

### 2.4 インストール確認

```bash
fistr1 --version   # または $FISTR1 --version
```

バージョン情報が表示されればOK。

## 3. 101_0を実行する(OpenFOAM CHT+輻射)

```bash
cd sample/101_0_openfoam_cht_radiation_box
source /usr/lib/openfoam/openfoam2512/etc/bashrc   # 1.2節、まだなら実行

./Allrun
```

`Allrun`は内部で`blockMesh → topoSet → splitMeshRegions → changeDictionary → chtMultiRegionFoam`を順に実行する(詳細は`Allrun.pre`, `Allrun`参照)。既定の`controlDict`は`endTime 600`(秒、ヒーター0-300秒ON→300-600秒OFFの1サイクル)。8000セル程度の背景メッシュであれば数分〜十数分程度で完了する。

途中経過は`postProcessing/solidSurfaceProbes/`, `postProcessing/heaterMatSurfaceProbes/`の温度時系列で確認できる。

掃除(生成物を消して初期状態に戻す)する場合:

```bash
./Allclean
```

## 4. 101_1を実行する(FrontISTR熱膨張、101_0の温度を引き継ぐ)

101_0の実行が完了したら、以下でFrontISTR熱膨張解析を実行する。

```bash
cd ../101_1_frontistr_cht_box_thermal_expansion
python3 -m venv .venv && source .venv/bin/activate   # 任意、Python仮想環境を使う場合
pip install -r requirements.txt

cd python

# 単一時刻(最新時刻)だけ計算する場合
python3 run_thermal_expansion.py --of-case ../../101_0_openfoam_cht_radiation_box --time latestTime

# 101_0の時刻歴すべてについて計算し、変位の時刻歴を得る場合
python3 run_thermal_expansion_timehistory.py --of-case ../../101_0_openfoam_cht_radiation_box
```

**事前に`export FISTR1=...`(2.3節)を済ませておくこと**。OpenFOAM環境(`source .../etc/bashrc`)も同じシェルでsource済みである必要がある(`postProcess`コマンドを内部で呼ぶため)。

出力は`data/summary.json`(単一時刻)または`data/timehistory.csv`, `data/timehistory.png`(時刻歴)。

## 5. よくあるトラブル

| 症状 | 原因・対処 |
|---|---|
| `blockMesh: command not found` | OpenFOAM環境をsourceし忘れている(1.2節) |
| `chtMultiRegionFoam`が`Foam::error()`で止まる | Foundation版(openfoam.org)を使っている可能性。`echo $WM_PROJECT_VERSION`が`v2512`のような形式か確認(`13`等の単純な整数ならFoundation版) |
| `fistr1: command not found` | 2.3節のパス設定を忘れている、またはビルドが失敗している(`build/`内のエラーログを確認) |
| `postProcess -func writeCellCentres`に失敗 | OpenFOAM環境がsourceされていない、または101_0がまだ`Allrun`で1度も実行されていない(時刻ディレクトリが無い) |
| FrontISTRのビルドで`gfortran: command not found` | `sudo apt-get install gfortran`を忘れている |
