# FrontISTR install notes

この環境では `fistr1` が見つかりませんでした。

また、Codex から `sudo apt update` を実行するとパスワード入力が必要になり、非対話実行では進められませんでした。そのため、以下の手順は WSL のターミナルでユーザー側が実行してください。

```bash
command -v fistr1
```

が何も返さない場合、FrontISTR をインストールして PATH を通してください。

## WSL / Ubuntu での推奨手順

### 1. apt で候補を確認する

まず apt のパッケージ候補を確認します。

```bash
sudo apt update
apt-cache search frontistr
apt-cache search hecmw
```

`frontistr` または `hecmw` 関連のパッケージが表示された場合は、そのパッケージを入れます。

例:

```bash
sudo apt install frontistr
```

または、検索結果に `frontistr` 以外のパッケージ名が出た場合:

```bash
sudo apt install <検索で出たFrontISTRパッケージ名>
```

インストール後に確認します。

```bash
which fistr1
fistr1 -v
```

### 2. apt に無い場合: 依存パッケージを入れる

apt に FrontISTR 本体が無い場合は、ソースビルド用の依存関係を入れます。

```bash
sudo apt install -y \
  build-essential \
  git \
  cmake \
  gfortran \
  openmpi-bin \
  libopenmpi-dev \
  libmetis-dev \
  libmumps-dev \
  libscalapack-openmpi-dev \
  libblas-dev \
  liblapack-dev
```

その後、FrontISTR を取得します。

```bash
mkdir -p ~/src
cd ~/src
git clone https://github.com/FrontISTR/FrontISTR.git
cd FrontISTR
```

ビルド方法は FrontISTR のバージョンで変わることがあるため、まず公式 README を確認してください。

```bash
ls
sed -n '1,220p' README.md
```

典型的には CMake ビルドになります。README に CMake 手順がある場合は、概ね次の形です。

```bash
mkdir -p build
cd build
cmake ..
make -j"$(nproc)"
sudo make install
```

インストール後に確認します。

```bash
which fistr1
fistr1 -v
```

`sudo make install` しない場合は、生成された `fistr1` の場所を PATH に追加します。

例:

```bash
export PATH="$HOME/src/FrontISTR/build:$PATH"
which fistr1
```

この設定を永続化する場合:

```bash
echo 'export PATH="$HOME/src/FrontISTR/build:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

## Codex に戻す前の確認

WSL ターミナルで次を実行し、結果を確認してください。

```bash
command -v fistr1
fistr1 -v
```

`command -v fistr1` がパスを返せば、Codex 側でも続きの実装確認に進めます。

## 今回のビルド結果

この環境では apt の `frontistr` パッケージは見つからなかったため、公式リポジトリから CMake でビルドした。

```bash
cd ~/src/FrontISTR
cmake -S . -B build-codex \
  -DCMAKE_BUILD_TYPE=RELEASE \
  -DCMAKE_INSTALL_PREFIX=$HOME/local/frontistr \
  -DWITH_MPI=OFF \
  -DWITH_TOOLS=OFF \
  -DWITH_MUMPS=OFF \
  -DWITH_METIS=OFF \
  -DWITH_REFINER=OFF \
  -DWITH_REVOCAP=OFF \
  -DWITH_ML=OFF \
  -DWITH_NETCDF=OFF
cmake --build build-codex -j2
cmake --install build-codex
```

生成された実行ファイル:

```bash
/home/kamakiri/local/frontistr/bin/fistr1
```

確認結果:

```bash
/home/kamakiri/local/frontistr/bin/fistr1 -v
```

FrontISTR 5.9、MPI disabled、OpenMP enabled、LAPACK enabled として起動した。

## このケースで必要な確認

インストール後、このリポジトリのルートで次を確認します。

```bash
command -v fistr1
fistr1 -v
```

`fistr1` が使えるようになったら、丸棒熱伝導ケースの `.msh` と `.cnt` を作り、OpenFOAM 版と同じ OI ループに接続します。
