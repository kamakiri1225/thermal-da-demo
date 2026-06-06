# Python 環境セットアップ

このディレクトリ以下のサンプルは Python 3.10 以上 (推奨: 3.12) を前提とします。

## 必要ライブラリ

| ライブラリ | 用途 |
|---|---|
| numpy | 数値計算 (行列演算・CSV入出力) |
| matplotlib | グラフ描画 |

## インストール手順

### 1. 仮想環境の作成（初回のみ）

```bash
cd "$(git rev-parse --show-toplevel)/sample"
python3 -m venv .venv
```

Windows (PowerShell) の場合:
```powershell
$repoRoot = git rev-parse --show-toplevel
Set-Location "$repoRoot\sample"
python -m venv .venv
```

### 2. 仮想環境の有効化

Linux / WSL:
```bash
source .venv/bin/activate
```

Windows (PowerShell):
```powershell
.venv\Scripts\Activate.ps1
```

### 3. ライブラリのインストール

```bash
pip install -r requirements.txt
```

### 4. インストール確認

```bash
python -c "import numpy, matplotlib; print('numpy:', numpy.__version__); print('matplotlib:', matplotlib.__version__)"
```

### 5. 仮想環境の無効化（使用後）

```bash
deactivate
```

## 各サンプルの実行方法

### OpenFOAM 丸棒 DA (002-1_laplacian_da_round_bar)

```bash
cd "$(git rev-parse --show-toplevel)/sample/002-1_laplacian_da_round_bar/oi"
python da_main.py          # データ同化実行
python plot_nodes_from_csv.py  # グラフ生成
```

### FrontISTR 丸棒 DA (002-1-FrontISTR_round_bar_da)

FrontISTR (`fistr1`) が別途必要です。
インストール方法は `002-1-FrontISTR_round_bar_da/INSTALL.md` を参照してください。

```bash
cd "$(git rev-parse --show-toplevel)/sample/002-1-FrontISTR_round_bar_da/oi"
python da_main.py          # データ同化実行
python plot_from_csv.py    # グラフ生成
```

## 注意事項

- OpenFOAM サンプルは WSL 上に OpenFOAM (laplacianFoam) がインストール済みであること
- FrontISTR サンプルは `fistr1` が PATH に通っているか `FISTR1` 環境変数で指定されていること
