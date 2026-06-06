# CalculiX install notes

この環境では `ccx` が見つかりませんでした。

```bash
command -v ccx
```

が何も返さない場合、CalculiX をインストールして PATH を通してください。

## WSL / Ubuntu での候補

### 1. apt で入れる

Ubuntu 系では次のパッケージ名で入ることがあります。

```bash
sudo apt update
sudo apt install calculix-ccx calculix-cgx
```

インストール後に確認します。

```bash
which ccx
ccx -v
which cgx
```

`ccx` はソルバ、`cgx` はプリポスト用です。このデータ同化ケースではまず `ccx` が必要です。

### 2. conda-forge を使う場合

conda 環境がある場合は、conda-forge に CalculiX パッケージがあるか確認します。

```bash
conda search -c conda-forge calculix
conda create -n calculix -c conda-forge calculix
conda activate calculix
which ccx
```

### 3. Windows バイナリを使う場合

Windows 版 CalculiX を使う場合は、`ccx.exe` のあるフォルダを Windows の PATH に追加します。

WSL から Windows 側の `ccx.exe` を呼ぶ場合は、パスや改行コード、出力ファイルの場所で問題が出やすいため、まずは WSL 内に Linux 版 `ccx` を入れる方が扱いやすいです。

## このケースで必要な確認

インストール後、このリポジトリのルートで次を確認します。

```bash
command -v ccx
ccx -v
```

`ccx` が使えるようになったら、丸棒熱伝導 `.inp` を作り、OpenFOAM 版と同じ OI ループに接続します。
