#!/bin/bash
# run.sh  (oi/)
# ==============
# 002-1_laplacian_da_round_bar/oi の実行スクリプト。
# 丸棒1D等価モデル + Optimal Interpolation によるデータ同化デモ。
#
# 実行方法:
#   cd sample/002-1_laplacian_da_round_bar/oi
#   chmod +x run.sh
#   ./run.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CASE_DIR="${SCRIPT_DIR}/case"
FOAM_BASHRC="/usr/lib/openfoam/openfoam2512/etc/bashrc"

echo "=========================================="
echo "  丸棒 laplacianFoam + OI データ同化デモ"
echo "=========================================="

if [ ! -f "$FOAM_BASHRC" ]; then
    echo "ERROR: OpenFOAM 環境が見つかりません: $FOAM_BASHRC"
    exit 1
fi

echo ""
echo "[1/2] blockMesh でメッシュを生成します..."
bash -c "source $FOAM_BASHRC && blockMesh -case $CASE_DIR"
echo "      メッシュ生成完了"

echo ""
echo "[2/2] Python DA スクリプト (OI) を実行します..."
echo "      (真値・DAなし・OIありの3本を実行します)"
python3 "${SCRIPT_DIR}/da_main.py"

echo ""
echo "=========================================="
echo "  完了。results/results_of_da.png を確認してください。"
echo "=========================================="
