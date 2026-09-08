#!/bin/sh
# 103_0 の OpenFOAM ベースケースを 102_0 から用意する。
# メッシュ(polyMesh)は 102_0 の Allrun.pre で生成済みのものを再利用する
# (「メッシュは今あるものを使う」方針)。重いのでリポジトリには含めない。
#
# 使い方(このフォルダ=openfoam/ をカレントにして):
#   sh setup_base_case.sh
#
# 前提: ../../102_0_openfoam_hollow_cylinder_heat_transfer が
#       Allrun.pre 済みで constant/{fluid,solid}/polyMesh と 0/ を持っていること。

set -e
HERE=$(cd "$(dirname "$0")" && pwd)
SRC="$HERE/../../102_0_openfoam_hollow_cylinder_heat_transfer"
DST="$HERE/base_case"

if [ ! -d "$SRC/constant/solid/polyMesh" ] || [ ! -d "$SRC/0/solid" ]; then
    echo "ERROR: 102_0 のメッシュ/初期場が見つかりません。" >&2
    echo "先に $SRC で ./Allrun.pre を実行してください。" >&2
    exit 1
fi

rm -rf "$DST"
mkdir -p "$DST"
cp -r "$SRC/constant" "$DST/constant"
cp -r "$SRC/system"   "$DST/system"
cp -r "$SRC/0"        "$DST/0"
echo "base_case を用意しました: $DST"
echo "固体セル数: $(grep -o 'nCells:[0-9]*' "$DST/constant/solid/polyMesh/owner" | head -1)"
