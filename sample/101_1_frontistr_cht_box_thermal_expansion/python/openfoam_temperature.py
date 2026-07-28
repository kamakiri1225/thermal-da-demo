"""
101_0(OpenFOAM chtMultiRegionFoam)の固体側の温度分布を読み込み、
FrontISTRメッシュの節点座標へ温度を写像する。

101_0では固体ブロックが "solid" と "heaterMat" の2リージョンに
分かれている(ヒートマットを独立リージョン化したため、
system/heaterMat/changeDictionaryDict参照)。101_1側は両方を1つの
連続した構造物として扱うため、"solid"だけでなく"heaterMat"の
セル中心・温度も合わせて読み込み、まとめて節点温度の写像に使う。

前提:
  101_0側で各リージョンについて
  `postProcess -func writeCellCentres -region <region> -time <T>` を
  実行済みで、`<T>/<region>/C` (セル中心座標) と `<T>/<region>/T` (温度) が
  同じセル順序で存在すること(OpenFOAMは同一メッシュ上のフィールドを
  常に同じセル順序で書き出すため、この2つを素直にzipすれば
  (座標, 温度)の対応が取れる)。

写像方法: 各FrontISTR節点について、最近傍k個のOpenFOAMセル中心の
逆距離加重平均(inverse distance weighting)を取る。101_1のメッシュ分割数を
101_0のcellSizeと揃えていれば、内部節点はほぼ隣接8セルの平均に一致する。

座標系: 101_0の固体は流体領域中央(既定ではx,y=0.4--0.6m)にある一方、
101_1の構造メッシュは原点始まり(x,y=0--0.2m)である。このため補間前に
両点群のバウンディングボックス中心を一致させる平行移動を行う。回転や
スケーリングは行わないため、形状寸法と軸方向は両モデルで一致している必要がある。
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np


def _read_of_scalar_field(path: Path) -> np.ndarray:
    """OpenFOAMのvolScalarField(uniform/nonuniform両対応)のinternalFieldを読む。"""
    text = path.read_text()
    m = re.search(r"internalField\s+uniform\s+([\-0-9.eE+]+)\s*;", text)
    if m:
        # uniform: セル数が分からないので呼び出し側でブロードキャストする
        return np.array([float(m.group(1))])

    # OpenFOAMはセル数が少ないと ``8(1 2 ...)`` の1行形式、通常は
    # ``8\n(\n1\n2\n...)`` の複数行形式で書くため、空白・改行を限定しない。
    m = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if not m:
        raise ValueError(f"internalFieldを解釈できません: {path}")
    n = int(m.group(1))
    vals = np.array([float(x) for x in m.group(2).split()], dtype=float)
    if len(vals) != n:
        raise ValueError(f"{path}: 宣言セル数{n}と実データ数{len(vals)}が不一致")
    return vals


def _read_of_vector_field(path: Path) -> np.ndarray:
    """OpenFOAMのvolVectorField(セル中心座標C等)のinternalFieldを(N,3)で読む。"""
    text = path.read_text()
    m = re.search(
        r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if not m:
        raise ValueError(f"internalField(vector)を解釈できません: {path}")
    n = int(m.group(1))
    body = m.group(2)
    triples = re.findall(r"\(([^()]+)\)", body)
    vecs = np.array([[float(x) for x in t.split()] for t in triples], dtype=float)
    if len(vecs) != n:
        raise ValueError(f"{path}: 宣言セル数{n}と実データ数{len(vecs)}が不一致")
    return vecs


SOLID_FAMILY_REGIONS = ("solid", "heaterMat")


def _load_region_cell_temperatures(of_case_dir: Path, time_dir: str, region: str) -> tuple[np.ndarray, np.ndarray]:
    """1リージョン分のセル中心座標(N,3)と温度(N,)を読む。"""
    t_path = Path(of_case_dir) / time_dir / region / "T"
    c_path = Path(of_case_dir) / time_dir / region / "C"
    if not c_path.exists():
        raise FileNotFoundError(
            f"{c_path} がありません。先に以下を実行してください:\n"
            f"  postProcess -func writeCellCentres -region {region} -time {time_dir} "
            f"-case {of_case_dir}"
        )
    centers = _read_of_vector_field(c_path)
    temps = _read_of_scalar_field(t_path)
    if len(temps) == 1 and len(centers) > 1:
        temps = np.full(len(centers), temps[0])
    if len(temps) != len(centers):
        raise ValueError(f"{region}: セル数不一致: T={len(temps)} C={len(centers)}")
    return centers, temps


def load_solid_cell_temperatures(
    of_case_dir: Path, time_dir: str, regions: tuple[str, ...] = SOLID_FAMILY_REGIONS
) -> tuple[np.ndarray, np.ndarray]:
    """101_0の指定時刻ディレクトリから、固体側(solid+heaterMat)全体の
    セル中心座標(N,3)と温度(N,)をまとめて読む。"""
    all_centers = []
    all_temps = []
    for region in regions:
        region_dir = Path(of_case_dir) / time_dir / region
        if not region_dir.is_dir():
            continue  # heaterMatが存在しない旧構成のケースにも対応
        centers, temps = _load_region_cell_temperatures(of_case_dir, time_dir, region)
        all_centers.append(centers)
        all_temps.append(temps)
    if not all_centers:
        raise FileNotFoundError(f"{regions} のいずれのリージョンも {of_case_dir}/{time_dir} に見つかりません。")
    return np.concatenate(all_centers, axis=0), np.concatenate(all_temps, axis=0)


def align_cell_centers_to_node_mesh(
    cell_centers: np.ndarray,
    node_coords: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """OpenFOAMセル中心をFrontISTR節点メッシュの座標系へ平行移動する。

    両モデルのバウンディングボックス中心の差を移動量とする。セル中心点群は
    境界面上に点を持たないが、対称な構造格子であれば中心位置は形状中心と一致
    する。戻り値は(移動後セル中心, 適用した移動ベクトル[m])。
    """
    cell_centers = np.asarray(cell_centers, dtype=float)
    node_coords = np.asarray(node_coords, dtype=float)
    if cell_centers.ndim != 2 or cell_centers.shape[1] != 3:
        raise ValueError("cell_centers must have shape (N, 3)")
    if node_coords.ndim != 2 or node_coords.shape[1] != 3:
        raise ValueError("node_coords must have shape (N, 3)")

    source_center = 0.5 * (cell_centers.min(axis=0) + cell_centers.max(axis=0))
    target_center = 0.5 * (node_coords.min(axis=0) + node_coords.max(axis=0))
    translation = target_center - source_center
    return cell_centers + translation, translation


def interpolate_to_nodes(
    node_coords: np.ndarray,
    cell_centers: np.ndarray,
    cell_temperatures: np.ndarray,
    k: int = 8,
) -> np.ndarray:
    """各節点座標へ、最近傍k個のセル中心からの逆距離加重平均で温度を写像する。"""
    node_coords = np.asarray(node_coords, dtype=float)
    n_nodes = node_coords.shape[0]
    result = np.zeros(n_nodes, dtype=float)

    k_eff = min(k, len(cell_centers))
    for idx in range(n_nodes):
        d = np.linalg.norm(cell_centers - node_coords[idx], axis=1)
        nearest = np.argsort(d)[:k_eff]
        dn = d[nearest]
        if dn[0] < 1e-9:
            result[idx] = cell_temperatures[nearest[0]]
            continue
        w = 1.0 / dn
        result[idx] = float(np.sum(w * cell_temperatures[nearest]) / np.sum(w))
    return result
