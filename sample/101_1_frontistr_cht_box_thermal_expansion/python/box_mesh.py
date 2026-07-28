"""
101_0の固体ブロック(既定: 200x200x400mm)と同じ寸法・同じ分割数の
構造六面体メッシュを作る。101_0のOpenFOAM固体メッシュ(cellSize基準の
一様グリッド)とノード位置を揃えることで、セル中心温度→節点温度の
写像を単純な最近傍平均で行えるようにしている。
"""

from __future__ import annotations

import numpy as np


def build_box_mesh(
    nx: int,
    ny: int,
    nz: int,
    lx: float,
    ly: float,
    lz: float,
) -> dict[str, object]:
    """原点(0,0,0)-(lx,ly,lz)の直方体を nx*ny*nz の構造六面体で分割する。

    ノードは (nx+1)*(ny+1)*(nz+1) 個、要素は nx*ny*nz 個(1次六面体、HEC-MW TYPE=361)。
    戻り値の "bottom_nodes" は z=0 面(底面)の節点ID一覧(FrontISTR側で
    「底面を固定」する境界条件に使う)。
    """
    node_ids = np.zeros((nx + 1, ny + 1, nz + 1), dtype=int)
    nodes: list[tuple[int, tuple[float, float, float]]] = []
    node_id = 1
    for k in range(nz + 1):
        z = lz * k / nz
        for j in range(ny + 1):
            y = ly * j / ny
            for i in range(nx + 1):
                x = lx * i / nx
                node_ids[i, j, k] = node_id
                nodes.append((node_id, (x, y, z)))
                node_id += 1

    elements: list[tuple[int, list[int]]] = []
    elem_id = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                n0 = int(node_ids[i, j, k])
                n1 = int(node_ids[i + 1, j, k])
                n2 = int(node_ids[i + 1, j + 1, k])
                n3 = int(node_ids[i, j + 1, k])
                n4 = int(node_ids[i, j, k + 1])
                n5 = int(node_ids[i + 1, j, k + 1])
                n6 = int(node_ids[i + 1, j + 1, k + 1])
                n7 = int(node_ids[i, j + 1, k + 1])
                elements.append((elem_id, [n0, n1, n2, n3, n4, n5, n6, n7]))
                elem_id += 1

    bottom_nodes = sorted(int(node_ids[i, j, 0]) for i in range(nx + 1) for j in range(ny + 1))
    top_nodes = sorted(int(node_ids[i, j, nz]) for i in range(nx + 1) for j in range(ny + 1))
    all_nodes = [nid for nid, _xyz in nodes]
    all_elements = [eid for eid, _conn in elements]

    return {
        "nodes": nodes,
        "elements": elements,
        "bottom_nodes": bottom_nodes,
        "top_nodes": top_nodes,
        "all_nodes": all_nodes,
        "all_elements": all_elements,
    }


def cell_centers(nx: int, ny: int, nz: int, lx: float, ly: float, lz: float) -> np.ndarray:
    """一様構造グリッドのセル中心座標を(nx*ny*nz, 3)で返す(参考・検証用)。
    実際の温度写像は101_0側のOpenFOAM出力(C, T)を直接読むため、
    本関数はメッシュ整合性の確認・単体テスト用途。"""
    dx, dy, dz = lx / nx, ly / ny, lz / nz
    pts = []
    for k in range(nz):
        z = (k + 0.5) * dz
        for j in range(ny):
            y = (j + 0.5) * dy
            for i in range(nx):
                x = (i + 0.5) * dx
                pts.append((x, y, z))
    return np.array(pts, dtype=float)
