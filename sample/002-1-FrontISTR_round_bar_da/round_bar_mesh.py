"""
Round-bar mesh generator for the FrontISTR data-assimilation demo.

The mesh is a simple structured hexahedral approximation of a cylindrical rod:
one square center block and four outer sector blocks, extruded along the axial
direction. The layout is intentionally similar in spirit to the OpenFOAM
blockMesh setup used in the paired laplacianFoam example.
"""

from __future__ import annotations

import math

import numpy as np


def build_round_bar_mesh(n_axial: int, length: float, radius: float) -> dict[str, object]:
    inner = radius / 2.5
    s_in = inner / math.sqrt(2.0)

    inner_pts = [
        (-s_in, -s_in),
        (s_in, -s_in),
        (s_in, s_in),
        (-s_in, s_in),
    ]

    blocks = [
        ("center", inner_pts[0], inner_pts[1], inner_pts[2], inner_pts[3], 2, 2, None, None),
        ("sector", inner_pts[0], inner_pts[1], None, None, 2, 1, -135.0, -45.0),
        ("sector", inner_pts[1], inner_pts[2], None, None, 2, 1, -45.0, 45.0),
        ("sector", inner_pts[2], inner_pts[3], None, None, 2, 1, 45.0, 135.0),
        ("sector", inner_pts[3], inner_pts[0], None, None, 2, 1, 135.0, 225.0),
    ]

    node_map: dict[tuple[float, float, float], int] = {}
    nodes: list[tuple[int, tuple[float, float, float]]] = []
    elements: list[tuple[int, list[int]]] = []

    def get_node(x: float, y: float, z: float) -> int:
        key = (round(x, 12), round(y, 12), round(z, 12))
        if key not in node_map:
            node_map[key] = len(node_map) + 1
            nodes.append((node_map[key], key))
        return node_map[key]

    elem_id = 1
    left_elements: list[int] = []
    for kind, a, b, c, d, ni, nj, theta0, theta1 in blocks:
        grid: list[list[list[int]]] = []
        for k in range(n_axial + 1):
            x = length * k / n_axial
            plane: list[list[int]] = []
            for j in range(nj + 1):
                row: list[int] = []
                v = j / nj
                for i in range(ni + 1):
                    u = i / ni
                    if kind == "center":
                        assert c is not None and d is not None
                        y = (1 - u) * (1 - v) * a[0] + u * (1 - v) * b[0] + u * v * c[0] + (1 - u) * v * d[0]
                        z = (1 - u) * (1 - v) * a[1] + u * (1 - v) * b[1] + u * v * c[1] + (1 - u) * v * d[1]
                    else:
                        assert theta0 is not None and theta1 is not None
                        inner_y = (1 - u) * a[0] + u * b[0]
                        inner_z = (1 - u) * a[1] + u * b[1]
                        theta = math.radians((1 - u) * theta0 + u * theta1)
                        outer_y = radius * math.cos(theta)
                        outer_z = radius * math.sin(theta)
                        y = (1 - v) * inner_y + v * outer_y
                        z = (1 - v) * inner_z + v * outer_z
                    row.append(get_node(x, y, z))
                plane.append(row)
            grid.append(plane)
        left_elems_this_block: list[int] = []
        for k in range(n_axial):
            for j in range(nj):
                for i in range(ni):
                    n0 = grid[k][j][i]
                    n1 = grid[k][j][i + 1]
                    n2 = grid[k][j + 1][i + 1]
                    n3 = grid[k][j + 1][i]
                    n4 = grid[k + 1][j][i]
                    n5 = grid[k + 1][j][i + 1]
                    n6 = grid[k + 1][j + 1][i + 1]
                    n7 = grid[k + 1][j + 1][i]
                    conn = [n0, n1, n2, n3, n4, n5, n6, n7]
                    if _hex_orientation(conn, nodes) < 0.0:
                        conn = [n0, n3, n2, n1, n4, n7, n6, n5]
                    elements.append((elem_id, conn))
                    if k == 0:
                        left_elems_this_block.append(elem_id)
                    elem_id += 1
        left_elements.extend(left_elems_this_block)

    left_nodes = sorted([node_id for node_id, (x, _y, _z) in nodes if abs(x) < 1e-12])
    right_nodes = sorted([node_id for node_id, (x, _y, _z) in nodes if abs(x - length) < 1e-12])
    return {
        "nodes": nodes,
        "elements": elements,
        "left_nodes": left_nodes,
        "right_nodes": right_nodes,
        "left_elements": left_elements,
    }


def _hex_orientation(conn: list[int], nodes: list[tuple[int, tuple[float, float, float]]]) -> float:
    coord = {node_id: np.array(xyz, dtype=float) for node_id, xyz in nodes}
    p0 = coord[conn[0]]
    return float(
        np.linalg.det(
            np.column_stack(
                [
                    coord[conn[1]] - p0,
                    coord[conn[3]] - p0,
                    coord[conn[4]] - p0,
                ]
            )
        )
    )
