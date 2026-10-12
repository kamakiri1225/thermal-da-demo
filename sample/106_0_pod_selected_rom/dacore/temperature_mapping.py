"""OpenFOAM cell temperatures -> FEM nodal temperatures (8-neighbour IDW).

Coordinates must use the same frame and units. The mapping is geometric;
Q-DEIM representative cells and gappy-POD reconstruction are separate steps.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

MAPPING_METHOD = "idw_k8_p1_v1"


def build_idw_mapping(cell_centres, fem_coords, k=8):
    """Return M (nodes x cells), with nonnegative weights and row sums 1.

    Weights are proportional to 1/distance. Coincident points use only the
    coincident cell(s), avoiding division by zero and preserving their value.
    """
    source = np.asarray(cell_centres, dtype=float)
    target = np.asarray(fem_coords, dtype=float)
    if source.ndim != 2 or target.ndim != 2 or source.shape[1] != 3 or target.shape[1] != 3:
        raise ValueError("cell_centres and fem_coords must have shape (n, 3)")
    if not len(source) or not np.isfinite(source).all() or not np.isfinite(target).all():
        raise ValueError("coordinates must be finite and the source nonempty")
    if k < 1:
        raise ValueError("k must be positive")
    neighbours = min(int(k), len(source))
    distances, indices = cKDTree(source).query(target, k=neighbours)
    distances = np.asarray(distances).reshape(len(target), neighbours)
    indices = np.asarray(indices).reshape(len(target), neighbours)
    # Scale each row before normalising to avoid overflow for small distances.
    exact = distances == 0
    weights = np.zeros_like(distances)
    coincident_rows = exact.any(axis=1)
    weights[coincident_rows] = exact[coincident_rows] / exact[coincident_rows].sum(axis=1, keepdims=True)
    other = ~coincident_rows
    weights[other] = distances[other, :1] / distances[other]
    weights[other] /= weights[other].sum(axis=1, keepdims=True)
    rows = np.repeat(np.arange(len(target)), neighbours)
    return csr_matrix((weights.ravel(), (rows, indices.ravel())), shape=(len(target), len(source)))


def map_temperature(mapping, values, axis=0):
    """Apply M along the cell axis; support fields, modes, and time histories."""
    moved = np.moveaxis(np.asarray(values), axis, 0)
    if moved.shape[0] != mapping.shape[1]:
        raise ValueError("temperature cell axis does not match mapping")
    mapped = mapping @ moved.reshape(moved.shape[0], -1)
    mapped = np.asarray(mapped).reshape((mapping.shape[0],) + moved.shape[1:])
    return np.moveaxis(mapped, 0, axis)


def mapping_cache_valid(path):
    """Never reuse an untagged nearest-neighbour cache as an IDW result."""
    if not Path(path).exists():
        return False
    with np.load(path) as data:
        return "mapping_method" in data and str(data["mapping_method"]) == MAPPING_METHOD


def mapping_fingerprint(*arrays):
    """Record the inputs used to form a mapping/operator for reproducibility."""
    digest = hashlib.sha256(MAPPING_METHOD.encode())
    for array in arrays:
        values = np.ascontiguousarray(array, dtype=np.float64)
        digest.update(str(values.shape).encode())
        digest.update(values.tobytes())
    return digest.hexdigest()
