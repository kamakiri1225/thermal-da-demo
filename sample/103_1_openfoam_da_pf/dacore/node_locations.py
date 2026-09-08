"""ROM 5ノードの3次元座標(102_0 の固体プローブ位置に対応)と中空円筒の寸法.

観測点・未観測点の可視化(PyVista)と、将来の実 OpenFOAM 連成で
プローブ位置を一元管理するための単一の情報源。
座標系・寸法は 102_0 の geometry.yaml / system/solid/probes と一致させている。
"""

from __future__ import annotations

# 中空円筒(試験体)寸法 [m] — 102_0 geometry.yaml より
OUTER_RADIUS = 0.0375      # 外径 75 mm
INNER_RADIUS = 0.020       # 内径 40 mm
HEIGHT = 0.1005            # 高さ 100.5 mm (z: 0 .. 0.1005)
MID_WALL_RADIUS = 0.5 * (OUTER_RADIUS + INNER_RADIUS)  # 肉厚中心 28.75 mm

# ヒータ(外周 +X 側)— 102_0 の heaterSelection と一致
HEATER_AXIAL_HEIGHT = 0.050        # 軸方向 50 mm
HEATER_CENTER_Z = 0.05025          # 中心高さ
HEATER_ARC_LENGTH = 0.100          # 外周に沿った周方向長さ 100 mm

# ROM ノードの座標 [m]。hot/mid/cold/top は 102_0 の probe と同一。
# core は肉厚中心の代表(潜在)ノードで、102_0 の probe には無い。
NODE_XYZ = {
    "hot":  (0.032, 0.0, 0.05025),          # +X, ヒータ直下側 (102_0 T_hot)
    "mid":  (0.0, 0.024, 0.05025),          # +Y, 周方向90度 (102_0 T_mid)
    "cold": (-0.032, 0.0, 0.05025),         # -X, 反ヒータ側 (102_0 T_cold)
    "top":  (0.025, 0.0, 0.095),            # 上面近傍 (102_0 T_top)
    "core": (0.0, -MID_WALL_RADIUS, 0.05025),  # 肉厚中心の代表(潜在)ノード
}
