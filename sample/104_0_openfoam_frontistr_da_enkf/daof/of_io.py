"""OpenFOAM フィールドファイルの読み書き(データ同化の状態⇄場の橋渡し).

field-space EnKF/PF では、各メンバーの固体温度場(全セル値)を状態ベクトルとして
読み出し、更新後に書き戻す必要がある。ここではその最小限の入出力を提供する:
    - read_internal_scalar : volScalarField の internalField(nonuniform)を配列で読む
    - read_internal_vector : volVectorField(セル中心 C など)を (N,3) で読む
    - write_solid_T        : 固体温度場ファイルを内部場+ヒータ発熱Qつきで書き出す
    - set_control_window   : controlDict の startTime/endTime/書き出しを1ウィンドウ用に設定
"""

from __future__ import annotations

import os
import re

import numpy as np

_HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| OpenFOAM field written by 103_0 field-space data assimilation driver         |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "{loc}";
    object      T;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 0 0 1 0 0 0 ];

"""

# 固体 T の boundaryField。ヒータ発熱 Q [W] は同化で推定するので都度書き換える。
_SOLID_BF = """boundaryField
{{
    solid_to_fluid
    {{
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        value           uniform 293.15;
        Tnbr            T;
        kappaMethod     solidThermo;
        kappa           none;
    }}
    heaterPower
    {{
        type            externalWallHeatFluxTemperature;
        mode            power;
        // ヒータ通電 0-300 s は推定 Q [W]、以降(300-600 s)は遮断(0 W)。
        // 102_0 / ROM 版と同じ加熱→冷却スケジュール。20 s 実行時は全区間 ON。
        Q               table
        (
            (0        {Q:.8g})
            (299.999  {Q:.8g})
            (300      0)
            (600      0)
        );
        kappaMethod     solidThermo;
        kappa           none;
        value           uniform 293.15;
    }}
}}


// ************************************************************************* //
"""


def read_internal_scalar(path, expect_n=None):
    """volScalarField の internalField を np.ndarray で返す(uniform/nonuniform 両対応)."""
    with open(path) as f:
        txt = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)\s*;",
                  txt, re.S)
    if m:
        n = int(m.group(1))
        vals = np.fromstring(m.group(2), sep=" ")
        if len(vals) != n:
            raise ValueError(f"{path}: expected {n} values, got {len(vals)}")
        return vals
    m = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)\s*;", txt)
    if m:
        if expect_n is None:
            raise ValueError(f"{path}: uniform field needs expect_n")
        return np.full(expect_n, float(m.group(1)))
    raise ValueError(f"{path}: could not parse internalField")


def read_internal_vector(path):
    """volVectorField(例: セル中心 C)の internalField を (N,3) で返す."""
    with open(path) as f:
        txt = f.read()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*(\d+)\s*\((.*?)\)\s*;",
                  txt, re.S)
    if not m:
        raise ValueError(f"{path}: could not parse vector internalField")
    n = int(m.group(1))
    body = m.group(2).replace("(", " ").replace(")", " ")
    vals = np.fromstring(body, sep=" ").reshape(-1, 3)
    if len(vals) != n:
        raise ValueError(f"{path}: expected {n} vectors, got {len(vals)}")
    return vals


def write_solid_T(path, values, Q, clip=(250.0, 400.0)):
    """固体温度場ファイルを書き出す.

    values : スカラー(uniform) または配列(nonuniform, セル数分)
    Q      : ヒータ発熱 [W]
    clip   : 物理的にありえない温度を丸める範囲 [K]
    """
    loc = os.path.basename(os.path.dirname(path)) or "solid"
    loc = os.path.join(os.path.basename(os.path.dirname(os.path.dirname(path))), "solid")
    out = _HEADER.format(loc=loc)
    if np.isscalar(values):
        out += f"internalField   uniform {float(values):.8g};\n\n"
    else:
        v = np.clip(np.asarray(values, dtype=float), *clip)
        body = "\n".join(f"{x:.8g}" for x in v)
        out += (f"internalField   nonuniform List<scalar>\n{len(v)}\n(\n"
                f"{body}\n)\n;\n\n")
    out += _SOLID_BF.format(Q=Q)
    with open(path, "w") as f:
        f.write(out)


def set_control_window(case_dir, t0, t1):
    """controlDict を1同化ウィンドウ [t0,t1] 用に設定する.

    startFrom startTime / startTime=t0 / endTime=t1 とし、t1 でのみ書き出す。
    """
    cd = os.path.join(case_dir, "system", "controlDict")
    with open(cd) as f:
        txt = f.read()
    txt = re.sub(r"startFrom\s+\w+;", "startFrom       startTime;", txt)
    txt = re.sub(r"startTime\s+[-\d.eE+]+;", f"startTime       {t0:g};", txt)
    txt = re.sub(r"stopAt\s+\w+;", "stopAt          endTime;", txt)
    txt = re.sub(r"endTime\s+[-\d.eE+]+;", f"endTime         {t1:g};", txt)
    # adjustTimeStep 使用時は adjustableRunTime にしないと書出/終了時刻を
    # 正確に踏めず、"2" ではなく "2.0016" のような時刻ディレクトリになる。
    txt = re.sub(r"writeControl\s+\w+;", "writeControl    adjustableRunTime;", txt)
    txt = re.sub(r"writeInterval\s+[-\d.eE+]+;",
                 f"writeInterval   {t1 - t0:g};", txt)
    txt = re.sub(r"purgeWrite\s+\d+;", "purgeWrite      0;", txt)
    with open(cd, "w") as f:
        f.write(txt)
