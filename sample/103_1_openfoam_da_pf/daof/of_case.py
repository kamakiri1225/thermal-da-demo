"""OpenFOAM メンバーケースのライフサイクル(生成・状態設定・実行・読み出し).

field-space データ同化の1メンバー = 1つの chtMultiRegionFoam ケース。
既存(102_0)メッシュを共有し、固体温度場の初期値とヒータ発熱 Q だけを変える。
"""

from __future__ import annotations

import os
import shutil
import subprocess

import numpy as np

from . import of_io

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BASE_CASE = os.path.join(ROOT, "openfoam", "base_case")


def solid_cell_centres():
    """固体セル中心 (N,3) を base_case の 0/solid/C から読む."""
    return of_io.read_internal_vector(os.path.join(BASE_CASE, "0", "solid", "C"))


def nearest_cells(coords, centres):
    """各観測座標に最も近い固体セルの index を返す."""
    idx = []
    for p in coords:
        d = np.linalg.norm(centres - np.asarray(p), axis=1)
        idx.append(int(d.argmin()))
    return np.array(idx, dtype=int)


def prepare_member(member_dir):
    """base_case を member_dir へ複製する(既存は削除)."""
    if os.path.exists(member_dir):
        shutil.rmtree(member_dir)
    shutil.copytree(BASE_CASE, member_dir, symlinks=True)


def set_solid_state(member_dir, time, T_values, Q):
    """指定時刻ディレクトリの 固体温度場 + ヒータ Q を書き込む."""
    tdir = os.path.join(member_dir, _tname(time), "solid")
    os.makedirs(tdir, exist_ok=True)
    of_io.write_solid_T(os.path.join(tdir, "T"), T_values, Q)


def read_solid_T(member_dir, time):
    """指定時刻の固体温度場を配列で読む."""
    return of_io.read_internal_scalar(
        os.path.join(member_dir, _tname(time), "solid", "T"))


def run_window(member_dir, t0, t1, logname=None):
    """chtMultiRegionFoam を [t0,t1] で直列実行する(ブロッキング)."""
    of_io.set_control_window(member_dir, t0, t1)
    log = os.path.join(member_dir, logname or f"log.run_{t0:g}_{t1:g}")
    with open(log, "w") as lf:
        p = subprocess.run(["chtMultiRegionFoam"], cwd=member_dir,
                           stdout=lf, stderr=subprocess.STDOUT)
    ok = p.returncode == 0 and os.path.isdir(os.path.join(member_dir, _tname(t1)))
    return ok, log


def _tname(t):
    """OpenFOAM 時刻ディレクトリ名(整数はそのまま, 小数は general 表記)."""
    if float(t) == int(t):
        return str(int(t))
    return ("%g" % float(t))


if __name__ == "__main__":
    # スモークテスト: 1メンバーを 0->2 s だけ回して場を読み戻す
    import time as _time
    C = solid_cell_centres()
    m = os.path.join(ROOT, "openfoam", "smoke_member")
    print("[smoke] prepare member ...")
    prepare_member(m)
    set_solid_state(m, 0.0, 300.0, Q=15.0)   # 一様 300 K から
    print("[smoke] run 0->2 s (radiation on, serial) ...")
    t0 = _time.time()
    ok, log = run_window(m, 0.0, 2.0)
    dt = _time.time() - t0
    print(f"[smoke] ok={ok} wall={dt:.0f}s log={os.path.relpath(log)}")
    if ok:
        T = read_solid_T(m, 2.0)
        print(f"[smoke] T@2s: min={T.min():.3f} max={T.max():.3f} mean={T.mean():.3f} K "
              f"(started uniform 300 K, heater warms +X side)")
