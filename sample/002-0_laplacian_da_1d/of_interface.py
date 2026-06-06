"""
of_interface.py
===============
OpenFOAM laplacianFoam ケースとの I/O インターフェース。

フィールドファイル（T）の読み書き、controlDict の更新、
ソルバーの起動を一元管理する。
外部ライブラリは使わず標準ライブラリと numpy のみで実装する。
"""

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

# T フィールドファイルのテンプレート
# Python の .format() で値を埋め込む（{n_cells}, {values}, {gradient}, {t_amb}）
_T_TEMPLATE = """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2512                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

internalField   nonuniform List<scalar>
{n_cells}
(
{values}
);

boundaryField
{{
    left
    {{
        type            fixedGradient;
        gradient        uniform {gradient};
    }}
    right
    {{
        type            fixedValue;
        value           uniform {t_amb};
    }}
    bottom
    {{
        type            empty;
    }}
    top
    {{
        type            empty;
    }}
    back
    {{
        type            empty;
    }}
    front
    {{
        type            empty;
    }}
}}

// ************************************************************************* //
"""

# OpenFOAM 環境設定スクリプト
_FOAM_BASHRC = "/usr/lib/openfoam/openfoam2512/etc/bashrc"


class OFInterface:
    """
    laplacianFoam ケースを1ステップずつ実行するインターフェース。

    DA ループから使う典型的な呼び出し順序:
        ofi = OFInterface(case_dir, n_cells=10, dt=10.0, t_amb=20.0)
        ofi.reset(initial_T=20.0)
        for k in range(N_STEPS):
            ofi.write_field(k*dt, x_current, gradient=g)
            ofi.run_step(k*dt, (k+1)*dt)
            x_pred = ofi.read_field((k+1)*dt)
            ...（KF更新）...
            x_current = x_kf

    Parameters
    ----------
    case_dir : Path
        OpenFOAM ケースのルートディレクトリ。
    n_cells : int
        x 方向のセル数（棒の節点数に対応）。
    dt : float
        時間刻み [s]。
    t_amb : float
        周囲温度（右端 fixedValue の値）[degC]。
    """

    def __init__(
        self,
        case_dir: Path,
        n_cells: int,
        dt: float,
        t_amb: float,
    ) -> None:
        self.case_dir = Path(case_dir)
        self.n_cells = n_cells
        self.dt = dt
        self.t_amb = t_amb

    def reset(self, initial_T: float) -> None:
        """0より大きい時刻ディレクトリをすべて削除して初期状態に戻す。"""
        for d in self.case_dir.iterdir():
            if d.is_dir() and _is_time_dir(d.name) and float(d.name) > 1e-12:
                shutil.rmtree(d)
        # 初期 T フィールドを書き直す
        self.write_field(0.0, np.ones(self.n_cells) * initial_T, gradient=0.0)

    def write_field(
        self,
        time: float,
        values: np.ndarray,
        gradient: float,
    ) -> None:
        """
        T フィールドを OpenFOAM ASCII 形式で書き込む。

        既存ファイルがあれば上書きする。
        """
        time_dir = self.case_dir / _format_time(time)
        time_dir.mkdir(exist_ok=True)

        values_str = "\n".join(f"{v:.10g}" for v in values)
        content = _T_TEMPLATE.format(
            n_cells=self.n_cells,
            values=values_str,
            gradient=f"{gradient:.10g}",
            t_amb=f"{self.t_amb:.10g}",
        )
        (time_dir / "T").write_text(content)

    def read_field(self, time: float) -> np.ndarray:
        """
        指定時刻の T フィールドの internalField を numpy 配列で返す。

        正確な時刻ディレクトリが見つからない場合は最も近い時刻を使う。
        """
        expected = self.case_dir / _format_time(time)
        if expected.is_dir():
            filepath = expected / "T"
        else:
            filepath = self._find_closest_field(time)

        text = filepath.read_text()
        return _parse_internal_field(text, self.n_cells)

    def run_step(self, start_time: float, end_time: float) -> None:
        """
        controlDict の endTime を更新してから laplacianFoam を実行する。

        startFrom は latestTime なので、write_field で書いた最新時刻から始まる。
        """
        self._update_end_time(end_time)
        self._run_solver()

    # ------------------------------------------------------------------ #
    # 内部メソッド                                                          #
    # ------------------------------------------------------------------ #

    def _update_end_time(self, end_time: float) -> None:
        """controlDict の endTime 行を書き換える。"""
        ctrl_path = self.case_dir / "system" / "controlDict"
        text = ctrl_path.read_text()
        # "endTime    <数値>;" を置換する
        text = re.sub(
            r"(endTime\s+)[0-9eE+\-.]+;",
            rf"\g<1>{end_time:.10g};",
            text,
        )
        ctrl_path.write_text(text)

    def _run_solver(self) -> None:
        """laplacianFoam をサブプロセスで実行する。"""
        log_path = self.case_dir / "log.laplacianFoam"
        cmd = (
            f"source {_FOAM_BASHRC} && "
            f"laplacianFoam -case {self.case_dir} "
            f"> {log_path} 2>&1"
        )
        result = subprocess.run(cmd, shell=True, executable="/bin/bash")
        if result.returncode != 0:
            raise RuntimeError(
                f"laplacianFoam が異常終了しました。"
                f"ログ: {log_path} を確認してください。"
            )

    def _find_closest_field(self, time: float) -> Path:
        """指定時刻に最も近い時刻ディレクトリの T ファイルを返す。"""
        candidates = [
            d for d in self.case_dir.iterdir()
            if d.is_dir() and _is_time_dir(d.name)
        ]
        if not candidates:
            raise FileNotFoundError(
                f"T フィールドが見つかりません (case={self.case_dir}, time={time})"
            )
        closest = min(candidates, key=lambda d: abs(float(d.name) - time))
        return closest / "T"


# ------------------------------------------------------------------ #
# モジュールレベルのユーティリティ関数                                   #
# ------------------------------------------------------------------ #

def _is_time_dir(name: str) -> bool:
    """ディレクトリ名が有効な数値（時刻）かを確認する。"""
    try:
        float(name)
        return True
    except ValueError:
        return False


def _format_time(time: float) -> str:
    """
    OpenFOAM の時刻ディレクトリ名形式に変換する。

    timePrecision 6 で general フォーマットを使うと、
    整数値に近ければ "10", "20" などの整数文字列になる。
    """
    if abs(time - round(time)) < 1e-8:
        return str(int(round(time)))
    return f"{time:.6g}"


def _parse_internal_field(text: str, n_cells: int) -> np.ndarray:
    """
    OpenFOAM ASCII フィールドファイルから internalField の値を読み取る。

    nonuniform List<scalar> と uniform value の両形式に対応する。
    """
    # nonuniform パターン: "internalField nonuniform List<scalar>\nN\n(\n値...\n);"
    m = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+\d+\s*\(\s*(.*?)\s*\)",
        text,
        re.DOTALL,
    )
    if m:
        nums = m.group(1).split()
        return np.array([float(x) for x in nums])

    # uniform パターン: "internalField uniform 20;"
    m = re.search(r"internalField\s+uniform\s+([0-9eE+\-.]+)", text)
    if m:
        return np.ones(n_cells) * float(m.group(1))

    raise ValueError(
        "internalField の解析に失敗しました。フィールドファイルの形式を確認してください。"
    )
