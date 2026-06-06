"""
丸棒 OI デモ用の FrontISTR 連携層。

このモジュールは、円柱の六面体メッシュを作成し、1 ステップ分の
FrontISTR 入力ファイルを書き出し、`fistr1` を実行し、節点温度を読む。
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

FISTR_BIN = Path(os.environ.get("FISTR1", "/home/kamakiri/local/frontistr/bin/fistr1"))
from round_bar_mesh import build_round_bar_mesh


class FrontISTRInterface:
    """OI ループから使う、丸棒 FrontISTR の入出力ラッパー。"""
    def __init__(
        self,
        case_dir: Path,
        n_axial: int,
        length_m: float,
        radius_m: float,
        t_amb: float,
        dt: float,
        conductivity: float = 160.0,
        density: float = 2700.0,
        specific_heat: float = 900.0,
    ) -> None:
        self.case_dir = Path(case_dir)
        self.n_axial = n_axial
        self.length_m = length_m
        self.radius_m = radius_m
        self.t_amb = t_amb
        self.dt = dt
        self.conductivity = conductivity
        self.density = density
        self.specific_heat = specific_heat
        self.node_ids: list[int] = []
        self.node_x: np.ndarray | None = None

    @property
    def n_nodes(self) -> int:
        if not self.node_ids:
            self._ensure_mesh()
        return len(self.node_ids)

    def reset(self) -> None:
        """case ディレクトリ、メッシュ、固定の制御ファイルを作り直す。"""
        if self.case_dir.exists():
            shutil.rmtree(self.case_dir)
        self.case_dir.mkdir(parents=True)
        self._write_mesh()
        self._write_hecmw_ctrl()

    def run_step(self, values: np.ndarray, left_flux: float, right_temp: float) -> np.ndarray:
        """`.cnt` を書き、`fistr1` を実行し、最新の節点温度を返す。"""
        self._ensure_mesh()
        values = np.asarray(values, dtype=float)
        if values.shape[0] != self.n_nodes:
            raise ValueError(f"values length must be {self.n_nodes}, got {values.shape[0]}")

        self._cleanup_previous_outputs()
        self._write_cnt(values, left_flux=left_flux, right_temp=right_temp)

        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "1")
        result = subprocess.run(
            [str(FISTR_BIN)],
            cwd=self.case_dir,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        (self.case_dir / "log.fistr1").write_text(result.stdout)
        if result.returncode != 0:
            raise RuntimeError(f"FrontISTR failed. See {self.case_dir / 'log.fistr1'}")
        return self.read_result()

    def read_result(self) -> np.ndarray:
        """最新の FrontISTR 結果ファイルを読んで節点温度へ戻す。"""
        candidates = sorted(self.case_dir.glob("round_bar.res.0.*"))
        if not candidates:
            raise FileNotFoundError(f"FrontISTR result not found in {self.case_dir}")
        text = candidates[-1].read_text()
        vals = _parse_frontistr_temperature(text)
        return np.array([vals[node_id] for node_id in self.node_ids], dtype=float)

    def axial_profile(self, values: np.ndarray) -> np.ndarray:
        """節点温度を軸方向平均に畳み込み、描画や RMSE に使える形にする。"""
        self._ensure_mesh()
        x = self.node_x
        assert x is not None
        values = np.asarray(values, dtype=float)
        profile = np.zeros(self.n_axial)
        counts = np.zeros(self.n_axial)
        dx = self.length_m / self.n_axial
        for xi, vi in zip(x, values):
            idx = int(round(xi / dx))
            idx = max(0, min(self.n_axial, idx))
            # 節点の位置を、描画用にセル中心に近い軸方向ビンへ丸める。
            if idx == self.n_axial:
                bin_idx = self.n_axial - 1
            else:
                bin_idx = min(self.n_axial - 1, idx)
            profile[bin_idx] += vi
            counts[bin_idx] += 1
        return profile / np.maximum(counts, 1)

    def axial_positions(self) -> np.ndarray:
        """各軸方向ビン中心の x 座標を返す。"""
        return (np.arange(self.n_axial) + 0.5) * (self.length_m / self.n_axial)

    def axial_representative_node_indices(self) -> list[int]:
        """各軸方向ビンから代表 1 節点を選び、1D の測定位置にする。"""
        self._ensure_mesh()
        x = self.node_x
        assert x is not None
        indices: list[int] = []
        coords = self._node_coordinates()
        radial = np.sqrt(coords[:, 1] ** 2 + coords[:, 2] ** 2)
        for pos in self.axial_positions():
            score = np.abs(x - pos) + radial
            indices.append(int(np.argmin(score)))
        return indices

    def _ensure_mesh(self) -> None:
        """未初期化なら、ディスク上のメッシュを読み込む。"""
        if not self.node_ids:
            if not (self.case_dir / "round_bar.msh").exists():
                self.reset()
            else:
                self._load_mesh_nodes()

    def _cleanup_previous_outputs(self) -> None:
        """前回の FrontISTR 実行で残った結果やログを消す。"""
        for pattern in ("round_bar.res.*", "round_bar.vis*", "FSTR.*", "*.log", "0.log"):
            for path in self.case_dir.glob(pattern):
                if path.is_file():
                    path.unlink()

    def _write_hecmw_ctrl(self) -> None:
        """FrontISTR が case 入力を読むための固定 HEC-MW 制御ファイルを書く。"""
        (self.case_dir / "hecmw_ctrl.dat").write_text(
            """\
##
## HEC-MW control file for FrontISTR round bar DA
##
!MESH, NAME=fstrMSH, TYPE=HECMW-ENTIRE
round_bar.msh
!CONTROL, NAME=fstrCNT
round_bar.cnt
!RESULT, NAME=fstrRES, IO=OUT
round_bar.res
!RESULT, NAME=vis_out, IO=OUT
round_bar_vis
"""
        )

    def _write_cnt(self, values: np.ndarray, left_flux: float, right_temp: float) -> None:
        # `round_bar.cnt` は読みやすいテンプレートとしてリポジトリに置いてある。
        # 各 DA ステップでは `fistr1` 実行前にここを書き換える。
        # 具体的には、節点初期温度全体、右端の固定温度、左端の熱流束を上書きする。
        lines = [
            "!STEP,INCMAX=10000",
            "!SOLUTION,TYPE=HEAT",
            "!HEAT",
            f" {self.dt:.10g}, {self.dt:.10g}",
            "!INITIAL_CONDITION,TYPE=TEMPERATURE",
        ]
        for node_id, value in zip(self.node_ids, values):
            lines.append(f" {node_id}, {value:.10g}")
        lines.extend(
            [
                "!FIXTEMP",
                f" RIGHT, {right_temp:.10g}",
                "!DFLUX",
                f" ELEFT, S1, {left_flux:.10g}",
                "!SOLVER,METHOD=1,PRECOND=2,ITERLOG=NO,TIMELOG=NO",
                " 2000, 1, 10, 10",
                " 1.0e-8, 1.0, 0.0",
                "!WRITE,RESULT,FREQUENCY=1",
                "!WRITE,VISUAL",
                "!VISUAL, method=PSR",
                "!surface_num=1",
                "!surface 1",
                "!output_type = VTK",
                "!END",
            ]
        )
        (self.case_dir / "round_bar.cnt").write_text("\n".join(lines) + "\n")

    def _write_mesh(self) -> None:
        """丸棒の六面体メッシュを組み立て、FrontISTR 形式で書き出す。"""
        mesh = build_round_bar_mesh(self.n_axial, self.length_m, self.radius_m)
        self.node_ids = [node_id for node_id, _coord in mesh["nodes"]]
        self.node_x = np.array([coord[0] for _node_id, coord in mesh["nodes"]], dtype=float)

        lines = ["!HEADER", " FrontISTR round bar heat model", "!NODE"]
        for node_id, (x, y, z) in mesh["nodes"]:
            lines.append(f"{node_id:8d}, {x:.12g}, {y:.12g}, {z:.12g}")
        lines.append("!ELEMENT, TYPE=361")
        for elem_id, conn in mesh["elements"]:
            lines.append(f"{elem_id:8d}, " + ", ".join(f"{n:8d}" for n in conn))
        lines.append("!SECTION, TYPE=SOLID, EGRP=EALL, MATERIAL=M1")
        lines.append(" 1.0")
        lines.append("!MATERIAL, NAME=M1, ITEM=3")
        lines.append("!ITEM=1")
        lines.append(f" {self.density:.10g}")
        lines.append("!ITEM=2")
        lines.append(f" {self.specific_heat:.10g}")
        lines.append("!ITEM=3")
        lines.append(f" {self.conductivity:.10g}, 0.0")
        lines.append(f" {self.conductivity:.10g}, 1000.0")
        lines.append("!NGROUP, NGRP=NALL")
        lines.extend(_format_group_ids(self.node_ids))
        lines.append("!NGROUP, NGRP=LEFT")
        lines.extend(_format_group_ids(mesh["left_nodes"]))
        lines.append("!NGROUP, NGRP=RIGHT")
        lines.extend(_format_group_ids(mesh["right_nodes"]))
        lines.append("!EGROUP, EGRP=EALL")
        lines.extend(_format_group_ids([elem_id for elem_id, _conn in mesh["elements"]]))
        lines.append("!EGROUP, EGRP=ELEFT")
        lines.extend(_format_group_ids(mesh["left_elements"]))
        lines.append("!END")
        (self.case_dir / "round_bar.msh").write_text("\n".join(lines) + "\n")

    def _load_mesh_nodes(self) -> None:
        """既存の `.msh` から節点 ID と x 座標を読む。"""
        text = (self.case_dir / "round_bar.msh").read_text()
        nodes: list[tuple[int, float]] = []
        in_nodes = False
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("!NODE"):
                in_nodes = True
                continue
            if in_nodes and stripped.startswith("!"):
                break
            if in_nodes and stripped:
                parts = [p.strip() for p in stripped.split(",")]
                nodes.append((int(parts[0]), float(parts[1])))
        self.node_ids = [node_id for node_id, _x in nodes]
        self.node_x = np.array([x for _node_id, x in nodes], dtype=float)

    def _node_coordinates(self) -> np.ndarray:
        """メッシュ全節点の xyz 座標配列を返す。"""
        text = (self.case_dir / "round_bar.msh").read_text()
        coords: list[tuple[float, float, float]] = []
        in_nodes = False
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith("!NODE"):
                in_nodes = True
                continue
            if in_nodes and stripped.startswith("!"):
                break
            if in_nodes and stripped:
                parts = [p.strip() for p in stripped.split(",")]
                coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
        return np.array(coords, dtype=float)


def _format_group_ids(ids: list[int], per_line: int = 8) -> list[str]:
    """節点 ID / 要素 ID を FrontISTR の group 行に整形する。"""
    lines = []
    for i in range(0, len(ids), per_line):
        lines.append(" " + ", ".join(str(x) for x in ids[i:i + per_line]))
    return lines


def _parse_frontistr_temperature(text: str) -> dict[int, float]:
    """FrontISTR の結果ファイルを node_id -> 温度 の対応に変換する。"""
    marker = re.search(r"^\s*TEMPERATURE\s*$", text, flags=re.MULTILINE)
    if not marker:
        raise ValueError("TEMPERATURE section not found in FrontISTR result")
    tokens = text[marker.end():].split()
    values: dict[int, float] = {}
    i = 0
    while i + 1 < len(tokens):
        try:
            node_id = int(tokens[i])
            value = float(tokens[i + 1])
        except ValueError:
            break
        values[node_id] = value
        i += 2
    if not values:
        raise ValueError("No temperature values parsed from FrontISTR result")
    return values
