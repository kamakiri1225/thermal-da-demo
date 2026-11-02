"""保存済みの配置を、方式ごとに大きな文字の平面位置図へ出力する。

解析・選定は再実行しない。円筒の上面投影の模式図であり、CAEコンタ図ではない。
実行: MPLCONFIGDIR=/tmp/mpl-fisher python3 run/make_blog012_position_figs.py
"""
from pathlib import Path
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from dacore import plots as _plots  # 日本語フォントとAgg
import matplotlib.pyplot as plt

ITEMS = [
    ("目的指向（本研究）", "goal", "反りを目的に選んだ配置", "#bd352b"),
    ("よく動く点に置く", "magnitude", "感度の大きさだけで選んだ配置", "#7753a0"),
    ("Q-DEIM（QR列ピボット）", "qdeim", "Q-DEIMで選んだ配置", "#2863aa"),
    ("A最適（温度の trace 最小）", "aopt", "A最適で選んだ配置", "#21845a"),
    ("D最適（温度の det 最小）", "dopt", "D最適で選んだ配置", "#b86a18"),
    ("A・O 自身（循環・参考）", "direct", "A・Oを直接観測する参考配置", "#606060"),
]


def main():
    result_path = ROOT / "results/placement_method_comparison.json"
    data = json.loads(result_path.read_text())
    methods = data["methods"]
    # 過去の結果ファイルの表記にも対応
    for alias in ("熱感度 |w| 最大（素朴）", "熱感度の大きさ最大（素朴）"):
        if alias in methods and "よく動く点に置く" not in methods:
            methods["よく動く点に置く"] = methods[alias]
    img = ROOT / "docs/img"
    plt.rcParams.update({"font.size": 26, "axes.titlesize": 28})
    theta = np.linspace(0, 2 * np.pi, 300)
    gallery = [
        "# blog_012：選定方式ごとの変位センサ位置",
        "",
        "[本文へ戻る](blog_012_fisher_sensor_placement.md)",
        "",
        "1方式1画像に分けた上面投影の位置図です。温度計P2・P0は全方式で固定です。",
        "S1・S2は各方式内での表示名で、方式をまたいで同じ点を表す名前ではありません。",
        "円筒の平面位置を示す模式図であり、CAEの分布図ではありません。",
        "矢印は変位の測定方向です。Uzは紙面から手前（上方）を測ります。",
        "",
    ]
    for key, slug, title, color in ITEMS:
        if key not in methods:
            continue
        record = methods[key]
        fig, ax = plt.subplots(figsize=(10, 10))
        fig.subplots_adjust(left=.02, right=.98, bottom=.10, top=.84)
        for radius in (37.5, 20.0):
            ax.plot(radius*np.cos(theta), radius*np.sin(theta),
                    color="#8696a7", lw=2.5)
        arc = np.linspace(-np.radians(76), np.radians(76), 80)
        ax.plot(38.6*np.cos(arc), 38.6*np.sin(arc), color="#ba3b2a", lw=7)
        ax.text(47, 0, "ヒータ側", color="#ba3b2a", rotation=90,
                fontsize=26, ha="center", va="center")
        for x, name in ((28.7, "A"), (-28.7, "O")):
            ax.plot(x, 0, "o", color="black", ms=11, zorder=4)
            ax.text(x, 7, name, fontsize=28, weight="bold", ha="center")
        for j, (xyz, comp) in enumerate(zip(record["xyz_mm"], record["comp"])):
            x, y, z = xyz
            ax.plot(x, y, "s", color=color, ms=16, mec="white", mew=2, zorder=5)
            if comp in ("Ux", "Uy"):
                dx, dy = (12, 0) if comp == "Ux" else (0, 12)
                ax.annotate("", (x+dx, y+dy), (x, y),
                            arrowprops=dict(arrowstyle="-|>", color=color, lw=3))
            # 呼び出し線とラベルは図形の外へ。左右で分けて重なりを避ける。
            label_pos = (-34 if x < 0 else 34, 54 if y >= 0 else -54)
            points = record["xyz_mm"]
            if len(points) == 2 and (points[0][1] >= 0) == (points[1][1] >= 0):
                left = int(points[1][0] < points[0][0])
                label_pos = (-34 if j == left else 34, 54 if y >= 0 else -54)
            ax.annotate(f"S{j+1}：{comp}", xy=(x,y), xytext=label_pos,
                        fontsize=28, weight="bold", color=color,
                        ha="center", va="center",
                        bbox=dict(fc="white", ec="none", pad=2),
                        arrowprops=dict(arrowstyle="-", color=color, lw=2.2))
        # 穴の領域に座標方向を示す（変位方向の矢印との混同を避け灰色）。
        for dx, dy, name in ((14, 0, "x"), (0, 14, "y")):
            ax.annotate("", (-8+dx, -6+dy), (-8, -6),
                        arrowprops=dict(arrowstyle="->", color="#536171", lw=2))
            ax.text(-8+dx+3, -6+dy, name, fontsize=24, color="#536171", va="center")
        ax.set(xlim=(-59,59), ylim=(-63,63), aspect="equal")
        ax.axis("off")
        fig.suptitle(title, fontsize=28, color=color, y=.98)
        fig.text(.5,.91, "温度計P2・P0は固定／上面から見た位置", fontsize=24,
                 ha="center")
        fig.text(.5,.055, "■ 選んだ変位センサ　● 評価点A・O", fontsize=24, ha="center")
        fig.text(.5,.015, "Ux：x方向　Uy：y方向　Uz：上下方向", fontsize=24, ha="center")
        path = img / f"blog012_positions_{slug}.png"
        fig.savefig(path, dpi=140, facecolor="white")
        plt.close(fig)
        gallery += [
            f"## {title}", "",
            f"![{title}](img/{path.name})", "",
        ]
        for j,(xyz,comp) in enumerate(zip(record["xyz_mm"],record["comp"])):
            gallery += [f"- S{j+1}：{comp}、位置 {tuple(xyz)} mm"]
        gallery += ["", "位置と方向は保存済みの選定結果から読み込んでいます。", ""]
        print(path)
    (ROOT/"docs/blog012_placement_gallery.md").write_text("\n".join(gallery)+"\n")


if __name__ == "__main__":
    main()
