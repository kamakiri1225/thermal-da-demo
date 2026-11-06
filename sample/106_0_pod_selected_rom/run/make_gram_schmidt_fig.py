"""blog_013 §3-4「重複する成分を除く」の説明図（2 次元の模式図）。
  python3 run/make_gram_schmidt_fig.py → docs/img/gram_schmidt.png
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 16})
fig, axs = plt.subplots(1, 3, figsize=(16, 5.4), gridspec_kw=dict(wspace=0.12))
v1 = np.array([1.0, 0.0]); w = np.array([0.97, 0.26])
arrow = lambda ax, p, q, c, lw=3: ax.annotate("", xy=q, xytext=p, arrowprops=dict(arrowstyle="-|>", color=c, lw=lw, mutation_scale=22))
for ax in axs:
    ax.set_xlim(-0.15, 1.25); ax.set_ylim(-0.15, 0.55); ax.set_aspect("equal"); ax.axis("off")
ax = axs[0]; ax.set_title("① 候補 2 本はほぼ同じ向き", fontsize=18)
arrow(ax, (0, 0), v1, "#1f77b4"); arrow(ax, (0, 0), w, "#d62728")
ax.text(1.02, -0.07, "$\\mathbf{v}_1$（1 本目）", color="#1f77b4", fontsize=16); ax.text(0.65, 0.3, "$\\mathbf{w}$（次の候補）", color="#d62728", fontsize=16)
ax = axs[1]; ax.set_title("② $\\mathbf{v}_1$ の向きの成分（影）を引く", fontsize=18)
arrow(ax, (0, 0), v1, "#1f77b4"); arrow(ax, (0, 0), w, "#d62728")
ax.plot([w[0], w[0]], [0, w[1]], ":", color="k", lw=1.5)
arrow(ax, (0, 0), (w[0], 0), "#9467bd", 5)
ax.text(0.15, -0.1, "影 $= \\mathbf{v}_1(\\mathbf{v}_1^{T}C\\mathbf{w})$", color="#9467bd", fontsize=16)
arrow(ax, (w[0], 0), w, "#2ca02c", 4); ax.text(1.0, 0.1, "残り", color="#2ca02c", fontsize=16)
ax = axs[2]; ax.set_title("③ 残りを長さ 1 にして 2 本目に", fontsize=18)
arrow(ax, (0, 0), v1, "#1f77b4"); arrow(ax, (0, 0), (0, 0.45), "#2ca02c")
ax.text(1.02, -0.07, "$\\mathbf{v}_1$", color="#1f77b4", fontsize=17); ax.text(0.04, 0.42, "$\\mathbf{v}_2$（新しい向き）", color="#2ca02c", fontsize=16)
fig.savefig(os.path.join(ROOT, "docs", "img", "gram_schmidt.png"), dpi=100, bbox_inches="tight", facecolor="white")
print("wrote gram_schmidt.png")
