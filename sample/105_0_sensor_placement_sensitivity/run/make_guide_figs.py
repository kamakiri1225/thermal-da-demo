"""docs/03_sensor_design_guide.md 用の解説図を描く.

出力:
  docs/img/guide_w_matrix.png  … W=K⁻¹H とは何か(行列の行/列の読み方)
  docs/img/guide_workflow.png  … 実験前センサ設計のワークフロー
  docs/img/guide_snr.png       … なぜ低感度点の観測は役に立たない(害になる)か
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
from dacore import plots as _p   # 日本語フォント登録
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
IMG=os.path.join(ROOT,"docs","img")


def fig_w_matrix():
    fig,ax=plt.subplots(figsize=(14,7.4)); ax.set_xlim(0,14); ax.set_ylim(0,7.4); ax.axis("off")
    ax.text(7,7.05,"感度行列 $W=K^{-1}H$ の読み方 —— 1つの行列に2つの使い道",
            ha="center",fontsize=17,weight="bold")

    # 式: u = W ΔT
    # u ベクトル
    ax.add_patch(Rectangle((1.0,1.6),0.8,3.6,fc="#fff2cc",ec="k"))
    ax.text(1.4,5.45,"変位\n$u$ [µm]",ha="center",fontsize=13)
    ax.text(1.4,1.2,"(全節点の変位)",ha="center",fontsize=10,color="gray")
    ax.text(2.35,3.4,"=",ha="center",fontsize=22)
    # W 行列
    ax.add_patch(Rectangle((2.9,1.6),4.6,3.6,fc="#eeeeee",ec="k"))
    ax.text(5.2,5.45,"感度行列 $W=K^{-1}H$",ha="center",fontsize=13)
    # 行ハイライト(赤)
    ax.add_patch(Rectangle((2.9,3.6),4.6,0.42,fc="#f4b6b6",ec="crimson",lw=2))
    ax.text(5.2,3.81,"行",ha="center",fontsize=12,color="crimson",weight="bold")
    # 列ハイライト(青)
    ax.add_patch(Rectangle((5.6,1.6),0.42,3.6,fc="#b6cdf4",ec="royalblue",lw=2,alpha=0.85))
    ax.text(5.81,1.32,"列",ha="center",fontsize=12,color="royalblue",weight="bold")
    ax.text(7.95,3.4,"×",ha="center",fontsize=20)
    # ΔT ベクトル
    ax.add_patch(Rectangle((8.45,1.6),0.8,3.6,fc="#dcefdc",ec="k"))
    ax.text(8.85,5.45,"温度変化\n$\\Delta T$ [K]",ha="center",fontsize=13)
    ax.text(8.85,1.2,"(全節点の温度)",ha="center",fontsize=10,color="gray")

    # 右側: 行/列の説明
    ax.add_patch(FancyBboxPatch((10.0,3.9),3.6,1.9,boxstyle="round,pad=0.12",
                                fc="#fdecec",ec="crimson",lw=1.5))
    ax.text(11.8,5.35,"行の読み方(赤)",ha="center",fontsize=13,color="crimson",weight="bold")
    ax.text(11.8,4.6,"「この点の変位は温度場に\nどれだけ敏感か」\n→ 変位計をどこに貼るか",
            ha="center",fontsize=11.5)
    ax.add_patch(FancyBboxPatch((10.0,1.3),3.6,1.9,boxstyle="round,pad=0.12",
                                fc="#eaf0fb",ec="royalblue",lw=1.5))
    ax.text(11.8,2.75,"列の読み方(青)",ha="center",fontsize=13,color="royalblue",weight="bold")
    ax.text(11.8,2.0,"「この場所の温度が動くと\n観測変位はどれだけ動くか」\n→ 温度を管理すべき場所",
            ha="center",fontsize=11.5)
    ax.annotate("",xy=(10.0,4.85),xytext=(7.5,3.81),
                arrowprops=dict(arrowstyle="->",color="crimson",lw=2))
    ax.annotate("",xy=(10.0,2.25),xytext=(5.81,1.6),
                arrowprops=dict(arrowstyle="->",color="royalblue",lw=2))

    ax.text(7,0.45,"Wは 形状・材料・拘束だけで決まる（ヒータ位置や運転条件に依存しない）"
            "→ 装置を作る前に計算できる",ha="center",fontsize=13.5,
            bbox=dict(boxstyle="round",fc="#fffbe6",ec="orange"))
    fig.tight_layout()
    fig.savefig(os.path.join(IMG,"guide_w_matrix.png"),dpi=140); plt.close(fig)
    print("[guide] guide_w_matrix.png")


def fig_workflow():
    fig,ax=plt.subplots(figsize=(14,5.6)); ax.set_xlim(0,14); ax.set_ylim(0,5.6); ax.axis("off")
    ax.text(7,5.25,"実験前センサ配置設計のワークフロー（すべて装置を作る前にできる）",
            ha="center",fontsize=16,weight="bold")
    steps=[
        ("① メッシュ準備","CAD→四面体341\n(六面体は6分割\nテトラ化)","#e8f0fe"),
        ("② DUMPW実行","パッチ版fistr1 1回\nW_diffとK,Hを出力\n(約15秒)","#e6f4ea"),
        ("③ W全体を構築","K,Hダンプから\n$W=K^{-1}H$\n(約2分)","#e6f4ea"),
        ("④ 選点","行感度マップから\n変位計の位置決定\n(固定近傍は除外)","#fdecec"),
        ("⑤ 事前リハーサル","データ同化(EnKF)を\n仮想実験で回して\n精度を確認","#fff4e5"),
        ("⑥ 装置製作・実験","最適配置で\nセンサを貼る","#f3e8fd"),
    ]
    x=0.4
    for i,(t,d,c) in enumerate(steps):
        ax.add_patch(FancyBboxPatch((x,1.6),1.95,2.6,boxstyle="round,pad=0.1",fc=c,ec="k"))
        ax.text(x+0.975,3.75,t,ha="center",fontsize=12.5,weight="bold")
        ax.text(x+0.975,2.6,d,ha="center",fontsize=10.5)
        if i<5:
            ax.annotate("",xy=(x+2.35,2.9),xytext=(x+1.98,2.9),
                        arrowprops=dict(arrowstyle="-|>",color="k",lw=2))
        x+=2.28
    ax.annotate("",xy=(2.0,1.25),xytext=(11.6,1.25),
                arrowprops=dict(arrowstyle="-",color="orange",lw=3))
    ax.text(6.8,0.7,"①〜⑤ = 計算だけ（数分）。実験のやり直し・貼り直しのコストがゼロになる",
            ha="center",fontsize=13.5,color="darkorange",weight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(IMG,"guide_workflow.png"),dpi=140); plt.close(fig)
    print("[guide] guide_workflow.png")


def fig_snr():
    kv=np.load(os.path.join(ROOT,"results","kinvh_sensitivity.npz"))
    hi_s=float(kv["row_sens"][kv["hi"][0]]); lo_s=float(kv["row_sens"][kv["lo"][0]])
    noise=0.3
    fig,axes=plt.subplots(1,2,figsize=(13.5,5.6))
    # 左: 信号とノイズの比較(棒, log)
    ax=axes[0]
    bars=ax.bar(["高感度点\n(HIGH1)","低感度点\n(LOW1)"],[hi_s,lo_s],
                color=["tab:red","tab:blue"],width=0.5)
    ax.axhline(noise,color="k",ls="--",lw=2)
    ax.text(1.32,noise*1.15,"観測ノイズ σ=0.3µm",fontsize=12,ha="right")
    ax.set_yscale("log"); ax.set_ylim(0.05,hi_s*6)
    ax.set_ylabel("温度場+1Kあたりの変位応答 [µm/K]",fontsize=12)
    ax.set_title("信号の大きさ vs ノイズ（実測値）",fontsize=14,pad=14)
    ax.text(0,hi_s*1.25,f"{hi_s:.1f} µm/K\n(ノイズの{hi_s/noise:.0f}倍)",ha="center",
            fontsize=12,color="tab:red",weight="bold")
    ax.text(1,lo_s*1.3,f"{lo_s:.2f} µm/K\n(ノイズと同等以下)",ha="center",
            fontsize=12,color="tab:blue",weight="bold")
    ax.grid(alpha=0.3,axis="y",which="both")
    # 右: 何が起きるか(説明ボックス)
    ax=axes[1]; ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis("off")
    ax.set_title("その結果、データ同化で何が起きるか",fontsize=14)
    ax.add_patch(FancyBboxPatch((0.3,5.6),9.4,3.6,boxstyle="round,pad=0.15",
                                fc="#fdecec",ec="crimson",lw=1.5))
    ax.text(5,8.6,"高感度点に貼った場合",fontsize=13,ha="center",color="crimson",weight="bold")
    ax.text(5,6.9,"信号がノイズを圧倒 → 変位1本で温度場が正しく縛られ、\n"
            "最初のサイクルから真値に追従（RMSE 0.0237K）",fontsize=12,ha="center")
    ax.add_patch(FancyBboxPatch((0.3,0.6),9.4,4.2,boxstyle="round,pad=0.15",
                                fc="#eaf0fb",ec="royalblue",lw=1.5))
    ax.text(5,4.2,"低感度点に貼った場合",fontsize=13,ha="center",color="royalblue",weight="bold")
    ax.text(5,2.2,"信号がノイズに埋もれる → EnKFがノイズを\n"
            "アンサンブルの偶然の相関で「本物」と誤解し、\n"
            "序盤に温度場を逆方向に補正（QoIが−1µmまで逆振れ、\n"
            "最終RMSEも0.0400Kと1.7倍悪い）",fontsize=12,ha="center")
    fig.tight_layout()
    fig.savefig(os.path.join(IMG,"guide_snr.png"),dpi=140); plt.close(fig)
    print("[guide] guide_snr.png")


def fig_equation():
    """変位計1本の読み値 = 未知の温度5040個の重み付き合計 = 方程式1本、の図解."""
    fig,ax=plt.subplots(figsize=(14,8.6)); ax.set_xlim(0,14); ax.set_ylim(0,8.6); ax.axis("off")
    ax.text(7,8.25,"変位計1本の読み値は「温度場についての方程式1本」",
            ha="center",fontsize=17,weight="bold")

    # --- 左: 円筒と変位計 ---
    ax.add_patch(Rectangle((0.9,3.3),2.2,3.4,fc="#eef3f8",ec="k"))
    rng=np.random.default_rng(3)
    for xx,zz in zip(rng.uniform(1.05,2.95,26),rng.uniform(3.45,6.55,26)):
        ax.plot(xx,zz,"o",ms=7,color="#f5b942",mec="orange")
        ax.text(xx,zz-0.02,"?",ha="center",va="center",fontsize=7,color="k")
    ax.text(2.0,2.85,"円筒の温度 $\\Delta T_1 \\cdots \\Delta T_{5040}$\n＝知りたい未知数(5040個)",
            ha="center",fontsize=12)
    # 変位計
    ax.add_patch(Rectangle((1.7,7.0),0.6,0.55,fc="#d9ead3",ec="k"))
    ax.annotate("",xy=(2.0,6.75),xytext=(2.0,7.0),
                arrowprops=dict(arrowstyle="-|>",color="k",lw=2))
    ax.text(2.0,7.85,"変位計(点i)\n読み値 $u_i$=4.2µm ←既知",ha="center",fontsize=12)

    # --- 右: 方程式の図解 ---
    y=5.0
    def box(x,w,txt,fc,ec,sub=None,subc="gray"):
        ax.add_patch(FancyBboxPatch((x,y-0.62),w,1.24,boxstyle="round,pad=0.07",fc=fc,ec=ec,lw=1.6))
        ax.text(x+w/2,y,txt,ha="center",va="center",fontsize=13)
        if sub: ax.text(x+w/2,y-1.05,sub,ha="center",fontsize=10.5,color=subc)
    box(4.0,1.5,"$u_i$\n読み値","#d9ead3","green","既知(測る)","green")
    ax.text(5.75,y,"=",ha="center",fontsize=20)
    box(6.0,1.35,"$W[i,1]$","#f4cccc","crimson","既知(事前計算)","crimson")
    ax.text(7.5,y,"×",ha="center",fontsize=15)
    box(7.65,1.2,"$\\Delta T_1$","#fff2cc","orange","未知","darkorange")
    ax.text(9.0,y,"+",ha="center",fontsize=15)
    box(9.2,1.35,"$W[i,2]$","#f4cccc","crimson","既知","crimson")
    ax.text(10.7,y,"×",ha="center",fontsize=15)
    box(10.85,1.2,"$\\Delta T_2$","#fff2cc","orange","未知","darkorange")
    ax.text(12.45,y,"+ … +",ha="center",fontsize=13)
    box(12.9,1.0,"5040\n項目","#f4cccc","crimson")
    ax.annotate("",xy=(4.0,y),xytext=(3.15,6.2),
                arrowprops=dict(arrowstyle="->",color="k",lw=1.5))

    # --- 下: 高感度 vs 低感度 ---
    ax.add_patch(FancyBboxPatch((3.9,1.9),4.7,1.6,boxstyle="round,pad=0.12",
                                fc="#fdecec",ec="crimson",lw=1.5))
    ax.text(6.25,3.15,"高感度点に貼ると(係数Wが大きい)",fontsize=12.5,ha="center",
            color="crimson",weight="bold")
    ax.text(6.25,2.4,"4.2µm = 大きな係数×温度たち\n→ 温度場の当たりが強く絞れる「使える方程式」",
            fontsize=11,ha="center")
    ax.add_patch(FancyBboxPatch((9.0,1.9),4.7,1.6,boxstyle="round,pad=0.12",
                                fc="#eaf0fb",ec="royalblue",lw=1.5))
    ax.text(11.35,3.15,"低感度点に貼ると(係数がほぼ全部0)",fontsize=12.5,ha="center",
            color="royalblue",weight="bold")
    ax.text(11.35,2.4,"読み値 ≈ 0×温度たち＋ノイズ\n→「0≈0」の式。温度場について何も教えない",
            fontsize=11,ha="center")
    ax.text(7,0.9,"変位計を1本貼る＝未知数5040個へ方程式を1本追加すること。"
            "係数(＝W行感度)が大きい場所ほど、その方程式は温度場を強く縛る",
            ha="center",fontsize=13,bbox=dict(boxstyle="round",fc="#fffbe6",ec="orange"))
    fig.tight_layout()
    fig.savefig(os.path.join(IMG,"guide_equation.png"),dpi=140); plt.close(fig)
    print("[guide] guide_equation.png")


if __name__=="__main__":
    fig_w_matrix(); fig_workflow(); fig_snr(); fig_equation()
