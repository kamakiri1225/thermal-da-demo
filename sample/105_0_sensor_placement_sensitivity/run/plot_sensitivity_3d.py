"""熱感度分布を3D円筒上に表示する(展開図より直感的な版).

results/sensitivity_uz.npz のパッチ感度を円筒表面に貼り、
観測点(変位2点=橙)とヒータ範囲を重ねて表示する。
出力: docs/img/sensitivity_3d.png
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
from dacore.node_locations import (OUTER_RADIUS, INNER_RADIUS, HEIGHT,
    HEATER_ARC_LENGTH, HEATER_AXIAL_HEIGHT, HEATER_CENTER_Z)
IMG=os.path.join(ROOT,"docs","img")

def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    d=np.load(os.path.join(ROOT,"results","sensitivity_uz.npz"))
    S=d["S_um"]; NTH_P=int(d["NTH_P"]); NZ_P=int(d["NZ_P"])

    # 円筒表面を細かい点群にして、各点に属するパッチの感度を割り当て
    # 連続した円筒サーフェス(StructuredGrid)にして面で塗る(点群は分かりにくいため)
    nth,nz=181,81
    th=np.linspace(0,2*np.pi,nth); zz=np.linspace(0,HEIGHT,nz)
    TH,ZZ=np.meshgrid(th,zz,indexing="ij")
    X=OUTER_RADIUS*np.cos(TH); Y=OUTER_RADIUS*np.sin(TH)
    # StructuredGrid はFortran順でフラット化されるため order='F' で対応付ける
    a=np.minimum((TH.ravel(order='F')%(2*np.pi)/(2*np.pi)*NTH_P).astype(int),NTH_P-1)
    b=np.minimum((ZZ.ravel(order='F')/HEIGHT*NZ_P).astype(int),NZ_P-1)

    vmax=abs(S).max()
    obs=[(0.028,0,HEIGHT,"uz_heater"),(-0.028,0,HEIGHT,"uz_opp")]
    pl=pv.Plotter(off_screen=True,shape=(1,2),window_size=(1500,760))
    for j,name in enumerate(["uz_heater(ヒータ側上面)","uz_opp(反対側上面)"]):
        pl.subplot(0,j)
        surf=pv.StructuredGrid(X,Y,ZZ)
        surf["sens"]=S[j,a,b]
        pl.add_mesh(surf,scalars="sens",cmap="coolwarm",clim=[-vmax,vmax],
                    smooth_shading=False,show_edges=False,
                    scalar_bar_args={"title":"dUz/dT [um/K]","title_font_size":20,"label_font_size":16})
        # 観測点マーカー
        for x,y,z,nm in obs:
            col="yellow" if nm in name else "darkorange"
            pl.add_mesh(pv.Sphere(radius=0.004,center=(x,y,z)),color=col)
        # このパネルの観測点(黄色)を明示
        tag = "obs: heater side (+X, yellow)" if j==0 else "obs: opposite side (-X, yellow)"
        pl.add_text(f"which cell's +1K moves {['uz_heater','uz_opp'][j]} ?\n{tag}",
                    font_size=15,color="black")
        pl.camera_position=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]
        pl.camera.zoom(1.2)
        pl.set_background("white")
    out=os.path.join(IMG,"sensitivity_3d.png")
    pl.screenshot(out); pl.close()
    print(f"wrote {os.path.relpath(out,ROOT)}")

if __name__=="__main__": main()
