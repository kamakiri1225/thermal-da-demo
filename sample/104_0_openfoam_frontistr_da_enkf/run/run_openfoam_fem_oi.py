"""実ソルバ(OpenFOAM+FrontISTR)で OI（最適内挿）データ同化を回す双子実験.

EnKF (of_fem_twin) との違い：アンサンブルを使わず **背景トラジェクトリ1本**だけ。
背景誤差共分散 B は事前に決め打ちで固定する（毎サイクルの再計算なし）:
  状態 x=[T固体場(Nc), Q]。 B = σ_T² diag(I,0) + σ_Q² v vᵀ,  v=[s_Q; 1].
  s_Q = ∂T/∂Q の場（cycle1の加熱場の形から作る）, duz/dQ は FrontISTR 2回で近似。
  固定ゲイン K = B Hᵀ (H B Hᵀ + R)⁻¹ ,  x_a = x_b + K(y - h(x_b)).
真値ランは 104 EnKF の既存 truth を再利用（背景10本だけ回せばよい＝EnKFの約1/5）。

出力: results/oi_fullsolver_history.csv/.npz, docs/img/oi_fullsolver_temp.png, oi_fullsolver_disp.png
再現(このフォルダをカレントに):
  OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 nohup python3 run/run_openfoam_fem_oi.py \
    > openfoam/run_fem_oi.log 2>&1 &
"""
from __future__ import annotations
import os, sys, time
import numpy as np, yaml
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from daof import of_case
from daof.of_case import BASE_CASE
from fem.fem_obs import displacement_obs
K=273.15
RESULTS=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
# --- 複数回（σ_Qを変えて）回すためのタグ。OI_TAG="" が1本目, "b" が2本目 等 ---
TAG=os.environ.get("OI_TAG","")
SUF=("_"+TAG) if TAG else ""
WORK=os.path.join(ROOT,"openfoam","run_fem_oi"+SUF)
TRUTH=os.path.join(ROOT,"openfoam","run_fem_enkf","truth")   # 既存の真値を再利用
# --- OIの決め打ち背景標準偏差（σ_Qは環境変数で切替＝調整の必要性を示すため） ---
SIGB_T=float(os.environ.get("OI_SIGB_T","1.5"))   # 背景温度場の標準偏差 [K]
SIGB_Q=float(os.environ.get("OI_SIGB_Q","6.0"))   # 背景発熱の標準偏差 [W]
FREE=os.environ.get("OI_FREE","0")=="1"           # 1: 同化なし(free run)=毎サイクル補正しない
# --- でたらめ初期（上側に外す：free-runが真値の上に来る向き） ---
T0_OFF=5.0     # 初期温度を +5K 高く
Q0_GUESS=22.0  # 発熱を過大（真値15W）


def _tname(t): return str(int(t)) if float(t)==int(t) else ("%g"%float(t))
def _log(m): print(m, flush=True)


def main():
    if not os.path.isdir(os.path.join(BASE_CASE,"constant","solid","polyMesh")):
        raise SystemExit("base_case のメッシュがありません。先に openfoam/setup_base_case.sh")
    os.makedirs(WORK, exist_ok=True); os.makedirs(RESULTS, exist_ok=True); os.makedirs(IMG, exist_ok=True)
    with open(os.path.join(ROOT,"openfoam","da_openfoam_config.yaml")) as f: cfg=yaml.safe_load(f)

    centres=of_case.solid_cell_centres(); Nc=len(centres)
    oc=cfg["observation"]; tnames=list(oc["temp_probes"])
    obs_cells=of_case.nearest_cells([oc["temp_probes"][n] for n in tnames], centres)
    n_t=len(obs_cells); n_d=len(oc["disp_names"]); n_obs=n_t+n_d
    R=np.diag([oc["temp_noise_C"]**2]*n_t + [oc["disp_noise_mm"]**2]*n_d)
    dt_obs=cfg["experiment"]["obs_interval_s"]
    times=list(np.round(np.arange(dt_obs, cfg["experiment"]["t_end_s"]+1e-9, dt_obs),6))
    Q_true=cfg["truth"]["Q_true_W"]; T_air=cfg["truth"]["T0_K"]
    seed=cfg["experiment"]["seed"]; rng_o=np.random.default_rng(seed)

    def member_obs(case_dir, t, fem_work):
        T=of_case.read_solid_T(case_dir, t)
        u=displacement_obs(case_dir, _tname(t), fem_work)
        return np.concatenate([T[obs_cells], u]), T

    # ---- 観測 y を既存 truth から作る（真値ランは再利用） ----
    _log(f"[oi] build observations from existing truth: {TRUTH}")
    y_obs={}; truth_T={0.0:np.full(Nc,T_air)}; truth_u={}
    for t in times:
        # 真値の場は TRUTH から読むだけ（読み取り専用）。FrontISTR作業は OI側 run_fem_oi に書く
        clean,Tf=member_obs(TRUTH, t, os.path.join(WORK, f"fem_truth_t{t:g}"))
        truth_T[t]=Tf; truth_u[t]=clean[n_t:].copy()
        y_obs[t]=clean+rng_o.normal(0.0, np.sqrt(np.diag(R)))
        _log(f"[oi]  truth t={t:g}s Tmax={Tf.max()-K:.2f}C uz_um={np.round(clean[n_t:]*1000,2)}")

    # ---- でたらめ初期の背景1本 ----
    of_case.prepare_member(WORK)
    of_case.set_solid_state(WORK, 0.0, T_air+T0_OFF, Q0_GUESS)
    Qb=float(Q0_GUESS)
    hist={"times":[0.0],"T_obs_analysis":[np.full(n_t,T_air+T0_OFF)],"T_obs_truth":[np.full(n_t,T_air)],
          "T_obs_forecast":[np.full(n_t,T_air+T0_OFF)],   # 補正前の予報値（温度）
          "u_analysis":[np.zeros(n_d)],"u_truth":[np.zeros(n_d)],"u_forecast":[np.zeros(n_d)],
          "q":[Qb],"free":bool(FREE),
          "rmse_field":[float(np.sqrt(np.mean((np.full(Nc,T_air+T0_OFF)-truth_T[0.0])**2)))]}
    gain=None; t0=0.0; t_start=time.time()

    for ci,t1 in enumerate(times,1):
        _log(f"[oi] === cycle {ci}/{len(times)}: forecast [{t0:g},{t1:g}] (1 background) ===")
        ok,logp=of_case.run_window(WORK, t0, t1)
        if not ok: raise RuntimeError(f"OpenFOAM failed [{t0},{t1}] see {logp}")
        yb, Tb = member_obs(WORK, t1, os.path.join(WORK, f"fem_t{t1:g}"))
        hist["T_obs_forecast"].append(yb[:n_t].copy()); hist["u_forecast"].append(yb[n_t:].copy())  # 補正前の予報値

        if (not FREE) and gain is None:
            # --- cycle1 の加熱場から B と 固定ゲインを構築 ---
            sQ=(Tb-T_air)/max(Qb,1e-6)                       # ∂T/∂Q の形 [K/W]
            # duz/dQ を FrontISTR 2回で近似（現場と +10%発熱相当のスケール場）
            hot=T_air+1.1*(Tb-T_air)
            of_case.set_solid_state(WORK, t1, hot, Qb)        # 一時的に書いて評価
            u_hot=displacement_obs(WORK, _tname(t1), os.path.join(WORK,"fem_dqdz"))
            of_case.set_solid_state(WORK, t1, Tb, Qb)         # 戻す
            duz_dQ=(u_hot - yb[n_t:])/(0.1*max(Qb,1e-6))      # [mm/W]
            sQc=sQ[obs_cells]                                 # 観測セルでの sQ
            # H B Hᵀ (n_obs×n_obs)
            HBHt=np.zeros((n_obs,n_obs))
            for a in range(n_t):
                for b in range(n_t): HBHt[a,b]=SIGB_T**2*(a==b)+SIGB_Q**2*sQc[a]*sQc[b]
                for b in range(n_d): HBHt[a,n_t+b]=SIGB_Q**2*sQc[a]*duz_dQ[b]
            for a in range(n_d):
                for b in range(n_t): HBHt[n_t+a,b]=SIGB_Q**2*duz_dQ[a]*sQc[b]
                for b in range(n_d): HBHt[n_t+a,n_t+b]=SIGB_Q**2*duz_dQ[a]*duz_dQ[b]
            S=HBHt+R; Sinv=np.linalg.inv(S)
            # B Hᵀ の field部(Nc×n_obs) と Q部(n_obs)
            BHt_f=np.zeros((Nc,n_obs)); BHt_q=np.zeros(n_obs)
            for i in range(n_t):
                BHt_f[:,i]=SIGB_Q**2*sQ*sQc[i]; BHt_f[obs_cells[i],i]+=SIGB_T**2
                BHt_q[i]=SIGB_Q**2*sQc[i]
            for j in range(n_d):
                BHt_f[:,n_t+j]=SIGB_Q**2*sQ*duz_dQ[j]; BHt_q[n_t+j]=SIGB_Q**2*duz_dQ[j]
            Kf=BHt_f@Sinv; Kq=BHt_q@Sinv                     # 固定ゲイン
            gain=(Kf,Kq)
            _log(f"[oi] built fixed B/gain: sQ range {sQ.min():.3f}..{sQ.max():.3f} K/W, duz/dQ={np.round(duz_dQ*1000,3)} um/W")

        if FREE:
            Ta=Tb; Qa=Qb                          # 同化なし＝補正しない（解析=予報）
        else:
            Kf,Kq=gain
            innov=y_obs[t1]-yb
            Ta=np.clip(Tb+Kf@innov, 250.0, 400.0)
            Qa=float(np.clip(Qb+Kq@innov, *cfg["filter"]["Q_bounds_W"]))
        of_case.set_solid_state(WORK, t1, Ta, Qa); Qb=Qa
        # 変位は「同化後の温度場」から解く。必ず set_solid_state で解析場を書いた後に評価する
        # （順序を誤ると予報場の変位＝u_forecast を記録してしまう）
        if FREE:
            u_ana=yb[n_t:]                                            # 同化なし＝予報のまま
        else:
            u_ana=displacement_obs(WORK,_tname(t1),os.path.join(WORK,f"fem_a{t1:g}"))

        hist["times"].append(float(t1))
        hist["T_obs_analysis"].append(Ta[obs_cells].copy()); hist["T_obs_truth"].append(truth_T[t1][obs_cells].copy())
        hist["u_analysis"].append(u_ana.copy())
        hist["u_truth"].append(truth_u[t1].copy()); hist["q"].append(Qa)
        hist["rmse_field"].append(float(np.sqrt(np.mean((Ta-truth_T[t1])**2))))
        el=time.time()-t_start
        _log(f"[oi]  analysis t={t1:g}s RMSE(field)={hist['rmse_field'][-1]:.3f}K Q={Qa:.2f}(true {Q_true}) elapsed {el/60:.1f}min")
        t0=t1

    for k in ["T_obs_analysis","T_obs_truth","T_obs_forecast","u_analysis","u_truth","u_forecast"]: hist[k]=np.array(hist[k])
    np.savez(os.path.join(RESULTS,f"oi_fullsolver{SUF}_history.npz"), **hist)
    total=time.time()-t_start
    yaml.safe_dump({"method":"OI full-solver (OpenFOAM+FrontISTR, single background, fixed B)",
        "n_background_per_cycle":1,"cycles":len(times),
        "rmse_field_initial_K":float(hist["rmse_field"][0]),"rmse_field_final_K":float(hist["rmse_field"][-1]),
        "Q_final_W":float(hist["q"][-1]),"Q_true_W":float(Q_true),
        "wall_seconds":float(total),"wall_hours":float(total/3600.0)},
        open(os.path.join(RESULTS,f"oi_fullsolver{SUF}_summary.yaml"),"w"), allow_unicode=True, sort_keys=False)
    _log(f"[oi] DONE. field RMSE {hist['rmse_field'][0]:.2f}->{hist['rmse_field'][-1]:.3f}K, "
         f"Q->{hist['q'][-1]:.2f}W, wall {total/3600:.2f}h")
    _make_figures(hist, tnames)


def _make_figures(hist, tnames):
    from dacore import plots as _p
    import matplotlib.pyplot as plt
    t=np.array(hist["times"])
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    cols=["tab:red","tab:blue"]
    for j in range(hist["T_obs_truth"].shape[1]):
        axes[0].plot(t, hist["T_obs_truth"][:,j]-K,"-",color=cols[j],lw=3,alpha=0.4,label=f"{tnames[j]} 真値")
        axes[0].plot(t, hist["T_obs_analysis"][:,j]-K,"--o",color=cols[j],ms=3,label=f"{tnames[j]} OI同化後")
    axes[0].set_xlabel("time [s]"); axes[0].set_ylabel("温度 [degC]"); axes[0].grid(alpha=.3); axes[0].legend(fontsize=9)
    axes[0].set_title("温度時刻歴（実ソルバOI：観測点）")
    for j in range(hist["u_truth"].shape[1]):
        axes[1].plot(t, hist["u_truth"][:,j]*1000,"-",color=cols[j],lw=3,alpha=0.4,label=f"変位{j} 真値")
        axes[1].plot(t, hist["u_analysis"][:,j]*1000,"--o",color=cols[j],ms=3,label=f"変位{j} OI同化後")
    axes[1].set_xlabel("time [s]"); axes[1].set_ylabel("上面 Uz [µm]"); axes[1].grid(alpha=.3); axes[1].legend(fontsize=9)
    axes[1].set_title("変位時刻歴（実ソルバOI：上面2点）")
    fig.suptitle("実ソルバ(OpenFOAM+FrontISTR) OI：背景1本・固定B（アンサンブル不要）",fontsize=13,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(os.path.join(IMG,f"oi_fullsolver{SUF}_temp.png"),dpi=140); plt.close(fig)
    print("[oi] wrote docs/img/oi_fullsolver_temp.png", flush=True)


if __name__=="__main__": main()
