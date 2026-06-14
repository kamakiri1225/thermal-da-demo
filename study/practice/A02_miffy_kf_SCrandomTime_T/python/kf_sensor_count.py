#!/usr/bin/env python3
"""
kf_sensor_count.py — センサー数 m の影響比較

q = 0: モデルは x_ss を知らない
Q = σ_q²·I: モデル誤差（P が収縮しないようにする）
30×30 グリッドで高速実行 (n=900)

出力:
  fig10_sensor_count.png     — 静止画（各 m × 主要ステップ）
  anim_kf_sensor_count.gif   — アニメーション（4パネル, センサー緑, 1ステップ刻み）
"""
import numpy as np
import scipy.sparse as sp
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os, time

NX, NY = 30, 30
n = NX * NY

T_BG = 0.; T_FUR = 100.; T_EYE = 5.; T_MOUTH = 5.
T_DRESS = 55.; T_COLLAR = 70.; T_OUTLINE = 8.

def make_miffy():
    f = np.full((NY, NX), T_BG)
    s = 30 / 60.0
    def fe(cy, cx, ry, rx, val):
        yy, xx = np.ogrid[:NY, :NX]
        f[((yy - int(cy*s)) / max(int(ry*s), 1))**2
          + ((xx - int(cx*s)) / max(int(rx*s), 1))**2 <= 1.0] = val
    def fr(y1, y2, x1, x2, val):
        f[max(0, int(y1*s)):min(NY-1, int(y2*s))+1,
          max(0, int(x1*s)):min(NX-1, int(x2*s))+1] = val
    fe(5,30,9,21,T_DRESS);  fr(0,11,7,52,T_DRESS); fr(10,12,24,36,T_COLLAR)
    fe(48,19,12,9,T_OUTLINE); fe(48,40,12,9,T_OUTLINE); fe(24,30,15,24,T_OUTLINE)
    fe(48,19,11,7,T_FUR);   fe(48,40,11,7,T_FUR)
    fr(35,40,12,27,T_FUR);  fr(35,40,33,48,T_FUR)
    fe(24,30,14,22,T_FUR);  fr(9,13,25,34,T_FUR)
    fe(24,21,2,3,T_EYE);    fe(24,39,2,3,T_EYE)
    for d in range(-1, 2):
        y_, x_ = min(max(int(18*s)+d, 0), NY-1), min(max(int(30*s)+d, 0), NX-1)
        f[y_, x_] = T_MOUTH
        f[min(max(int(18*s)-d, 0), NY-1), x_] = T_MOUTH
    return f

def build_laplacian(NY, NX):
    n = NY * NX
    row, col, data = [], [], []
    for iy in range(NY):
        for ix in range(NX):
            k = iy * NX + ix
            row.append(k); col.append(k); data.append(-4.)
            for diy, dix in [(-1,0),(1,0),(0,-1),(0,1)]:
                jy, jx = iy+diy, ix+dix
                if 0 <= jy < NY and 0 <= jx < NX:
                    row.append(k); col.append(jy*NX+jx); data.append(1.)
    return sp.csr_matrix((data,(row,col)), shape=(n,n), dtype=np.float64)

x_ss  = make_miffy().flatten()
L     = build_laplacian(NY, NX)
ALPHA = 0.5; DT = 0.2
M_sp  = sp.eye(n, format='csr') + DT * ALPHA * L

SIGMA_R  = 5.0
SIGMA_Q  = 8.0
N_STEPS  = 100
SIGMA_B0 = 30.0
L_CORR   = 5.0
M_LIST   = [1, 5, 10, 30, 50]

iy_all = np.arange(n) // NX
ix_all = np.arange(n) % NX

def run_kf(m_obs, seed=42):
    rng = np.random.default_rng(seed)
    CHUNK = 100
    P = np.zeros((n, n), dtype=np.float64)
    for i in range(0, n, CHUNK):
        j  = min(i + CHUNK, n)
        dy = (iy_all[i:j, None] - iy_all[None, :]).astype(np.float32)
        dx = (ix_all[i:j, None] - ix_all[None, :]).astype(np.float32)
        P[i:j, :] = (SIGMA_B0**2 * np.exp(-0.5*(dy**2+dx**2)/L_CORR**2)).astype(np.float64)

    x_kf = np.zeros(n)
    history = {0: {'x': x_kf.copy(), 'obs': None,
                   'rmse': np.sqrt(np.mean(x_kf**2 + x_ss**2 - 2*x_kf*x_ss))}}

    for step in range(1, N_STEPS + 1):
        x_f = M_sp @ x_kf
        MP  = M_sp.dot(P); P_f = (M_sp.dot(MP.T)).T; P_f = 0.5*(P_f+P_f.T)
        np.fill_diagonal(P_f, P_f.diagonal() + SIGMA_Q**2)

        obs_idx = np.sort(rng.choice(n, m_obs, replace=False))
        y = x_ss[obs_idx] + rng.normal(0., SIGMA_R, m_obs)

        PHt  = P_f[:, obs_idx]
        S    = P_f[np.ix_(obs_idx, obs_idx)] + SIGMA_R**2 * np.eye(m_obs)
        K    = PHt @ np.linalg.inv(S)
        x_kf = x_f + K @ (y - x_f[obs_idx])
        HPf  = P_f[obs_idx, :]; P = P_f - K @ HPf; P = 0.5*(P+P.T)

        rmse = np.sqrt(np.mean((x_kf - x_ss)**2))
        history[step] = {'x': x_kf.copy(), 'obs': obs_idx.copy(), 'rmse': rmse}

    return history

print("KF 計算中（4パターン, q=0, Q=σ_q²·I）...")
all_hist = {}
for m in M_LIST:
    print(f"  m={m:2d}", end=' ', flush=True)
    t0 = time.time()
    all_hist[m] = run_kf(m)
    print(f"-> {time.time()-t0:.1f}s  RMSE_final={all_hist[m][N_STEPS]['rmse']:.2f} degC")

script_dir = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(script_dir, '..', 'img')
os.makedirs(OUT, exist_ok=True)

# ─────────────────────────────────────────
#  静止画: 各 m × 主要ステップ
# ─────────────────────────────────────────
SHOW_STEPS = [0, 5, 20, 100]
CMAP = 'RdYlBu_r'

fig_s, axes_s = plt.subplots(len(M_LIST), len(SHOW_STEPS),
                              figsize=(len(SHOW_STEPS)*2.4+0.4, len(M_LIST)*2.6+1.0))
fig_s.patch.set_facecolor('#111111')

for row, m in enumerate(M_LIST):
    for col, s in enumerate(SHOW_STEPS):
        ax = axes_s[row, col]
        ax.set_facecolor('#111111'); ax.set_xticks([]); ax.set_yticks([])
        for sp_ in ax.spines.values(): sp_.set_edgecolor('#444444')

        d = all_hist[m][s]
        ax.imshow(d['x'].reshape(NY, NX), vmin=0, vmax=T_FUR,
                  cmap=CMAP, origin='lower', interpolation='nearest')

        if d['obs'] is not None:
            sz = max(60 // m, 8)
            iy_o = d['obs'] // NX; ix_o = d['obs'] % NX
            ax.scatter(ix_o, iy_o, s=sz, c='#00E676', marker='+',
                       linewidths=max(1.5-m*0.02, 0.8), zorder=5, alpha=0.95)

        ax.text(0.5, -0.07, f'RMSE={d["rmse"]:.1f}',
                transform=ax.transAxes, ha='center', color='cyan', fontsize=8)
        if col == 0:
            ax.set_ylabel(f'm = {m} sensors/step',
                          color='#64B5F6', fontsize=9, fontweight='bold')
        if row == 0:
            ax.set_title(f'Step {s}', color='white', fontsize=9, fontweight='bold')

fig_s.suptitle(
    f'KF sensor count comparison  (q=0, σ_q={SIGMA_Q}, σ_r={SIGMA_R} degC, 30×30 grid)\n'
    f'green + = sensor position  |  DA recovers Miffy from observations only',
    color='white', fontsize=10, fontweight='bold')
plt.tight_layout(rect=[0, 0.0, 1.0, 0.93])
out_s = os.path.join(OUT, 'fig10_sensor_count.png')
fig_s.savefig(out_s, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig_s)
print(f"保存: {out_s}")

# ─────────────────────────────────────────
#  GIF アニメーション（4パネル, 1ステップ刻り）
# ─────────────────────────────────────────
print(f"\nGIF 作成中（{N_STEPS+1} フレーム × {len(M_LIST)} パネル）...")
SENSOR_COLOR = '#00E676'
ANIM_STEPS   = list(range(0, N_STEPS + 1))

fig_a = plt.figure(figsize=(len(M_LIST)*3.4 + 0.8, 4.6))
fig_a.patch.set_facecolor('#111111')
gs = fig_a.add_gridspec(1, len(M_LIST)+1,
                        left=0.03, right=0.93, bottom=0.14, top=0.84, wspace=0.08,
                        width_ratios=[1]*len(M_LIST) + [0.06])

ims_a   = []
scats_a = []
rmse_txts = []
step_texts = []

for col, m in enumerate(M_LIST):
    ax = fig_a.add_subplot(gs[0, col])
    ax.set_facecolor('#111111'); ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values(): sp_.set_edgecolor('#555555')

    d0 = all_hist[m][0]
    im = ax.imshow(d0['x'].reshape(NY, NX), vmin=0, vmax=T_FUR,
                   cmap=CMAP, origin='lower', interpolation='nearest')
    ims_a.append(im)

    sc = ax.scatter([], [], s=35, c=SENSOR_COLOR, marker='+',
                    linewidths=1.5, zorder=5, alpha=0.95)
    scats_a.append(sc)

    ax.set_title(f'm = {m} sensors/step', color='#64B5F6',
                 fontsize=9.5, fontweight='bold', pad=4)

    if col == 0:
        st = ax.text(0.5, 0.95, f'Step   0 / {N_STEPS}',
                     transform=ax.transAxes, ha='center', va='top',
                     color='yellow', fontsize=12, fontweight='bold',
                     bbox=dict(facecolor='#000000', alpha=0.75,
                               edgecolor='yellow', linewidth=1.2, pad=3))
    else:
        st = ax.text(0., 0., '', transform=ax.transAxes)
    step_texts.append(st)

    rt = ax.text(0.5, 0.04, f'RMSE={d0["rmse"]:.1f} degC',
                 transform=ax.transAxes, ha='center', va='bottom',
                 color='cyan', fontsize=9, fontweight='bold',
                 bbox=dict(facecolor='#111111', alpha=0.7, edgecolor='none', pad=1))
    rmse_txts.append(rt)

cbar_ax = fig_a.add_subplot(gs[0, -1])
sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(0, T_FUR))
cbar = fig_a.colorbar(sm, cax=cbar_ax)
cbar.ax.tick_params(colors='white', labelsize=8)
cbar.set_label('Temp [degC]', color='white', fontsize=8, labelpad=6)

fig_a.text(0.47, 0.97,
           f'KF Data Assimilation — Sensor count comparison  '
           f'(q=0, σ_q={SIGMA_Q} degC, n={n} cells)',
           ha='center', va='top', color='white', fontsize=10.5, fontweight='bold')
fig_a.text(0.47, 0.02,
           f'green + = sensor position (re-randomized every step)',
           ha='center', va='bottom', color=SENSOR_COLOR, fontsize=9)

def update_anim(frame_idx):
    step = ANIM_STEPS[frame_idx]
    artists = []
    for col, m in enumerate(M_LIST):
        d = all_hist[m][step]
        ims_a[col].set_data(d['x'].reshape(NY, NX))
        artists.append(ims_a[col])

        if d['obs'] is not None:
            iy_o = d['obs'] // NX; ix_o = d['obs'] % NX
            scats_a[col].set_offsets(np.column_stack([ix_o, iy_o]))
        else:
            scats_a[col].set_offsets(np.empty((0, 2)))
        artists.append(scats_a[col])

        rmse_txts[col].set_text(f'RMSE={d["rmse"]:.1f} degC')
        artists.append(rmse_txts[col])

        if col == 0:
            step_texts[0].set_text(f'Step {step:3d} / {N_STEPS}')
            artists.append(step_texts[0])

    return artists

t0 = time.time()
ani = animation.FuncAnimation(fig_a, update_anim, frames=len(ANIM_STEPS),
                               interval=150, blit=False)
out_gif = os.path.join(OUT, 'anim_kf_sensor_count.gif')
ani.save(out_gif, writer=animation.PillowWriter(fps=8), dpi=120)
plt.close(fig_a)
print(f"保存: {out_gif}  ({time.time()-t0:.1f}s)")
print("\n完了 →", OUT)
