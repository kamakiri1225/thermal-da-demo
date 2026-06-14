#!/usr/bin/env python3
"""
enkf_fixsensor.py -- Case 1: EnKF with fixed sensor positions, sensor count comparison

Goal: recover static Miffy temperature field from sparse fixed sensors.
      Compare m = [1, 10, 30, 50] fixed sensor positions.

Key choices:
  - No physical forecast (truth is static): X_forecast = X_analysis
  - Local multiplicative inflation: cells far from sensors get NO inflation
    => prevents ensemble divergence for uncorrected far-field cells
  - Gaussian localization radius r_loc = 8 cells

Grid: 30x30 (n=900), N=50 ensemble, 50 DA cycles

Outputs:
  fig01_rmse.png          -- RMSE convergence (4 sensor counts)
  fig02_final_fields.png  -- Final ensemble mean fields (2x2)
  anim_enkf_fixsensor.gif -- 4-panel animation
"""
import numpy as np
import scipy.sparse as sp
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.ndimage import gaussian_filter
import os, time

# ── Grid ──────────────────────────────────────────────────
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
    fe( 5,30, 9,21,T_DRESS);  fr( 0,11, 7,52,T_DRESS); fr(10,12,24,36,T_COLLAR)
    fe(48,19,12, 9,T_OUTLINE); fe(48,40,12, 9,T_OUTLINE); fe(24,30,15,24,T_OUTLINE)
    fe(48,19,11, 7,T_FUR);    fe(48,40,11, 7,T_FUR)
    fr(35,40,12,27,T_FUR);   fr(35,40,33,48,T_FUR)
    fe(24,30,14,22,T_FUR);   fr( 9,13,25,34,T_FUR)
    fe(24,21, 2, 3,T_EYE);   fe(24,39, 2, 3,T_EYE)
    for d in range(-1, 2):
        y_ = min(max(int(18*s)+d, 0), NY-1)
        x_ = min(max(int(30*s)+d, 0), NX-1)
        f[y_, x_] = T_MOUTH
        f[min(max(int(18*s)-d,0),NY-1), x_] = T_MOUTH
    return f

x_ss   = make_miffy().flatten()
iy_all = np.arange(n) // NX
ix_all = np.arange(n) % NX

# ── EnKF parameters ───────────────────────────────────────
N_ENS   = 50
SIGMA_R = 5.0    # observation noise std [degC]
INFL    = 1.10   # max local inflation factor
R_LOC   = 8.0    # Gaussian localization radius [cells]
N_CYC   = 50     # DA cycles (static truth: 50 cycles sufficient)
SIGMA_B = 30.0   # initial ensemble std
CORR_L  = 5.0    # initial spatial correlation length [cells]

M_LIST = [1, 10, 30, 50]

# ── EnKF runner (fixed sensors, local inflation) ──────────
def run_enkf(m_obs, seed=42):
    rng = np.random.default_rng(seed)

    # Fixed sensor positions (chosen once)
    obs_idx = np.sort(rng.choice(n, m_obs, replace=False))
    obs_iy  = obs_idx // NX
    obs_ix  = obs_idx % NX

    # Localization weights
    dy  = iy_all[:, None] - obs_iy[None, :]   # (n, m)
    dx  = ix_all[:, None] - obs_ix[None, :]
    rho = np.exp(-0.5 * (dy**2 + dx**2) / R_LOC**2)

    # Local inflation factor per cell: max localization across all sensors
    # => cells near sensors get full INFL; far cells get 1.0 (no inflation)
    coverage   = np.max(rho, axis=1)           # (n,): max rho over all sensors
    infl_cell  = 1.0 + (INFL - 1.0) * coverage  # (n,): between 1.0 and INFL

    # Initial ensemble: spatially correlated noise, zero mean
    X = np.zeros((n, N_ENS))
    for i in range(N_ENS):
        noise     = rng.normal(0., 1., (NY, NX))
        noise_cor = gaussian_filter(noise, sigma=CORR_L)
        noise_cor = noise_cor / (noise_cor.std() + 1e-10) * SIGMA_B
        X[:, i]   = noise_cor.flatten()

    rmse0   = np.sqrt(np.mean((X.mean(axis=1) - x_ss)**2))
    history = {0: {'x': X.mean(axis=1).copy(), 'rmse': rmse0}}

    for cyc in range(1, N_CYC + 1):
        # No physical forecast (static truth): X_f = X
        X_f = X.copy()

        # Ensemble anomalies
        xbar = X_f.mean(axis=1, keepdims=True)
        Xp   = X_f - xbar
        HXp  = Xp[obs_idx, :]      # (m, N)

        # Innovation covariance & localized Kalman gain
        S     = (HXp @ HXp.T) / (N_ENS - 1) + SIGMA_R**2 * np.eye(m_obs)
        PHt   = (Xp @ HXp.T) / (N_ENS - 1)
        K     = np.linalg.solve(S.T, PHt.T).T   # (n, m)
        K_loc = K * rho                           # Schur-product localization

        # Perturbed-observation stochastic EnKF
        y     = x_ss[obs_idx] + rng.normal(0., SIGMA_R, m_obs)
        eps   = rng.normal(0., SIGMA_R, (m_obs, N_ENS))
        innov = y[:, None] + eps - X_f[obs_idx, :]

        X_a = X_f + K_loc @ innov

        # Local multiplicative inflation (only near sensors)
        xbar_a = X_a.mean(axis=1, keepdims=True)
        X      = xbar_a + infl_cell[:, None] * (X_a - xbar_a)

        rmse = np.sqrt(np.mean((X.mean(axis=1) - x_ss)**2))
        history[cyc] = {'x': X.mean(axis=1).copy(), 'rmse': rmse}

    return history, obs_idx


# ── Run all m values ──────────────────────────────────────
print(f"EnKF fixed sensors (static, local-INFL={INFL}, r_loc={R_LOC}, N={N_ENS}) ...")
all_hist = {}
all_obs  = {}
for m in M_LIST:
    print(f"  m={m:2d}", end=' ', flush=True)
    t0 = time.time()
    hist, obs = run_enkf(m)
    all_hist[m] = hist
    all_obs[m]  = obs
    print(f"-> {time.time()-t0:.1f}s  RMSE_final={hist[N_CYC]['rmse']:.2f} degC")

script_dir = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(script_dir, '..', 'img')
os.makedirs(OUT, exist_ok=True)

CMAP   = 'RdYlBu_r'
COLORS = ['#FF6B6B', '#FFD93D', '#6BCB77', '#4D96FF']
FS_T   = 15
FS_L   = 13
FS_TK  = 11
FS_LEG = 12

# ── fig01: RMSE convergence ───────────────────────────────
fig1, ax1 = plt.subplots(figsize=(11, 5.5))
fig1.patch.set_facecolor('#111111')
ax1.set_facecolor('#111111')

cyc_arr = np.arange(N_CYC + 1)
for m, col in zip(M_LIST, COLORS):
    rmse_arr = [all_hist[m][c]['rmse'] for c in cyc_arr]
    ax1.plot(cyc_arr, rmse_arr, color=col, lw=2.5, label=f'm = {m} sensors')

ax1.set_xlabel('DA Cycle', color='white', fontsize=FS_L)
ax1.set_ylabel('RMSE [degC]', color='white', fontsize=FS_L)
ax1.set_title(
    f'EnKF Fixed Sensors -- RMSE Convergence  (N={N_ENS}, r_loc={R_LOC} cells)\n'
    f'Local inflation only near sensors -- cells without coverage keep prior',
    color='white', fontsize=FS_T, fontweight='bold')
ax1.tick_params(colors='white', labelsize=FS_TK)
ax1.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=FS_LEG)
ax1.grid(alpha=0.25, color='gray', lw=0.5)
for sp_ in ax1.spines.values(): sp_.set_edgecolor('#444444')

plt.tight_layout()
out1 = os.path.join(OUT, 'fig01_rmse.png')
fig1.savefig(out1, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig1)
print(f"\nSaved: {out1}")

# ── fig02: Final fields (2x2) ─────────────────────────────
fig2, axes2 = plt.subplots(2, 2, figsize=(11, 11))
fig2.patch.set_facecolor('#111111')

for ax, m, col in zip(axes2.flatten(), M_LIST, COLORS):
    ax.set_facecolor('#111111')
    ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values(): sp_.set_edgecolor('#333333')

    xf     = all_hist[m][N_CYC]['x'].reshape(NY, NX)
    rmse_f = all_hist[m][N_CYC]['rmse']
    ax.imshow(xf, vmin=0, vmax=T_FUR, cmap=CMAP, origin='lower', interpolation='nearest')

    obs = all_obs[m]
    sz  = max(300 // max(m, 1), 20)
    ax.scatter(obs % NX, obs // NX, s=sz, c='#00E676', marker='+',
               linewidths=max(2.0 - m*0.03, 0.9), zorder=5, alpha=0.95)

    ax.set_title(f'm = {m} sensors   RMSE = {rmse_f:.1f} degC',
                 color=col, fontsize=FS_T, fontweight='bold')

fig2.suptitle(
    f'EnKF Fixed Sensors -- Final State (Cycle {N_CYC})\n'
    'green + = fixed sensor positions',
    color='white', fontsize=FS_T, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.96])
out2 = os.path.join(OUT, 'fig02_final_fields.png')
fig2.savefig(out2, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig2)
print(f"Saved: {out2}")

# ── GIF animation (4 panels) ──────────────────────────────
print(f"\nMaking GIF ({N_CYC+1} frames x {len(M_LIST)} panels) ...")

fig_a = plt.figure(figsize=(len(M_LIST)*3.8 + 0.8, 5.5))
fig_a.patch.set_facecolor('#111111')
gs = fig_a.add_gridspec(1, len(M_LIST)+1,
                        left=0.03, right=0.93, bottom=0.10, top=0.84,
                        wspace=0.06,
                        width_ratios=[1]*len(M_LIST) + [0.06])

ims_a     = []
rmse_txts = []
step_txts = []

for col_i, (m, mcol) in enumerate(zip(M_LIST, COLORS)):
    ax = fig_a.add_subplot(gs[0, col_i])
    ax.set_facecolor('#111111'); ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values(): sp_.set_edgecolor('#333333')

    d0 = all_hist[m][0]
    im = ax.imshow(d0['x'].reshape(NY, NX), vmin=0, vmax=T_FUR,
                   cmap=CMAP, origin='lower', interpolation='nearest')
    ims_a.append(im)

    obs = all_obs[m]
    sz  = max(300 // max(m, 1), 20)
    ax.scatter(obs % NX, obs // NX, s=sz, c='#00E676', marker='+',
               linewidths=max(2.0 - m*0.03, 0.9), zorder=5, alpha=0.95)

    ax.set_title(f'm = {m} sensors', color=mcol, fontsize=FS_T,
                 fontweight='bold', pad=4)

    if col_i == 0:
        st = ax.text(0.5, 0.95, f'Cycle   0 / {N_CYC}',
                     transform=ax.transAxes, ha='center', va='top',
                     color='yellow', fontsize=13, fontweight='bold',
                     bbox=dict(facecolor='#000000', alpha=0.75,
                               edgecolor='yellow', linewidth=1.2, pad=3))
    else:
        st = ax.text(0., 0., '', transform=ax.transAxes)
    step_txts.append(st)

    rt = ax.text(0.5, 0.04, f'RMSE={d0["rmse"]:.1f} degC',
                 transform=ax.transAxes, ha='center', va='bottom',
                 color='cyan', fontsize=FS_TK, fontweight='bold',
                 bbox=dict(facecolor='#111111', alpha=0.7, edgecolor='none', pad=1))
    rmse_txts.append(rt)

cbar_ax = fig_a.add_subplot(gs[0, -1])
sm      = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(0, T_FUR))
cbar    = fig_a.colorbar(sm, cax=cbar_ax)
cbar.ax.tick_params(colors='white', labelsize=FS_TK)
cbar.set_label('Temp [degC]', color='white', fontsize=FS_TK, labelpad=4)

fig_a.text(0.47, 0.97,
           f'EnKF Fixed Sensors -- local inflation only  (N={N_ENS}, r_loc={R_LOC} cells)',
           ha='center', va='top', color='white', fontsize=FS_T, fontweight='bold')
fig_a.text(0.47, 0.02,
           'green + = fixed sensor position  (same every cycle)',
           ha='center', va='bottom', color='#00E676', fontsize=FS_L)


def update_anim(frame):
    artists = []
    for col_i, m in enumerate(M_LIST):
        d = all_hist[m][frame]
        ims_a[col_i].set_data(d['x'].reshape(NY, NX))
        artists.append(ims_a[col_i])
        rmse_txts[col_i].set_text(f'RMSE={d["rmse"]:.1f} degC')
        artists.append(rmse_txts[col_i])
    step_txts[0].set_text(f'Cycle {frame:3d} / {N_CYC}')
    artists.append(step_txts[0])
    return artists


t0  = time.time()
ani = animation.FuncAnimation(fig_a, update_anim, frames=N_CYC + 1,
                               interval=150, blit=False)
out_gif = os.path.join(OUT, 'anim_enkf_fixsensor.gif')
ani.save(out_gif, writer=animation.PillowWriter(fps=5), dpi=120)
plt.close(fig_a)
print(f"Saved: {out_gif}  ({time.time()-t0:.1f}s)")
print("\nDone ->", OUT)
