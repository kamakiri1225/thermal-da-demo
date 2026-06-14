"""
EnKF data assimilation demo: Miffy pixel art emerges from random noise.

Key design choices (after testing):
  1. Random sensor locations EACH CYCLE  (not fixed)
     Reason: fixed sensors + smooth prior → boundary cells get wrong-direction
     updates that accumulate. Random sensors average out these errors.
  2. Spatially correlated prior (Gaussian-filtered noise, L=5 cells)
     Reason: iid uniform noise → zero off-diagonal covariance → K[i,j]≈0
     for non-sensor cells → no spatial propagation of information.
  3. Post-update multiplicative inflation (alpha=1.05)
     Reason: without inflation, ensemble spread collapses → K→0 → filter
     freezes. Post-update inflation keeps the filter sensitive.
  4. Small localization radius (r_loc=8 cells)
     Reason: larger radius (12 cells) + smooth prior → boundary confusion.

State    : n = 2400  (40x60 temperature field, flattened)
Ensemble : N = 50 members, initialized from spatially correlated noise (mean=50°C)
Obs      : m = 30 RANDOM sensors per cycle (resampled every cycle)
DA method: Stochastic EnKF + Gaussian covariance localization (r=8 cells)
           + post-update multiplicative inflation (alpha=1.05)
Cycles   : 50

Outputs:
  img/fig01_enkf_snapshots.png   snapshot grid + RMSE panel
  img/fig02_rmse_convergence.png RMSE convergence (standalone)
  img/anim_enkf_miffy.gif        animation of Miffy emerging
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
import os
from scipy.ndimage import gaussian_filter

matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.use('Agg')  # headless rendering

# ============================================================
# Grid constants
# ============================================================
NX, NY = 60, 60
n = NX * NY   # 3600

# ============================================================
# Miffy true field  (identical to A01_miffy)
# ============================================================
T_BG = 0.0; T_OUTLINE = 8.0; T_FUR = 100.0
T_EYE = 5.0; T_MOUTH = 5.0; T_DRESS = 55.0; T_COLLAR = 70.0

def fill_rect(f, y1, y2, x1, x2, val):
    f[max(0, y1):min(NY-1, y2)+1, max(0, x1):min(NX-1, x2)+1] = val

def fill_ellipse(f, cy, cx, ry, rx, val):
    yy, xx = np.ogrid[:NY, :NX]
    f[((yy - cy) / ry)**2 + ((xx - cx) / rx)**2 <= 1.0] = val

def make_miffy():
    f = np.full((NY, NX), T_BG)
    fill_ellipse(f,  5, 30,  9, 21, T_DRESS)
    fill_rect(f,  0, 11,  7, 52, T_DRESS)
    fill_rect(f, 10, 12, 24, 36, T_COLLAR)
    fill_ellipse(f, 48, 19, 12,  9, T_OUTLINE)
    fill_ellipse(f, 48, 40, 12,  9, T_OUTLINE)
    fill_ellipse(f, 24, 30, 15, 24, T_OUTLINE)
    fill_ellipse(f, 48, 19, 11,  7, T_FUR)
    fill_ellipse(f, 48, 40, 11,  7, T_FUR)
    fill_rect(f, 35, 40, 12, 27, T_FUR)
    fill_rect(f, 35, 40, 33, 48, T_FUR)
    fill_ellipse(f, 24, 30, 14, 22, T_FUR)
    fill_rect(f,  9, 13, 25, 34, T_FUR)
    fill_ellipse(f, 24, 21,  2,  3, T_EYE)
    fill_ellipse(f, 24, 39,  2,  3, T_EYE)
    for d in range(-2, 3):
        f[min(max(18+d, 0), NY-1), min(max(30+d, 0), NX-1)] = T_MOUTH
        f[min(max(18-d, 0), NY-1), min(max(30+d, 0), NX-1)] = T_MOUTH
    return f

x_true_2d = make_miffy()
x_true    = x_true_2d.flatten()

# ============================================================
# EnKF parameters
# ============================================================
rng     = np.random.default_rng(42)
N_ENS   = 50     # ensemble size
M_OBS   = 30     # sensors per cycle (resampled randomly each cycle)
SIGMA_R = 5.0    # observation error std [°C]
N_CYC   = 50     # DA cycles
R_LOC   = 8.0    # Gaussian localization radius [cells]
SIGMA_B = 30.0   # background uncertainty [°C]
CORR_L  = 5.0    # prior spatial correlation length [cells]
INFL    = 1.05   # post-update multiplicative inflation

iy_all = np.arange(n) // NX
ix_all = np.arange(n) % NX

print(f"EnKF: n={n}, N={N_ENS}, m={M_OBS}/cycle (random), "
      f"sigma_r={SIGMA_R}, r_loc={R_LOC}, inflation={INFL}")
print(f"Prior: mean=50°C ± {SIGMA_B}°C (corr_L={CORR_L} cells)")

# ============================================================
# Initialize ensemble with spatially correlated noise
# ============================================================
X = np.zeros((n, N_ENS))
for i in range(N_ENS):
    noise     = rng.normal(0.0, 1.0, size=(NY, NX))
    noise_cor = gaussian_filter(noise, sigma=CORR_L)
    noise_cor = noise_cor / (noise_cor.std() + 1e-8) * SIGMA_B
    X[:, i]   = (50.0 + noise_cor).flatten()

# ============================================================
# One EnKF cycle
# ============================================================
def enkf_step(X):
    N = X.shape[1]

    # Random sensor locations this cycle
    obs_idx = np.sort(rng.choice(n, size=M_OBS, replace=False))
    obs_iy  = obs_idx // NX
    obs_ix  = obs_idx % NX
    R_diag  = np.full(M_OBS, SIGMA_R**2)

    # Localization for this cycle's sensors
    dy  = iy_all[:, None] - obs_iy[None, :]   # (n, m)
    dx  = ix_all[:, None] - obs_ix[None, :]
    rho = np.exp(-0.5 * (dy**2 + dx**2) / R_LOC**2)

    # Observations from true field + noise
    y = x_true[obs_idx] + rng.normal(0.0, SIGMA_R, size=M_OBS)

    # Ensemble anomalies
    xbar = X.mean(axis=1, keepdims=True)
    Xp   = X - xbar                          # (n, N)

    # Project to obs space
    HXp = Xp[obs_idx, :]                     # (m, N)

    # Innovation covariance  (m x m)
    S = (HXp @ HXp.T) / (N - 1) + np.diag(R_diag)

    # Localized Kalman gain  (n x m)
    K_loc = (Xp @ HXp.T) / (N - 1) @ np.linalg.inv(S) * rho

    # Stochastic EnKF update
    eps   = rng.normal(0.0, SIGMA_R, size=(M_OBS, N))
    innov = y[:, None] + eps - X[obs_idx, :]
    X_new = X + K_loc @ innov

    # Post-update multiplicative inflation (prevents ensemble collapse)
    xbar_new = X_new.mean(axis=1, keepdims=True)
    X_new    = xbar_new + INFL * (X_new - xbar_new)

    return X_new, obs_idx

# ============================================================
# Run DA cycles
# ============================================================
all_means    = []
rmse_hist    = []
obs_per_cycle = []  # sensor locations for animation

def snapshot(X, obs_idx=None):
    xm   = X.mean(axis=1)
    rmse = np.sqrt(np.mean((xm - x_true)**2))
    all_means.append(xm.reshape(NY, NX).copy())
    rmse_hist.append(rmse)
    obs_per_cycle.append(obs_idx)

snapshot(X, obs_idx=None)
print(f"\nCycle  0  RMSE = {rmse_hist[0]:.2f}°C  (correlated random initial state)")

for cyc in range(1, N_CYC + 1):
    X, obs_idx = enkf_step(X)
    snapshot(X, obs_idx)
    if cyc in (1, 5, 10, 20, 30, 50):
        print(f"Cycle {cyc:2d}  RMSE = {rmse_hist[cyc]:.2f}°C")

print(f"\nConvergence: {rmse_hist[0]:.1f} -> {rmse_hist[-1]:.1f}°C  "
      f"(reduction {(1-rmse_hist[-1]/rmse_hist[0])*100:.0f}%)")

# ============================================================
# Plotting helpers
# ============================================================
CMAP = 'RdYlBu_r'
VMIN, VMAX = 0, 100
FS_TITLE = 17
FS_LABEL = 15
FS_TICK  = 13
FS_LEGEND = 13

def style_ax(ax, ticks=False):
    ax.set_facecolor('#111111')
    if not ticks:
        ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_edgecolor('gray')

script_dir = os.path.dirname(os.path.abspath(__file__))
img_dir    = os.path.join(script_dir, '..', 'img')
os.makedirs(img_dir, exist_ok=True)

# ============================================================
# Figure 1: Snapshot grid  (8 cycles + true field + RMSE)
# ============================================================
SNAP_CYCLES = [0, 1, 5, 10, 20, 30, 40, 50]

fig1, axes1 = plt.subplots(2, 5, figsize=(22, 9))
fig1.patch.set_facecolor('#111111')
axes1 = axes1.flatten()

last_im = None
for ai, cyc in enumerate(SNAP_CYCLES):
    ax = axes1[ai]
    style_ax(ax)
    im = ax.imshow(all_means[cyc], origin='lower', cmap=CMAP,
                   vmin=VMIN, vmax=VMAX, extent=[0, NX, 0, NY],
                   aspect='equal', interpolation='nearest')
    if obs_per_cycle[cyc] is not None:
        oi = obs_per_cycle[cyc]
        ax.scatter(oi % NX + 0.5, oi // NX + 0.5, c='lime', s=14,
                   marker='x', linewidths=1.2, alpha=0.9)
    ax.set_title(f'Cycle {cyc}   RMSE={rmse_hist[cyc]:.1f}°C',
                 fontsize=FS_TITLE, color='white', fontweight='bold')
    last_im = im

# Panel 9: True field
ax9 = axes1[8]
style_ax(ax9)
ax9.imshow(x_true_2d, origin='lower', cmap=CMAP,
           vmin=VMIN, vmax=VMAX, extent=[0, NX, 0, NY],
           aspect='equal', interpolation='nearest')
ax9.set_title('True Field (Miffy)', fontsize=FS_TITLE, color='white', fontweight='bold')

# Panel 10: RMSE
ax10 = axes1[9]
ax10.set_facecolor('#111111')
cyc_arr = np.arange(N_CYC + 1)
ax10.fill_between(cyc_arr, rmse_hist, alpha=0.15, color='white')
ax10.plot(cyc_arr, rmse_hist, 'w-', lw=2)
for cyc in SNAP_CYCLES:
    ax10.axvline(cyc, color='cyan', lw=0.8, alpha=0.5)
ax10.set_xlabel('DA Cycle', color='white', fontsize=FS_LABEL)
ax10.set_ylabel('RMSE [°C]', color='white', fontsize=FS_LABEL)
ax10.set_title('RMSE Convergence', fontsize=FS_TITLE, color='white', fontweight='bold')
ax10.tick_params(colors='white', labelsize=FS_TICK)
ax10.grid(alpha=0.3, color='gray', lw=0.5)
for sp in ax10.spines.values():
    sp.set_edgecolor('gray')

cbar1 = fig1.colorbar(last_im, ax=axes1[:9].tolist(),
                      orientation='vertical', fraction=0.012, pad=0.01)
cbar1.set_label('Temperature [°C]', color='white', fontsize=FS_LABEL)
cbar1.ax.tick_params(colors='white', labelsize=FS_TICK)
cbar1.ax.yaxis.label.set_color('white')
cbar1.outline.set_edgecolor('gray')

plt.suptitle(
    f'EnKF: Miffy Emerges from Random Noise  '
    f'(N={N_ENS}, m={M_OBS}/cycle random, σ_r={SIGMA_R}, r_loc={R_LOC} cells)',
    fontsize=19, fontweight='bold', color='white')
plt.tight_layout()
out1 = os.path.join(img_dir, 'fig01_enkf_snapshots.png')
fig1.savefig(out1, dpi=150, bbox_inches='tight', facecolor=fig1.get_facecolor())
print(f"\nSaved: {out1}")
plt.close(fig1)

# ============================================================
# Figure 2: RMSE convergence (standalone)
# ============================================================
fig2, ax2 = plt.subplots(figsize=(11, 4.8))
fig2.patch.set_facecolor('#111111')
style_ax(ax2, ticks=True)
ax2.fill_between(cyc_arr, rmse_hist, alpha=0.15, color='white')
ax2.plot(cyc_arr, rmse_hist, 'w-', lw=2.5, label='Ensemble mean RMSE')
for cyc in SNAP_CYCLES:
    ax2.axvline(cyc, color='cyan', lw=0.8, alpha=0.6,
                label='Snapshot' if cyc == SNAP_CYCLES[0] else '')
ax2.set_xlabel('DA Cycle', color='white', fontsize=FS_LABEL)
ax2.set_ylabel('RMSE [°C]', color='white', fontsize=FS_LABEL)
ax2.set_title(
    f'EnKF Convergence  (N={N_ENS}, m={M_OBS}/cycle, σ_r={SIGMA_R}, r_loc={R_LOC})',
    color='white', fontsize=FS_TITLE, fontweight='bold')
ax2.tick_params(colors='white', labelsize=FS_TICK)
ax2.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=FS_LEGEND)
ax2.grid(alpha=0.3, color='gray', lw=0.5)
plt.tight_layout()
out2 = os.path.join(img_dir, 'fig02_rmse_convergence.png')
fig2.savefig(out2, dpi=150, bbox_inches='tight', facecolor=fig2.get_facecolor())
print(f"Saved: {out2}")
plt.close(fig2)

# ============================================================
# Animation: Miffy emerging side-by-side with true field
# ============================================================
fig_a, (ax_da, ax_tr) = plt.subplots(1, 2, figsize=(14, 9.5))
fig_a.patch.set_facecolor('#111111')
style_ax(ax_da); style_ax(ax_tr)

im_da = ax_da.imshow(all_means[0], origin='lower', cmap=CMAP,
                     vmin=VMIN, vmax=VMAX, extent=[0, NX, 0, NY],
                     aspect='equal', interpolation='nearest')
sensor_scatter = ax_da.scatter([], [], c='lime', s=22,
                               marker='x', linewidths=1.5, alpha=0.9,
                               label='Sensors (30/cycle)')
ax_da.legend(facecolor='black', edgecolor='gray', labelcolor='white',
             fontsize=FS_LEGEND, loc='lower right')
ax_da.set_title('EnKF Ensemble Mean', color='white', fontweight='bold', fontsize=FS_TITLE)

ax_tr.imshow(x_true_2d, origin='lower', cmap=CMAP,
             vmin=VMIN, vmax=VMAX, extent=[0, NX, 0, NY],
             aspect='equal', interpolation='nearest')
ax_tr.set_title('True State (Miffy)', color='white', fontweight='bold', fontsize=FS_TITLE)

cbar_a = fig_a.colorbar(im_da, ax=[ax_da, ax_tr],
                        orientation='vertical', fraction=0.03, pad=0.02)
cbar_a.set_label('Temperature [°C]', color='white', fontsize=FS_LABEL)
cbar_a.ax.tick_params(colors='white', labelsize=FS_TICK)
cbar_a.ax.yaxis.label.set_color('white')
cbar_a.outline.set_edgecolor('gray')

sup = fig_a.suptitle(
    f'Cycle  0  |  RMSE = {rmse_hist[0]:.1f}°C',
    fontsize=19, color='white', fontweight='bold')
plt.tight_layout()

def update(frame):
    im_da.set_data(all_means[frame])
    oi = obs_per_cycle[frame]
    if oi is not None:
        sensor_scatter.set_offsets(
            np.c_[oi % NX + 0.5, oi // NX + 0.5])
    sup.set_text(
        f'Cycle {frame:2d}  |  RMSE = {rmse_hist[frame]:.1f}°C  '
        f'|  N={N_ENS}, m={M_OBS}/cycle, σ_r={SIGMA_R}')
    return [im_da, sensor_scatter, sup]

ani = animation.FuncAnimation(
    fig_a, update, frames=N_CYC + 1, interval=300, blit=False)
out_gif = os.path.join(img_dir, 'anim_enkf_miffy.gif')
ani.save(out_gif, writer='pillow', fps=4, dpi=100)
print(f"Saved animation: {out_gif}")
plt.close(fig_a)

# Save RMSE history for fig05 comparison in oi_miffy.py
np.save(os.path.join(img_dir, 'enkf_rmse.npy'), np.array(rmse_hist))
print("\nDone.")
