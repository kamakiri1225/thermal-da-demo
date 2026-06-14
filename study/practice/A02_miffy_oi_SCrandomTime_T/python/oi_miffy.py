"""
OI (Optimal Interpolation) data assimilation demo: Miffy emerges from noise.

OI uses an ANALYTIC background covariance B (Gaussian spatial correlation),
unlike EnKF which estimates B from the ensemble. This makes OI:
  + Guaranteed convergent (no sampling error)
  + Simpler (no ensemble needed)
  - Fixed covariance structure (cannot adapt to actual field)
  - No uncertainty quantification (no spread estimate)

Algorithm per cycle:
  1. Sample m sensors (random locations each cycle)
  2. BHt (n x m) = sigma_b^2 * exp(-d^2 / 2L^2)  [grid-to-sensor covariance]
  3. HBHt (m x m) = sigma_b^2 * exp(-d_obs^2 / 2L^2)  [sensor-sensor covariance]
  4. K = BHt @ inv(HBHt + R)                          [Kalman gain]
  5. x^a = x^b + K @ (y - Hx^b)                       [OI analysis]

State    : n = 2400  (40x60), single estimate (no ensemble)
Obs      : m = 30 RANDOM sensors per cycle
Method   : Iterative OI, sigma_b decays each cycle
Cycles   : 50
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
import os
from scipy.ndimage import gaussian_filter

matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.use('Agg')

NX, NY = 60, 60
n = NX * NY

T_BG=0.0; T_OUTLINE=8.0; T_FUR=100.0; T_EYE=5.0; T_MOUTH=5.0; T_DRESS=55.0; T_COLLAR=70.0

def fill_rect(f, y1, y2, x1, x2, val):
    f[max(0,y1):min(NY-1,y2)+1, max(0,x1):min(NX-1,x2)+1] = val
def fill_ellipse(f, cy, cx, ry, rx, val):
    yy,xx = np.ogrid[:NY,:NX]; f[((yy-cy)/ry)**2+((xx-cx)/rx)**2<=1.0] = val
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

rng = np.random.default_rng(42)

# OI parameters
M_OBS    = 60     # sensors per cycle (random)
SIGMA_R  = 5.0    # observation error [°C]
N_CYC    = 50     # cycles
SIGMA_B0 = 30.0   # initial background error [°C]
SIGMA_B_MIN = 3.0 # minimum background error
DECAY    = 0.88   # sigma_b decay rate per cycle
L_CORR   = 5.0    # spatial correlation length [cells]

iy_all = np.arange(n) // NX
ix_all = np.arange(n) % NX

# Initial state: smooth noise centered at 50°C
noise = rng.normal(0,1,(NY,NX))
noise_cor = gaussian_filter(noise, sigma=5.0)
x_est = (50.0 + noise_cor/(noise_cor.std()+1e-8)*SIGMA_B0).flatten()

print(f"OI: n={n}, m={M_OBS}/cycle (random), sigma_r={SIGMA_R}, L={L_CORR} cells")
print(f"sigma_b: {SIGMA_B0}°C → {SIGMA_B_MIN}°C (decay={DECAY})")

all_states   = [x_est.reshape(NY,NX).copy()]
rmse_hist    = [np.sqrt(np.mean((x_est - x_true)**2))]
obs_per_cycle = [None]
print(f"\nCycle  0  RMSE = {rmse_hist[0]:.2f}°C")

sigma_b = SIGMA_B0

for cyc in range(1, N_CYC+1):
    # Random sensors
    obs_idx = np.sort(rng.choice(n, size=M_OBS, replace=False))
    obs_iy  = obs_idx // NX
    obs_ix  = obs_idx % NX

    # Current background error
    sb = max(SIGMA_B0 * DECAY**cyc, SIGMA_B_MIN)

    # BHt (n x m): grid-to-sensor Gaussian covariance
    dy_all = iy_all[:, None] - obs_iy[None, :]
    dx_all = ix_all[:, None] - obs_ix[None, :]
    BHt = sb**2 * np.exp(-0.5 * (dy_all**2 + dx_all**2) / L_CORR**2)  # (n, m)

    # HBHt (m x m): sensor-sensor covariance + obs error
    dy_obs = obs_iy[:, None] - obs_iy[None, :]
    dx_obs = obs_ix[:, None] - obs_ix[None, :]
    HBHt = sb**2 * np.exp(-0.5 * (dy_obs**2 + dx_obs**2) / L_CORR**2)  # (m, m)
    HBHt += np.diag(np.full(M_OBS, SIGMA_R**2))                          # + R

    # Kalman gain  (n x m)
    K = BHt @ np.linalg.inv(HBHt)

    # Observations + OI update
    y     = x_true[obs_idx] + rng.normal(0.0, SIGMA_R, M_OBS)
    innov = y - x_est[obs_idx]
    x_est = x_est + K @ innov

    rmse = np.sqrt(np.mean((x_est - x_true)**2))
    all_states.append(x_est.reshape(NY, NX).copy())
    rmse_hist.append(rmse)
    obs_per_cycle.append(obs_idx)

    if cyc in (1, 5, 10, 20, 30, 50):
        print(f"Cycle {cyc:2d}  RMSE = {rmse:.2f}°C  (sigma_b={sb:.1f}°C)")

print(f"\nConvergence: {rmse_hist[0]:.1f} -> {rmse_hist[-1]:.1f}°C  "
      f"(reduction {(1-rmse_hist[-1]/rmse_hist[0])*100:.0f}%)")

# ============================================================
# Plotting
# ============================================================
CMAP = 'RdYlBu_r'; VMIN, VMAX = 0, 100
FS_TITLE = 17
FS_LABEL = 15
FS_TICK  = 13
FS_LEGEND = 13

def style_ax(ax, ticks=False):
    ax.set_facecolor('#111111')
    if not ticks: ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_edgecolor('gray')

script_dir = os.path.dirname(os.path.abspath(__file__))
img_dir    = os.path.join(script_dir, '..', 'img')
os.makedirs(img_dir, exist_ok=True)

SNAP_CYCLES = [0, 1, 5, 10, 20, 30, 40, 50]

# Figure 1: Snapshots
fig1, axes1 = plt.subplots(2, 5, figsize=(22, 9))
fig1.patch.set_facecolor('#111111')
axes1 = axes1.flatten()
last_im = None
for ai, cyc in enumerate(SNAP_CYCLES):
    ax = axes1[ai]; style_ax(ax)
    im = ax.imshow(all_states[cyc], origin='lower', cmap=CMAP,
                   vmin=VMIN, vmax=VMAX, extent=[0,NX,0,NY],
                   aspect='equal', interpolation='nearest')
    oi = obs_per_cycle[cyc]
    if oi is not None:
        ax.scatter(oi%NX+0.5, oi//NX+0.5, c='lime', s=14,
                   marker='x', linewidths=1.2, alpha=0.9)
    ax.set_title(f'Cycle {cyc}   RMSE={rmse_hist[cyc]:.1f}°C',
                 fontsize=FS_TITLE, color='white', fontweight='bold')
    last_im = im

ax9 = axes1[8]; style_ax(ax9)
ax9.imshow(x_true_2d, origin='lower', cmap=CMAP, vmin=VMIN, vmax=VMAX,
           extent=[0,NX,0,NY], aspect='equal', interpolation='nearest')
ax9.set_title('True Field (Miffy)', fontsize=FS_TITLE, color='white', fontweight='bold')

ax10 = axes1[9]; ax10.set_facecolor('#111111')
cyc_arr = np.arange(N_CYC+1)
ax10.fill_between(cyc_arr, rmse_hist, alpha=0.15, color='white')
ax10.plot(cyc_arr, rmse_hist, 'w-', lw=2)
for cyc in SNAP_CYCLES: ax10.axvline(cyc, color='cyan', lw=0.8, alpha=0.5)
ax10.set_xlabel('DA Cycle', color='white', fontsize=FS_LABEL)
ax10.set_ylabel('RMSE [°C]', color='white', fontsize=FS_LABEL)
ax10.set_title('RMSE Convergence', fontsize=FS_TITLE, color='white', fontweight='bold')
ax10.tick_params(colors='white', labelsize=FS_TICK)
ax10.grid(alpha=0.3, color='gray', lw=0.5)
for sp in ax10.spines.values(): sp.set_edgecolor('gray')

cbar1 = fig1.colorbar(last_im, ax=axes1[:9].tolist(),
                      orientation='vertical', fraction=0.012, pad=0.01)
cbar1.set_label('Temperature [°C]', color='white', fontsize=FS_LABEL)
cbar1.ax.tick_params(colors='white', labelsize=FS_TICK); cbar1.ax.yaxis.label.set_color('white')
cbar1.outline.set_edgecolor('gray')

plt.suptitle(
    f'OI: Miffy Emerges from Random Noise  '
    f'(m={M_OBS}/cycle random, σ_r={SIGMA_R}, L={L_CORR} cells)',
    fontsize=19, fontweight='bold', color='white')
plt.tight_layout()
out1 = os.path.join(img_dir, 'fig03_oi_snapshots.png')
fig1.savefig(out1, dpi=150, bbox_inches='tight', facecolor=fig1.get_facecolor())
print(f"\nSaved: {out1}")
plt.close(fig1)

# Figure 2: RMSE
fig2, ax2 = plt.subplots(figsize=(11, 4.8))
fig2.patch.set_facecolor('#111111')
style_ax(ax2, ticks=True)
ax2.fill_between(cyc_arr, rmse_hist, alpha=0.15, color='white')
ax2.plot(cyc_arr, rmse_hist, 'w-', lw=2.5, label='OI analysis RMSE')
for cyc in SNAP_CYCLES:
    ax2.axvline(cyc, color='cyan', lw=0.8, alpha=0.6,
                label='Snapshot' if cyc==SNAP_CYCLES[0] else '')
ax2.set_xlabel('DA Cycle', color='white', fontsize=FS_LABEL)
ax2.set_ylabel('RMSE [°C]', color='white', fontsize=FS_LABEL)
ax2.set_title(f'OI Convergence  (m={M_OBS}/cycle, σ_r={SIGMA_R}, L={L_CORR})',
              color='white', fontsize=FS_TITLE, fontweight='bold')
ax2.tick_params(colors='white', labelsize=FS_TICK)
ax2.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=FS_LEGEND)
ax2.grid(alpha=0.3, color='gray', lw=0.5)
plt.tight_layout()
out2 = os.path.join(img_dir, 'fig04_oi_rmse.png')
fig2.savefig(out2, dpi=150, bbox_inches='tight', facecolor=fig2.get_facecolor())
print(f"Saved: {out2}")
plt.close(fig2)

# Animation
fig_a, (ax_da, ax_tr) = plt.subplots(1, 2, figsize=(14, 9.5))
fig_a.patch.set_facecolor('#111111')
style_ax(ax_da); style_ax(ax_tr)
im_da = ax_da.imshow(all_states[0], origin='lower', cmap=CMAP,
                     vmin=VMIN, vmax=VMAX, extent=[0,NX,0,NY],
                     aspect='equal', interpolation='nearest')
sc = ax_da.scatter([], [], c='lime', s=22, marker='x',
                   linewidths=1.5, alpha=0.9, label=f'Sensors ({M_OBS}/cycle)')
ax_da.legend(facecolor='black', edgecolor='gray', labelcolor='white', fontsize=FS_LEGEND, loc='lower right')
ax_da.set_title('OI Analysis', color='white', fontweight='bold', fontsize=FS_TITLE)
ax_tr.imshow(x_true_2d, origin='lower', cmap=CMAP, vmin=VMIN, vmax=VMAX,
             extent=[0,NX,0,NY], aspect='equal', interpolation='nearest')
ax_tr.set_title('True State (Miffy)', color='white', fontweight='bold', fontsize=FS_TITLE)
cbar_a = fig_a.colorbar(im_da, ax=[ax_da,ax_tr], orientation='vertical',
                        fraction=0.03, pad=0.02)
cbar_a.set_label('Temperature [°C]', color='white', fontsize=FS_LABEL)
cbar_a.ax.tick_params(colors='white', labelsize=FS_TICK); cbar_a.ax.yaxis.label.set_color('white')
cbar_a.outline.set_edgecolor('gray')
sup = fig_a.suptitle(f'Cycle  0  |  RMSE = {rmse_hist[0]:.1f}°C',
                     fontsize=19, color='white', fontweight='bold')
plt.tight_layout()

def update(frame):
    im_da.set_data(all_states[frame])
    oi = obs_per_cycle[frame]
    if oi is not None:
        sc.set_offsets(np.c_[oi%NX+0.5, oi//NX+0.5])
    sup.set_text(f'Cycle {frame:2d}  |  RMSE = {rmse_hist[frame]:.1f}°C  '
                 f'|  OI, m={M_OBS}/cycle, L={L_CORR} cells')
    return [im_da, sc, sup]

ani = animation.FuncAnimation(fig_a, update, frames=N_CYC+1, interval=300, blit=False)
out_gif = os.path.join(img_dir, 'anim_oi_miffy.gif')
ani.save(out_gif, writer='pillow', fps=4, dpi=100)
print(f"Saved animation: {out_gif}")
plt.close(fig_a)

# ── fig05: EnKF vs OI RMSE comparison ─────────────────────────────────────
enkf_rmse_path = os.path.join(script_dir, '..', '..', 'A02_miffy_enkf_T', 'img', 'enkf_rmse.npy')
if os.path.exists(enkf_rmse_path):
    enkf_rmse = np.load(enkf_rmse_path)
    cyc_arr_e = np.arange(len(enkf_rmse))
    cyc_arr_o = np.arange(len(rmse_hist))

    fig5, ax5 = plt.subplots(figsize=(12, 6))
    fig5.patch.set_facecolor('#111111')
    style_ax(ax5, ticks=True)

    ax5.plot(cyc_arr_e, enkf_rmse, color='#4FC3F7', lw=2.5,
             label=f'EnKF  (N=50, m=30/cycle, r_loc=8, α=1.05)')
    ax5.plot(cyc_arr_o, rmse_hist, color='#FF8A65', lw=2.5, ls='--',
             label=f'OI  (m={M_OBS}/cycle, L={L_CORR} cells, σ_b decay={DECAY})')
    ax5.fill_between(cyc_arr_e, enkf_rmse, alpha=0.08, color='#4FC3F7')
    ax5.fill_between(cyc_arr_o, rmse_hist, alpha=0.08, color='#FF8A65')

    ax5.set_xlabel('DA Cycle', color='white', fontsize=FS_LABEL)
    ax5.set_ylabel('RMSE  [°C]', color='white', fontsize=FS_LABEL)
    ax5.set_title('EnKF vs OI — RMSE Convergence', color='white',
                  fontsize=FS_TITLE, fontweight='bold')
    ax5.tick_params(colors='white', labelsize=FS_TICK)
    ax5.legend(facecolor='#222222', edgecolor='gray', labelcolor='white',
               fontsize=FS_LEGEND)
    ax5.grid(alpha=0.3, color='gray', lw=0.5)
    ax5.annotate(f'EnKF Cycle 50: {enkf_rmse[-1]:.1f}°C',
                 xy=(50, enkf_rmse[-1]), xytext=(38, enkf_rmse[-1] + 5),
                 color='#4FC3F7', fontsize=FS_TICK,
                 arrowprops=dict(arrowstyle='->', color='#4FC3F7'))
    ax5.annotate(f'OI Cycle 50: {rmse_hist[-1]:.1f}°C',
                 xy=(50, rmse_hist[-1]), xytext=(35, rmse_hist[-1] + 8),
                 color='#FF8A65', fontsize=FS_TICK,
                 arrowprops=dict(arrowstyle='->', color='#FF8A65'))
    plt.tight_layout()
    out5 = os.path.join(img_dir, 'fig05_enkf_vs_oi.png')
    fig5.savefig(out5, dpi=150, bbox_inches='tight', facecolor=fig5.get_facecolor())
    print(f"Saved: {out5}")
    plt.close(fig5)
else:
    print("Skip fig05: run enkf_miffy.py first to generate enkf_rmse.npy")

print("\nDone.")
