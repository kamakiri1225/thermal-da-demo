#!/usr/bin/env python3
"""
enkf_heatid.py -- Case 2: Heat identification via ESMDA

Algorithm: deterministic ESMDA / optimal interpolation (Gaussian-process regression).

This is a static linear direct-observation problem, so the Gaussian posterior mean
can be evaluated analytically. The animation follows the ESMDA likelihood-tempering
path from 0% to 100% observation information. This removes finite-ensemble sampling
noise and cannot diverge.

State  : x (temperature field, n=900 cells)
Forward: H = I_obs (direct temperature observation)
Prior  : uniform mean + shape-free stationary Matern-3/2 covariance
Sensors: geometry-only space-filling layout (does not use the Miffy truth)
Output : top row = x_pred, bottom row = q_derived = -alpha * L * x_mean
Grid   : 30x30 (n=900), N_ITER=30
"""
import numpy as np
import scipy.sparse as sp
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os, time

script_dir = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(script_dir, '..', 'img')
os.makedirs(OUT, exist_ok=True)

OUTPUT_FILES = [
    'fig00_sensor_layout.png',
    'fig01_rmse.png',
    'fig02_xpred_final.png',
    'fig03_qest_final.png',
    'fig04_esmda_oi_concept.png',
    'anim_enkf_heatid.gif',
]

for filename in OUTPUT_FILES:
    path = os.path.join(OUT, filename)
    if os.path.isfile(path):
        os.remove(path)
print(f"Initialized output directory: {os.path.abspath(OUT)}")

# ── Grid & physics ──────────────────────────────────────
NX, NY = 30, 30
n      = NX * NY
ALPHA  = 0.5

def make_miffy():
    f = np.full((NY, NX), 0.)
    s = 30 / 60.0
    def fe(cy, cx, ry, rx, v):
        yy, xx = np.ogrid[:NY, :NX]
        f[((yy-int(cy*s))/max(int(ry*s),1))**2
          +((xx-int(cx*s))/max(int(rx*s),1))**2 <= 1.0] = v
    def fr(y1, y2, x1, x2, v):
        f[max(0,int(y1*s)):min(NY-1,int(y2*s))+1,
          max(0,int(x1*s)):min(NX-1,int(x2*s))+1] = v
    fe( 5,30, 9,21, 55.); fr( 0,11, 7,52, 55.); fr(10,12,24,36, 70.)
    fe(48,19,12, 9,  8.); fe(48,40,12, 9,  8.); fe(24,30,15,24,  8.)
    fe(48,19,11, 7,100.); fe(48,40,11, 7,100.)
    fr(35,40,12,27,100.); fr(35,40,33,48,100.)
    fe(24,30,14,22,100.); fr( 9,13,25,34,100.)
    fe(24,21, 2, 3,  5.); fe(24,39, 2, 3,  5.)
    for d in range(-1, 2):
        y_ = min(max(int(18*s)+d,0),NY-1)
        x_ = min(max(int(30*s)+d,0),NX-1)
        f[y_,x_] = 5.
        f[min(max(int(18*s)-d,0),NY-1),x_] = 5.
    return f

def build_laplacian(ny, nx):
    nn = ny*nx
    row, col, data = [], [], []
    for iy in range(ny):
        for ix in range(nx):
            k = iy*nx+ix
            row.append(k); col.append(k); data.append(-4.)
            for diy, dix in [(-1,0),(1,0),(0,-1),(0,1)]:
                jy, jx = iy+diy, ix+dix
                if 0 <= jy < ny and 0 <= jx < nx:
                    row.append(k); col.append(jy*nx+jx); data.append(1.)
    return sp.csr_matrix((data,(row,col)), shape=(nn,nn), dtype=np.float64)

x_ss   = make_miffy().flatten()
L      = build_laplacian(NY, NX)
q_true = -ALPHA * (L @ x_ss)

iy_all = np.arange(n) // NX
ix_all = np.arange(n) % NX

print(f"Grid: {NY}x{NX}={n}")
print(f"q_true: std={q_true.std():.1f}  max={q_true.max():.1f}  min={q_true.min():.1f}")

# ── Truth-independent, nested sensor positions ──────────
def make_space_filling_sensor_order():
    """Greedy farthest-point order based only on grid coordinates."""
    coords = np.column_stack((iy_all, ix_all))
    center = (NY//2)*NX + NX//2
    selected = [center]
    min_dist2 = np.sum((coords - coords[center])**2, axis=1).astype(float)
    min_dist2[center] = -1.0

    for _ in range(1, n):
        idx = int(np.argmax(min_dist2))
        selected.append(idx)
        dist2 = np.sum((coords - coords[idx])**2, axis=1)
        min_dist2 = np.minimum(min_dist2, dist2)
        min_dist2[selected] = -1.0
    return np.asarray(selected, dtype=int)

SENSOR_ORDER = make_space_filling_sensor_order()

def sensor_subsets(m):
    return SENSOR_ORDER[:m].copy()

M_LIST = [30, 100, 300, 600]
COLORS = ['#FF6B6B', '#FFD93D', '#6BCB77', '#4D96FF']

# Print sensor temperatures
print("\nSensor temperatures (x_ss at each position):")
for m in M_LIST:
    obs = sensor_subsets(m)
    vals = x_ss[obs]
    print(f"  m={m:3d}: min={vals.min():.0f}, max={vals.max():.0f}, "
          f"unique={np.unique(vals)}")

# ── Deterministic ESMDA / optimal interpolation ──────────
PRIOR_MEAN = 20.0  # uniform initial estimate; contains no Miffy shape
SIGMA_B = 40.0     # prior standard deviation [degC]
SIGMA_R = 3.0      # observation noise [degC]
CORR_L = 1.4       # Matern-3/2 correlation length [cells]
N_ITER = 30        # likelihood-tempering frames

coords = np.column_stack((iy_all, ix_all))
dy_all = coords[:, None, 0] - coords[None, :, 0]
dx_all = coords[:, None, 1] - coords[None, :, 1]
distance = np.sqrt(dy_all**2 + dx_all**2)
matern_arg = np.sqrt(3.0) * distance / CORR_L
PRIOR_COV = SIGMA_B**2 * (1.0 + matern_arg) * np.exp(-matern_arg)

print(f"\nDeterministic ESMDA/OI params: N_ITER={N_ITER}, "
      f"SIGMA_B={SIGMA_B}, SIGMA_R={SIGMA_R}, Matern L={CORR_L}")

def run_esmda(obs_idx):
    prior = np.full(n, PRIOR_MEAN)
    innovation = x_ss[obs_idx] - PRIOR_MEAN
    cov_xo = PRIOR_COV[:, obs_idx]
    cov_oo = PRIOR_COV[np.ix_(obs_idx, obs_idx)]

    # One eigendecomposition supports the full likelihood-tempering path:
    # x(t) = mu + C_xo [C_oo + R/t]^-1 (y - H mu), 0 < t <= 1.
    eigval, eigvec = np.linalg.eigh(cov_oo)
    projected_innovation = eigvec.T @ innovation

    rmse0 = np.sqrt(np.mean((prior - x_ss)**2))
    history = {0: {
        'x':      prior.copy(),
        'q':      -ALPHA * (L @ prior),
        'rmse_x': rmse0,
    }}

    for it in range(1, N_ITER+1):
        information_fraction = it / N_ITER
        denominator = eigval + SIGMA_R**2 / information_fraction
        weights = eigvec @ (projected_innovation / denominator)
        xmean = prior + cov_xo @ weights
        history[it] = {
            'x':      xmean.copy(),
            'q':      -ALPHA * (L @ xmean),
            'rmse_x': np.sqrt(np.mean((xmean - x_ss)**2)),
        }

    return history

# ── Run ─────────────────────────────────────────────────
print(f"\nRunning ESMDA ...")
all_hist = {}
all_obs  = {}
for m in M_LIST:
    obs_idx = sensor_subsets(m)
    print(f"  m={m:2d}", end=' ', flush=True)
    t0 = time.time()
    hist = run_esmda(obs_idx)
    all_hist[m] = hist
    all_obs[m]  = obs_idx
    print(f"-> {time.time()-t0:.1f}s  RMSE: {hist[0]['rmse_x']:.1f} -> {hist[N_ITER]['rmse_x']:.2f}°C")

CMAP_T = 'RdYlBu_r'
CMAP_Q = 'seismic'
FS_T   = 15; FS_L = 13; FS_TK = 11; FS_LEG = 12
cyc_arr = np.arange(N_ITER+1)

# Robust common color range: ignore the most extreme 1% so weak heat-source
# structure remains visible. Actual min/max values are still printed per panel.
final_q = np.concatenate([all_hist[m][N_ITER]['q'] for m in M_LIST])
Q_LIM = max(10.0, float(np.percentile(np.abs(final_q), 99.0)))
print(f"Heat visualization range: {-Q_LIM:.1f} to {Q_LIM:.1f} "
      f"(99th percentile; actual {final_q.min():.1f} to {final_q.max():.1f})")

# ── fig04: deterministic ESMDA/OI concept ───────────────
# One-dimensional example using the same Matern-3/2 covariance as the 2-D case.
x_line = np.linspace(0.0, 29.0, 300)
sensor_x = np.array([3.0, 10.0, 18.0, 26.0])
sensor_y = np.array([0.0, 85.0, 35.0, 100.0])

def matern32_cov(a, b):
    dist_ = np.abs(a[:, None] - b[None, :])
    z_ = np.sqrt(3.0) * dist_ / CORR_L
    return SIGMA_B**2 * (1.0 + z_) * np.exp(-z_)

cov_line_obs = matern32_cov(x_line, sensor_x)
cov_obs_obs = matern32_cov(sensor_x, sensor_x)
concept_curves = []
for fraction in [0.0, 0.15, 0.45, 1.0]:
    if fraction == 0.0:
        estimate = np.full_like(x_line, PRIOR_MEAN)
    else:
        weights = np.linalg.solve(
            cov_obs_obs + (SIGMA_R**2/fraction)*np.eye(len(sensor_x)),
            sensor_y - PRIOR_MEAN,
        )
        estimate = PRIOR_MEAN + cov_line_obs @ weights
    concept_curves.append((fraction, estimate))

fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(13, 5.2))
fig4.patch.set_facecolor('#111111')
for ax_ in (ax4a, ax4b):
    ax_.set_facecolor('#111111')
    ax_.tick_params(colors='white')
    ax_.grid(alpha=0.2, color='gray')
    for sp in ax_.spines.values():
        sp.set_edgecolor('#555555')

curve_colors = ['#74ADD1', '#FFD166', '#F8961E', '#EF476F']
for (fraction, estimate), color in zip(concept_curves, curve_colors):
    label = 'prior (0%)' if fraction == 0.0 else f'observation information {fraction:.0%}'
    ax4a.plot(x_line, estimate, color=color, lw=2.5, label=label)
ax4a.scatter(sensor_x, sensor_y, s=80, c='#00E676', marker='o',
             edgecolors='white', linewidths=1.2, zorder=5, label='sensors')
ax4a.set_ylim(-8, 110)
ax4a.set_xlabel('Position', color='white')
ax4a.set_ylabel('Temperature [degC]', color='white')
ax4a.set_title('Likelihood tempering: prior -> posterior',
               color='white', fontweight='bold')
ax4a.legend(facecolor='#222222', edgecolor='gray', labelcolor='white',
            fontsize=9, loc='upper left')

distance_line = np.linspace(0.0, 8.0, 300)
z_line = np.sqrt(3.0) * distance_line / CORR_L
corr_line = (1.0 + z_line) * np.exp(-z_line)
ax4b.plot(distance_line, corr_line, color='#4D96FF', lw=3)
ax4b.axvline(CORR_L, color='#FFD166', ls='--', lw=2,
             label=f'correlation length L={CORR_L}')
ax4b.fill_between(distance_line, 0, corr_line, color='#4D96FF', alpha=0.2)
ax4b.set_ylim(0, 1.05)
ax4b.set_xlabel('Distance from a sensor [cells]', color='white')
ax4b.set_ylabel('Correlation / influence', color='white')
ax4b.set_title('Matern-3/2 spreads sensor information locally',
               color='white', fontweight='bold')
ax4b.legend(facecolor='#222222', edgecolor='gray', labelcolor='white')

fig4.suptitle('Deterministic ESMDA / Optimal Interpolation (OI)',
              color='white', fontsize=16, fontweight='bold')
fig4.tight_layout(rect=[0, 0, 1, 0.93])
out4 = os.path.join(OUT, 'fig04_esmda_oi_concept.png')
fig4.savefig(out4, dpi=160, bbox_inches='tight', facecolor='#111111')
plt.close(fig4)
print(f"Saved: {out4}")

# ── fig00: sensor layout ─────────────────────────────────
fig0, ax0 = plt.subplots(figsize=(7, 7.5))
fig0.patch.set_facecolor('#111111'); ax0.set_facecolor('#111111')
im0 = ax0.imshow(x_ss.reshape(NY,NX), vmin=0, vmax=100, cmap=CMAP_T,
                 origin='lower', interpolation='nearest')
for m, col in zip(M_LIST, COLORS):
    obs = all_obs[m]
    ax0.scatter(obs%NX, obs//NX, s=250, c=col, marker='o',
                edgecolors='white', linewidths=1.2, zorder=5,
                label=f'm={m}', alpha=0.85)
ax0.legend(facecolor='#222222', edgecolor='gray', labelcolor='white',
           fontsize=FS_LEG, loc='lower right')
ax0.set_title('Sensor layout on Miffy truth', color='white',
              fontsize=FS_T, fontweight='bold')
ax0.set_xticks([]); ax0.set_yticks([])
for sp in ax0.spines.values(): sp.set_edgecolor('#333333')
cb0 = fig0.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)
cb0.ax.tick_params(colors='white', labelsize=9)
cb0.set_label('Temperature [degC]', color='white', fontsize=9)
out0 = os.path.join(OUT, 'fig00_sensor_layout.png')
fig0.savefig(out0, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig0)
print(f"\nSaved: {out0}")

# ── fig01: RMSE ──────────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(10, 5.5))
fig1.patch.set_facecolor('#111111'); ax1.set_facecolor('#111111')
for m, col in zip(M_LIST, COLORS):
    arr = [all_hist[m][c]['rmse_x'] for c in cyc_arr]
    ax1.plot(cyc_arr, arr, color=col, lw=2.5, label=f'm = {m} sensors')
ax1.set_xlabel('ESMDA Iteration', color='white', fontsize=FS_L)
ax1.set_ylabel('x RMSE [degC]', color='white', fontsize=FS_L)
ax1.set_title(f'RMSE Convergence (ESMDA, N_ITER={N_ITER})',
              color='white', fontsize=FS_T, fontweight='bold')
ax1.tick_params(colors='white', labelsize=FS_TK)
ax1.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=FS_LEG)
ax1.grid(alpha=0.25, color='gray', lw=0.5)
for sp in ax1.spines.values(): sp.set_edgecolor('#444444')
plt.tight_layout()
out1 = os.path.join(OUT, 'fig01_rmse.png')
fig1.savefig(out1, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig1)
print(f"Saved: {out1}")

# ── fig02/03: final fields ───────────────────────────────
for tag, key, cmap, vmin_, vmax_, ttl in [
    ('fig02_xpred_final', 'x', CMAP_T, 0,      100,   'x_pred (temperature)'),
    ('fig03_qest_final',  'q', CMAP_Q, -Q_LIM, Q_LIM,
     f'q_derived = -alpha*L*x  (display: {-Q_LIM:.1f} to {Q_LIM:.1f})'),
]:
    fig_, axes_ = plt.subplots(2, 2, figsize=(11, 11))
    fig_.patch.set_facecolor('#111111')
    im_ref = None
    for ax_, m, col in zip(axes_.flatten(), M_LIST, COLORS):
        ax_.set_facecolor('#111111'); ax_.set_xticks([]); ax_.set_yticks([])
        for sp in ax_.spines.values(): sp.set_edgecolor('#333333')
        field = all_hist[m][N_ITER][key].reshape(NY, NX)
        im_ref = ax_.imshow(field, vmin=vmin_, vmax=vmax_, cmap=cmap,
                            origin='lower', interpolation='nearest')
        obs = all_obs[m]
        sz = max(160.0/np.sqrt(m), 5.0)
        lw = max(1.4-m/1000.0, 0.4)
        ax_.scatter(obs%NX, obs//NX, s=sz, c='#00E676', marker='+',
                    linewidths=lw, zorder=5, alpha=0.65)
        rmse = all_hist[m][N_ITER]['rmse_x']
        detail = f'xRMSE={rmse:.1f} degC'
        if key == 'q':
            detail += f'\nq min={field.min():.1f}, max={field.max():.1f}'
        ax_.set_title(f'm={m}  {detail}', color=col,
                      fontsize=FS_T, fontweight='bold')
    fig_.suptitle(f'Final {ttl} (Iter {N_ITER})',
                  color='white', fontsize=FS_T, fontweight='bold')
    cb_ = fig_.colorbar(im_ref, ax=axes_.ravel().tolist(),
                        fraction=0.035, pad=0.03)
    cb_.ax.tick_params(colors='white', labelsize=FS_TK)
    cb_.set_label('Temperature [degC]' if key == 'x' else 'Heat q [deg/step]',
                  color='white', fontsize=FS_L)
    fig_.subplots_adjust(left=0.03, right=0.90, bottom=0.04, top=0.90,
                         wspace=0.08, hspace=0.18)
    out_ = os.path.join(OUT, f'{tag}.png')
    fig_.savefig(out_, dpi=150, bbox_inches='tight', facecolor='#111111')
    plt.close(fig_)
    print(f"Saved: {out_}")

# ── GIF ──────────────────────────────────────────────────
print(f"\nMaking GIF ({N_ITER+1} frames) ...")
nM    = len(M_LIST)
fig_a = plt.figure(figsize=(nM*3.8+0.9, 9.5))
fig_a.patch.set_facecolor('#111111')
gs = fig_a.add_gridspec(2, nM+1,
                        left=0.02, right=0.93, bottom=0.06, top=0.88,
                        hspace=0.08, wspace=0.06,
                        width_ratios=[1]*nM + [0.06])

ims_x = []; ims_q = []; rmse_txts = []; step_txt = None

for ci, (m, mcol) in enumerate(zip(M_LIST, COLORS)):
    d0  = all_hist[m][0]
    obs = all_obs[m]
    sz  = max(140.0/np.sqrt(m), 4.0)
    lw  = max(1.3-m/1000.0, 0.35)

    ax_x = fig_a.add_subplot(gs[0, ci])
    ax_x.set_facecolor('#111111'); ax_x.set_xticks([]); ax_x.set_yticks([])
    for sp in ax_x.spines.values(): sp.set_edgecolor('#333333')
    im_x = ax_x.imshow(d0['x'].reshape(NY,NX), vmin=0, vmax=100,
                        cmap=CMAP_T, origin='lower', interpolation='nearest')
    ax_x.scatter(obs%NX, obs//NX, s=sz, c='#00E676', marker='+',
                 linewidths=lw, zorder=5, alpha=0.6)
    ax_x.set_title(f'm={m} sensors', color=mcol, fontsize=FS_T,
                   fontweight='bold', pad=4)
    rt = ax_x.text(0.5, 0.04, f'xRMSE={d0["rmse_x"]:.1f}',
                   transform=ax_x.transAxes, ha='center', va='bottom',
                   color='cyan', fontsize=FS_TK, fontweight='bold',
                   bbox=dict(facecolor='#111111', alpha=0.7, edgecolor='none', pad=1))
    if ci == 0:
        step_txt = ax_x.text(0.5, 0.95, f'Iter   0 / {N_ITER}',
                             transform=ax_x.transAxes, ha='center', va='top',
                             color='yellow', fontsize=12, fontweight='bold',
                             bbox=dict(facecolor='#000000', alpha=0.75,
                                       edgecolor='yellow', lw=1.2, pad=3))
    ims_x.append(im_x); rmse_txts.append(rt)

    ax_q = fig_a.add_subplot(gs[1, ci])
    ax_q.set_facecolor('#111111'); ax_q.set_xticks([]); ax_q.set_yticks([])
    for sp in ax_q.spines.values(): sp.set_edgecolor('#333333')
    im_q = ax_q.imshow(d0['q'].reshape(NY,NX), vmin=-Q_LIM, vmax=Q_LIM,
                        cmap=CMAP_Q, origin='lower', interpolation='nearest')
    ax_q.scatter(obs%NX, obs//NX, s=sz, c='#00E676', marker='+',
                 linewidths=lw, zorder=5, alpha=0.6)
    if ci == 0:
        ax_q.set_ylabel(f'q derived\n(+-{Q_LIM:.0f})',
                        color='#FF9999', fontsize=FS_TK, labelpad=4)
    ims_q.append(im_q)

cb_x = fig_a.colorbar(plt.cm.ScalarMappable(cmap=CMAP_T,
       norm=plt.Normalize(0,100)), cax=fig_a.add_subplot(gs[0,-1]))
cb_x.ax.tick_params(colors='white', labelsize=9)
cb_x.set_label('Temp [degC]', color='white', fontsize=9)

cb_q = fig_a.colorbar(plt.cm.ScalarMappable(cmap=CMAP_Q,
       norm=plt.Normalize(-Q_LIM, Q_LIM)), cax=fig_a.add_subplot(gs[1,-1]))
cb_q.ax.tick_params(colors='white', labelsize=9)
cb_q.set_label('Heat [deg/step]', color='white', fontsize=9)

fig_a.text(0.47, 0.96,
              f'Deterministic ESMDA/OI from uniform prior -- Matern L={CORR_L}, N_ITER={N_ITER}',
           ha='center', va='top', color='white', fontsize=FS_T, fontweight='bold')
fig_a.text(0.47, 0.02,
           'top: temperature  |  bottom: q_derived=-alpha*L*x (red=source, blue=sink) | green+=sensor',
           ha='center', va='bottom', color='#aaaaaa', fontsize=FS_L)


def update_anim(frame):
    arts = []
    for ci, m in enumerate(M_LIST):
        d = all_hist[m][frame]
        ims_x[ci].set_data(d['x'].reshape(NY,NX)); arts.append(ims_x[ci])
        ims_q[ci].set_data(d['q'].reshape(NY,NX)); arts.append(ims_q[ci])
        rmse_txts[ci].set_text(f'xRMSE={d["rmse_x"]:.1f}'); arts.append(rmse_txts[ci])
    step_txt.set_text(f'Iter {frame:3d} / {N_ITER}'); arts.append(step_txt)
    return arts


t0  = time.time()
ani = animation.FuncAnimation(fig_a, update_anim, frames=N_ITER+1,
                               interval=150, blit=False)
out_gif = os.path.join(OUT, 'anim_enkf_heatid.gif')
ani.save(out_gif, writer=animation.PillowWriter(fps=6), dpi=110)
plt.close(fig_a)
print(f"Saved: {out_gif}  ({time.time()-t0:.1f}s)")
print("\nDone ->", OUT)
