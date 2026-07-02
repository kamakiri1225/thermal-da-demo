#!/usr/bin/env python3
"""
esmda_qfull_ossan.py -- Ossan heat identification via Q-state ESMDA

Same full-inverse formulation as A02_miffy_esmda_SCfixTimefix_Q/python/
esmda_qfull.py, applied to the 60x60 ossan field (PNG pixel art + anatomy
temperature overrides, shared with A03_ossan_enkf_SCrandomTime_Q).

State  : q (heat-source field, n=3600 cells)
Forward: x(q) = L^-1 (-q/alpha), observed at fixed sensors
Prior  : zero-mean stationary Gaussian random fields (no ossan shape)
Sensors: geometry-only space-filling layout
Output : top row = x(q_mean), bottom row = q_mean vs q_true
Grid   : 60x60 (n=3600), N=300, N_ITER=30, R_LOC=20 (domain is 2x miffy's,
         so the nonlocal q-T sensitivity reaches farther in cell units)
"""
import os
# Cap BLAS threads before importing numpy: the matrices here are small
# (<= 3600x1200) and uncapped OpenBLAS threading is slower than serial.
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '4')
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.ndimage import gaussian_filter
from PIL import Image
import csv
import time

script_dir = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(script_dir, '..', 'img')
os.makedirs(OUT, exist_ok=True)

OUTPUT_FILES = [
    'fig00_sensor_layout.png',
    'fig05_qfull_rmse.png',
    'fig06_qfull_qest_final.png',
    'fig07_qfull_xrecon_final.png',
    'fig08_qfull_vs_derived.png',
    'fig09_qfull_sensor_count_sweep.png',
    'anim_esmda_qfull.gif',
]

for filename in OUTPUT_FILES:
    path = os.path.join(OUT, filename)
    if os.path.isfile(path):
        os.remove(path)
print(f"Initialized output directory: {os.path.abspath(OUT)}")

# ── Grid & physics ──────────────────────────────────────
NX, NY = 60, 60
n      = NX * NY
ALPHA  = 0.5

IMG_PATH = os.path.join(script_dir, '..', '..', 'A01_ossan', 'img', '001_ossan.png')
T_BG, T_SHIRT, T_COLLAR, T_SKIN = 0.0, 50.0, 65.0, 100.0
T_HAT, T_BROW, T_EYE, T_MOUTH = 25.0, 15.0, 5.0, 10.0

def make_ossan():
    """PNG pixel art -> 36x36 mask -> 60x60 field with anatomy overrides.

    Identical to A03_ossan_enkf_SCrandomTime_Q so results are comparable.
    """
    img_arr = np.array(Image.open(IMG_PATH).convert('RGB'))
    ART_START, ART_SIZE, GRID_N = 24, 1206, 36
    cell = ART_SIZE / GRID_N
    grid_36 = np.zeros((GRID_N, GRID_N), dtype=np.uint8)
    for row in range(GRID_N):
        for col in range(GRID_N):
            py = int(ART_START + (row + 0.5) * cell)
            px = int(ART_START + (col + 0.5) * cell)
            grid_36[row, col] = 255 if img_arr[py, px, 0] > 150 else 0
    pil_60 = Image.fromarray(grid_36, mode='L').resize((NX, NY), Image.NEAREST)
    mask = np.flipud(np.array(pil_60)) > 127
    field = np.where(mask, T_SKIN, T_BG).astype(float)
    yy = np.arange(NY)[:, None]
    xx = np.arange(NX)[None, :]
    def override(region, val):
        field[region & mask] = val
    override(yy >= 48, T_HAT)
    override(yy <= 13, T_SHIRT)
    override((yy >= 14) & (yy <= 21) & (xx >= 20) & (xx <= 40), T_COLLAR)
    override((yy >= 42) & (yy <= 45) & (xx >= 18) & (xx <= 43), T_BROW)
    override((yy >= 36) & (yy <= 41) & (xx >= 17) & (xx <= 27), T_EYE)
    override((yy >= 36) & (yy <= 41) & (xx >= 33) & (xx <= 43), T_EYE)
    override((yy >= 22) & (yy <= 27) & (xx >= 14) & (xx <= 46), T_MOUTH)
    return field

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

x_ss   = make_ossan().flatten()
L      = build_laplacian(NY, NX)
q_true = -ALPHA * (L @ x_ss)
LU     = spl.splu(L.tocsc())

def forward_x(Q):
    """Forward model: q -> steady temperature. Works on (n,) or (n, N)."""
    return LU.solve(-Q / ALPHA)

assert np.abs(forward_x(q_true) - x_ss).max() < 1e-9

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

M_LIST = [100, 300, 600, 1200]
SENSOR_SWEEP = [50, 100, 200, 300, 400, 600, 800, 1200]
SWEEP_SEEDS = [42, 43, 44]
COLORS = ['#FF6B6B', '#FFD93D', '#6BCB77', '#4D96FF']

# ── Stochastic ESMDA parameters ──────────────────────────
N_ENS = 300
SIGMA_QB = float(np.round(q_true.std(), 1))  # prior spread matches q_true scale
SIGMA_R = 3.0      # observation noise [degC]
CORR_L = 1.5       # Gaussian-filter correlation scale [cells]
# The 60x60 domain doubles the nonlocal reach of L^-1 in cell units compared
# with the 30x30 miffy case (R_LOC=10 there): 10 still shows a transient
# blow-up here, 20 is clean with the same final accuracy.
R_LOC = 20.0       # covariance-localization radius [cells]
R_LOC_TSTATE = 5.0 # baseline keeps the direct-observation value
N_ITER = 30
SIGMA_R_EFF = SIGMA_R * np.sqrt(N_ITER)

print(f"\nQ-state ESMDA params: N={N_ENS}, N_ITER={N_ITER}, "
      f"SIGMA_QB={SIGMA_QB}, SIGMA_R={SIGMA_R}, "
      f"CORR_L={CORR_L}, R_LOC={R_LOC}")

def initial_ensemble(rng, spread):
    X = np.empty((n, N_ENS))
    for i in range(N_ENS):
        noise = rng.normal(0.0, 1.0, (NY, NX))
        X[:, i] = gaussian_filter(noise, sigma=CORR_L).ravel()
    X -= X.mean(axis=1, keepdims=True)
    X *= spread / (X.std() + 1e-10)
    return X

def localization(obs_idx, radius):
    obs_iy = obs_idx // NX
    obs_ix = obs_idx % NX
    dy = iy_all[:, None] - obs_iy[None, :]
    dx = ix_all[:, None] - obs_ix[None, :]
    return np.exp(-0.5*(dy**2 + dx**2)/radius**2)

def run_esmda_qfull(obs_idx, seed=42, keep_history=False):
    """ESMDA with state = q; sensors observe x(q) = L^-1(-q/alpha)."""
    m_obs = len(obs_idx)
    rho = localization(obs_idx, R_LOC)

    rng = np.random.default_rng(seed)
    Q = initial_ensemble(rng, SIGMA_QB)

    def snapshot(it):
        qmean = Q.mean(axis=1)
        xmean = forward_x(qmean)
        return {
            'q':      qmean if keep_history else None,
            'x':      xmean if keep_history else None,
            'rmse_q': np.sqrt(np.mean((qmean - q_true)**2)),
            'rmse_x': np.sqrt(np.mean((xmean - x_ss)**2)),
        }

    history = {0: snapshot(0)}

    for it in range(1, N_ITER+1):
        X_temp = forward_x(Q)
        Qp = Q - Q.mean(axis=1, keepdims=True)
        Tp = X_temp - X_temp.mean(axis=1, keepdims=True)
        HTp = Tp[obs_idx, :]
        S = (HTp @ HTp.T)/(N_ENS-1) + SIGMA_R_EFF**2*np.eye(m_obs)
        PHt = (Qp @ HTp.T)/(N_ENS-1)
        K = np.linalg.solve(S.T, PHt.T).T
        eps = rng.normal(0.0, SIGMA_R_EFF, (m_obs, N_ENS))
        innovation = x_ss[obs_idx, None] + eps - X_temp[obs_idx, :]
        Q += (K*rho) @ innovation
        history[it] = snapshot(it)

    return history

def run_esmda_tstate_qderived(obs_idx, seed=42):
    """Baseline: T-state ESMDA; returns final derived-q field and RMSE."""
    m_obs = len(obs_idx)
    rho = localization(obs_idx, R_LOC_TSTATE)
    rng = np.random.default_rng(seed)
    X = initial_ensemble(rng, 40.0) + 20.0   # SIGMA_B=40, PRIOR_MEAN=20

    for it in range(N_ITER):
        Xp = X - X.mean(axis=1, keepdims=True)
        HXp = Xp[obs_idx, :]
        S = (HXp @ HXp.T)/(N_ENS-1) + SIGMA_R_EFF**2*np.eye(m_obs)
        PHt = (Xp @ HXp.T)/(N_ENS-1)
        K = np.linalg.solve(S.T, PHt.T).T
        eps = rng.normal(0.0, SIGMA_R_EFF, (m_obs, N_ENS))
        innovation = x_ss[obs_idx, None] + eps - X[obs_idx, :]
        X += (K*rho) @ innovation

    xmean = X.mean(axis=1)
    q_derived = -ALPHA * (L @ xmean)
    return q_derived, np.sqrt(np.mean((q_derived - q_true)**2))

# ── Run ─────────────────────────────────────────────────
print(f"\nRunning Q-state ESMDA ...")
all_hist = {}
all_obs  = {}
sweep_rmse_q = {}
sweep_rmse_x = {}
for m in SENSOR_SWEEP:
    obs_idx = sensor_subsets(m)
    print(f"  m={m:4d}", end=' ', flush=True)
    start = time.perf_counter()
    final_q = []
    final_x = []
    for seed in SWEEP_SEEDS:
        keep = (seed == SWEEP_SEEDS[0] and m in M_LIST)
        hist = run_esmda_qfull(obs_idx, seed=seed, keep_history=keep)
        final_q.append(hist[N_ITER]['rmse_q'])
        final_x.append(hist[N_ITER]['rmse_x'])
        if keep:
            all_hist[m] = hist
            all_obs[m] = obs_idx
    sweep_rmse_q[m] = np.asarray(final_q)
    sweep_rmse_x[m] = np.asarray(final_x)
    elapsed = time.perf_counter() - start
    print(f"-> {elapsed:.1f}s  qRMSE={np.mean(final_q):.2f} "
          f"+- {np.std(final_q, ddof=1):.2f}  xRMSE={np.mean(final_x):.2f} degC")

print(f"\nRunning T-state baseline (derived q) for comparison ...")
baseline_qderived = {}
baseline_rmse_q = {}
for m in M_LIST:
    obs_idx = sensor_subsets(m)
    q_derived, rq = run_esmda_tstate_qderived(obs_idx, seed=SWEEP_SEEDS[0])
    baseline_qderived[m] = q_derived
    baseline_rmse_q[m] = rq
    print(f"  m={m:4d}: derived-q RMSE={rq:.2f}  "
          f"(Q-state: {all_hist[m][N_ITER]['rmse_q']:.2f})")

# Point temperatures cannot constrain the cell-scale edges of q_true, so
# also score both estimates against a smoothed truth: the recoverable scale.
Q_SMOOTH_SIGMA = 1.0
q_true_sm = gaussian_filter(q_true.reshape(NY, NX), sigma=Q_SMOOTH_SIGMA).ravel()

def corr(a, b):
    return float(np.corrcoef(a, b)[0, 1])

print(f"\nCorrelation with smoothed q_true (sigma={Q_SMOOTH_SIGMA}):")
for m in M_LIST:
    c_qfull = corr(all_hist[m][N_ITER]['q'], q_true_sm)
    c_deriv = corr(baseline_qderived[m], q_true_sm)
    print(f"  m={m:4d}: Q-state={c_qfull:.3f}  derived={c_deriv:.3f}")

sweep_m = np.asarray(SENSOR_SWEEP)
sweep_q_mean = np.asarray([sweep_rmse_q[m].mean() for m in SENSOR_SWEEP])
sweep_q_std = np.asarray([sweep_rmse_q[m].std(ddof=1) for m in SENSOR_SWEEP])
sweep_x_mean = np.asarray([sweep_rmse_x[m].mean() for m in SENSOR_SWEEP])
sweep_x_std = np.asarray([sweep_rmse_x[m].std(ddof=1) for m in SENSOR_SWEEP])

docs_dir = os.path.abspath(os.path.join(script_dir, '..', 'docs'))
os.makedirs(docs_dir, exist_ok=True)
csv_path = os.path.join(docs_dir, 'sensor_count_sweep_qfull.csv')
with open(csv_path, 'w', newline='', encoding='utf-8') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['sensor_count',
                     'q_rmse_mean', 'q_rmse_std',
                     'x_rmse_mean_degC', 'x_rmse_std_degC',
                     *[f'q_rmse_seed_{seed}' for seed in SWEEP_SEEDS]])
    for i, m in enumerate(sweep_m):
        writer.writerow([m,
                         f'{sweep_q_mean[i]:.6f}', f'{sweep_q_std[i]:.6f}',
                         f'{sweep_x_mean[i]:.6f}', f'{sweep_x_std[i]:.6f}',
                         *[f'{v:.6f}' for v in sweep_rmse_q[int(m)]]])
print(f"Saved: {csv_path}")

CMAP_T = 'RdYlBu_r'
CMAP_Q = 'seismic'
FS_T   = 15; FS_L = 13; FS_TK = 11; FS_LEG = 12
cyc_arr = np.arange(N_ITER+1)

# Common q color range from the estimates so their structure stays visible
# (q_true edge amplitudes are several times larger and would wash them out).
final_q = np.concatenate([all_hist[m][N_ITER]['q'] for m in M_LIST])
Q_LIM = max(10.0, float(np.percentile(np.abs(final_q), 99.0)))
print(f"Heat visualization range: {-Q_LIM:.1f} to {Q_LIM:.1f} "
      f"(99th percentile of estimates; q_true is clipped on this scale)")

def dark_ax(ax):
    ax.set_facecolor('#111111')
    ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values():
        sp_.set_edgecolor('#333333')

# ── fig00: sensor layout ─────────────────────────────────
fig0, ax0 = plt.subplots(figsize=(7, 7.5))
fig0.patch.set_facecolor('#111111'); ax0.set_facecolor('#111111')
im0 = ax0.imshow(x_ss.reshape(NY,NX), vmin=0, vmax=100, cmap=CMAP_T,
                 origin='lower', interpolation='nearest')
for m, col in zip(M_LIST, COLORS):
    obs = all_obs[m]
    ax0.scatter(obs%NX, obs//NX, s=max(220.0/np.sqrt(m/30.0), 8.0), c=col,
                marker='o', edgecolors='white', linewidths=0.8, zorder=5,
                label=f'm={m}', alpha=0.85)
ax0.legend(facecolor='#222222', edgecolor='gray', labelcolor='white',
           fontsize=FS_LEG, loc='lower right')
ax0.set_title('Sensor layout on Ossan truth', color='white',
              fontsize=FS_T, fontweight='bold')
ax0.set_xticks([]); ax0.set_yticks([])
for sp_ in ax0.spines.values(): sp_.set_edgecolor('#333333')
cb0 = fig0.colorbar(im0, ax=ax0, fraction=0.046, pad=0.04)
cb0.ax.tick_params(colors='white', labelsize=9)
cb0.set_label('Temperature [degC]', color='white', fontsize=9)
out0 = os.path.join(OUT, 'fig00_sensor_layout.png')
fig0.savefig(out0, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig0)
print(f"\nSaved: {out0}")

# ── fig05: RMSE convergence ──────────────────────────────
fig5, (ax5a, ax5b) = plt.subplots(1, 2, figsize=(14, 5.5))
fig5.patch.set_facecolor('#111111')
for ax_, key, ylabel, title in [
    (ax5a, 'rmse_q', 'q RMSE', 'Heat-source RMSE convergence'),
    (ax5b, 'rmse_x', 'x RMSE [degC]', 'Temperature RMSE of x(q_mean)'),
]:
    ax_.set_facecolor('#111111')
    for m, col in zip(M_LIST, COLORS):
        arr = [all_hist[m][c][key] for c in cyc_arr]
        ax_.plot(cyc_arr, arr, color=col, lw=2.5, label=f'm = {m} sensors')
    if key == 'rmse_q':
        for m, col in zip(M_LIST, COLORS):
            ax_.axhline(baseline_rmse_q[m], color=col, ls=':', lw=1.3, alpha=0.7)
        ax_.plot([], [], color='gray', ls=':', lw=1.3,
                 label='derived q = -aLx (T-state)')
    ax_.set_xlabel('ESMDA Iteration', color='white', fontsize=FS_L)
    ax_.set_ylabel(ylabel, color='white', fontsize=FS_L)
    ax_.set_title(title, color='white', fontsize=FS_T, fontweight='bold')
    ax_.tick_params(colors='white', labelsize=FS_TK)
    ax_.legend(facecolor='#222222', edgecolor='gray', labelcolor='white',
               fontsize=10)
    ax_.grid(alpha=0.25, color='gray', lw=0.5)
    for sp_ in ax_.spines.values(): sp_.set_edgecolor('#444444')
fig5.suptitle(f'Q-state ESMDA on Ossan (full inverse problem, N_ITER={N_ITER})',
              color='white', fontsize=FS_T, fontweight='bold')
plt.tight_layout(rect=(0, 0, 1, 0.94))
out5 = os.path.join(OUT, 'fig05_qfull_rmse.png')
fig5.savefig(out5, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig5)
print(f"Saved: {out5}")

# ── fig06/07: final fields ───────────────────────────────
for tag, key, cmap, vmin_, vmax_, ttl, rmse_key, unit in [
    ('fig06_qfull_qest_final', 'q', CMAP_Q, -Q_LIM, Q_LIM,
     'q_est (state variable)', 'rmse_q', ''),
    ('fig07_qfull_xrecon_final', 'x', CMAP_T, 0, 100,
     'x(q_est) via forward solve', 'rmse_x', ' degC'),
]:
    fig_, axes_ = plt.subplots(2, 2, figsize=(11, 11))
    fig_.patch.set_facecolor('#111111')
    im_ref = None
    for ax_, m, col in zip(axes_.flatten(), M_LIST, COLORS):
        dark_ax(ax_)
        field = all_hist[m][N_ITER][key].reshape(NY, NX)
        im_ref = ax_.imshow(field, vmin=vmin_, vmax=vmax_, cmap=cmap,
                            origin='lower', interpolation='nearest')
        obs = all_obs[m]
        sz = max(160.0/np.sqrt(m), 3.0)
        lw = max(1.4-m/1500.0, 0.3)
        ax_.scatter(obs%NX, obs//NX, s=sz, c='#00E676', marker='+',
                    linewidths=lw, zorder=5, alpha=0.65)
        rmse = all_hist[m][N_ITER][rmse_key]
        ax_.set_title(f'm={m}  RMSE={rmse:.1f}{unit}', color=col,
                      fontsize=FS_T, fontweight='bold')
    fig_.suptitle(f'Final {ttl} (Iter {N_ITER})',
                  color='white', fontsize=FS_T, fontweight='bold')
    cb_ = fig_.colorbar(im_ref, ax=axes_.ravel().tolist(),
                        fraction=0.035, pad=0.03)
    cb_.ax.tick_params(colors='white', labelsize=FS_TK)
    cb_.set_label('Heat q [deg/step]' if key == 'q' else 'Temperature [degC]',
                  color='white', fontsize=FS_L)
    fig_.subplots_adjust(left=0.03, right=0.90, bottom=0.04, top=0.90,
                         wspace=0.08, hspace=0.18)
    out_ = os.path.join(OUT, f'{tag}.png')
    fig_.savefig(out_, dpi=150, bbox_inches='tight', facecolor='#111111')
    plt.close(fig_)
    print(f"Saved: {out_}")

# ── fig08: truth vs Q-state vs derived baseline ──────────
# Per-panel robust scale: the fields differ in amplitude by design (edges vs
# smoothed estimate), so structure is what the panels compare.
m_show = M_LIST[-1]
q_show = all_hist[m_show][N_ITER]['q']
fig8, axes8 = plt.subplots(1, 4, figsize=(21, 6))
fig8.patch.set_facecolor('#111111')
panels8 = [
    (q_true, 'q_true = -alpha*L*x_ss\n(cell-scale edge field)'),
    (q_true_sm, f'q_true smoothed (sigma={Q_SMOOTH_SIGMA})\n= recoverable scale'),
    (q_show,
     f'Q-state ESMDA  (m={m_show})\ncorr vs smoothed = '
     f'{corr(q_show, q_true_sm):.2f}'),
    (baseline_qderived[m_show],
     f'derived q = -alpha*L*x  (m={m_show})\ncorr vs smoothed = '
     f'{corr(baseline_qderived[m_show], q_true_sm):.2f}'),
]
for ax_, (field, title) in zip(axes8, panels8):
    dark_ax(ax_)
    lim = max(1.0, float(np.percentile(np.abs(field), 99.0)))
    im8 = ax_.imshow(field.reshape(NY, NX), vmin=-lim, vmax=lim,
                     cmap=CMAP_Q, origin='lower', interpolation='nearest')
    ax_.set_title(title, color='white', fontsize=FS_T, fontweight='bold')
    cb8 = fig8.colorbar(im8, ax=ax_, fraction=0.046, pad=0.03)
    cb8.ax.tick_params(colors='white', labelsize=9)
fig8.suptitle('Ossan heat-source identification: only the smoothed scale is '
              'recoverable (each panel on its own robust scale)',
              color='white', fontsize=FS_T+1, fontweight='bold')
fig8.tight_layout(rect=(0, 0, 1, 0.93))
out8 = os.path.join(OUT, 'fig08_qfull_vs_derived.png')
fig8.savefig(out8, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig8)
print(f"Saved: {out8}")

# ── fig09: sensor-count sweep ────────────────────────────
fig9, (ax9a, ax9b) = plt.subplots(1, 2, figsize=(14, 5.5))
fig9.patch.set_facecolor('#111111')
for ax_ in (ax9a, ax9b):
    ax_.set_facecolor('#111111')
    ax_.tick_params(colors='white', labelsize=FS_TK)
    ax_.grid(alpha=0.25, color='gray', lw=0.5)
    for sp_ in ax_.spines.values():
        sp_.set_edgecolor('#444444')

ax9a.errorbar(sweep_m, sweep_q_mean, yerr=sweep_q_std, color='#EF476F',
              marker='o', ms=7, lw=2.5, capsize=4,
              label=f'q RMSE mean +- std ({len(SWEEP_SEEDS)} seeds)')
ax9a.axhline(q_true.std(), color='#FFD166', ls='--', lw=1.5,
             label=f'q_true std = {q_true.std():.1f} (zero-info level)')
ax9a.set_xlabel('Number of fixed sensors', color='white', fontsize=FS_L)
ax9a.set_ylabel('Final q RMSE', color='white', fontsize=FS_L)
ax9a.set_title('Heat-source accuracy vs sensor count', color='white',
               fontsize=FS_T, fontweight='bold')
ax9a.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=9)

ax9b.errorbar(sweep_m, sweep_x_mean, yerr=sweep_x_std, color='#4D96FF',
              marker='o', ms=7, lw=2.5, capsize=4,
              label=f'x RMSE mean +- std ({len(SWEEP_SEEDS)} seeds)')
ax9b.set_xlabel('Number of fixed sensors', color='white', fontsize=FS_L)
ax9b.set_ylabel('Final temperature RMSE [degC]', color='white', fontsize=FS_L)
ax9b.set_title('Temperature accuracy of x(q_mean)', color='white',
               fontsize=FS_T, fontweight='bold')
ax9b.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=9)

fig9.suptitle('Q-state ESMDA on Ossan: sensor-count study',
              color='white', fontsize=FS_T, fontweight='bold')
fig9.tight_layout(rect=[0, 0, 1, 0.92])
out9 = os.path.join(OUT, 'fig09_qfull_sensor_count_sweep.png')
fig9.savefig(out9, dpi=160, bbox_inches='tight', facecolor='#111111')
plt.close(fig9)
print(f"Saved: {out9}")

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
    sz  = max(140.0/np.sqrt(m), 2.5)
    lw  = max(1.3-m/1500.0, 0.3)

    ax_x = fig_a.add_subplot(gs[0, ci])
    dark_ax(ax_x)
    im_x = ax_x.imshow(d0['x'].reshape(NY,NX), vmin=0, vmax=100,
                        cmap=CMAP_T, origin='lower', interpolation='nearest')
    ax_x.scatter(obs%NX, obs//NX, s=sz, c='#00E676', marker='+',
                 linewidths=lw, zorder=5, alpha=0.6)
    ax_x.set_title(f'm={m} sensors', color=mcol, fontsize=FS_T,
                   fontweight='bold', pad=4)
    rt = ax_x.text(0.5, 0.04, f'qRMSE={d0["rmse_q"]:.1f}',
                   transform=ax_x.transAxes, ha='center', va='bottom',
                   color='cyan', fontsize=FS_TK, fontweight='bold',
                   bbox=dict(facecolor='#111111', alpha=0.7, edgecolor='none', pad=1))
    if ci == 0:
        ax_x.set_ylabel('x(q_mean)\nforward solve', color='#99CCFF',
                        fontsize=FS_TK, labelpad=4)
        step_txt = ax_x.text(0.5, 0.95, f'Iter   0 / {N_ITER}',
                             transform=ax_x.transAxes, ha='center', va='top',
                             color='yellow', fontsize=12, fontweight='bold',
                             bbox=dict(facecolor='#000000', alpha=0.75,
                                       edgecolor='yellow', lw=1.2, pad=3))
    ims_x.append(im_x); rmse_txts.append(rt)

    ax_q = fig_a.add_subplot(gs[1, ci])
    dark_ax(ax_q)
    im_q = ax_q.imshow(d0['q'].reshape(NY,NX), vmin=-Q_LIM, vmax=Q_LIM,
                        cmap=CMAP_Q, origin='lower', interpolation='nearest')
    ax_q.scatter(obs%NX, obs//NX, s=sz, c='#00E676', marker='+',
                 linewidths=lw, zorder=5, alpha=0.6)
    if ci == 0:
        ax_q.set_ylabel(f'q_mean (state)\n(+-{Q_LIM:.0f})',
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
           f'Q-state ESMDA on Ossan (full inverse) -- N={N_ENS}, CORR_L={CORR_L}, R_LOC={R_LOC}',
           ha='center', va='top', color='white', fontsize=FS_T, fontweight='bold')
fig_a.text(0.47, 0.02,
           'top: x(q_mean) forward solve  |  bottom: q_mean state '
           '(red=source, blue=sink) | green+=sensor',
           ha='center', va='bottom', color='#aaaaaa', fontsize=FS_L)


def update_anim(frame):
    arts = []
    for ci, m in enumerate(M_LIST):
        d = all_hist[m][frame]
        ims_x[ci].set_data(d['x'].reshape(NY,NX)); arts.append(ims_x[ci])
        ims_q[ci].set_data(d['q'].reshape(NY,NX)); arts.append(ims_q[ci])
        rmse_txts[ci].set_text(f'qRMSE={d["rmse_q"]:.1f}'); arts.append(rmse_txts[ci])
    step_txt.set_text(f'Iter {frame:3d} / {N_ITER}'); arts.append(step_txt)
    return arts


t0  = time.time()
ani = animation.FuncAnimation(fig_a, update_anim, frames=N_ITER+1,
                               interval=150, blit=False)
out_gif = os.path.join(OUT, 'anim_esmda_qfull.gif')
ani.save(out_gif, writer=animation.PillowWriter(fps=6), dpi=110)
plt.close(fig_a)
print(f"Saved: {out_gif}  ({time.time()-t0:.1f}s)")
print("\nDone ->", OUT)
