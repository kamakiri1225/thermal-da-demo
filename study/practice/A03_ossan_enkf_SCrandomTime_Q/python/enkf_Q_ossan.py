"""
A003_ossan_enkf_Q — Physical DA: Estimate Q(x,y) so that T(Q) → T_ossan.

Purpose:
  Use EnKF with -∇²T = Q as the forward model to estimate Q(x,y) from sparse
  T observations.

T_true strategy:
  1. Load the PNG as a binary mask (salmon=100°C, background=0°C) at the
     natural 36×36 cell resolution, then nearest-neighbour upscale to 60×60.
     This preserves the exact pixel-art outline from the image.
  2. Within the ossan silhouette, override anatomical sub-regions with
     distinct temperatures (hat, shirt, collar, eyebrows, eyes, mouth).
     The background/outline shape is never changed — only interior pixels.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
import os, warnings
from PIL import Image
from scipy.ndimage import gaussian_filter
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import factorized

matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
warnings.filterwarnings('ignore')

NX, NY = 60, 60
n = NX * NY

# ─── おっさん温度分布 ─────────────────────────────────────────────────
IMG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        '..', '..', 'A01_ossan', 'img', '001_ossan.png')

T_BG     =  0.0   # background
T_SHIRT  = 50.0   # shirt / body
T_COLLAR = 65.0   # jacket collar
T_SKIN   = 100.0  # face, neck, ears
T_HAT    = 25.0   # cap
T_BROW   = 15.0   # eyebrows
T_EYE    =  5.0   # eyes
T_MOUTH  = 10.0   # smile

def make_ossan():
    """
    Load PNG → exact 36x36 pixel-art mask → scale to 60x60 (nearest-neighbour)
    → flip y → override anatomy regions with distinct temperatures.
    Shape comes entirely from the PNG; temperatures are added inside.
    """
    img_arr = np.array(Image.open(IMG_PATH).convert('RGB'))

    # ── Step 1: sample 36×36 pixel art (row 0 = top of image) ──────
    ART_START, ART_SIZE, GRID_N = 24, 1206, 36
    cell = ART_SIZE / GRID_N
    grid_36 = np.zeros((GRID_N, GRID_N), dtype=np.uint8)
    for row in range(GRID_N):
        for col in range(GRID_N):
            py = int(ART_START + (row + 0.5) * cell)
            px = int(ART_START + (col + 0.5) * cell)
            grid_36[row, col] = 255 if img_arr[py, px, 0] > 150 else 0

    # ── Step 2: scale 36→60 (nearest neighbour) + flip y ──────────
    pil_60 = Image.fromarray(grid_36, mode='L').resize((NX, NY), Image.NEAREST)
    mask = np.flipud(np.array(pil_60)) > 127   # True = ossan pixel, y=0 at bottom

    # ── Step 3: base field: skin everywhere in mask ─────────────────
    field = np.where(mask, T_SKIN, T_BG).astype(float)

    # ── Step 4: anatomy overrides (only inside mask) ─────────────────
    yy = np.arange(NY)[:, None]   # (NY, 1)
    xx = np.arange(NX)[None, :]   # (1, NX)

    def override(region, val):
        field[region & mask] = val

    # Hat (rows 1-6 in 36-grid → y=48-57 in 60-grid)
    override(yy >= 48, T_HAT)

    # Shirt / body (rows 28-35 in 36-grid → y=0-13 in 60-grid)
    override(yy <= 13, T_SHIRT)

    # Collar (narrow region above shirt, centred)
    override((yy >= 14) & (yy <= 21) & (xx >= 20) & (xx <= 40), T_COLLAR)

    # Eyebrows (row 9 in 36-grid → y=43-44)
    override((yy >= 42) & (yy <= 45) & (xx >= 18) & (xx <= 43), T_BROW)

    # Eyes (rows 12-13 → y=37-41)
    override((yy >= 36) & (yy <= 41) & (xx >= 17) & (xx <= 27), T_EYE)   # left
    override((yy >= 36) & (yy <= 41) & (xx >= 33) & (xx <= 43), T_EYE)   # right

    # Mouth / smile (rows 20-22 → y=22-27)
    override((yy >= 22) & (yy <= 27) & (xx >= 14) & (xx <= 46), T_MOUTH)

    return field

T_true_2d = make_ossan()         # (NY, NX)
T_true    = T_true_2d.flatten()  # (n,)

print(f"T_true (Ossan PNG+anatomy): {T_true.min():.1f} ~ {T_true.max():.1f} °C")
_unique = np.unique(T_true)
for t in _unique:
    print(f"    T={t:5.1f}°C: {(T_true==t).sum():5d} cells")

# ─── 2D ラプラシアン行列（ディリクレ境界 T=0） ──────────────────────
def build_laplacian():
    A = lil_matrix((n, n))
    for iy in range(NY):
        for ix in range(NX):
            idx = iy * NX + ix
            if ix == 0 or ix == NX-1 or iy == 0 or iy == NY-1:
                A[idx, idx] = 1.0
            else:
                A[idx, idx]    = -4.0
                A[idx, idx-1]  =  1.0
                A[idx, idx+1]  =  1.0
                A[idx, idx-NX] =  1.0
                A[idx, idx+NX] =  1.0
    return csr_matrix(A)

print("Laplacian: build & LU factorize...")
A_lap     = build_laplacian()
solve_lap = factorized(A_lap)
print("Done")

iy_all = np.arange(n) // NX
ix_all = np.arange(n) % NX
bnd    = (ix_all == 0) | (ix_all == NX-1) | (iy_all == 0) | (iy_all == NY-1)
interior = ~bnd

def forward_T(Q_flat):
    rhs = -Q_flat.copy()
    rhs[bnd] = 0.0
    return solve_lap(rhs)

# 理論 Q_true: おっさん温度をポアソン方程式に代入して逆算
Q_true_th    = -(A_lap @ T_true)
Q_true_th_2d = Q_true_th.reshape(NY, NX)

print(f"Q_true (理論値): {Q_true_th.min():.1f} ~ {Q_true_th.max():.1f}")

rng = np.random.default_rng(42)

# ─── EnKF パラメータ ─────────────────────────────────────────────────
N_ENS   = 50
M_OBS   = 60
SIGMA_R = 3.0
N_CYC   = 50
Q_MEAN  = 0.0
SIGMA_Q = 30.0
CORR_L  = 5.0
R_LOC   = 15.0
INFL    = 1.05

print(f"\nEnKF: N={N_ENS}, m={M_OBS}/cycle, σ_r={SIGMA_R}, R_loc={R_LOC}, α={INFL}")

interior_idx = np.where(interior)[0]

# ─── 初期アンサンブル ────────────────────────────────────────────────
X = np.zeros((n, N_ENS))
for i in range(N_ENS):
    noise = gaussian_filter(rng.normal(0, 1, (NY, NX)), sigma=CORR_L)
    X[:, i] = Q_MEAN + noise.flatten() / (noise.std() + 1e-8) * SIGMA_Q

T_ens = np.zeros((n, N_ENS))
for i in range(N_ENS):
    T_ens[:, i] = forward_T(X[:, i])

def rmse(a, b): return np.sqrt(np.mean((a - b)**2))

rmse_T_hist = [rmse(T_ens.mean(1), T_true)]
T_snaps     = [T_ens.mean(1).reshape(NY, NX).copy()]
Q_snaps     = [X.mean(1).reshape(NY, NX).copy()]
obs_hist    = [None]

print(f"\nCycle  0  RMSE_T={rmse_T_hist[0]:.2f}°C")

R_diag = np.full(M_OBS, SIGMA_R**2)

# ─── EnKF DA サイクル ────────────────────────────────────────────────
for cyc in range(1, N_CYC + 1):
    obs_idx = np.sort(rng.choice(interior_idx, M_OBS, replace=False))
    obs_iy  = obs_idx // NX
    obs_ix  = obs_idx % NX
    obs_hist.append(obs_idx)

    dy  = iy_all[:, None] - obs_iy[None, :]
    dx  = ix_all[:, None] - obs_ix[None, :]
    rho = np.exp(-0.5 * (dy**2 + dx**2) / R_LOC**2)

    y = T_true[obs_idx] + rng.normal(0, SIGMA_R, M_OBS)

    Xp  = X     - X.mean(1, keepdims=True)
    Tp  = T_ens - T_ens.mean(1, keepdims=True)
    HXp = Tp[obs_idx, :]

    S     = HXp @ HXp.T / (N_ENS - 1) + np.diag(R_diag)
    K_raw = Xp @ HXp.T / (N_ENS - 1) @ np.linalg.inv(S)
    K_loc = K_raw * rho

    eps   = rng.normal(0, SIGMA_R, (M_OBS, N_ENS))
    innov = y[:, None] + eps - T_ens[obs_idx, :]
    X_new = X + K_loc @ innov

    xbar_n = X_new.mean(1, keepdims=True)
    X = xbar_n + INFL * (X_new - xbar_n)

    for i in range(N_ENS):
        T_ens[:, i] = forward_T(X[:, i])

    rt = rmse(T_ens.mean(1), T_true)
    rmse_T_hist.append(rt)
    T_snaps.append(T_ens.mean(1).reshape(NY, NX).copy())
    Q_snaps.append(X.mean(1).reshape(NY, NX).copy())

    if cyc in (1, 5, 10, 20, 30, 50):
        print(f"Cycle {cyc:2d}  RMSE_T={rt:.2f}°C")

pct = (1 - rmse_T_hist[-1] / rmse_T_hist[0]) * 100
print(f"\nT 収束: {rmse_T_hist[0]:.2f} → {rmse_T_hist[-1]:.2f}°C  ({pct:.0f}% 改善)")

# ─── 可視化 ──────────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
img_dir    = os.path.join(script_dir, '..', 'img')
os.makedirs(img_dir, exist_ok=True)

CMAP_T = 'RdYlBu_r'
CMAP_Q = 'bwr'
SNAP_CYCLES = [0, 1, 5, 10, 20, 30, 40, 50]

FS_TITLE = 17; FS_LABEL = 15; FS_TICK = 13

def dark_ax(ax):
    ax.set_facecolor('#111111')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_edgecolor('gray')

cyc_arr = np.arange(N_CYC + 1)

# ── fig01: T スナップショット ─────────────────────────────────────────
fig1, axes1 = plt.subplots(2, 5, figsize=(20, 9))
fig1.patch.set_facecolor('#111111')
axes1 = axes1.flatten()
fig1.subplots_adjust(left=0.03, right=0.88, top=0.90, bottom=0.06,
                     hspace=0.38, wspace=0.12)

for ai, cyc in enumerate(SNAP_CYCLES):
    ax = axes1[ai]; dark_ax(ax)
    im = ax.imshow(T_snaps[cyc], origin='lower', cmap=CMAP_T,
                   vmin=0, vmax=100, extent=[0, NX, 0, NY], aspect='equal')
    oi = obs_hist[cyc]
    if oi is not None:
        ax.scatter(oi % NX + 0.5, oi // NX + 0.5, c='lime', s=14,
                   marker='x', linewidths=1.2, alpha=0.8)
    ax.set_title(f'Cycle {cyc}   RMSE={rmse_T_hist[cyc]:.1f}°C',
                 fontsize=FS_TITLE, color='white', fontweight='bold')
last_im = im

ax_tr = axes1[8]; dark_ax(ax_tr)
ax_tr.imshow(T_true_2d, origin='lower', cmap=CMAP_T,
             vmin=0, vmax=100, extent=[0, NX, 0, NY], aspect='equal')
ax_tr.set_title('True T (Ossan)', fontsize=FS_TITLE, color='white', fontweight='bold')

ax_rm = axes1[9]; ax_rm.set_facecolor('#111111')
ax_rm.fill_between(cyc_arr, rmse_T_hist, alpha=0.15, color='white')
ax_rm.plot(cyc_arr, rmse_T_hist, 'w-', lw=2.5)
for c in SNAP_CYCLES: ax_rm.axvline(c, color='cyan', lw=0.8, alpha=0.5)
ax_rm.set_xlabel('DA Cycle', color='white', fontsize=FS_LABEL)
ax_rm.set_ylabel('RMSE_T [°C]', color='white', fontsize=FS_LABEL)
ax_rm.set_title('T Convergence', fontsize=FS_TITLE, color='white', fontweight='bold')
ax_rm.tick_params(colors='white', labelsize=FS_TICK)
ax_rm.grid(alpha=0.3, color='gray', lw=0.6)
for sp in ax_rm.spines.values(): sp.set_edgecolor('gray')

cax_t = fig1.add_axes([0.895, 0.08, 0.018, 0.80])
cbar1 = fig1.colorbar(last_im, cax=cax_t)
cbar1.set_label('Temperature [°C]', color='white', fontsize=FS_LABEL)
cbar1.ax.tick_params(colors='white', labelsize=FS_TICK)
cbar1.outline.set_edgecolor('gray')
fig1.suptitle(
    f'A003_ossan_enkf_Q: T(Q) → T_ossan  (m={M_OBS}/cycle, σ_r={SIGMA_R}°C, R_loc={R_LOC} cells)',
    fontsize=19, color='white', fontweight='bold')
out1 = os.path.join(img_dir, 'fig01_T_snapshots.png')
fig1.savefig(out1, dpi=150, bbox_inches='tight', facecolor=fig1.get_facecolor())
print(f"\nSaved: {out1}")
plt.close(fig1)

# ── fig02: Q スナップショット ─────────────────────────────────────────
# Use one shared Q range across all Q panels so the colors are directly comparable.
Q_shared_max = max(
    np.percentile(np.abs(Q_snaps[-1]), 90) * 1.1,
    1e-6,
)

fig2, axes2 = plt.subplots(2, 5, figsize=(26, 12))
fig2.patch.set_facecolor('#111111')
axes2 = axes2.flatten()
fig2.subplots_adjust(left=0.03, right=0.88, top=0.90, bottom=0.06,
                     hspace=0.38, wspace=0.12)

for ai, cyc in enumerate(SNAP_CYCLES):
    ax = axes2[ai]; dark_ax(ax)
    imq = ax.imshow(Q_snaps[cyc], origin='lower', cmap=CMAP_Q,
                    vmin=-Q_shared_max, vmax=Q_shared_max,
                    extent=[0, NX, 0, NY], aspect='equal')
    ax.set_title(f'Cycle {cyc}   Q analysis', fontsize=FS_TITLE,
                 color='white', fontweight='bold')
last_imq = imq

ax_qt = axes2[8]; dark_ax(ax_qt)
imqt = ax_qt.imshow(Q_true_th_2d, origin='lower', cmap=CMAP_Q,
                    vmin=-Q_shared_max, vmax=Q_shared_max,
                    extent=[0, NX, 0, NY], aspect='equal')
ax_qt.set_title('Q_true = -∇²T_ossan\n(theoretical target)',
                fontsize=FS_TITLE, color='white', fontweight='bold')

ax_rq = axes2[9]; ax_rq.set_facecolor('#111111')
ax_rq.fill_between(cyc_arr, rmse_T_hist, alpha=0.15, color='white')
ax_rq.plot(cyc_arr, rmse_T_hist, 'w-', lw=2.5)
for c in SNAP_CYCLES: ax_rq.axvline(c, color='cyan', lw=0.8, alpha=0.5)
ax_rq.set_xlabel('DA Cycle', color='white', fontsize=FS_LABEL)
ax_rq.set_ylabel('RMSE_T [°C]', color='white', fontsize=FS_LABEL)
ax_rq.set_title('T Convergence', fontsize=FS_TITLE, color='white', fontweight='bold')
ax_rq.tick_params(colors='white', labelsize=FS_TICK)
ax_rq.grid(alpha=0.3, color='gray', lw=0.6)
for sp in ax_rq.spines.values(): sp.set_edgecolor('gray')

cax1 = fig2.add_axes([0.895, 0.52, 0.018, 0.36])
cb1  = fig2.colorbar(last_imq, cax=cax1)
cb1.set_label('Q analysis [red: source / blue: sink]',
              color='white', fontsize=FS_LABEL)
cb1.ax.tick_params(colors='white', labelsize=FS_TICK)
cb1.outline.set_edgecolor('gray')

cax2 = fig2.add_axes([0.895, 0.07, 0.018, 0.36])
cb2  = fig2.colorbar(imqt, cax=cax2)
cb2.set_label('Q_true [red: source / blue: sink]',
              color='white', fontsize=FS_LABEL)
cb2.ax.tick_params(colors='white', labelsize=FS_TICK)
cb2.outline.set_edgecolor('gray')

fig2.suptitle(
    f'A003_ossan_enkf_Q: Estimated Q(x,y)  (m={M_OBS}/cycle, σ_r={SIGMA_R}°C)',
    fontsize=19, color='white', fontweight='bold')
out2 = os.path.join(img_dir, 'fig02_Q_snapshots.png')
fig2.savefig(out2, dpi=150, bbox_inches='tight', facecolor=fig2.get_facecolor())
print(f"Saved: {out2}")
plt.close(fig2)

# ── fig03: 最終比較 ───────────────────────────────────────────────────
fig3, axes3 = plt.subplots(2, 3, figsize=(18, 10))
fig3.patch.set_facecolor('#111111')

items = [
    (T_snaps[0],   'Initial T',              CMAP_T,  0,           100),
    (T_snaps[-1],  'Analysis T  (Cycle 50)', CMAP_T,  0,           100),
    (T_true_2d,    'True T  (Ossan)',         CMAP_T,  0,           100),
    (Q_snaps[0],   'Initial Q',              CMAP_Q, -Q_shared_max, Q_shared_max),
    (Q_snaps[-1],  'Analysis Q  (Cycle 50)', CMAP_Q, -Q_shared_max, Q_shared_max),
    (Q_true_th_2d, 'True Q = -∇²T_ossan',    CMAP_Q, -Q_shared_max, Q_shared_max),
]
for ax, (data, title, cmap, vmin, vmax) in zip(axes3.flatten(), items):
    dark_ax(ax)
    im3 = ax.imshow(data, origin='lower', cmap=cmap,
                    vmin=vmin, vmax=vmax, extent=[0, NX, 0, NY], aspect='equal')
    ax.set_title(title, fontsize=FS_TITLE, color='white', fontweight='bold')
    cb = plt.colorbar(im3, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.tick_params(colors='white', labelsize=FS_TICK)
    cb.outline.set_edgecolor('gray')

plt.suptitle(
    f'A003_ossan_enkf_Q: Initial / Analysis / True  '
    f'(RMSE_T: {rmse_T_hist[0]:.1f} → {rmse_T_hist[-1]:.1f}°C,  {pct:.0f}% improvement)',
    fontsize=19, color='white', fontweight='bold')
plt.tight_layout()
out3 = os.path.join(img_dir, 'fig03_comparison.png')
fig3.savefig(out3, dpi=150, bbox_inches='tight', facecolor=fig3.get_facecolor())
print(f"Saved: {out3}")
plt.close(fig3)

# ── Animation（4 パネル） ─────────────────────────────────────────────
fig_a, axes_a = plt.subplots(2, 2, figsize=(16, 10.5))
fig_a.patch.set_facecolor('#111111')
(ax_Ta, ax_Tt), (ax_Qa, ax_Qt) = axes_a
for ax in axes_a.flatten(): dark_ax(ax)

im_T = ax_Ta.imshow(T_snaps[0], origin='lower', cmap=CMAP_T,
                    vmin=0, vmax=100, extent=[0, NX, 0, NY], aspect='equal')
sc   = ax_Ta.scatter([], [], c='lime', s=20, marker='x',
                     linewidths=1.5, alpha=0.9, label='Sensors')
ax_Ta.legend(facecolor='black', edgecolor='gray', labelcolor='white',
             fontsize=FS_TICK, loc='lower right')
ax_Ta.set_title('T analysis (estimated)', color='white', fontweight='bold', fontsize=FS_TITLE)
cb = fig_a.colorbar(im_T, ax=ax_Ta, fraction=0.046, pad=0.04)
cb.set_label('T [°C]', color='white', fontsize=FS_LABEL); cb.ax.tick_params(colors='white', labelsize=FS_TICK)

ax_Tt.imshow(T_true_2d, origin='lower', cmap=CMAP_T,
             vmin=0, vmax=100, extent=[0, NX, 0, NY], aspect='equal')
ax_Tt.set_title('True T  (target: Ossan)', color='white', fontweight='bold', fontsize=FS_TITLE)
cb2 = fig_a.colorbar(
    plt.cm.ScalarMappable(cmap=CMAP_T, norm=plt.Normalize(0, 100)),
    ax=ax_Tt, fraction=0.046, pad=0.04)
cb2.set_label('T [°C]', color='white', fontsize=FS_LABEL); cb2.ax.tick_params(colors='white', labelsize=FS_TICK)

im_Q = ax_Qa.imshow(Q_snaps[0], origin='lower', cmap=CMAP_Q,
                    vmin=-Q_shared_max, vmax=Q_shared_max,
                    extent=[0, NX, 0, NY], aspect='equal')
ax_Qa.set_title('Q analysis (estimated heat source)', color='white', fontweight='bold', fontsize=FS_TITLE)
cb3 = fig_a.colorbar(im_Q, ax=ax_Qa, fraction=0.046, pad=0.04)
cb3.set_label('Q [red: source / blue: sink]', color='white', fontsize=FS_LABEL)
cb3.ax.tick_params(colors='white', labelsize=FS_TICK)

ax_Qt.imshow(Q_true_th_2d, origin='lower', cmap=CMAP_Q,
             vmin=-Q_shared_max, vmax=Q_shared_max,
             extent=[0, NX, 0, NY], aspect='equal')
ax_Qt.set_title('True Q  (ideal: $-\\nabla^2 T_{\\rm ossan}$)',
                color='white', fontweight='bold', fontsize=FS_TITLE)
cb4 = fig_a.colorbar(
    plt.cm.ScalarMappable(cmap=CMAP_Q,
                          norm=plt.Normalize(-Q_shared_max, Q_shared_max)),
    ax=ax_Qt, fraction=0.046, pad=0.04)
cb4.set_label('Q [red: source / blue: sink]', color='white', fontsize=FS_LABEL)
cb4.ax.tick_params(colors='white', labelsize=FS_TICK)

sup = fig_a.suptitle(
    f'Cycle  0  |  RMSE_T={rmse_T_hist[0]:.1f}°C',
    fontsize=19, color='white', fontweight='bold')
plt.tight_layout()

def update(frame):
    im_T.set_data(T_snaps[frame])
    im_Q.set_data(Q_snaps[frame])
    oi = obs_hist[frame]
    if oi is not None:
        sc.set_offsets(np.c_[oi % NX + 0.5, oi // NX + 0.5])
    sup.set_text(f'Cycle {frame:2d}  |  RMSE_T={rmse_T_hist[frame]:.1f}°C'
                 f'  |  m={M_OBS}/cycle, σ_r={SIGMA_R}°C')
    return [im_T, im_Q, sc, sup]

ani = animation.FuncAnimation(fig_a, update, frames=N_CYC + 1, interval=300, blit=False)
out_gif = os.path.join(img_dir, 'anim_enkf_Q_ossan.gif')
ani.save(out_gif, writer='pillow', fps=4, dpi=100)
print(f"Saved animation: {out_gif}")
plt.close(fig_a)

# ── fig04: Final Q analysis → forward heat solve → T comparison ─────────────
Q_final  = Q_snaps[-1]
T_from_Q = forward_T(Q_final.flatten()).reshape(NY, NX)
err_map  = T_from_Q - T_true_2d

rmse_final  = rmse(T_from_Q.flatten(), T_true)
err_abs_max = np.percentile(np.abs(err_map), 97) * 1.05
Q4_max = np.percentile(np.abs(Q_final), 90) * 1.1

fig4, axes4 = plt.subplots(1, 4, figsize=(24, 7))
fig4.patch.set_facecolor('#111111')
fig4.subplots_adjust(left=0.04, right=0.97, top=0.88, bottom=0.05, wspace=0.35)

panels4 = [
    (Q_final,   f'Q analysis  (Cycle {N_CYC})\n[red=heat source  blue=cooling]',
     CMAP_Q, -Q4_max, Q4_max, 'Q  [a.u.]'),
    (T_from_Q,  'T from Q\n(forward: solve −∇²T = Q)',
     CMAP_T, 0, 100, 'T [degC]'),
    (T_true_2d, 'True T  (target: Ossan)',
     CMAP_T, 0, 100, 'T [degC]'),
    (err_map,   f'Error  T_from_Q − T_true\nRMSE = {rmse_final:.1f} degC',
     'bwr', -err_abs_max, err_abs_max, 'dT [degC]'),
]

for ax, (data, title, cmap, vmin, vmax, clabel) in zip(axes4, panels4):
    dark_ax(ax)
    im4 = ax.imshow(data, origin='lower', cmap=cmap,
                    vmin=vmin, vmax=vmax, extent=[0, NX, 0, NY], aspect='equal')
    ax.set_title(title, fontsize=FS_TITLE, color='white', fontweight='bold', linespacing=1.5)
    cb4 = fig4.colorbar(im4, ax=ax, fraction=0.046, pad=0.04)
    cb4.set_label(clabel, color='white', fontsize=FS_LABEL)
    cb4.ax.tick_params(colors='white', labelsize=FS_TICK)
    cb4.outline.set_edgecolor('gray')

fig4.suptitle(
    f'A003_ossan_enkf_Q:  Final Q analysis  →  forward heat solve  →  T comparison'
    f'  (RMSE_T = {rmse_final:.1f} degC)',
    fontsize=19, color='white', fontweight='bold')
out4 = os.path.join(img_dir, 'fig04_forward_heat_analysis.png')
fig4.savefig(out4, dpi=150, bbox_inches='tight', facecolor=fig4.get_facecolor())
print(f"Saved: {out4}")
plt.close(fig4)

print("\nDone.")
