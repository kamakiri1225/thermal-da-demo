"""
A03 - Physical DA: Estimate heat source Q(x,y) so that T(Q) → T_miffy.

目的:
  ミッフィーの温度分布 T_true を正解として与え、
  -∇²T = Q を前向きモデルとして、スパースな T 観測から
  Q(x,y) をデータ同化（EnKF）で逆推定する。

  → Q_analysis を前向きモデルに入れたとき T_analysis ≈ T_miffy になることを目指す。

A02 との違い:
  A02: 状態変数 = T（直接推定、前向きモデルなし）
  A03: 状態変数 = Q（間接推定、T = solve(-∇²T=Q) という前向きモデルあり）

前向きモデル: 2D 定常熱伝導（ポアソン方程式）
    -∇²T = Q    (Dirichlet T=0 on all boundaries)
    → scipy.sparse.linalg.factorized で LU 分解して高速ソルブ
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
import os, warnings
from scipy.ndimage import gaussian_filter
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import factorized

matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
warnings.filterwarnings('ignore')

NX, NY = 60, 60
n = NX * NY

# ─── ミッフィー温度分布 ───────────────────────────────────────────
# NX=60, NY=60 に合わせて x 座標を 1.5 倍（旧 NX=40 から拡大）
T_BG=0.0; T_OUTLINE=8.0; T_FUR=100.0; T_EYE=5.0
T_MOUTH=5.0; T_DRESS=55.0; T_COLLAR=70.0

def fill_rect(f, y1, y2, x1, x2, val):
    f[max(0,y1):min(NY,y2+1), max(0,x1):min(NX,x2+1)] = val

def fill_ellipse(f, cy, cx, ry, rx, val):
    yy, xx = np.ogrid[:NY, :NX]
    f[((yy-cy)/ry)**2 + ((xx-cx)/rx)**2 <= 1.0] = val

def make_miffy():
    """ミッフィー温度分布 (NY×NX の 2次元配列)  NX=60, NY=60"""
    f = np.full((NY, NX), T_BG)
    fill_ellipse(f,  5, 30,  9, 21, T_DRESS)
    fill_rect(f,  0, 11,  7, 52, T_DRESS)
    fill_rect(f, 10, 12, 24, 36, T_COLLAR)
    fill_ellipse(f, 48, 19, 12,  9, T_OUTLINE)
    fill_ellipse(f, 48, 40, 12,  9, T_OUTLINE)
    fill_ellipse(f, 48, 19, 11,  7, T_FUR)
    fill_ellipse(f, 48, 40, 11,  7, T_FUR)
    fill_rect(f, 35, 40, 12, 27, T_FUR)
    fill_rect(f, 35, 40, 33, 48, T_FUR)
    fill_ellipse(f, 24, 30, 15, 24, T_OUTLINE)
    fill_ellipse(f, 24, 30, 14, 22, T_FUR)
    fill_rect(f,  9, 13, 25, 34, T_FUR)
    fill_ellipse(f, 24, 21,  2,  3, T_EYE)
    fill_ellipse(f, 24, 39,  2,  3, T_EYE)
    for d in range(-2, 3):
        f[min(max(18+d,0),NY-1), min(max(30+d,0),NX-1)] = T_MOUTH
        f[min(max(18-d,0),NY-1), min(max(30+d,0),NX-1)] = T_MOUTH
    return f

# ─── 2D ラプラシアン行列（ディリクレ境界 T=0） ──────────────────
def build_laplacian():
    A = lil_matrix((n, n))
    for iy in range(NY):
        for ix in range(NX):
            idx = iy * NX + ix
            if ix==0 or ix==NX-1 or iy==0 or iy==NY-1:
                A[idx, idx] = 1.0            # 境界セル: T=0 (恒等行)
            else:
                A[idx, idx]   = -4.0
                A[idx, idx-1]  = 1.0
                A[idx, idx+1]  = 1.0
                A[idx, idx-NX] = 1.0
                A[idx, idx+NX] = 1.0
    return csr_matrix(A)

print("Laplacian 構築・LU 分解中...")
A_lap     = build_laplacian()
solve_lap = factorized(A_lap)       # LU 分解を一度だけ（以降はバック代入のみ）
print("完了")

iy_all = np.arange(n) // NX
ix_all = np.arange(n) % NX
bnd    = (ix_all==0) | (ix_all==NX-1) | (iy_all==0) | (iy_all==NY-1)
interior = ~bnd   # 内部セルのマスク

def forward_T(Q_flat):
    """前向きモデル: -∇²T = Q を解いて T を返す（境界 T=0）"""
    rhs = -Q_flat.copy()
    rhs[bnd] = 0.0          # 境界の寄与をゼロにする
    return solve_lap(rhs)

# ─── 正解: ミッフィーの温度分布 ─────────────────────────────────
T_true_2d = make_miffy()            # (NY, NX)  ← A02 と同じ温度場
T_true    = T_true_2d.flatten()     # (n,)      ← EnKF 計算用 1次元ベクトル

# 理論的な Q_true: ミッフィー温度を Poisson 方程式に代入して逆算
# Q_true_th = -(A_lap @ T_true)  → 内部では -∇²T_miffy, 境界では -T_miffy
# これは正負の値を持つ（発熱＋冷却源の分布）
Q_true_th   = -(A_lap @ T_true)
Q_true_th_2d = Q_true_th.reshape(NY, NX)

print(f"\nT_true (Miffy): {T_true.min():.1f} ~ {T_true.max():.1f} °C")
print(f"Q_true (理論値): {Q_true_th.min():.1f} ~ {Q_true_th.max():.1f}")

rng = np.random.default_rng(42)

# ─── EnKF パラメータ ─────────────────────────────────────────────
N_ENS   = 50
M_OBS   = 60      # 1サイクルあたりセンサ数（内部セルのみ）
SIGMA_R = 3.0     # 観測誤差 [°C]
N_CYC   = 50
Q_MEAN  = 0.0     # 初期アンサンブル平均 Q
SIGMA_Q = 30.0    # 初期アンサンブル偏差
CORR_L  = 5.0     # 初期アンサンブルの空間相関長 [cells]
R_LOC   = 15.0    # ローカライゼーション半径 [cells]
INFL    = 1.05    # 乗法インフレーション

print(f"\nEnKF: N={N_ENS}, m={M_OBS}/cycle, σ_r={SIGMA_R}, R_loc={R_LOC}, α={INFL}")

# 内部セルのインデックス（センサは内部のみに配置）
interior_idx = np.where(interior)[0]   # shape: (~2204,)

# ─── 初期アンサンブル（Q） ────────────────────────────────────────
X = np.zeros((n, N_ENS))    # 状態ベクトル = Q
for i in range(N_ENS):
    noise = gaussian_filter(rng.normal(0, 1, (NY, NX)), sigma=CORR_L)
    X[:, i] = Q_MEAN + noise.flatten() / (noise.std()+1e-8) * SIGMA_Q

# 初期アンサンブルに対して前向きモデルを実行
T_ens = np.zeros((n, N_ENS))
for i in range(N_ENS):
    T_ens[:, i] = forward_T(X[:, i])

def rmse(a, b): return np.sqrt(np.mean((a-b)**2))

rmse_T_hist = [rmse(T_ens.mean(1), T_true)]
T_snaps     = [T_ens.mean(1).reshape(NY, NX).copy()]
Q_snaps     = [X.mean(1).reshape(NY, NX).copy()]
obs_hist    = [None]

print(f"\nCycle  0  RMSE_T={rmse_T_hist[0]:.2f}°C")

R_diag = np.full(M_OBS, SIGMA_R**2)

# ─── EnKF DA サイクル ─────────────────────────────────────────────
for cyc in range(1, N_CYC+1):
    # 内部セルからランダムにセンサ選択
    obs_idx = np.sort(rng.choice(interior_idx, M_OBS, replace=False))
    obs_iy  = obs_idx // NX
    obs_ix  = obs_idx % NX
    obs_hist.append(obs_idx)

    # ローカライゼーション行列: Q セル i とセンサ j の距離
    dy  = iy_all[:, None] - obs_iy[None, :]
    dx  = ix_all[:, None] - obs_ix[None, :]
    rho = np.exp(-0.5 * (dy**2 + dx**2) / R_LOC**2)   # (n, m)

    # 観測（T_true のスパースサンプリング）
    y = T_true[obs_idx] + rng.normal(0, SIGMA_R, M_OBS)

    # アンサンブル統計量
    Xp  = X    - X.mean(1, keepdims=True)       # Q 偏差  (n, N)
    Tp  = T_ens - T_ens.mean(1, keepdims=True)  # T 偏差  (n, N)
    HXp = Tp[obs_idx, :]                         # 観測空間の T 偏差  (m, N)

    # イノベーション共分散 (m, m)
    S = HXp @ HXp.T / (N_ENS-1) + np.diag(R_diag)

    # カルマンゲイン: Cov(Q, T_obs) @ S^{-1}  →  (n, m)
    K_raw = Xp @ HXp.T / (N_ENS-1) @ np.linalg.inv(S)
    K_loc = K_raw * rho     # ローカライゼーション（シュール積）

    # 確率的 EnKF 更新
    eps   = rng.normal(0, SIGMA_R, (M_OBS, N_ENS))
    innov = y[:, None] + eps - T_ens[obs_idx, :]   # (m, N)
    X_new = X + K_loc @ innov

    # 乗法インフレーション
    xbar_n = X_new.mean(1, keepdims=True)
    X = xbar_n + INFL * (X_new - xbar_n)

    # 前向きモデルを再実行
    for i in range(N_ENS):
        T_ens[:, i] = forward_T(X[:, i])

    rt = rmse(T_ens.mean(1), T_true)
    rmse_T_hist.append(rt)
    T_snaps.append(T_ens.mean(1).reshape(NY, NX).copy())
    Q_snaps.append(X.mean(1).reshape(NY, NX).copy())

    if cyc in (1, 5, 10, 20, 30, 50):
        print(f"Cycle {cyc:2d}  RMSE_T={rt:.2f}°C")

pct = (1 - rmse_T_hist[-1]/rmse_T_hist[0]) * 100
print(f"\nT 収束: {rmse_T_hist[0]:.2f} → {rmse_T_hist[-1]:.2f}°C  ({pct:.0f}% 改善)")

# ─── 可視化 ──────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
img_dir    = os.path.join(script_dir, '..', 'img')
os.makedirs(img_dir, exist_ok=True)

CMAP_T  = 'RdYlBu_r'
CMAP_Q  = 'bwr'          # 正負の Q を表現（青=冷却源、赤=熱源）
SNAP_CYCLES = [0, 1, 5, 10, 20, 30, 40, 50]

def dark_ax(ax, ticks=False):
    ax.set_facecolor('#111111')
    if not ticks: ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_edgecolor('gray')

# ── fig01: T スナップショット ─────────────────────────────────────
fig1, axes1 = plt.subplots(2, 5, figsize=(20, 9))
fig1.patch.set_facecolor('#111111')
axes1 = axes1.flatten()

FS_TITLE = 17; FS_LABEL = 15; FS_TICK = 13

fig1.subplots_adjust(left=0.03, right=0.88, top=0.90, bottom=0.06,
                     hspace=0.38, wspace=0.12)

for ai, cyc in enumerate(SNAP_CYCLES):
    ax = axes1[ai]; dark_ax(ax)
    im = ax.imshow(T_snaps[cyc], origin='lower', cmap=CMAP_T,
                   vmin=0, vmax=100, extent=[0,NX,0,NY], aspect='equal')
    oi = obs_hist[cyc]
    if oi is not None:
        ax.scatter(oi%NX+0.5, oi//NX+0.5, c='lime', s=14,
                   marker='x', linewidths=1.2, alpha=0.8)
    ax.set_title(f'Cycle {cyc}   RMSE={rmse_T_hist[cyc]:.1f}°C',
                 fontsize=FS_TITLE, color='white', fontweight='bold')
last_im = im

ax_tr = axes1[8]; dark_ax(ax_tr)
ax_tr.imshow(T_true_2d, origin='lower', cmap=CMAP_T,
             vmin=0, vmax=100, extent=[0,NX,0,NY], aspect='equal')
ax_tr.set_title('True T (Miffy)', fontsize=FS_TITLE, color='white', fontweight='bold')

ax_rm = axes1[9]; ax_rm.set_facecolor('#111111')
cyc_arr = np.arange(N_CYC+1)
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
    f'A02_miffy_enkf_Q: T(Q) → T_miffy  (m={M_OBS}/cycle, σ_r={SIGMA_R}°C, R_loc={R_LOC} cells)',
    fontsize=19, color='white', fontweight='bold')
out1 = os.path.join(img_dir, 'fig01_T_snapshots.png')
fig1.savefig(out1, dpi=150, bbox_inches='tight', facecolor=fig1.get_facecolor())
print(f"\nSaved: {out1}")
plt.close(fig1)

# ── fig02: Q スナップショット ─────────────────────────────────────
# Use one shared Q range across all Q panels so the colors are directly comparable.
Q_shared_max = max(
    np.percentile(np.abs(Q_snaps[-1]), 90) * 1.1,
    1e-6,
)

FS_TITLE = 17
FS_LABEL = 15
FS_TICK  = 13

fig2, axes2 = plt.subplots(2, 5, figsize=(26, 12))
fig2.patch.set_facecolor('#111111')
axes2 = axes2.flatten()
fig2.subplots_adjust(left=0.03, right=0.88, top=0.90, bottom=0.06,
                     hspace=0.38, wspace=0.12)

for ai, cyc in enumerate(SNAP_CYCLES):
    ax = axes2[ai]; dark_ax(ax)
    imq = ax.imshow(Q_snaps[cyc], origin='lower', cmap=CMAP_Q,
                    vmin=-Q_shared_max, vmax=Q_shared_max,
                    extent=[0,NX,0,NY], aspect='equal')
    ax.set_title(f'Cycle {cyc}   Q analysis', fontsize=FS_TITLE,
                 color='white', fontweight='bold')
last_imq = imq

# Q_true: use the same range as the analysis panels
ax_qt = axes2[8]; dark_ax(ax_qt)
imqt = ax_qt.imshow(Q_true_th_2d, origin='lower', cmap=CMAP_Q,
                    vmin=-Q_shared_max, vmax=Q_shared_max,
                    extent=[0,NX,0,NY], aspect='equal')
ax_qt.set_title('Q_true = -∇²T_miffy\n(theoretical target)',
                fontsize=FS_TITLE, color='white', fontweight='bold')

# RMSE グラフ（カラーバーと重ならないよう独立）
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

# カラーバー①: Q analysis（右上）
cax1 = fig2.add_axes([0.895, 0.52, 0.018, 0.36])
cb1  = fig2.colorbar(last_imq, cax=cax1)
cb1.set_label('Q analysis [red: source / blue: sink]',
              color='white', fontsize=FS_LABEL)
cb1.ax.tick_params(colors='white', labelsize=FS_TICK)
cb1.outline.set_edgecolor('gray')

# カラーバー②: Q_true（右下）
cax2 = fig2.add_axes([0.895, 0.07, 0.018, 0.36])
cb2  = fig2.colorbar(imqt, cax=cax2)
cb2.set_label('Q_true [red: source / blue: sink]',
              color='white', fontsize=FS_LABEL)
cb2.ax.tick_params(colors='white', labelsize=FS_TICK)
cb2.outline.set_edgecolor('gray')

fig2.suptitle(
    f'A02_miffy_enkf_Q: Estimated Q(x,y) — heat source distribution'
    f'  (m={M_OBS}/cycle, σ_r={SIGMA_R}°C)',
    fontsize=19, color='white', fontweight='bold')

out2 = os.path.join(img_dir, 'fig02_Q_snapshots.png')
fig2.savefig(out2, dpi=150, bbox_inches='tight', facecolor=fig2.get_facecolor())
print(f"Saved: {out2}")
plt.close(fig2)

# ── fig03: 最終比較（Initial / Analysis / True） ─────────────────
fig3, axes3 = plt.subplots(2, 3, figsize=(18, 10))
fig3.patch.set_facecolor('#111111')

items = [
    (T_snaps[0],     'Initial T',              CMAP_T,  0,            100),
    (T_snaps[-1],    'Analysis T  (Cycle 50)', CMAP_T,  0,            100),
    (T_true_2d,      'True T  (Miffy)',         CMAP_T,  0,            100),
    (Q_snaps[0],     'Initial Q',              CMAP_Q, -Q_shared_max, Q_shared_max),
    (Q_snaps[-1],    'Analysis Q  (Cycle 50)', CMAP_Q, -Q_shared_max, Q_shared_max),
    (Q_true_th_2d,   'True Q = -∇²T_miffy',    CMAP_Q, -Q_shared_max, Q_shared_max),
]
for ax, (data, title, cmap, vmin, vmax) in zip(axes3.flatten(), items):
    dark_ax(ax)
    im3 = ax.imshow(data, origin='lower', cmap=cmap,
                    vmin=vmin, vmax=vmax, extent=[0,NX,0,NY], aspect='equal')
    ax.set_title(title, fontsize=FS_TITLE, color='white', fontweight='bold')
    cb = plt.colorbar(im3, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.tick_params(colors='white', labelsize=FS_TICK)
    cb.outline.set_edgecolor('gray')

plt.suptitle(
    f'A02_miffy_enkf_Q: Initial / Analysis / True  '
    f'(RMSE_T: {rmse_T_hist[0]:.1f} → {rmse_T_hist[-1]:.1f}°C,  {pct:.0f}% improvement)',
    fontsize=19, color='white', fontweight='bold')
plt.tight_layout()
out3 = os.path.join(img_dir, 'fig03_comparison.png')
fig3.savefig(out3, dpi=150, bbox_inches='tight', facecolor=fig3.get_facecolor())
print(f"Saved: {out3}")
plt.close(fig3)

# ── Animation（4 パネル）────────────────────────────────────────
#  左上: T analysis（推定温度場）   右上: True T（目標ミッフィー温度）
#  左下: Q analysis（推定発熱量）   右下: True Q（理想発熱量 = -∇²T_miffy）
fig_a, axes_a = plt.subplots(2, 2, figsize=(16, 10.5))
fig_a.patch.set_facecolor('#111111')
(ax_Ta, ax_Tt), (ax_Qa, ax_Qt) = axes_a
for ax in axes_a.flatten(): dark_ax(ax)

im_T = ax_Ta.imshow(T_snaps[0], origin='lower', cmap=CMAP_T,
                    vmin=0, vmax=100, extent=[0,NX,0,NY], aspect='equal')
sc   = ax_Ta.scatter([], [], c='lime', s=20, marker='x',
                     linewidths=1.5, alpha=0.9, label='Sensors')
ax_Ta.legend(facecolor='black', edgecolor='gray', labelcolor='white',
             fontsize=FS_TICK, loc='lower right')
ax_Ta.set_title('T analysis (estimated)', color='white', fontweight='bold', fontsize=FS_TITLE)
cb = fig_a.colorbar(im_T, ax=ax_Ta, fraction=0.046, pad=0.04)
cb.set_label('T [°C]', color='white', fontsize=FS_LABEL); cb.ax.tick_params(colors='white', labelsize=FS_TICK)

ax_Tt.imshow(T_true_2d, origin='lower', cmap=CMAP_T,
             vmin=0, vmax=100, extent=[0,NX,0,NY], aspect='equal')
ax_Tt.set_title('True T  (target: Miffy)', color='white',
                fontweight='bold', fontsize=FS_TITLE)
cb2 = fig_a.colorbar(
    plt.cm.ScalarMappable(cmap=CMAP_T, norm=plt.Normalize(0, 100)),
    ax=ax_Tt, fraction=0.046, pad=0.04)
cb2.set_label('T [°C]', color='white', fontsize=FS_LABEL); cb2.ax.tick_params(colors='white', labelsize=FS_TICK)

im_Q = ax_Qa.imshow(Q_snaps[0], origin='lower', cmap=CMAP_Q,
                    vmin=-Q_shared_max, vmax=Q_shared_max,
                    extent=[0,NX,0,NY], aspect='equal')
ax_Qa.set_title('Q analysis (estimated heat source)', color='white', fontweight='bold', fontsize=FS_TITLE)
cb3 = fig_a.colorbar(im_Q, ax=ax_Qa, fraction=0.046, pad=0.04)
cb3.set_label('Q [red: source / blue: sink]', color='white', fontsize=FS_LABEL)
cb3.ax.tick_params(colors='white', labelsize=FS_TICK)

ax_Qt.imshow(Q_true_th_2d, origin='lower', cmap=CMAP_Q,
             vmin=-Q_shared_max, vmax=Q_shared_max,
             extent=[0,NX,0,NY], aspect='equal')
ax_Qt.set_title('True Q  (ideal: $-\\nabla^2 T_{\\rm miffy}$)',
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
        sc.set_offsets(np.c_[oi%NX+0.5, oi//NX+0.5])
    sup.set_text(f'Cycle {frame:2d}  |  RMSE_T={rmse_T_hist[frame]:.1f}°C'
                 f'  |  m={M_OBS}/cycle, σ_r={SIGMA_R}°C')
    return [im_T, im_Q, sc, sup]

ani = animation.FuncAnimation(fig_a, update, frames=N_CYC+1, interval=300, blit=False)
out_gif = os.path.join(img_dir, 'anim_enkf_Q_miffy.gif')
ani.save(out_gif, writer='pillow', fps=4, dpi=100)
print(f"Saved animation: {out_gif}")
plt.close(fig_a)

# ── fig04: Final Q analysis → forward heat solve → T comparison ─────
Q_final  = Q_snaps[-1]
T_from_Q = forward_T(Q_final.flatten()).reshape(NY, NX)
err_map  = T_from_Q - T_true_2d

rmse_final  = rmse(T_from_Q.flatten(), T_true)
err_abs_max = np.percentile(np.abs(err_map), 97) * 1.05

fig4, axes4 = plt.subplots(1, 4, figsize=(24, 7))
fig4.patch.set_facecolor('#111111')
fig4.subplots_adjust(left=0.04, right=0.97, top=0.88, bottom=0.05,
                     wspace=0.35)

panels4 = [
    (Q_final,   f'Q analysis  (Cycle {N_CYC})\n[red=heat source  blue=cooling]',
     CMAP_Q, -Q_shared_max, Q_shared_max, 'Q  [a.u.]'),
    (T_from_Q,  'T from Q\n(forward: solve −∇²T = Q)',
     CMAP_T, 0, 100, 'T [degC]'),
    (T_true_2d, 'True T  (target: Miffy)',
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
    f'A02_miffy_enkf_Q:  Final Q analysis  →  forward heat solve  →  T comparison'
    f'  (RMSE_T = {rmse_final:.1f} degC)',
    fontsize=19, color='white', fontweight='bold')
out4 = os.path.join(img_dir, 'fig04_forward_heat_analysis.png')
fig4.savefig(out4, dpi=150, bbox_inches='tight', facecolor=fig4.get_facecolor())
print(f"Saved: {out4}")
plt.close(fig4)

print("\nDone.")
