#!/usr/bin/env python3
"""
kf_da_off.py — DA 停止後の挙動デモ

q = 0: モデルは x_ss を知らない
Q = σ_q²·I: モデル誤差
60×60 グリッド

Phase 1 (step 0-100):   KF DA → ミッフィーに近づく
Phase 2 (step 101-300): DA off → 純拡散でミッフィーが消える

出力:
  fig11_da_off.png     — 静止画（主要ステップ比較）
  anim_kf_da_off.gif   — アニメーション（DA on → off の推移）
"""
import numpy as np
import scipy.sparse as sp
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os, time

NX, NY = 60, 60
n = NX * NY
T_BG = 0.; T_FUR = 100.; T_EYE = 5.; T_MOUTH = 5.
T_DRESS = 55.; T_COLLAR = 70.; T_OUTLINE = 8.

def fill_ellipse(f, cy, cx, ry, rx, val):
    yy, xx = np.ogrid[:NY, :NX]
    f[((yy-cy)/ry)**2 + ((xx-cx)/rx)**2 <= 1.0] = val

def fill_rect(f, y1, y2, x1, x2, val):
    f[max(0,y1):min(NY-1,y2)+1, max(0,x1):min(NX-1,x2)+1] = val

def make_miffy():
    f = np.full((NY, NX), T_BG)
    fill_ellipse(f,5,30,9,21,T_DRESS);  fill_rect(f,0,11,7,52,T_DRESS)
    fill_rect(f,10,12,24,36,T_COLLAR)
    fill_ellipse(f,48,19,12,9,T_OUTLINE); fill_ellipse(f,48,40,12,9,T_OUTLINE)
    fill_ellipse(f,24,30,15,24,T_OUTLINE)
    fill_ellipse(f,48,19,11,7,T_FUR);  fill_ellipse(f,48,40,11,7,T_FUR)
    fill_rect(f,35,40,12,27,T_FUR);   fill_rect(f,35,40,33,48,T_FUR)
    fill_ellipse(f,24,30,14,22,T_FUR); fill_rect(f,9,13,25,34,T_FUR)
    fill_ellipse(f,24,21,2,3,T_EYE);  fill_ellipse(f,24,39,2,3,T_EYE)
    for d in range(-2, 3):
        f[min(max(18+d,0),NY-1), min(max(30+d,0),NX-1)] = T_MOUTH
        f[min(max(18-d,0),NY-1), min(max(30+d,0),NX-1)] = T_MOUTH
    return f

def build_laplacian(NY, NX):
    n = NY*NX; row, col, data = [], [], []
    for iy in range(NY):
        for ix in range(NX):
            k = iy*NX+ix; row.append(k); col.append(k); data.append(-4.)
            for diy, dix in [(-1,0),(1,0),(0,-1),(0,1)]:
                jy, jx = iy+diy, ix+dix
                if 0 <= jy < NY and 0 <= jx < NX:
                    row.append(k); col.append(jy*NX+jx); data.append(1.)
    return sp.csr_matrix((data,(row,col)), shape=(n,n), dtype=np.float64)

x_ss  = make_miffy().flatten()
L     = build_laplacian(NY, NX)
ALPHA = 0.5; DT = 0.2
M_sp  = sp.eye(n, format='csr') + DT * ALPHA * L

M_OBS    = 50
SIGMA_R  = 5.0
SIGMA_Q  = 8.0
N_DA     = 100
N_FREE   = 200
N_TOTAL  = N_DA + N_FREE
SIGMA_B0 = 30.0
L_CORR   = 5.0

rng    = np.random.default_rng(42)
iy_all = np.arange(n) // NX
ix_all = np.arange(n) % NX

# ── 初期 P ──────────────────────────────────────────────────────────
print("初期 P 行列を構築中...")
CHUNK = 300
P = np.zeros((n,n), dtype=np.float64)
for i in range(0, n, CHUNK):
    j  = min(i+CHUNK, n)
    dy = (iy_all[i:j,None] - iy_all[None,:]).astype(np.float32)
    dx = (ix_all[i:j,None] - ix_all[None,:]).astype(np.float32)
    P[i:j,:] = (SIGMA_B0**2 * np.exp(-0.5*(dy**2+dx**2)/L_CORR**2)).astype(np.float64)

x_kf = np.zeros(n)
hist = {0: {'x': x_kf.copy(),
            'rmse': np.sqrt(np.mean(x_kf**2 + x_ss**2 - 2*x_kf*x_ss))}}

# ── Phase 1: KF DA ──────────────────────────────────────────────────
print(f"Phase 1: KF DA (step 1-{N_DA}, q=0, σ_q={SIGMA_Q}) ...")
t0 = time.time()
for step in range(1, N_DA + 1):
    x_f = M_sp @ x_kf
    MP  = M_sp.dot(P); P_f = (M_sp.dot(MP.T)).T; P_f = 0.5*(P_f+P_f.T)
    np.fill_diagonal(P_f, P_f.diagonal() + SIGMA_Q**2)

    obs = np.sort(rng.choice(n, M_OBS, replace=False))
    y   = x_ss[obs] + rng.normal(0., SIGMA_R, M_OBS)
    PHt = P_f[:, obs]; S = P_f[np.ix_(obs,obs)] + SIGMA_R**2*np.eye(M_OBS)
    K   = PHt @ np.linalg.inv(S)
    x_kf = x_f + K @ (y - x_f[obs])
    HPf = P_f[obs,:]; P = P_f - K@HPf; P = 0.5*(P+P.T)

    rmse = np.sqrt(np.mean((x_kf - x_ss)**2))
    hist[step] = {'x': x_kf.copy(), 'rmse': rmse}

print(f"  完了: RMSE={hist[N_DA]['rmse']:.2f} degC  ({time.time()-t0:.1f}s)")

# ── Phase 2: DA off（純拡散）────────────────────────────────────────
print(f"Phase 2: DA off 純拡散 (step {N_DA+1}-{N_TOTAL}) ...")
x_free = x_kf.copy()
t0 = time.time()
for step in range(N_DA + 1, N_TOTAL + 1):
    x_free = M_sp @ x_free
    rmse   = np.sqrt(np.mean((x_free - x_ss)**2))
    hist[step] = {'x': x_free.copy(), 'rmse': rmse}

print(f"  完了: RMSE={hist[N_TOTAL]['rmse']:.2f} degC  ({time.time()-t0:.1f}s)")

script_dir = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(script_dir, '..', 'img')
os.makedirs(OUT, exist_ok=True)

# ─────────────────────────────────────────
#  静止画
# ─────────────────────────────────────────
CMAP = 'RdYlBu_r'
SNAP_STEPS = [0, 20, 50, 100, 120, 150, 200, 300]

fig_s, axes_s = plt.subplots(2, len(SNAP_STEPS), figsize=(len(SNAP_STEPS)*2.3+0.5, 6.0))
fig_s.patch.set_facecolor('#111111')

for col, s in enumerate(SNAP_STEPS):
    tc = '#64B5F6' if s <= N_DA else '#FF8A65'
    label = 'DA on' if s <= N_DA else 'DA off'

    ax_est = axes_s[0, col]
    ax_est.set_facecolor('#111111'); ax_est.set_xticks([]); ax_est.set_yticks([])
    for sp_ in ax_est.spines.values():
        sp_.set_edgecolor('yellow' if s == N_DA else '#444444')
        sp_.set_linewidth(2.0 if s == N_DA else 1.0)
    ax_est.imshow(hist[s]['x'].reshape(NY,NX), vmin=0, vmax=T_FUR,
                  cmap=CMAP, origin='lower', interpolation='nearest')
    ax_est.set_title(f'Step {s}\n({label})', color=tc, fontsize=8.5, fontweight='bold')
    ax_est.text(0.5, 0.03, f'RMSE={hist[s]["rmse"]:.1f}',
                transform=ax_est.transAxes, ha='center', va='bottom',
                color='cyan', fontsize=8,
                bbox=dict(facecolor='#111111', alpha=0.6, edgecolor='none', pad=1))

    ax_truth = axes_s[1, col]
    ax_truth.set_facecolor('#111111'); ax_truth.set_xticks([]); ax_truth.set_yticks([])
    for sp_ in ax_truth.spines.values(): sp_.set_edgecolor('#444444')
    ax_truth.imshow(x_ss.reshape(NY,NX), vmin=0, vmax=T_FUR,
                    cmap=CMAP, origin='lower', interpolation='nearest')

axes_s[0, 0].set_ylabel('KF estimate', color='white', fontsize=9)
axes_s[1, 0].set_ylabel('Truth (static)', color='white', fontsize=9)

fig_s.suptitle(
    f'KF DA on/off  (q=0, σ_q={SIGMA_Q}, σ_r={SIGMA_R} degC, m={M_OBS})\n'
    f'Step 0-{N_DA}: DA on  →  Step {N_DA+1}-{N_TOTAL}: DA off (pure diffusion)',
    color='white', fontsize=10, fontweight='bold')
plt.tight_layout(rect=[0, 0, 1, 0.91])
out_s = os.path.join(OUT, 'fig11_da_off.png')
fig_s.savefig(out_s, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig_s)
print(f"保存: {out_s}")

# ─────────────────────────────────────────
#  GIF アニメーション（5ステップ刻り）
# ─────────────────────────────────────────
print(f"\nGIF 作成中 ...")
ANIM_STEPS = list(range(0, N_TOTAL + 1, 5))

fig_a = plt.figure(figsize=(9.5, 5.2))
fig_a.patch.set_facecolor('#111111')
gs = fig_a.add_gridspec(1, 2, left=0.05, right=0.88, bottom=0.10, top=0.82, wspace=0.06)
ax_est   = fig_a.add_subplot(gs[0, 0])
ax_truth = fig_a.add_subplot(gs[0, 1])

for ax, title, tc in [
    (ax_est,   'KF estimate', '#4FC3F7'),
    (ax_truth, 'Truth (Miffy)', 'white'),
]:
    ax.set_facecolor('#111111'); ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values(): sp_.set_edgecolor('#555555')
    ax.set_title(title, color=tc, fontsize=11, fontweight='bold', pad=5)

im_est   = ax_est.imshow(hist[0]['x'].reshape(NY,NX), vmin=0, vmax=T_FUR,
                         cmap=CMAP, origin='lower', interpolation='nearest')
im_truth = ax_truth.imshow(x_ss.reshape(NY,NX), vmin=0, vmax=T_FUR,
                           cmap=CMAP, origin='lower', interpolation='nearest')

step_txt = ax_est.text(
    0.5, 0.96, f'Step   0 / {N_TOTAL}',
    transform=ax_est.transAxes, ha='center', va='top',
    color='yellow', fontsize=13, fontweight='bold',
    bbox=dict(facecolor='#000000', alpha=0.75, edgecolor='yellow', linewidth=1.5, pad=3))

phase_txt = ax_truth.text(
    0.5, 0.96, '  DA on  ',
    transform=ax_truth.transAxes, ha='center', va='top',
    color='#64B5F6', fontsize=11, fontweight='bold',
    bbox=dict(facecolor='#000000', alpha=0.75, edgecolor='#64B5F6', linewidth=1.5, pad=3))

rmse_txt = ax_est.text(
    0.5, 0.04, f'RMSE={hist[0]["rmse"]:.2f} degC',
    transform=ax_est.transAxes, ha='center', va='bottom',
    color='cyan', fontsize=10, fontweight='bold',
    bbox=dict(facecolor='#111111', alpha=0.7, edgecolor='none', pad=1))

cbar_ax = fig_a.add_axes([0.90, 0.15, 0.022, 0.60])
sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(0, T_FUR))
cbar = fig_a.colorbar(sm, cax=cbar_ax)
cbar.ax.tick_params(colors='white', labelsize=9)
cbar.set_label('Temp [degC]', color='white', fontsize=9, labelpad=6)

fig_a.text(0.46, 0.92,
           f'KF DA on/off  (q=0, σ_q={SIGMA_Q}, σ_r={SIGMA_R} degC, m={M_OBS})',
           ha='center', va='bottom', color='white', fontsize=10.5, fontweight='bold')
fig_a.text(0.46, 0.02,
           f'Step 0-{N_DA}: DA on  |  Step {N_DA+1}-{N_TOTAL}: DA off — pure diffusion  '
           f'(yellow border)',
           ha='center', va='bottom', color='#aaaaaa', fontsize=8.5)

def update_anim(frame_idx):
    step = ANIM_STEPS[frame_idx]
    d = hist[step]
    im_est.set_data(d['x'].reshape(NY, NX))
    step_txt.set_text(f'Step {step:3d} / {N_TOTAL}')
    rmse_txt.set_text(f'RMSE={d["rmse"]:.2f} degC')

    if step <= N_DA:
        phase_txt.set_text('  DA on  ')
        phase_txt.set_color('#64B5F6')
        phase_txt.get_bbox_patch().set_edgecolor('#64B5F6')
        for sp_ in ax_est.spines.values():
            sp_.set_edgecolor('#555555'); sp_.set_linewidth(1.0)
    else:
        phase_txt.set_text('  DA off  ')
        phase_txt.set_color('#FF8A65')
        phase_txt.get_bbox_patch().set_edgecolor('#FF8A65')
        for sp_ in ax_est.spines.values():
            sp_.set_edgecolor('yellow'); sp_.set_linewidth(2.0)

    return [im_est, im_truth, step_txt, phase_txt, rmse_txt]

t0 = time.time()
ani = animation.FuncAnimation(fig_a, update_anim, frames=len(ANIM_STEPS),
                               interval=120, blit=False)
out_gif = os.path.join(OUT, 'anim_kf_da_off.gif')
ani.save(out_gif, writer=animation.PillowWriter(fps=8), dpi=110)
plt.close(fig_a)
print(f"保存: {out_gif}  ({time.time()-t0:.1f}s)")
print("\n完了 →", OUT)
