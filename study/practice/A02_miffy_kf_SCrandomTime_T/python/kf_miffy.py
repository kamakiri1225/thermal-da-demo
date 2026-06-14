#!/usr/bin/env python3
"""
kf_miffy.py — 熱伝導方程式を予測ステップに持つカルマンフィルタ (KF) データ同化デモ

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
問題設定（正直な設定）:
  真の場     : ミッフィー定常温度場 x_ss（静的, DA では使わない）
  予測モデル : 純粋な熱拡散  x^f = M x^a   (q = 0)
               モデルは x_ss を一切知らない
  モデル誤差 : Q = σ_q²·I  (真の場は拡散しないが, モデルは拡散を予測するため)
  観測       : 毎ステップ m 点をランダムに選び, 真値 + ノイズ を取得

  → DA なし  : x → 0 (熱拡散で一様になる)
  → DA あり  : スパース観測だけで x → x_ss (ミッフィー) に収束

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
KF アルゴリズム (1ステップ):
  予測:
    x^f = M x^a                              (熱拡散で 1 ステップ先を予測)
    P^f = M P^a M^T + σ_q²·I                (モデル誤差 Q で P の縮退を防ぐ)
  解析:
    K = P^f H^T (H P^f H^T + R)^{-1}
    x^a = x^f + K (y - H x^f)
    P^a = P^f - K H P^f

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
主要パラメータ:
  n       = 3600  (60×60 格子)
  m       = 30    センサー数/ステップ
  ALPHA   = 0.5   熱拡散率 [cells²/step]
  DT      = 0.2   時間刻み  (dt·α = 0.1 < 0.125 → 安定)
  SIGMA_Q = 8.0   モデル誤差標準偏差 [°C] (Q = σ_q²·I)
  N_STEPS = 100
"""

import numpy as np
import scipy.sparse as sp
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os, time

# ─────────────────────────────────────────
#  グリッド & ミッフィー定義
# ─────────────────────────────────────────
NX, NY = 60, 60
n = NX * NY

T_BG = 0.0; T_FUR = 100.0; T_EYE = 5.0; T_MOUTH = 5.0
T_DRESS = 55.0; T_COLLAR = 70.0; T_OUTLINE = 8.0

def fill_ellipse(f, cy, cx, ry, rx, val):
    yy, xx = np.ogrid[:NY, :NX]
    f[((yy - cy) / ry)**2 + ((xx - cx) / rx)**2 <= 1.0] = val

def fill_rect(f, y1, y2, x1, x2, val):
    f[max(0, y1):min(NY-1, y2)+1, max(0, x1):min(NX-1, x2)+1] = val

def make_miffy():
    f = np.full((NY, NX), T_BG)
    fill_ellipse(f,  5, 30,  9, 21, T_DRESS)
    fill_rect   (f,  0, 11,  7, 52, T_DRESS)
    fill_rect   (f, 10, 12, 24, 36, T_COLLAR)
    fill_ellipse(f, 48, 19, 12,  9, T_OUTLINE)
    fill_ellipse(f, 48, 40, 12,  9, T_OUTLINE)
    fill_ellipse(f, 24, 30, 15, 24, T_OUTLINE)
    fill_ellipse(f, 48, 19, 11,  7, T_FUR)
    fill_ellipse(f, 48, 40, 11,  7, T_FUR)
    fill_rect   (f, 35, 40, 12, 27, T_FUR)
    fill_rect   (f, 35, 40, 33, 48, T_FUR)
    fill_ellipse(f, 24, 30, 14, 22, T_FUR)
    fill_rect   (f,  9, 13, 25, 34, T_FUR)
    fill_ellipse(f, 24, 21,  2,  3, T_EYE)
    fill_ellipse(f, 24, 39,  2,  3, T_EYE)
    for d in range(-2, 3):
        f[min(max(18+d, 0), NY-1), min(max(30+d, 0), NX-1)] = T_MOUTH
        f[min(max(18-d, 0), NY-1), min(max(30+d, 0), NX-1)] = T_MOUTH
    return f

# ─────────────────────────────────────────
#  2D ラプラシアン行列（5点差分）
# ─────────────────────────────────────────
def build_laplacian(NY, NX):
    n = NY * NX
    row, col, data = [], [], []
    for iy in range(NY):
        for ix in range(NX):
            k = iy * NX + ix
            row.append(k); col.append(k); data.append(-4.0)
            for diy, dix in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                jy, jx = iy + diy, ix + dix
                if 0 <= jy < NY and 0 <= jx < NX:
                    row.append(k)
                    col.append(jy * NX + jx)
                    data.append(1.0)
    return sp.csr_matrix((data, (row, col)), shape=(n, n), dtype=np.float64)

# ─────────────────────────────────────────
#  物理モデル設定
# ─────────────────────────────────────────
x_ss = make_miffy().flatten()    # ミッフィー定常解（観測生成と評価にのみ使用）
L    = build_laplacian(NY, NX)

ALPHA = 0.5   # 熱拡散率 [cells²/step]
DT    = 0.2   # 時間刻み  dt·α = 0.1 ≤ 0.125 → 安定

# q = 0: モデルは熱源を持たない（x_ss の情報をモデルに与えない）
M_sp = sp.eye(n, format='csr') + DT * ALPHA * L   # 状態遷移行列 M（疎）

# ─────────────────────────────────────────
#  KF パラメータ
# ─────────────────────────────────────────
M_OBS   = 50      # センサー数/ステップ
SIGMA_R = 5.0     # 観測誤差 [°C]
SIGMA_Q = 8.0     # モデル誤差 [°C]  (Q = σ_q²·I)
N_STEPS = 100     # DA ステップ数
SIGMA_B0= 30.0    # 初期不確かさ σ_b [°C]
L_CORR  = 5.0     # 初期ガウス相関長 [cells]

rng    = np.random.default_rng(42)
iy_all = np.arange(n) // NX
ix_all = np.arange(n) % NX

# ─────────────────────────────────────────
#  初期誤差共分散 P（ガウス型, チャンク構築）
# ─────────────────────────────────────────
print("\n初期 P 行列を構築中（ガウス共分散）...")
t0 = time.time()
CHUNK = 300
P = np.zeros((n, n), dtype=np.float64)
for i in range(0, n, CHUNK):
    j  = min(i + CHUNK, n)
    dy = (iy_all[i:j, None] - iy_all[None, :]).astype(np.float32)
    dx = (ix_all[i:j, None] - ix_all[None, :]).astype(np.float32)
    P[i:j, :] = (SIGMA_B0**2 * np.exp(-0.5 * (dy**2 + dx**2) / L_CORR**2)).astype(np.float64)
print(f"  完了: {time.time()-t0:.1f}s, {P.nbytes/1e6:.0f} MB")

# ─────────────────────────────────────────
#  初期状態・ベースライン
# ─────────────────────────────────────────
x_kf   = np.zeros(n)   # KF 推定値（コールドスタート）
x_phys = np.zeros(n)   # 物理モデルのみ（センサーなし, q=0 で拡散）

rmse_kf   = [np.sqrt(np.mean((x_kf   - x_ss)**2))]
rmse_phys = [np.sqrt(np.mean((x_phys - x_ss)**2))]
p_std     = [np.sqrt(P.diagonal().mean())]

# アニメーション用: 5ステップ刻み (0, 5, 10, ..., 100) → 21フレーム
SNAP_STEPS = set(range(0, N_STEPS + 1, 5))
snaps = {0: (x_kf.copy(), x_phys.copy(), P.diagonal().copy())}

print(f"\nStep   0: KF RMSE = {rmse_kf[0]:.2f}°C,  Physics-only RMSE = {rmse_phys[0]:.2f}°C")
print(f"KF ループ開始 (N_STEPS={N_STEPS}, m={M_OBS}, σ_r={SIGMA_R}°C, σ_q={SIGMA_Q}°C) ...")
t_start = time.time()

# ─────────────────────────────────────────
#  KF DA ループ
# ─────────────────────────────────────────
for step in range(1, N_STEPS + 1):

    # ── 予測ステップ（熱拡散, q=0）────────────────────────────────
    x_f = M_sp @ x_kf                           # x^f = M x^a  (q=0)

    MP  = M_sp.dot(P)
    P_f = (M_sp.dot(MP.T)).T
    P_f = 0.5 * (P_f + P_f.T)
    # モデル誤差 Q = σ_q²·I を加算: P^f ← M P M^T + σ_q²·I
    np.fill_diagonal(P_f, P_f.diagonal() + SIGMA_Q**2)

    # ── センサー観測 ──────────────────────────────────────────────
    obs_idx = np.sort(rng.choice(n, M_OBS, replace=False))
    y = x_ss[obs_idx] + rng.normal(0.0, SIGMA_R, M_OBS)

    # ── 解析ステップ ──────────────────────────────────────────────
    PHt  = P_f[:, obs_idx]
    HPHt = P_f[np.ix_(obs_idx, obs_idx)]
    S    = HPHt + SIGMA_R**2 * np.eye(M_OBS)

    K     = PHt @ np.linalg.inv(S)
    innov = y - x_f[obs_idx]
    x_kf  = x_f + K @ innov

    HPf = P_f[obs_idx, :]
    P   = P_f - K @ HPf
    P   = 0.5 * (P + P.T)

    # ── 物理ベースライン（q=0, DA なし → 拡散して消える）────────
    x_phys = M_sp @ x_phys

    # ── 記録 ──────────────────────────────────────────────────────
    rmse_kf.append(  np.sqrt(np.mean((x_kf   - x_ss)**2)))
    rmse_phys.append(np.sqrt(np.mean((x_phys - x_ss)**2)))
    p_std.append(    np.sqrt(P.diagonal().mean()))

    if step in SNAP_STEPS:
        snaps[step] = (x_kf.copy(), x_phys.copy(), P.diagonal().copy())

    if step in (1, 5, 10, 20, 50, 100):
        elapsed = time.time() - t_start
        print(f"Step {step:3d}: KF={rmse_kf[-1]:.2f}°C,  Phys={rmse_phys[-1]:.2f}°C,  "
              f"√P̄_ii={p_std[-1]:.2f}°C  ({elapsed:.1f}s)")

print(f"\n合計時間: {time.time()-t_start:.1f}s")
print(f"収束: KF {rmse_kf[0]:.1f}→{rmse_kf[-1]:.1f}°C  "
      f"(物理のみ {rmse_phys[0]:.1f}→{rmse_phys[-1]:.1f}°C)")

# ─────────────────────────────────────────
#  RMSE を npy に保存
# ─────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(script_dir, '..', 'img')
os.makedirs(OUT, exist_ok=True)

np.save(os.path.join(OUT, 'kf_rmse.npy'),   np.array(rmse_kf))
np.save(os.path.join(OUT, 'phys_rmse.npy'), np.array(rmse_phys))

# ─────────────────────────────────────────
#  可視化
# ─────────────────────────────────────────
CMAP  = 'RdYlBu_r'
VMAX  = T_FUR

def dark_ax(ax):
    ax.set_facecolor('#111111')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_edgecolor('gray')

snap_steps = sorted(snaps.keys())
n_snap = len(snap_steps)

# ── 図1: スナップショット比較 ─────────────────────────────────────
fig1, axes1 = plt.subplots(3, n_snap, figsize=(n_snap * 2.5, 8))
fig1.patch.set_facecolor('#111111')
for col, s in enumerate(snap_steps):
    xk, xp, pd = snaps[s]
    ax_kf  = axes1[0, col]
    ax_ph  = axes1[1, col]
    ax_std = axes1[2, col]

    dark_ax(ax_kf)
    ax_kf.imshow(xk.reshape(NY, NX), vmin=0, vmax=VMAX, cmap=CMAP,
                 origin='lower', interpolation='nearest')
    ax_kf.set_title(f'Step {s}\nRMSE={rmse_kf[s]:.1f}', color='white', fontsize=9, fontweight='bold')

    dark_ax(ax_ph)
    ax_ph.imshow(xp.reshape(NY, NX), vmin=0, vmax=VMAX, cmap=CMAP,
                 origin='lower', interpolation='nearest')
    ax_ph.set_title(f'Phys={rmse_phys[s]:.1f}', color='#aaaaaa', fontsize=8)

    dark_ax(ax_std)
    ax_std.imshow(pd.reshape(NY, NX), vmin=0, vmax=SIGMA_B0**2, cmap='Blues',
                  origin='lower', interpolation='nearest')
    ax_std.set_title(f'sqrtP={p_std[s]:.1f}', color='#aaaaff', fontsize=8)

axes1[0, 0].set_ylabel('KF estimate', color='white', fontsize=9)
axes1[1, 0].set_ylabel('Physics only\n(q=0, diffuses)', color='#aaaaaa', fontsize=9)
axes1[2, 0].set_ylabel('Uncertainty\n√P_ii', color='#aaaaff', fontsize=9)

fig1.suptitle(
    f'Kalman Filter + Heat Diffusion  (m={M_OBS} sensors/step, σ_r={SIGMA_R}, σ_q={SIGMA_Q} degC)\n'
    f'q=0: model has NO knowledge of x_ss — DA recovers Miffy from sparse observations only',
    color='white', fontsize=10, fontweight='bold')
plt.tight_layout()
out1 = os.path.join(OUT, 'fig01_kf_snapshots.png')
fig1.savefig(out1, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig1)
print(f"保存: {out1}")

# ── 図2: RMSE 収束曲線 ────────────────────────────────────────────
steps = np.arange(N_STEPS + 1)

fig2, ax2 = plt.subplots(figsize=(10, 5))
fig2.patch.set_facecolor('#111111')
ax2.set_facecolor('#111111')
ax2.semilogy(steps, rmse_kf,   color='#4FC3F7', lw=2.5,
             label=f'KF  (heat diffusion + sensors, m={M_OBS}/step)')
ax2.semilogy(steps, rmse_phys, color='#FF8A65', lw=2, ls='--',
             label='Physics only (q=0, no sensors) — diffuses to 0')
ax2.axhline(SIGMA_R, color='gray', ls=':', lw=1.5, label=f'Obs. error σ_r={SIGMA_R} degC')

ax2.set_xlabel('Step', color='white', fontsize=12)
ax2.set_ylabel('RMSE [degC]', color='white', fontsize=12)
ax2.set_title('KF (heat diffusion forecast, q=0) vs Physics-only  —  RMSE convergence',
              color='white', fontsize=13, fontweight='bold')
ax2.tick_params(colors='white', labelsize=10)
ax2.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=9)
ax2.grid(alpha=0.3, color='gray', lw=0.5)
for sp in ax2.spines.values(): sp.set_edgecolor('gray')

plt.tight_layout()
out2 = os.path.join(OUT, 'fig02_kf_rmse.png')
fig2.savefig(out2, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig2)
print(f"保存: {out2}")

# ── 図3: P 対角の空間分布（不確かさマップ） ──────────────────────
snap_for_p = [s for s in snap_steps if s > 0]
n_p = len(snap_for_p)
fig3, axes3 = plt.subplots(1, n_p, figsize=(n_p * 2.5 + 1.2, 3.2))
fig3.patch.set_facecolor('#111111')
fig3.subplots_adjust(left=0.02, right=0.85, top=0.84, bottom=0.02, wspace=0.08)
for col, s in enumerate(snap_for_p):
    ax = axes3[col]
    ax.set_facecolor('#111111')
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_edgecolor('gray')
    _, _, pd = snaps[s]
    im = ax.imshow(pd.reshape(NY, NX), vmin=0, vmax=SIGMA_B0**2,
                   cmap='hot_r', origin='lower', interpolation='nearest')
    ax.set_title(f'Step {s}\nsqrtP={p_std[s]:.1f}', color='white', fontsize=8)

cbar_ax3 = fig3.add_axes([0.87, 0.12, 0.025, 0.7])
cbar3 = fig3.colorbar(im, cax=cbar_ax3)
cbar3.set_label('Variance P_ii [degC^2]', color='white', fontsize=8)
cbar3.ax.tick_params(colors='white', labelsize=7)
cbar3.ax.yaxis.label.set_color('white')
cbar3.outline.set_edgecolor('gray')

fig3.suptitle('Uncertainty map: P diagonal — kept non-zero by model error Q = σ_q²·I',
              color='white', fontsize=9, y=0.98)
out3 = os.path.join(OUT, 'fig03_kf_pdiag.png')
fig3.savefig(out3, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig3)
print(f"保存: {out3}")

# ── 図4: 最終比較（正解 / KF 推定 / 残差） ──────────────────────
fig4, axes4 = plt.subplots(1, 3, figsize=(13, 5))
fig4.patch.set_facecolor('#111111')
titles4 = ['Truth (Miffy)',
           f'KF estimate (Step {N_STEPS})\nRMSE={rmse_kf[-1]:.2f} degC',
           'Residual  KF - Truth']
cmaps4  = [CMAP, CMAP, 'RdBu']
vmins4  = [0, 0, -30]
vmaxs4  = [VMAX, VMAX, 30]
data4   = [x_ss, x_kf, x_kf - x_ss]
for ax, d, t, cm, vmin, vmax in zip(axes4, data4, titles4, cmaps4, vmins4, vmaxs4):
    dark_ax(ax)
    ax.imshow(d.reshape(NY, NX), vmin=vmin, vmax=vmax, cmap=cm,
              origin='lower', interpolation='nearest')
    ax.set_title(t, color='white', fontsize=10, fontweight='bold')

plt.tight_layout()
out4 = os.path.join(OUT, 'fig04_kf_final.png')
fig4.savefig(out4, dpi=150, bbox_inches='tight', facecolor='#111111')
plt.close(fig4)
print(f"保存: {out4}")

# ── GIF アニメーション（FuncAnimation, 5ステップ刻み）────────────
print("GIF アニメーション作成中 ...")

ANIM_STEPS = sorted(SNAP_STEPS)   # [0, 5, 10, ..., 100]

fig_g = plt.figure(figsize=(11, 5.5))
fig_g.patch.set_facecolor('#111111')
gs = fig_g.add_gridspec(1, 2, left=0.04, right=0.97, bottom=0.06, top=0.82, wspace=0.06)
ax_kf = fig_g.add_subplot(gs[0, 0])
ax_pd = fig_g.add_subplot(gs[0, 1])
dark_ax(ax_kf); dark_ax(ax_pd)

ax_kf.set_title('KF estimate  (q=0: model knows nothing)', color='#4FC3F7', fontsize=11)
ax_pd.set_title('Uncertainty  √P_ii', color='#aaaaff', fontsize=11)

fig_g.text(0.5, 0.95,
           f'KF Data Assimilation  (m={M_OBS} sensors/step, σ_r={SIGMA_R}, σ_q={SIGMA_Q} degC)',
           ha='center', va='top', color='white', fontsize=11, fontweight='bold')

x0, _, pd0 = snaps[0]
im_kf = ax_kf.imshow(x0.reshape(NY, NX), vmin=0, vmax=VMAX,
                     cmap=CMAP, origin='lower', interpolation='nearest')
im_pd = ax_pd.imshow(pd0.reshape(NY, NX), vmin=0, vmax=SIGMA_B0**2,
                     cmap='hot_r', origin='lower', interpolation='nearest')

# テキストはサブプロット内（タイトルと重ならない位置）
step_txt = ax_kf.text(0.5, 0.96, 'Step   0 / 100',
                      transform=ax_kf.transAxes, ha='center', va='top',
                      color='yellow', fontsize=12, fontweight='bold',
                      bbox=dict(facecolor='#000000', alpha=0.75,
                                edgecolor='yellow', linewidth=1.5, pad=3))
rmse_txt = ax_kf.text(0.5, 0.04, f'RMSE={rmse_kf[0]:.2f} degC',
                      transform=ax_kf.transAxes, ha='center', va='bottom',
                      color='cyan', fontsize=10, fontweight='bold',
                      bbox=dict(facecolor='#111111', alpha=0.7, edgecolor='none', pad=1))

def update(frame_idx):
    s = ANIM_STEPS[frame_idx]
    xk, _, pd = snaps[s]
    im_kf.set_data(xk.reshape(NY, NX))
    im_pd.set_data(pd.reshape(NY, NX))
    step_txt.set_text(f'Step {s:3d} / {N_STEPS}')
    rmse_txt.set_text(f'RMSE={rmse_kf[s]:.2f} degC')
    return [im_kf, im_pd, step_txt, rmse_txt]

ani = animation.FuncAnimation(fig_g, update, frames=len(ANIM_STEPS),
                               interval=500, blit=False)
out_gif = os.path.join(OUT, 'anim_kf_miffy.gif')
ani.save(out_gif, writer=animation.PillowWriter(fps=2), dpi=120)
plt.close(fig_g)
print(f"保存: {out_gif}")

print("\n全出力完了 →", OUT)
