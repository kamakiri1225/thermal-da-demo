"""
1次元OI法：背景場 + 観測 → 解析場
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'

# ===== 背景誤差共分散行列（ガウス型）=====
def make_B(x_grid, sigma_b, L):
    d = x_grid[:, None] - x_grid[None, :]          # (n, n) 距離行列
    return sigma_b**2 * np.exp(-0.5 * (d / L)**2)

# ===== OI法 =====
def optimal_interpolation(x_grid, xb, obs_locs, obs_vals, sigma_b, L, sigma_r):
    n = len(x_grid)
    m = len(obs_locs)

    B = make_B(x_grid, sigma_b, L)                  # (n, n)

    # 観測演算子 H：最近傍グリッド点から取得
    H = np.zeros((m, n))
    for i, loc in enumerate(obs_locs):
        H[i, np.argmin(np.abs(x_grid - loc))] = 1.0

    R   = sigma_r**2 * np.eye(m)
    BHT = B @ H.T                                    # (n, m)
    S   = H @ BHT + R                                # (m, m) = HBH^T + R
    K   = BHT @ np.linalg.inv(S)                    # (n, m) カルマンゲイン

    innovation = obs_vals - H @ xb                   # (m,)
    xa         = xb + K @ innovation                 # (n,)
    Pa         = (np.eye(n) - K @ H) @ B             # (n, n)
    sigma_a    = np.sqrt(np.diag(Pa))                # (n,)

    return xa, sigma_a

# ===== 設定 =====
np.random.seed(42)
x_grid = np.linspace(0, 100, 200)

x_true = 20 + 0.05 * x_grid + 2.0 * np.sin(x_grid / 15)   # 真値
xb     = 20 + 0.05 * x_grid                                  # 背景場

obs_locs = np.array([15.0, 35.0, 55.0, 75.0, 90.0])
obs_vals = np.interp(obs_locs, x_grid, x_true) + np.random.normal(0, 0.8, len(obs_locs))

sigma_b = 2.0
L       = 20.0
sigma_r = 0.8

xa, sigma_a = optimal_interpolation(x_grid, xb, obs_locs, obs_vals, sigma_b, L, sigma_r)

# ===== プロット1：解析場 =====
fig, axes = plt.subplots(2, 1, figsize=(13, 8.5))

ax = axes[0]
ax.plot(x_grid, x_true, 'k--', lw=1.5,       label='True state')
ax.plot(x_grid, xb,     'b-',  lw=1.5, alpha=0.6, label='Background $x^b$')
ax.plot(x_grid, xa,     'r-',  lw=2,          label='Analysis $x^a$ (OI)')
ax.fill_between(x_grid, xa - 2*sigma_a, xa + 2*sigma_a,
                color='red', alpha=0.15, label='±2σ (analysis error)')
ax.scatter(obs_locs, obs_vals, s=80, color='green', zorder=5, label='Observations')
ax.set_ylabel('Temperature [deg C]', fontsize=13)
ax.set_title(f'OI analysis  (sigma_b={sigma_b}, L={L} km, sigma_r={sigma_r})', fontsize=14)
ax.tick_params(labelsize=12)
ax.legend(loc='upper left', fontsize=11)
ax.grid(alpha=0.3)

# ===== プロット2：sigma_b 感度 =====
ax = axes[1]
for sb, lbl, col in [(0.5, 'sigma_b=0.5 (trust model)',  'navy'),
                     (2.0, 'sigma_b=2.0 (standard)',      'crimson'),
                     (5.0, 'sigma_b=5.0 (trust obs)',     'darkorange')]:
    xa_s, _ = optimal_interpolation(x_grid, xb, obs_locs, obs_vals, sb, L, sigma_r)
    ax.plot(x_grid, xa_s, color=col, lw=1.8, label=lbl)

ax.plot(x_grid, x_true, 'k--', lw=1.5, label='True state')
ax.scatter(obs_locs, obs_vals, s=80, color='green', zorder=5)
ax.set_xlabel('Distance [km]', fontsize=13)
ax.set_ylabel('Temperature [deg C]', fontsize=13)
ax.set_title('Sensitivity to background error sigma_b', fontsize=14)
ax.tick_params(labelsize=12)
ax.legend(fontsize=11)
ax.grid(alpha=0.3)

plt.tight_layout()
import os
img_dir = os.path.join(os.path.dirname(__file__), '..', 'img')
os.makedirs(img_dir, exist_ok=True)
out = os.path.join(img_dir, 'fig02_oi_1d.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved: {out}")
plt.show()
