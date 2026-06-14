"""
スカラー版OI法：カルマンゲインの理論値とモンテカルロの比較
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'DejaVu Sans'

np.random.seed(42)
N_trial = 100_000

sigma_b = 2.0
sigma_r = 1.0
x_true  = 5.0

# --- 理論値 ---
K_theory       = sigma_b**2 / (sigma_b**2 + sigma_r**2)
sigma_a_theory = np.sqrt((1 - K_theory) * sigma_b**2)

print("=" * 40)
print(f"  sigma_b = {sigma_b},  sigma_r = {sigma_r}")
print(f"  Kalman gain K   (theory) : {K_theory:.4f}")
print(f"  Analysis error  (theory) : {sigma_a_theory:.4f}")

# --- モンテカルロ ---
eps_b = np.random.normal(0, sigma_b, N_trial)
eps_r = np.random.normal(0, sigma_r, N_trial)
x_b   = x_true + eps_b
y     = x_true + eps_r
x_a   = x_b + K_theory * (y - x_b)
sigma_a_mc = np.std(x_a - x_true)

print(f"  Analysis error  (MC)     : {sigma_a_mc:.4f}")
print("=" * 40)

# --- B/R 比を変えたときの K の変化 ---
sigma_b_range = np.linspace(0.1, 5.0, 200)
K_range = sigma_b_range**2 / (sigma_b_range**2 + sigma_r**2)

# --- プロット ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# (1) 誤差ヒストグラム比較
bins = np.linspace(-8, 8, 80)
axes[0].hist(eps_b,        bins=bins, density=True, color='steelblue',  alpha=0.7, label=f'Background  sigma={sigma_b}')
axes[0].hist(eps_r,        bins=bins, density=True, color='darkorange', alpha=0.7, label=f'Observation sigma={sigma_r}')
axes[0].hist(x_a - x_true, bins=bins, density=True, color='crimson',    alpha=0.8, label=f'Analysis    sigma={sigma_a_mc:.3f}')
axes[0].set_title('Error distributions', fontsize=14)
axes[0].set_xlabel('Error', fontsize=13)
axes[0].set_ylabel('Probability density', fontsize=13)
axes[0].tick_params(labelsize=12)
axes[0].legend(fontsize=11)
axes[0].axvline(0, color='k', lw=1)
axes[0].grid(alpha=0.3)

# (2) K の sigma_b 依存性
axes[1].plot(sigma_b_range, K_range, 'b-', lw=2)
axes[1].axvline(sigma_b, color='r', lw=1.5, ls='--', label=f'sigma_b={sigma_b}')
axes[1].axhline(K_theory, color='r', lw=1.5, ls='--', label=f'K={K_theory:.3f}')
axes[1].set_xlabel('Background error sigma_b', fontsize=13)
axes[1].set_ylabel('Kalman gain K', fontsize=13)
axes[1].set_title('K = sigma_b^2 / (sigma_b^2 + sigma_r^2)', fontsize=14)
axes[1].tick_params(labelsize=12)
axes[1].legend(fontsize=11)
axes[1].grid(alpha=0.3)

# (3) 解析誤差分散 sigma_a の sigma_b 依存性
sigma_a_range = np.sqrt((1 - K_range) * sigma_b_range**2)
axes[2].plot(sigma_b_range, sigma_b_range, 'b--', lw=1.5, label='sigma_b (background)')
axes[2].axhline(sigma_r, color='darkorange', lw=1.5, ls='--', label=f'sigma_r={sigma_r} (obs)')
axes[2].plot(sigma_b_range, sigma_a_range, 'r-', lw=2, label='sigma_a (analysis)')
axes[2].scatter([sigma_b], [sigma_a_theory], s=80, color='crimson', zorder=5)
axes[2].set_xlabel('Background error sigma_b', fontsize=13)
axes[2].set_ylabel('Error std', fontsize=13)
axes[2].set_title('Analysis error < min(sigma_b, sigma_r)', fontsize=14)
axes[2].tick_params(labelsize=12)
axes[2].legend(fontsize=11)
axes[2].grid(alpha=0.3)

plt.suptitle(f'Scalar OI: K={K_theory:.3f}, sigma_a theory={sigma_a_theory:.4f}, MC={sigma_a_mc:.4f}',
             fontsize=16)
plt.tight_layout()
import os
img_dir = os.path.join(os.path.dirname(__file__), '..', 'img')
os.makedirs(img_dir, exist_ok=True)
out = os.path.join(img_dir, 'fig01_oi_scalar.png')
plt.savefig(out, dpi=150, bbox_inches='tight')
print(f"Saved: {out}")
plt.show()
