"""
Spatial explanation figures for OI and EnKF.

Generates:
  img/fig_ex01_sensor_and_domain.png  センサ配置と計算領域
  img/fig_ex02_oi_covariance.png      OI の B 行列：センサからの影響範囲
  img/fig_ex03_enkf_localization.png  EnKF のローカライゼーション ρ
  img/fig_ex04_kalman_gain.png        カルマンゲイン K の空間分布（OI vs EnKF）
  img/fig_ex05_innovation.png         イノベーション（観測−予測）と補正の流れ
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
import os
from scipy.ndimage import gaussian_filter

matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'IPAexGothic'

NX, NY = 60, 60
n = NX * NY
rng = np.random.default_rng(42)

# ── ミッフィー正解場 ──────────────────────────────────────────
T_BG=0.0; T_OUTLINE=8.0; T_FUR=100.0; T_EYE=5.0; T_MOUTH=5.0; T_DRESS=55.0; T_COLLAR=70.0
def fill_rect(f, y1, y2, x1, x2, val):
    f[max(0,y1):min(NY-1,y2)+1, max(0,x1):min(NX-1,x2)+1] = val
def fill_ellipse(f, cy, cx, ry, rx, val):
    yy, xx = np.ogrid[:NY,:NX]
    f[((yy-cy)/ry)**2+((xx-cx)/rx)**2<=1.0] = val
def make_miffy():
    f = np.full((NY,NX), T_BG)
    fill_ellipse(f,5,30,9,21,T_DRESS); fill_rect(f,0,11,7,52,T_DRESS)
    fill_rect(f,10,12,24,36,T_COLLAR)
    fill_ellipse(f,48,19,12,9,T_OUTLINE); fill_ellipse(f,48,40,12,9,T_OUTLINE)
    fill_ellipse(f,24,30,15,24,T_OUTLINE)
    fill_ellipse(f,48,19,11,7,T_FUR); fill_ellipse(f,48,40,11,7,T_FUR)
    fill_rect(f,35,40,12,27,T_FUR); fill_rect(f,35,40,33,48,T_FUR)
    fill_ellipse(f,24,30,14,22,T_FUR); fill_rect(f,9,13,25,34,T_FUR)
    fill_ellipse(f,24,21,2,3,T_EYE); fill_ellipse(f,24,39,2,3,T_EYE)
    for d in range(-2,3):
        f[min(max(18+d,0),NY-1), min(max(30+d,0),NX-1)] = T_MOUTH
        f[min(max(18-d,0),NY-1), min(max(30+d,0),NX-1)] = T_MOUTH
    return f

x_true_2d = make_miffy()
x_true    = x_true_2d.flatten()
iy_all    = np.arange(n) // NX
ix_all    = np.arange(n) % NX

BG   = '#111111'
CMAP = 'RdYlBu_r'
FS_TITLE = 16; FS_LABEL = 13; FS_TICK = 11; FS_ANNO = 12

def dark_ax(ax, ticks=False):
    ax.set_facecolor(BG)
    if not ticks:
        ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_edgecolor('gray')

script_dir = os.path.dirname(os.path.abspath(__file__))
img_dir    = os.path.join(script_dir, '..', 'img')
os.makedirs(img_dir, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# fig_ex01: センサ配置と計算領域の説明
# ─────────────────────────────────────────────────────────────────────────────
obs_idx_30 = np.sort(rng.choice(n, size=30,  replace=False))
obs_idx_60 = np.sort(rng.choice(n, size=60,  replace=False))

fig, axes = plt.subplots(1, 3, figsize=(18, 7))
fig.patch.set_facecolor(BG)

# パネル1: 正解場（ミッフィー）
ax = axes[0]; dark_ax(ax)
im = ax.imshow(x_true_2d, origin='lower', cmap=CMAP, vmin=0, vmax=100,
               extent=[0, NX, 0, NY], aspect='equal')
ax.set_title('正解場 x_true（ミッフィー）\nアルゴリズムは直接見えない',
             color='white', fontsize=FS_TITLE, fontweight='bold')
ax.text(30, -5, 'n = 60×60 = 3600 セル', color='#aaaaaa',
        fontsize=FS_ANNO, ha='center', transform=ax.transData,
        clip_on=False)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('温度 [°C]', color='white', fontsize=FS_LABEL)
cbar.ax.tick_params(colors='white', labelsize=FS_TICK)
cbar.ax.yaxis.label.set_color('white'); cbar.outline.set_edgecolor('gray')

# パネル2: EnKF センサ（30点）
ax = axes[1]; dark_ax(ax)
ax.imshow(x_true_2d, origin='lower', cmap=CMAP, vmin=0, vmax=100,
          extent=[0, NX, 0, NY], aspect='equal', alpha=0.35)
ax.scatter(obs_idx_30 % NX + 0.5, obs_idx_30 // NX + 0.5,
           c='lime', s=60, marker='x', linewidths=2, zorder=5, label='センサ位置')
# センサ1点から影響範囲を円で示す
cx, cy = obs_idx_30[0] % NX + 0.5, obs_idx_30[0] // NX + 0.5
circle = plt.Circle((cx, cy), 8, color='cyan', fill=False, lw=1.5,
                     linestyle='--', label='ローカライゼーション半径 r=8')
ax.add_patch(circle)
ax.scatter([cx], [cy], c='cyan', s=100, marker='*', zorder=6)
ax.legend(facecolor='#222222', edgecolor='gray', labelcolor='white',
          fontsize=FS_TICK, loc='upper right')
ax.set_title('EnKF：センサ 30 点/サイクル\nローカライゼーション半径 r=8 cells',
             color='white', fontsize=FS_TITLE, fontweight='bold')

# パネル3: OI センサ（60点）
ax = axes[2]; dark_ax(ax)
ax.imshow(x_true_2d, origin='lower', cmap=CMAP, vmin=0, vmax=100,
          extent=[0, NX, 0, NY], aspect='equal', alpha=0.35)
ax.scatter(obs_idx_60 % NX + 0.5, obs_idx_60 // NX + 0.5,
           c='orange', s=40, marker='x', linewidths=1.5, zorder=5, label='センサ位置')
cx2, cy2 = obs_idx_60[0] % NX + 0.5, obs_idx_60[0] // NX + 0.5
circle2 = plt.Circle((cx2, cy2), 5, color='yellow', fill=False, lw=1.5,
                      linestyle='--', label='ガウス相関長 L=5 cells')
ax.add_patch(circle2)
ax.scatter([cx2], [cy2], c='yellow', s=100, marker='*', zorder=6)
ax.legend(facecolor='#222222', edgecolor='gray', labelcolor='white',
          fontsize=FS_TICK, loc='upper right')
ax.set_title('OI：センサ 60 点/サイクル\nガウス相関長 L=5 cells',
             color='white', fontsize=FS_TITLE, fontweight='bold')

fig.suptitle('計算領域（60×60 グリッド）とセンサ配置',
             color='white', fontsize=FS_TITLE+2, fontweight='bold')
plt.tight_layout()
out = os.path.join(img_dir, 'fig_ex01_sensor_and_domain.png')
fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=BG)
print(f'Saved: {out}'); plt.close(fig)

# ─────────────────────────────────────────────────────────────────────────────
# fig_ex02: OI の B 行列——センサ 1 点から格子全体への影響範囲
# ─────────────────────────────────────────────────────────────────────────────
# センサをグリッド中央に置いたときの BHt の 1 列（その 1 センサ分）
sensor_y, sensor_x = 30, 30  # グリッド中央
dy = iy_all - sensor_y
dx = ix_all - sensor_x
L = 5.0
sigma_b = 20.0  # 例：サイクル5付近の値

BHt_one = sigma_b**2 * np.exp(-0.5 * (dy**2 + dx**2) / L**2)
BHt_map = BHt_one.reshape(NY, NX)

# 距離に対する BHt の断面（y=30 の行）
dist_arr = np.arange(-NX//2, NX//2)
BHt_slice = sigma_b**2 * np.exp(-0.5 * dist_arr**2 / L**2)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.patch.set_facecolor(BG)

# 左: B の空間マップ
ax = axes[0]; dark_ax(ax)
im = ax.imshow(BHt_map, origin='lower', cmap='hot', vmin=0,
               extent=[0, NX, 0, NY], aspect='equal')
ax.scatter([sensor_x+0.5], [sensor_y+0.5], c='cyan', s=150, marker='*',
           zorder=5, label=f'センサ位置 ({sensor_x},{sensor_y})')
for r in [5, 10, 15]:
    c = plt.Circle((sensor_x+0.5, sensor_y+0.5), r, color='white',
                   fill=False, lw=1, linestyle='--', alpha=0.5)
    ax.add_patch(c)
    ax.text(sensor_x+0.5+r, sensor_y+0.5, f'd={r}', color='white',
            fontsize=9, va='center')
ax.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=FS_TICK)
ax.set_title('BHt（セル i ─ センサ j 間の共分散）\nσ_b²×exp(−d²/2L²)',
             color='white', fontsize=FS_TITLE, fontweight='bold')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('[°C²]', color='white', fontsize=FS_LABEL)
cbar.ax.tick_params(colors='white'); cbar.ax.yaxis.label.set_color('white')
cbar.outline.set_edgecolor('gray')

# 中: 断面プロット（距離 vs B 値）
ax = axes[1]; ax.set_facecolor(BG)
ax.plot(dist_arr, BHt_slice, color='orange', lw=2.5)
ax.axvline( L, color='white', lw=1, ls='--', alpha=0.7)
ax.axvline(-L, color='white', lw=1, ls='--', alpha=0.7)
ax.axhline(sigma_b**2 * np.exp(-0.5), color='gray', lw=0.8, ls=':')
ax.text( L+1, sigma_b**2*0.05, f'L={int(L)} cells', color='white', fontsize=FS_ANNO)
ax.text(-NX//2+1, sigma_b**2*np.exp(-0.5)+2, f'exp(−0.5)×σ_b² ≈ 61%', color='gray', fontsize=10)
ax.set_xlabel('センサからの距離 [cells]', color='white', fontsize=FS_LABEL)
ax.set_ylabel('共分散値 [°C²]', color='white', fontsize=FS_LABEL)
ax.set_title('距離と共分散の関係\n（センサを通る水平断面）',
             color='white', fontsize=FS_TITLE, fontweight='bold')
ax.tick_params(colors='white', labelsize=FS_TICK)
ax.grid(alpha=0.3, color='gray'); ax.set_facecolor(BG)
for sp in ax.spines.values(): sp.set_edgecolor('gray')

# 右: センサ → 格子への情報伝播の矢印図
ax = axes[2]; dark_ax(ax)
ax.imshow(BHt_map, origin='lower', cmap='hot', vmin=0,
          extent=[0, NX, 0, NY], aspect='equal', alpha=0.7)
ax.scatter([sensor_x+0.5], [sensor_y+0.5], c='cyan', s=200, marker='*', zorder=6)
# 代表的な格子点への矢印（距離ごと）
for (ty, tx, label) in [(30,35,'d=5\n61%'), (30,40,'d=10\n14%'), (30,45,'d=15\n≈1%')]:
    val = np.exp(-0.5*((ty-sensor_y)**2+(tx-sensor_x)**2)/L**2)*100
    ax.annotate('', xy=(tx+0.5, ty+0.5), xytext=(sensor_x+0.5, sensor_y+0.5),
                arrowprops=dict(arrowstyle='->', color='white', lw=1.5))
    ax.scatter([tx+0.5], [ty+0.5], c='white', s=40, zorder=5)
    ax.text(tx+1, ty+1, f'd={tx-sensor_x}\n{val:.0f}%', color='white', fontsize=9)
ax.set_title('センサ 1 点から格子への情報伝播\n（OI）近いほど大きく修正される',
             color='white', fontsize=FS_TITLE, fontweight='bold')

fig.suptitle('OI の背景誤差共分散 B：センサからの影響範囲（σ_b=20°C, L=5 cells）',
             color='white', fontsize=FS_TITLE+1, fontweight='bold')
plt.tight_layout()
out = os.path.join(img_dir, 'fig_ex02_oi_covariance.png')
fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=BG)
print(f'Saved: {out}'); plt.close(fig)

# ─────────────────────────────────────────────────────────────────────────────
# fig_ex03: EnKF のローカライゼーション ρ vs OI の B
# ─────────────────────────────────────────────────────────────────────────────
R_LOC = 8.0
rho_one = np.exp(-0.5 * (dy**2 + dx**2) / R_LOC**2)
rho_map = rho_one.reshape(NY, NX)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.patch.set_facecolor(BG)

# 左: ρ の空間マップ
ax = axes[0]; dark_ax(ax)
im = ax.imshow(rho_map, origin='lower', cmap='Blues', vmin=0, vmax=1,
               extent=[0, NX, 0, NY], aspect='equal')
ax.scatter([sensor_x+0.5], [sensor_y+0.5], c='red', s=200, marker='*', zorder=6,
           label=f'センサ位置 ({sensor_x},{sensor_y})')
for r in [8, 16, 24]:
    c = plt.Circle((sensor_x+0.5, sensor_y+0.5), r, color='white',
                   fill=False, lw=1.2, linestyle='--', alpha=0.7)
    ax.add_patch(c)
    ax.text(sensor_x+0.5+r, sensor_y+0.5, f'r={r}', color='white',
            fontsize=9, va='center')
ax.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=FS_TICK)
ax.set_title('EnKF ローカライゼーション ρ\nexp(−d²/2r²), r=8 cells',
             color='white', fontsize=FS_TITLE, fontweight='bold')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('ρ の値（0〜1）', color='white', fontsize=FS_LABEL)
cbar.ax.tick_params(colors='white'); cbar.ax.yaxis.label.set_color('white')
cbar.outline.set_edgecolor('gray')

# 中: OI の B（L=5）と EnKF の ρ（r=8）の比較断面
ax = axes[1]; ax.set_facecolor(BG)
d_arr = np.arange(0, 41)
oi_b   = np.exp(-0.5 * d_arr**2 / L**2)
enkf_r = np.exp(-0.5 * d_arr**2 / R_LOC**2)
ax.plot(d_arr, oi_b,   color='orange', lw=2.5, label=f'OI : exp(−d²/2L²) L=5')
ax.plot(d_arr, enkf_r, color='cyan',   lw=2.5, label=f'EnKF: exp(−d²/2r²) r=8')
ax.axhline(0.05, color='gray', lw=0.8, ls=':', label='≈5%（影響ほぼなし）')
ax.set_xlabel('センサからの距離 [cells]', color='white', fontsize=FS_LABEL)
ax.set_ylabel('重みの値（0〜1）', color='white', fontsize=FS_LABEL)
ax.set_title('OI の空間相関 vs EnKF のローカライゼーション\n（断面比較）',
             color='white', fontsize=FS_TITLE, fontweight='bold')
ax.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=FS_TICK)
ax.tick_params(colors='white', labelsize=FS_TICK)
ax.grid(alpha=0.3, color='gray')
for sp in ax.spines.values(): sp.set_edgecolor('gray')

# 右: 2つの役割の説明（テキスト図）
ax = axes[2]; ax.set_facecolor(BG)
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis('off')
props_oi   = dict(boxstyle='round', facecolor='#442200', edgecolor='orange', alpha=0.9)
props_enkf = dict(boxstyle='round', facecolor='#002244', edgecolor='cyan',   alpha=0.9)
ax.text(5, 9.0, 'OI の空間相関（BHt）', color='orange', fontsize=FS_LABEL,
        fontweight='bold', ha='center', bbox=props_oi)
ax.text(5, 7.5, '「近いセルほど同じ温度になる\nはず」という仮定を\nσ_b²×exp(−d²/2L²) で表現\n→ B 行列そのものを計算する',
        color='white', fontsize=10, ha='center', va='top',
        bbox=dict(boxstyle='round', facecolor='#111111', edgecolor='orange', alpha=0.7))
ax.text(5, 4.5, 'EnKF のローカライゼーション（ρ）', color='cyan', fontsize=FS_LABEL,
        fontweight='bold', ha='center', bbox=props_enkf)
ax.text(5, 3.0, 'アンサンブル推定の誤差（N=50 の\nサンプリング誤差）による\n偽の遠距離相関を除去するため\n→ K にかけてフィルターする',
        color='white', fontsize=10, ha='center', va='top',
        bbox=dict(boxstyle='round', facecolor='#111111', edgecolor='cyan', alpha=0.7))
ax.text(5, 0.3, '役割は似て見えるが目的が異なる',
        color='#aaaaaa', fontsize=11, ha='center', style='italic')
ax.set_title('2つの「距離に基づく重み」の役割の違い',
             color='white', fontsize=FS_TITLE, fontweight='bold')

fig.suptitle('センサからの影響範囲：OI（B の形状）vs EnKF（ローカライゼーション ρ）',
             color='white', fontsize=FS_TITLE+1, fontweight='bold')
plt.tight_layout()
out = os.path.join(img_dir, 'fig_ex03_enkf_localization.png')
fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=BG)
print(f'Saved: {out}'); plt.close(fig)

# ─────────────────────────────────────────────────────────────────────────────
# fig_ex04: カルマンゲイン K の空間分布——センサ 1 点が全グリッドをどう更新するか
# ─────────────────────────────────────────────────────────────────────────────
# OI の K の 1 列（センサ j=中央）
M_OBS_OI = 60; SIGMA_R = 5.0; sigma_b_oi = 20.0
# センサを 1 点だけ配置して K の 1 列を計算
obs_y1, obs_x1 = 30, 30
dy1 = iy_all - obs_y1; dx1 = ix_all - obs_x1
BHt1 = sigma_b_oi**2 * np.exp(-0.5*(dy1**2+dx1**2)/L**2)  # (n,)
HBHt1 = sigma_b_oi**2  + SIGMA_R**2  # センサ 1 点のとき: 自分自身の共分散 + R
K_oi_col = BHt1 / HBHt1
K_oi_map = K_oi_col.reshape(NY, NX)

# EnKF の K の 1 列：アンサンブルから模擬（ランダムサンプル）
# 初期アンサンブルを生成
X_init = np.zeros((n, 50))
for i in range(50):
    noise = rng.normal(0, 1, (NY, NX))
    nc = gaussian_filter(noise, sigma=5.0)
    nc = nc / (nc.std()+1e-8) * 30.0
    X_init[:, i] = (50.0 + nc).flatten()
xbar_init = X_init.mean(axis=1, keepdims=True)
Xp_init   = X_init - xbar_init
HXp1 = Xp_init[[obs_y1*NX+obs_x1], :]  # (1, 50)
S1   = (HXp1 @ HXp1.T) / 49 + SIGMA_R**2  # スカラー
K_enkf_col = (Xp_init @ HXp1.T) / 49 / S1  # (n, 1)
rho1 = np.exp(-0.5*(dy1**2+dx1**2)/8.0**2)
K_enkf_col_loc = K_enkf_col.flatten() * rho1
K_enkf_map = K_enkf_col_loc.reshape(NY, NX)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.patch.set_facecolor(BG)

# 左: OI の K（1 センサ分）
ax = axes[0]; dark_ax(ax)
vmax = max(np.abs(K_oi_map).max(), np.abs(K_enkf_map).max())
im = ax.imshow(K_oi_map, origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax,
               extent=[0,NX,0,NY], aspect='equal')
ax.scatter([obs_x1+0.5], [obs_y1+0.5], c='black', s=200, marker='*', zorder=6,
           label=f'センサ ({obs_x1},{obs_y1})')
ax.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=FS_TICK)
ax.set_title('OI のカルマンゲイン K の 1 列\n（センサ 1 点が各セルを修正する強さ）',
             color='white', fontsize=FS_TITLE, fontweight='bold')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('K の値', color='white', fontsize=FS_LABEL)
cbar.ax.tick_params(colors='white'); cbar.ax.yaxis.label.set_color('white')
cbar.outline.set_edgecolor('gray')

# 中: EnKF の K（1 センサ分、ローカライゼーション後）
ax = axes[1]; dark_ax(ax)
im = ax.imshow(K_enkf_map, origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax,
               extent=[0,NX,0,NY], aspect='equal')
ax.scatter([obs_x1+0.5], [obs_y1+0.5], c='black', s=200, marker='*', zorder=6,
           label=f'センサ ({obs_x1},{obs_y1})')
circle_loc = plt.Circle((obs_x1+0.5, obs_y1+0.5), 8, color='white',
                         fill=False, lw=1.5, linestyle='--', alpha=0.7,
                         label='r_loc=8 cells')
ax.add_patch(circle_loc)
ax.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=FS_TICK)
ax.set_title('EnKF のカルマンゲイン K の 1 列\n（ローカライゼーション適用済み）',
             color='white', fontsize=FS_TITLE, fontweight='bold')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('K の値', color='white', fontsize=FS_LABEL)
cbar.ax.tick_params(colors='white'); cbar.ax.yaxis.label.set_color('white')
cbar.outline.set_edgecolor('gray')

# 右: 違いの説明テキスト
ax = axes[2]; ax.set_facecolor(BG)
ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis('off')
ax.text(5, 9.5, 'K[i,j] の意味', color='white', fontsize=FS_LABEL+1,
        fontweight='bold', ha='center')
ax.text(5, 8.5,
        'セルi の温度を修正する量 =\n  K[i,j]  ×  イノベーション[j]\n\nイノベーション[j] =\n  観測値[j] − 推定値[j]',
        color='white', fontsize=11, ha='center', va='top',
        bbox=dict(boxstyle='round', facecolor='#1a1a1a', edgecolor='gray'))
ax.text(5, 5.0, 'OI のK', color='orange', fontsize=FS_LABEL, fontweight='bold', ha='center')
ax.text(5, 4.2,
        '常にガウス型の滑らかな分布\n→形は毎サイクル変わらない\n（σ_b の大きさだけが変わる）',
        color='#cccccc', fontsize=10, ha='center', va='top',
        bbox=dict(boxstyle='round', facecolor='#1a1a1a', edgecolor='orange', alpha=0.7))
ax.text(5, 2.0, 'EnKF のK', color='cyan', fontsize=FS_LABEL, fontweight='bold', ha='center')
ax.text(5, 1.2,
        'アンサンブルの実際の構造を反映\n→ ミッフィーの形に応じて変化\n（境界では急激に変わる）',
        color='#cccccc', fontsize=10, ha='center', va='top',
        bbox=dict(boxstyle='round', facecolor='#1a1a1a', edgecolor='cyan', alpha=0.7))
ax.set_title('OI vs EnKF：\nK の形状の違いと意味',
             color='white', fontsize=FS_TITLE, fontweight='bold')

fig.suptitle('カルマンゲイン K の空間分布：センサ 1 点が全格子セルをどう修正するか',
             color='white', fontsize=FS_TITLE+1, fontweight='bold')
plt.tight_layout()
out = os.path.join(img_dir, 'fig_ex04_kalman_gain.png')
fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=BG)
print(f'Saved: {out}'); plt.close(fig)

# ─────────────────────────────────────────────────────────────────────────────
# fig_ex05: イノベーション——観測と予測のずれが格子全体に伝わる仕組み
# ─────────────────────────────────────────────────────────────────────────────
# 初期推定値（まだミッフィーを知らない状態）
x_init_mean = X_init.mean(axis=1).reshape(NY, NX)
innov_map = np.zeros((NY, NX))

# センサ 5 点を代表例として示す
sample_sensors = obs_idx_30[:5]
for sidx in sample_sensors:
    sy, sx = sidx // NX, sidx % NX
    obs_val  = x_true[sidx]
    pred_val = x_init_mean[sy, sx]
    innov_val = obs_val - pred_val
    # このセンサのカルマンゲイン（OI 簡易版）
    dy_s = iy_all - sy; dx_s = ix_all - sx
    k_col = 25.0 * np.exp(-0.5*(dy_s**2+dx_s**2)/L**2) / (25.0+25.0)
    innov_map += (k_col * innov_val).reshape(NY, NX)

fig, axes = plt.subplots(1, 4, figsize=(22, 6))
fig.patch.set_facecolor(BG)

# パネル1: 正解場
ax = axes[0]; dark_ax(ax)
ax.imshow(x_true_2d, origin='lower', cmap=CMAP, vmin=0, vmax=100,
          extent=[0,NX,0,NY], aspect='equal')
ax.scatter(sample_sensors % NX + 0.5, sample_sensors // NX + 0.5,
           c='lime', s=120, marker='*', zorder=5)
ax.set_title('正解場 x_true\n（センサが測る値）', color='white',
             fontsize=FS_TITLE, fontweight='bold')
for sidx in sample_sensors:
    ax.text(sidx%NX+1, sidx//NX+1, f'{x_true[sidx]:.0f}°C',
            color='lime', fontsize=9, fontweight='bold')

# パネル2: 初期推定値
ax = axes[1]; dark_ax(ax)
ax.imshow(x_init_mean, origin='lower', cmap=CMAP, vmin=0, vmax=100,
          extent=[0,NX,0,NY], aspect='equal')
ax.scatter(sample_sensors % NX + 0.5, sample_sensors // NX + 0.5,
           c='lime', s=120, marker='*', zorder=5)
ax.set_title('初期推定値 x_est\n（まだミッフィーを知らない状態）', color='white',
             fontsize=FS_TITLE, fontweight='bold')
for sidx in sample_sensors:
    ax.text(sidx%NX+1, sidx//NX+1, f'{x_init_mean[sidx//NX, sidx%NX]:.0f}°C',
            color='yellow', fontsize=9, fontweight='bold')

# パネル3: イノベーション = 観測 - 予測
innov_vals = x_true[sample_sensors] - x_init_mean[sample_sensors//NX, sample_sensors%NX]
ax = axes[2]; dark_ax(ax)
ax.imshow(np.zeros((NY,NX)), origin='lower', cmap='gray', vmin=0, vmax=1,
          extent=[0,NX,0,NY], aspect='equal')
for i, sidx in enumerate(sample_sensors):
    sy, sx = sidx // NX, sidx % NX
    iv = innov_vals[i]
    color = '#ff4444' if iv > 0 else '#4444ff'
    ax.scatter([sx+0.5], [sy+0.5], c=color, s=200, marker='*', zorder=5)
    ax.text(sx+1, sy+1, f'Δ={iv:+.0f}°C', color=color, fontsize=10, fontweight='bold')
ax.set_title('イノベーション\n観測値 y − 推定値 Hx\n（赤:+、青:−）', color='white',
             fontsize=FS_TITLE, fontweight='bold')
ax.text(30, 3, '← この差が「修正の指令」になる', color='white',
        fontsize=10, ha='center')

# パネル4: 補正量が格子全体に伝播した結果
vmax_inn = np.abs(innov_map).max()
ax = axes[3]; dark_ax(ax)
im = ax.imshow(innov_map, origin='lower', cmap='RdBu_r',
               vmin=-vmax_inn, vmax=vmax_inn,
               extent=[0,NX,0,NY], aspect='equal')
ax.scatter(sample_sensors % NX + 0.5, sample_sensors // NX + 0.5,
           c='lime', s=120, marker='*', zorder=5, label='センサ位置')
ax.legend(facecolor='#222222', edgecolor='gray', labelcolor='white', fontsize=FS_TICK)
ax.set_title('補正量 K @ innov\nセンサから格子全体へ伝播した\n修正量（赤:+、青:−）', color='white',
             fontsize=FS_TITLE, fontweight='bold')
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('[°C]', color='white', fontsize=FS_LABEL)
cbar.ax.tick_params(colors='white'); cbar.ax.yaxis.label.set_color('white')
cbar.outline.set_edgecolor('gray')

fig.suptitle('イノベーション（観測−予測）がカルマンゲイン K を通じて格子全体へ伝播する仕組み',
             color='white', fontsize=FS_TITLE+1, fontweight='bold')
plt.tight_layout()
out = os.path.join(img_dir, 'fig_ex05_innovation.png')
fig.savefig(out, dpi=150, bbox_inches='tight', facecolor=BG)
print(f'Saved: {out}'); plt.close(fig)

print('\nAll explanation figures saved.')
