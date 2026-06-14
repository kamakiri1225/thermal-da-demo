"""
Pixel-art ossan temperature field — matched to reference image.

Reference features:
  - Pure black background
  - Skin face fills most of grid
  - Short sparse hair marks at crown (not full coverage)
  - Thin drooping eyebrows above large glasses
  - LARGE rectangular glasses with bright cyan lenses
  - Tiny pupils inside glasses
  - Small nose and neutral mouth
  - Dark collar at bottom

Grid: NX=40 (x), NY=60 (y, bottom->top)
OpenFOAM 2D: z=1 cell. Cell index = ix + iy*NX
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

matplotlib.rcParams['font.family'] = 'DejaVu Sans'

NX, NY = 40, 60

# ===== Temperature palette =====
T_BG    =  0.0   # background (black)
T_CLOTH = 10.0   # dark clothing/collar
T_HAIR  = 18.0   # dark hair marks
T_FRAME = 12.0   # glasses frame (dark)
T_EYE   =  6.0   # pupils
T_BROW  = 14.0   # eyebrows
T_MOUTH =  8.0   # mouth
T_LENS  = 48.0   # glasses lens (bright cyan)
T_SKIN  = 100.0  # skin
T_NOSE  =  92.0  # nose

field = np.full((NY, NX), T_BG)

def fill_rect(f, y1, y2, x1, x2, val):
    y1, y2 = max(0, y1), min(NY-1, y2)
    x1, x2 = max(0, x1), min(NX-1, x2)
    f[y1:y2+1, x1:x2+1] = val

def fill_ellipse(f, cy, cx, ry, rx, val):
    y, x = np.ogrid[:NY, :NX]
    mask = ((y - cy) / ry)**2 + ((x - cx) / rx)**2 <= 1.0
    f[mask] = val

# =====================================================
# Ossan — matched to reference pixel art
# Face: y=17..53, x=5..35  (cx=20, cy=35, ry=18, rx=15)
# =====================================================

# --- Dark collar / clothing (bottom) ---
fill_rect(field, 0, 15, 0, 39, T_CLOTH)
fill_ellipse(field, cy=12, cx=20, ry=7, rx=20, val=T_CLOTH)

# --- Neck ---
fill_rect(field, 14, 19, 16, 24, T_SKIN)

# --- Face ---
fill_ellipse(field, cy=35, cx=20, ry=18, rx=15, val=T_SKIN)

# --- Ears ---
fill_ellipse(field, cy=36, cx=4,  ry=4, rx=3, val=T_SKIN)
fill_ellipse(field, cy=36, cx=36, ry=4, rx=3, val=T_SKIN)

# --- Glasses (LARGE — dominant feature) ---
# Left frame (outer border)
fill_rect(field, 32, 47, 7, 18, T_FRAME)
# Left lens (bright cyan fill)
fill_rect(field, 33, 46, 8, 17, T_LENS)
# Right frame
fill_rect(field, 32, 47, 22, 33, T_FRAME)
# Right lens
fill_rect(field, 33, 46, 23, 32, T_LENS)
# Bridge between lenses
fill_rect(field, 38, 41, 19, 21, T_FRAME)
# Ear pieces
fill_rect(field, 39, 40,  4,  7, T_FRAME)   # left ear piece
fill_rect(field, 39, 40, 33, 36, T_FRAME)   # right ear piece

# --- Eyes (tiny pupils inside lenses) ---
fill_rect(field, 39, 40, 11, 12, T_EYE)
fill_rect(field, 39, 40, 27, 28, T_EYE)

# --- Eyebrows (drooping: inner=high, outer=low) ---
# Left brow: /  inner(right) is HIGH, outer(left) is LOW
fill_rect(field, 51, 51, 13, 15, T_BROW)   # inner (high)
fill_rect(field, 50, 50, 11, 13, T_BROW)   # middle
fill_rect(field, 49, 49,  9, 11, T_BROW)   # outer (low)
# Right brow: \  inner(left) is HIGH, outer(right) is LOW
fill_rect(field, 51, 51, 25, 27, T_BROW)   # inner (high)
fill_rect(field, 50, 50, 27, 29, T_BROW)   # middle
fill_rect(field, 49, 49, 29, 31, T_BROW)   # outer (low)

# --- Hair: short sparse marks at crown (buzz-cut style) ---
# ~6 small dark strokes staggered across top of head
hair_marks = [
    (52, 11), (53, 14), (52, 17),
    (53, 20), (52, 23), (53, 27), (52, 30),
]
for hy, hx in hair_marks:
    fill_rect(field, hy, hy+1, hx, hx+1, T_HAIR)

# --- Nose (simple, just below glasses) ---
fill_ellipse(field, cy=28, cx=20, ry=2, rx=2, val=T_NOSE)
# Nostrils
fill_rect(field, 26, 27, 17, 18, T_BROW)
fill_rect(field, 26, 27, 22, 23, T_BROW)

# --- Mouth (thin neutral line) ---
fill_rect(field, 21, 21, 15, 25, T_MOUTH)

# ===== Stats =====
palette = {
    T_BG:'BG', T_CLOTH:'CLOTH', T_HAIR:'HAIR', T_FRAME:'FRAME',
    T_EYE:'EYE', T_BROW:'BROW', T_MOUTH:'MOUTH',
    T_LENS:'LENS', T_SKIN:'SKIN', T_NOSE:'NOSE'
}
print("Temperature distribution:")
for v in sorted(np.unique(field)):
    lbl = palette.get(v, '?')
    print(f"  T={v:6.1f}  ({lbl:6s}): {np.sum(field==v):5d} cells")

# ===== Custom colormap (matches pixel-art reference) =====
color_map = {
    T_BG:    [0.00, 0.00, 0.00],   # pure black background
    T_MOUTH: [0.10, 0.05, 0.05],   # dark mouth
    T_EYE:   [0.05, 0.05, 0.05],   # near black pupils
    T_CLOTH: [0.08, 0.08, 0.10],   # near black clothing
    T_FRAME: [0.08, 0.08, 0.08],   # dark glasses frame
    T_BROW:  [0.15, 0.08, 0.02],   # dark brown brows/nostrils
    T_HAIR:  [0.12, 0.07, 0.02],   # dark brown hair marks
    T_LENS:  [0.05, 0.85, 1.00],   # bright CYAN lenses
    T_NOSE:  [0.96, 0.75, 0.58],   # nose skin
    T_SKIN:  [0.97, 0.79, 0.63],   # skin tone (peach)
}

rgb = np.zeros((NY, NX, 3))
for t_val, color in color_map.items():
    mask = (np.abs(field - t_val) < 0.5)
    rgb[mask] = color

# ===== Plot =====
fig, axes = plt.subplots(1, 2, figsize=(16, 8.5))

ax = axes[0]
ax.imshow(rgb, origin='lower', extent=[0, NX, 0, NY],
          aspect='equal', interpolation='nearest')
ax.set_title('Ossan Pixel Art  (reference match)', fontsize=12, fontweight='bold',
             color='white')
ax.set_facecolor('black')
fig.patch.set_facecolor('#111111')
ax.set_xlabel('x [cells]', color='white', fontsize=13)
ax.set_ylabel('y [cells]', color='white', fontsize=13)
ax.tick_params(colors='white', labelsize=12)
ax.title.set_fontsize(15)
ax.set_xticks(np.arange(0, NX+1, 5))
ax.set_yticks(np.arange(0, NY+1, 5))
ax.grid(alpha=0.2, color='gray', lw=0.5)
for spine in ax.spines.values():
    spine.set_edgecolor('gray')

ax2 = axes[1]
im2 = ax2.imshow(field, origin='lower', cmap='RdYlBu_r',
                 vmin=0, vmax=100, extent=[0, NX, 0, NY],
                 aspect='equal', interpolation='nearest')
cbar = plt.colorbar(im2, ax=ax2, label='Temperature [deg C]')
cbar.ax.tick_params(labelsize=12)
cbar.ax.yaxis.label.set_size(13)
ax2.set_title('Temperature Field (DA target)', fontsize=15, fontweight='bold')
ax2.set_xlabel('x [cells]', fontsize=13)
ax2.set_ylabel('y [cells]', fontsize=13)
ax2.tick_params(labelsize=12)
ax2.set_xticks(np.arange(0, NX+1, 5))
ax2.set_yticks(np.arange(0, NY+1, 5))
ax2.grid(alpha=0.2, color='gray', lw=0.5)

plt.suptitle('True State: Ossan with Glasses  (40x60 grid)',
             fontsize=18, fontweight='bold')
plt.tight_layout()

script_dir = os.path.dirname(os.path.abspath(__file__))
img_dir = os.path.join(script_dir, '..', 'img')
os.makedirs(img_dir, exist_ok=True)
out_png = os.path.join(img_dir, 'fig01_ossan_true_field.png')
plt.savefig(out_png, dpi=150, bbox_inches='tight', facecolor=fig.get_facecolor())
print(f"\nSaved PNG: {out_png}")
plt.show()

# ===== Write OpenFOAM 0/T =====
values  = field.flatten(order='C')
n_cells = len(values)
assert n_cells == NX * NY

of_dir = os.path.join(script_dir, '..', '0')
os.makedirs(of_dir, exist_ok=True)
out_T = os.path.join(of_dir, 'T')

header = """\
/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     |
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      T;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];

"""

internal = f"internalField   nonuniform List<scalar>\n{n_cells}\n(\n"
internal += "\n".join(f"{v:.4g}" for v in values)
internal += "\n);\n"

boundary = """
boundaryField
{
    top          { type zeroGradient; }
    bottom       { type zeroGradient; }
    left         { type zeroGradient; }
    right        { type zeroGradient; }
    frontAndBack { type empty; }
}

// ************************************************************************* //
"""

with open(out_T, 'w') as f:
    f.write(header + internal + boundary)

print(f"Saved OpenFOAM 0/T: {out_T}")
print(f"  n_cells={n_cells}  T_range=[{values.min():.0f}, {values.max():.0f}]")
