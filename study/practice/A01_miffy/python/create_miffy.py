"""
Miffy (Nijntje) pixel-art temperature field.
Grid: NX=40 (x), NY=60 (y, bottom->top)
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

matplotlib.rcParams['font.family'] = 'DejaVu Sans'

NX, NY = 40, 60

T_BG      =  0.0
T_OUTLINE =  8.0
T_FUR     = 100.0
T_EYE     =  5.0
T_MOUTH   =  5.0
T_DRESS   = 55.0
T_COLLAR  = 70.0

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
# Miffy layout (y=0 = bottom):
#   y= 0-10: dress
#   y=10-37: face  (cy=24, ry=14, rx=15) <- wider/fluffier
#   y=37-57: ears  (cy=48, ry=11, rx=5)  <- shorter than before
# =====================================================

# --- Dress / body ---
fill_ellipse(field, 5, 20, 9, 14, T_DRESS)
fill_rect(field, 0, 11, 5, 35, T_DRESS)
# Collar bow
fill_rect(field, 10, 12, 16, 24, T_COLLAR)

# --- Ear outlines (dark border) ---
fill_ellipse(field, 48, 13, 12, 6, T_OUTLINE)   # left ear outline
fill_ellipse(field, 48, 27, 12, 6, T_OUTLINE)   # right ear outline

# --- Face outline (wider = ふっくら) ---
fill_ellipse(field, 24, 20, 15, 16, T_OUTLINE)  # face outline

# --- Ear fills (white) ---
fill_ellipse(field, 48, 13, 11, 5, T_FUR)        # left ear
fill_ellipse(field, 48, 27, 11, 5, T_FUR)        # right ear

# --- Ear-face junction fills (connect ears to head) ---
fill_rect(field, 35, 40, 8,  18, T_FUR)          # left ear base connects to head
fill_rect(field, 35, 40, 22, 32, T_FUR)          # right ear base

# --- Face fill (draw last so it covers junction neatly) ---
fill_ellipse(field, 24, 20, 14, 15, T_FUR)       # face (wide, fluffy)

# --- Neck ---
fill_rect(field, 9, 13, 17, 23, T_FUR)

# --- Eyes (small dark ovals, centered in face) ---
fill_ellipse(field, 24, 14, 2, 2, T_EYE)
fill_ellipse(field, 24, 26, 2, 2, T_EYE)

# --- Mouth: x cross centered at (y=18, x=20) ---
cy_m, cx_m = 18, 20
for d in range(-2, 3):
    field[min(max(cy_m + d, 0), NY-1), min(max(cx_m + d, 0), NX-1)] = T_MOUTH
    field[min(max(cy_m - d, 0), NY-1), min(max(cx_m + d, 0), NX-1)] = T_MOUTH

# ===== Stats =====
palette = {T_BG:'BG', T_OUTLINE:'OUTLINE', T_FUR:'FUR',
           T_EYE:'EYE', T_MOUTH:'MOUTH', T_DRESS:'DRESS', T_COLLAR:'COLLAR'}
print("Temperature distribution:")
for v in sorted(np.unique(field)):
    lbl = palette.get(v, '?')
    print(f"  T={v:6.1f}  ({lbl:7s}): {np.sum(field==v):5d} cells")

# ===== Custom colormap =====
color_map = {
    T_BG:      [0.00, 0.00, 0.00],
    T_OUTLINE: [0.06, 0.04, 0.02],
    T_MOUTH:   [0.05, 0.05, 0.05],
    T_EYE:     [0.05, 0.05, 0.05],
    T_DRESS:   [0.30, 0.60, 1.00],
    T_COLLAR:  [1.00, 0.85, 0.10],
    T_FUR:     [1.00, 0.99, 0.96],
}
rgb = np.zeros((NY, NX, 3))
for t_val, color in color_map.items():
    mask = (np.abs(field - t_val) < 0.5)
    rgb[mask] = color

# ===== Plot =====
fig, axes = plt.subplots(1, 2, figsize=(16, 8.5))
fig.patch.set_facecolor('#111111')

# ---- Left: pixel art (black background, white text) ----
ax = axes[0]
ax.set_facecolor('black')
ax.imshow(rgb, origin='lower', extent=[0, NX, 0, NY],
          aspect='equal', interpolation='nearest')
ax.set_title('Miffy Pixel Art', fontsize=15, fontweight='bold', color='white')
ax.set_xlabel('x [cells]', color='white', fontsize=13)
ax.set_ylabel('y [cells]', color='white', fontsize=13)
ax.tick_params(colors='white', labelsize=12)
for spine in ax.spines.values():
    spine.set_edgecolor('gray')
ax.set_xticks(np.arange(0, NX+1, 5))
ax.set_yticks(np.arange(0, NY+1, 5))
ax.grid(alpha=0.2, color='gray', lw=0.5)

# White labels on black bg
for lbl, y, x in [
    ('Ear L', 48, 13), ('Ear R', 48, 27),
    ('Eye L', 24, 14), ('Eye R', 24, 26),
    ('x mouth', 18, 20), ('Dress', 5, 20),
]:
    ax.text(x, y, lbl, ha='center', va='center', fontsize=10,
            color='white', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.1', fc='black', alpha=0.6))

# ---- Right: heatmap ----
ax2 = axes[1]
ax2.set_facecolor('#111111')
im2 = ax2.imshow(field, origin='lower', cmap='RdYlBu_r',
                 vmin=0, vmax=100, extent=[0, NX, 0, NY],
                 aspect='equal', interpolation='nearest')
cbar = plt.colorbar(im2, ax=ax2, label='Temperature [deg C]')
cbar.ax.yaxis.label.set_color('white')
cbar.ax.tick_params(colors='white', labelsize=12)
cbar.ax.yaxis.label.set_size(13)
cbar.outline.set_edgecolor('gray')
ax2.set_title('Temperature Field (DA target)', fontsize=15, fontweight='bold', color='white')
ax2.set_xlabel('x [cells]', color='white', fontsize=13)
ax2.set_ylabel('y [cells]', color='white', fontsize=13)
ax2.tick_params(colors='white', labelsize=12)
for spine in ax2.spines.values():
    spine.set_edgecolor('gray')
ax2.set_xticks(np.arange(0, NX+1, 5))
ax2.set_yticks(np.arange(0, NY+1, 5))
ax2.grid(alpha=0.2, color='gray', lw=0.5)

plt.suptitle('True State: Miffy (Nijntje)  (40x60 grid)',
             fontsize=18, fontweight='bold', color='white')
plt.tight_layout()

script_dir = os.path.dirname(os.path.abspath(__file__))
img_dir = os.path.join(script_dir, '..', 'img')
os.makedirs(img_dir, exist_ok=True)
out_png = os.path.join(img_dir, 'fig01_miffy_true_field.png')
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

print(f"Saved: {out_T}  ({n_cells} cells, T=[{values.min():.0f},{values.max():.0f}])")
