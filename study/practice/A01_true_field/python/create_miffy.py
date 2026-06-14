"""
Miffy (Nijntje) pixel-art temperature field for data assimilation demo.

Grid: NX=40 (x), NY=60 (y, bottom->top)
OpenFOAM 2D: z=1 cell (frontAndBack=empty)
Cell index = ix + iy*NX
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

matplotlib.rcParams['font.family'] = 'DejaVu Sans'

NX, NY = 40, 60

# ===== Temperature palette =====
T_BG    =  0.0   # background
T_FUR   = 100.0  # white fur (face, ears, body)
T_EYE   =  5.0   # black eyes
T_MOUTH =  5.0   # black mouth (x cross)
T_DRESS = 55.0   # colored dress
T_RIBBON= 75.0   # ribbon / accessory (optional)

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
# Miffy  (Nijntje by Dick Bruna)
# Simple, iconic: round face, long ears, x mouth
# =====================================================

# --- Dress / body ---
fill_ellipse(field, cy=8,  cx=20, ry=10, rx=13, val=T_DRESS)
fill_rect(field, 0, 16, 7, 33, T_DRESS)

# --- Collar ribbon ---
fill_rect(field, 15, 17, 15, 25, T_RIBBON)
fill_ellipse(field, cy=16, cx=20, ry=2, rx=3, val=T_RIBBON)

# --- Neck ---
fill_rect(field, 17, 20, 17, 23, T_FUR)

# --- Ears (long upright ovals — Miffy's signature!) ---
fill_ellipse(field, cy=53, cx=13, ry=9, rx=4, val=T_FUR)   # left ear
fill_ellipse(field, cy=53, cx=27, ry=9, rx=4, val=T_FUR)   # right ear

# --- Face (large round oval) ---
fill_ellipse(field, cy=33, cx=20, ry=15, rx=14, val=T_FUR)

# --- Eyes (small oval dots) ---
fill_ellipse(field, cy=38, cx=15, ry=2, rx=2, val=T_EYE)
fill_ellipse(field, cy=38, cx=25, ry=2, rx=2, val=T_EYE)

# --- Mouth: iconic x (cross) shape ---
# / diagonal: bottom-left -> top-right
# \ diagonal: top-left -> bottom-right
cy_m, cx_m = 29, 20
for d in range(-2, 3):
    yd = min(max(cy_m + d, 0), NY-1)
    xd = min(max(cx_m + d, 0), NX-1)
    field[yd, xd] = T_MOUTH                           # / line
    yd2 = min(max(cy_m - d, 0), NY-1)
    field[yd2, xd] = T_MOUTH                          # \ line

# Thicken the cross slightly (1 cell wide becomes 2)
for d in range(-1, 2):
    xd = min(max(cx_m + d, 0), NX-1)
    yd_up   = min(max(cy_m + abs(d), 0), NY-1)
    yd_down = min(max(cy_m - abs(d), 0), NY-1)
    field[yd_up,   xd] = T_MOUTH
    field[yd_down, xd] = T_MOUTH

# ===== Stats =====
palette = {T_BG:'BG', T_FUR:'FUR', T_EYE:'EYE',
           T_MOUTH:'MOUTH', T_DRESS:'DRESS', T_RIBBON:'RIBBON'}
print("Temperature distribution:")
for v in sorted(np.unique(field)):
    lbl = palette.get(v, '?')
    print(f"  T={v:6.1f}  ({lbl:6s}): {np.sum(field==v):5d} cells")

# ===== Custom colormap =====
color_map = {
    T_BG:    [0.75, 0.88, 1.00],   # sky blue background
    T_MOUTH: [0.05, 0.05, 0.05],   # black
    T_EYE:   [0.05, 0.05, 0.05],   # black
    T_DRESS: [0.20, 0.55, 0.95],   # blue dress
    T_RIBBON:[0.95, 0.40, 0.40],   # red/pink ribbon
    T_FUR:   [1.00, 1.00, 0.98],   # near white fur
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
ax.set_title('Miffy Pixel Art (true colors)', fontsize=15, fontweight='bold')
ax.set_xlabel('x [cells]', fontsize=13)
ax.set_ylabel('y [cells]', fontsize=13)
ax.tick_params(labelsize=12)
ax.set_xticks(np.arange(0, NX+1, 5))
ax.set_yticks(np.arange(0, NY+1, 5))
ax.grid(alpha=0.2, color='gray', lw=0.5)
for lbl, y, x in [
    ('Ear L', 53, 13), ('Ear R', 53, 27),
    ('Eye L', 38, 15), ('Eye R', 38, 25),
    ('x mouth', 29, 20),
    ('Dress',   8, 20), ('Ribbon', 16, 20),
]:
    ax.text(x, y, lbl, ha='center', va='center', fontsize=10,
            color='black', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.6))

ax2 = axes[1]
im2 = ax2.imshow(field, origin='lower', cmap='RdYlBu_r',
                 vmin=0, vmax=100, extent=[0, NX, 0, NY],
                 aspect='equal', interpolation='nearest')
cbar = plt.colorbar(im2, ax=ax2, label='Temperature [deg C]')
cbar.ax.tick_params(labelsize=12)
cbar.ax.yaxis.label.set_size(13)
ax2.set_title('Temperature Field (for data assimilation)', fontsize=15, fontweight='bold')
ax2.set_xlabel('x [cells]', fontsize=13)
ax2.set_ylabel('y [cells]', fontsize=13)
ax2.tick_params(labelsize=12)
ax2.set_xticks(np.arange(0, NX+1, 5))
ax2.set_yticks(np.arange(0, NY+1, 5))
ax2.grid(alpha=0.2, color='gray', lw=0.5)

plt.suptitle('True State: Miffy (Nijntje)  (40x60 grid)', fontsize=18, fontweight='bold')
plt.tight_layout()

script_dir = os.path.dirname(os.path.abspath(__file__))
img_dir = os.path.join(script_dir, '..', 'img')
os.makedirs(img_dir, exist_ok=True)
out_png = os.path.join(img_dir, 'fig02_miffy_true_field.png')
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\nSaved PNG: {out_png}")
plt.show()

# ===== Write OpenFOAM 0/T_miffy =====
values  = field.flatten(order='C')
n_cells = len(values)
assert n_cells == NX * NY

of_dir = os.path.join(script_dir, '..', '0')
os.makedirs(of_dir, exist_ok=True)
out_T = os.path.join(of_dir, 'T_miffy')   # separate from ossan 0/T

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

print(f"Saved OpenFOAM 0/T_miffy: {out_T}")
print(f"  n_cells={n_cells}  T_range=[{values.min():.0f}, {values.max():.0f}]")
print("\nTo use as target: cp 0/T_miffy 0/T  (in practice/A01_true_field/)")
