"""ani/ の GIF を、ヒータ（+X 側）が右手前に見える向きで描き直す（pvpython で実行）。

  左：OpenFOAM の温度（固体）と、円筒の軸を通る断面での空気の流速の矢印（矢印の色＝空気の温度）
  右：FrontISTR の変位（2000 倍に変形して表示、色＝変位の大きさ）
  元の ani/ani.gif と ani/ani.*.png はそのまま残し、ani_heater_front/ に出力する。

  実行（ParaView 6 系の pvpython）：
    python3 ../107_machining_center_rom/common/pvpython_launch.py python/make_ani_heater_front.py
  時刻の文字は TIME_POS=upper（左上、既定）か TIME_POS=lower（左下）。コマは ani_heater_front/frames_{upper,lower}/ に出る。
  そのあと python3 ani_heater_front/make_gif.py で GIF にする。
"""
import math
import os
from paraview.simple import *  # noqa: F401,F403

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)
FISTR = os.path.join(os.path.dirname(CASE), '102_1_frontistr_hollow_cylinder_thermal_expansion', 'case')
POS = os.environ.get('TIME_POS', 'upper')            # 時刻の文字の場所：upper（左上）／lower（左下）
OUT = os.path.join(CASE, 'ani_heater_front', f'frames_{POS}')
os.makedirs(OUT, exist_ok=True)
SCALE = 2000.0
TMIN, TMAX = 20.0, 26.0
UMAX = 6.5   # µm

# 106 の図と同じ向き：+X（ヒータ）が右手前、-Y 側から斜め上に見る
FOCAL = [0.0, 0.0, 0.05]
DIR = [190.0, -230.0, 210.0]
n = math.sqrt(sum(d * d for d in DIR))
DIR = [d / n for d in DIR]

times = [5 * k for k in range(121)]
only = os.environ.get('ONLY')
if only:
    times = [int(only)]

# ---- 左：OpenFOAM ----
foam = OpenFOAMReader(FileName=os.path.join(CASE, 'hollowCylinder.foam'))
foam.MeshRegions = ['/fluid/internalMesh', '/solid/internalMesh']; foam.SkipZeroTime = 0
foam.CellArrays = ['T', 'U']
foam.UpdatePipelineInformation()
fluid = ExtractBlock(Input=foam); fluid.Selectors = ['/Root/fluid']
solid = ExtractBlock(Input=foam); solid.Selectors = ['/Root/solid']
solidC = Calculator(Input=CellDatatoPointData(Input=MergeBlocks(Input=solid)))
solidC.ResultArrayName = 'T_C'; solidC.Function = 'T-273.15'
fl = CellDatatoPointData(Input=MergeBlocks(Input=fluid))
flC = Calculator(Input=fl); flC.ResultArrayName = 'T_C'; flC.Function = 'T-273.15'
# 円筒の軸を通り、視線に向いた断面（円筒の奥に見える）
nx, ny = DIR[0], DIR[1]; nn = math.hypot(nx, ny)
sl = Slice(Input=flC); sl.SliceType = 'Plane'
sl.SliceType.Origin = [0, 0, 0]; sl.SliceType.Normal = [nx / nn, ny / nn, 0]
mask = MaskPoints(Input=sl); mask.OnRatio = 6; mask.MaximumNumberofPoints = 4000
mask.GenerateVertices = 0
gl = Glyph(Input=mask, GlyphType='Arrow'); gl.OrientationArray = ['POINTS', 'U']
gl.ScaleArray = ['POINTS', 'U']; gl.ScaleFactor = 0.25; gl.GlyphMode = 'All Points'

layout = CreateLayout('two')
vL = CreateView('RenderView'); AssignViewToLayout(view=vL, layout=layout, hint=0)
layout.SplitHorizontal(0, 0.5)
vR = CreateView('RenderView'); AssignViewToLayout(view=vR, layout=layout, hint=2)
for v in (vL, vR):
    v.Background = [1, 1, 1]; v.OrientationAxesVisibility = 1
    try:
        v.UseColorPaletteForBackground = 0
    except AttributeError:
        pass

lutT = GetColorTransferFunction('T_C'); lutT.ApplyPreset('Turbo', True)
lutT.RescaleTransferFunction(TMIN, TMAX); lutT.AutomaticRescaleRangeMode = 'Never'
dS = Show(solidC, vL); ColorBy(dS, ('POINTS', 'T_C')); dS.LookupTable = lutT; dS.Ambient = 0.35; dS.Diffuse = 0.8
# 背景を白にするため、断面の面は描かず、断面上の流速の矢印を空気の温度で色付けする
dG = Show(gl, vL); ColorBy(dG, ('POINTS', 'T_C')); dG.LookupTable = lutT
bar = GetScalarBar(lutT, vL); bar.Title = 'T [degC]'; bar.ComponentTitle = ''
bar.Visibility = 1; bar.TitleFontSize = 22; bar.LabelFontSize = 20; bar.TitleColor = bar.LabelColor = [0, 0, 0]
bar.UseCustomLabels = 1; bar.CustomLabels = [20, 21, 22, 23, 24, 25, 26]; bar.AddRangeLabels = 0; bar.LabelFormat = '%.0f'; bar.AutomaticLabelFormat = 0
bar.WindowLocation = 'Any Location'; bar.Position = [0.86, 0.25]; bar.ScalarBarLength = 0.5
txt = Text(); dT = Show(txt, vL); dT.FontSize = 46; dT.Color = [0, 0, 0]; dT.WindowLocation = 'Upper Left Corner' if POS == 'upper' else 'Lower Left Corner'; dT.Bold = 1

# ---- 右：FrontISTR ----
def fistr_file(t):
    return os.path.join(FISTR, f't_{t}', 'hollow_cylinder_thermal_expansion_vis_psf.0001.pvtu')


fr = XMLPartitionedUnstructuredGridReader(FileName=[fistr_file(times[0])])
warp = WarpByVector(Input=fr); warp.Vectors = ['POINTS', 'DISPLACEMENT']; warp.ScaleFactor = SCALE
um = Calculator(Input=warp); um.ResultArrayName = 'u_um'; um.Function = 'mag(DISPLACEMENT)*1e6'
dW = Show(um, vR); ColorBy(dW, ('POINTS', 'u_um'))
lutU = GetColorTransferFunction('u_um'); lutU.ApplyPreset('Jet', True)
lutU.RescaleTransferFunction(0.0, UMAX); lutU.AutomaticRescaleRangeMode = 'Never'
dW.LookupTable = lutU; dW.Ambient = 0.35; dW.Diffuse = 0.8
fe = FeatureEdges(Input=ExtractSurface(Input=fr)); fe.FeatureAngle = 30
dO = Show(fe, vR); dO.ColorArrayName = ['POINTS', '']
dO.AmbientColor = dO.DiffuseColor = [0.55, 0.55, 0.55]
barU = GetScalarBar(lutU, vR); barU.Title = '|u| [um]'; barU.ComponentTitle = ''
barU.Visibility = 1; barU.TitleFontSize = 22; barU.LabelFontSize = 20; barU.TitleColor = barU.LabelColor = [0, 0, 0]
barU.WindowLocation = 'Any Location'; barU.Position = [0.86, 0.25]; barU.ScalarBarLength = 0.5
barU.UseCustomLabels = 1; barU.CustomLabels = [0, 1, 2, 3, 4, 5, 6]; barU.AddRangeLabels = 0; barU.LabelFormat = '%.0f'; barU.AutomaticLabelFormat = 0
txtR = Text(); txtR.Text = f'deformation x{SCALE:.0f} (gray = undeformed)'
dTR = Show(txtR, vR); dTR.FontSize = 22; dTR.Color = [0, 0, 0]; dTR.WindowLocation = 'Upper Center'

def set_cam(v, dist, zoom_focal):
    c = v.GetActiveCamera()
    v.CameraFocalPoint = zoom_focal
    v.CameraPosition = [zoom_focal[i] + dist * DIR[i] for i in range(3)]
    v.CameraViewUp = [0, 0, 1]; v.CameraParallelProjection = 0


set_cam(vL, 0.45, FOCAL)
set_cam(vR, 0.30, FOCAL)
layout.SetSize(1496, 764)

for t in times:
    foam.UpdatePipeline(float(t))
    for v in (vL,):
        v.ViewTime = float(t)
    txt.Text = f'Time: {t:.1f} s' + ('  (heater ON 15 W)' if t <= 300 else '  (heater OFF)')
    fr.FileName = [fistr_file(t)]
    lutT.RescaleTransferFunction(TMIN, TMAX); lutU.RescaleTransferFunction(0.0, UMAX)   # 全コマ同じ色の範囲
    SaveScreenshot(os.path.join(OUT, f'ani.{t // 5:04d}.png'), layout, ImageResolution=[1496, 764])
    print('frame', t, flush=True)
