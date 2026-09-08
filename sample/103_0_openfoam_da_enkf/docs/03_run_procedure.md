# 実行手順(他人が確実に再現できるように)

このフォルダには **2つの実装**がある。まず軽い ROM 版で挙動を確認し、
その後 OpenFOAM 連成版を回すのがおすすめ。

- **ROM 版(高速・数秒)**: 集中定数モデルによる高速プロトタイプ。`dacore/`
- **OpenFOAM 連成版(本命・数十分)**: 各メンバーが実際に chtMultiRegionFoam を回す。`daof/`

---

## 0. 依存関係

```
pip install --user -r requirements.txt
# numpy / scipy / matplotlib / PyYAML、可視化に pyvista
```

OpenFOAM 連成版は **OpenFOAM v2512** 環境(`chtMultiRegionFoam` が PATH にあること)が必要。

```
source /usr/lib/openfoam/openfoam2512/etc/bashrc   # 環境に合わせて
which chtMultiRegionFoam
```

WSL2 では numpy/matplotlib のスレッド暴走を防ぐため、実行時に
`OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4` を付ける。

---

## 1. ROM 版(まずこちらで挙動確認)

```
cd sample/103_0_openfoam_da_enkf
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 dacore/calibrate.py   # 102_0 にROMを校正(初回のみ)
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/run_enkf.py       # EnKF 双子実験
```

出力:
- 図: `docs/img/enkf_state_convergence.png`(温度の補正時刻歴)、`enkf_rmse.png`、`enkf_params.png`
- 数値: `results/enkf_history.csv`、`results/enkf_summary.yaml`

観測点・未観測点の3D図:

```
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/plot_probe_locations.py
# -> docs/img/probe_locations_3d.png, probe_locations_top.png
```

---

## 2. OpenFOAM 連成版(本命)

### 2.1 ベースケースを用意(既存メッシュを再利用)

事前に `../102_0_openfoam_hollow_cylinder_heat_transfer` を `./Allrun.pre` 済みにしておく
(polyMesh と 0/ が必要)。その上で:

```
cd sample/103_0_openfoam_da_enkf/openfoam
sh setup_base_case.sh          # 102_0 の constant/system/0 を base_case へコピー
cd ..
```

### 2.2 実行(計算が長いのでバックグラウンド推奨)

```
cd sample/103_0_openfoam_da_enkf
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
  nohup python3 run/run_openfoam_enkf.py > openfoam/run_enkf.log 2>&1 &
```

進捗の確認:

```
tail -f openfoam/run_enkf.log                       # 各サイクルの RMSE / Q
grep '^Time = ' openfoam/run_enkf/member_00/log.run_* | tail   # ソルバの到達時刻
```

出力:
- 図: `docs/img/openfoam_enkf_rmse.png`、`openfoam_enkf_Q.png`、`openfoam_enkf_field.png`
- 数値: `results/openfoam_enkf_history.csv`、`results/openfoam_enkf_summary.yaml`
- 中間場: `openfoam/run_enkf/fields/*.npy`(可視化用、Git 管理外)

### 2.3 規模の調整

`openfoam/da_openfoam_config.yaml` を編集:

- `ensemble.n_members`: メンバー数(多いほど場の再現性↑・計算↑)
- `experiment.t_end_s` / `obs_interval_s`: 時間範囲と観測間隔
- `observation.probes`: 観測点の追加・変更

> **目安(既存67kメッシュ・輻射あり・直列)**: 起動直後は約 30–40 秒/(実時間1秒)。
> 最小実行(N=5, 0→20 s, 2サイクル)で truth + アンサンブル合わせて数十分程度。
> フル解像度・多メンバーにすると数時間規模になるので、まず最小実行で確認すること。

---

## 3. 粒子フィルタ(PF)との比較

PF 版は `../103_1_openfoam_da_pf` にある。同一の真値・観測で EnKF と PF を比較できる。
詳細は `103_1_openfoam_da_pf/docs/04_enkf_vs_pf.md` を参照。
