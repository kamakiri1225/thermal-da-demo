# 実行手順(粒子フィルタ版)

構成は 103_0(EnKF)と同じ。まず軽い ROM 版で挙動を確認し、その後 OpenFOAM 連成版を回す。

## 0. 依存関係

```
pip install --user -r requirements.txt
```

OpenFOAM 連成版は OpenFOAM v2512(`chtMultiRegionFoam` が PATH にあること)が必要。
WSL2 では `OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4` を付けて実行する。

## 1. ROM 版(高速)

```
cd sample/103_1_openfoam_da_pf
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 dacore/calibrate.py   # 初回のみ
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/run_pf.py         # PF 双子実験
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/compare_enkf_pf.py # EnKF vs PF 比較図
```

出力: `docs/img/pf_state_convergence.png` ほか、`docs/img/enkf_vs_pf_rmse.png`、
`results/pf_history.csv`、`results/compare_summary.yaml`。

> `run/compare_enkf_pf.py` は隣の `103_0_openfoam_da_enkf/config/da_config.yaml` を
> 参照するため、両フォルダが並んで存在している必要がある。

## 2. OpenFOAM 連成版

```
cd sample/103_1_openfoam_da_pf/openfoam
sh setup_base_case.sh          # 102_0 のメッシュ入りケースを複製
cd ..
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
  nohup python3 run/run_openfoam_pf.py > openfoam/run_pf.log 2>&1 &
tail -f openfoam/run_pf.log
```

出力: `docs/img/openfoam_pf_rmse.png`、`openfoam_pf_Q.png`、`openfoam_pf_field.png`、
`results/openfoam_pf_history.csv` / `openfoam_pf_summary.yaml`。

規模は `openfoam/da_openfoam_config.yaml`(`ensemble.n_members` など)で調整。
PF は次元の呪いに弱いので、安定させたい場合は粒子数を EnKF より多めにする。

## 3. 注意

- 既存67kメッシュ・輻射あり・直列では、1メンバーの短ウィンドウでも数分かかる。
  まず最小実行(N=5, 0→20 s, 2サイクル)で確認すること。
- 真値・観測は seed で決まるので、EnKF(103_0)と PF(103_1)は同一観測を同化する。
