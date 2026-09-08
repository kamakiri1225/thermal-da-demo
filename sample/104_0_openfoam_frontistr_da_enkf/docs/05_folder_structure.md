# フォルダ構成ガイド ― どこに何があるか

104(OpenFOAM×FrontISTR データ同化)フォルダの中身を、フォルダごとに整理する。
「どこを見ればいいか分からない」を解消するための地図。

```
104_0_openfoam_frontistr_da_enkf/
├── README.md              … このケースの概要(まずここ)
├── requirements.txt       … 必要な Python ライブラリ
├── config/                … 設定ファイルと校正結果
├── dacore/                … 【ROM版】高速な代用モデル + データ同化の中身
├── daof/                  … 【OpenFOAM版】実ソルバを回す連成コード
├── fem/                   … FrontISTR 変位の観測演算子
├── run/                   … 実行スクリプト(ここを python3 で叩く)
├── docs/                  … 解説書(.md)と図(img/)
├── results/               … 数値結果(CSV/YAML)
├── openfoam/              … OpenFOAM ケース(計算の重い中身。Git管理外)
└── paraview/              … ParaView 用 VTU/PVD(Git管理外, 再生成可能)
```

---

## 各フォルダの役割

### `config/` … 設定と校正結果
| ファイル | 中身 |
|----------|------|
| `da_config.yaml` | データ同化の設定(時間範囲・観測点・ノイズ・アンサンブル数など) |
| `rom_calibrated.yaml` | ROM を 102_0 に校正した結果(熱容量 C・コンダクタンス K・放熱 h) |
| `displacement_operator.yaml` | 温度→変位の線形写像 D(102_1 FrontISTR に校正) |
| `openfoam_102_0_temperature_history.csv` | 校正に使う 102_0 の温度履歴 |

### `dacore/` … 【ROM版】高速モデル + データ同化アルゴリズム
数秒で回る縮約モデルと、フィルタ本体。**仕組みの理解・お試しはこちら**。
| ファイル | 中身 |
|----------|------|
| `cht_rom.py` | 5ノード集中定数モデル(前進計算) |
| `calibrate.py` | ROM を 102_0 に最小二乗校正 |
| `displacement.py` | 温度→上面変位の線形オペレータ(FrontISTR校正) |
| `enkf.py` | アンサンブルカルマンフィルタ(解析更新) |
| `pf.py` | 粒子フィルタ(参考) |
| `ensemble.py` | でたらめ初期アンサンブルの生成・前進 |
| `observations.py` | 真値ラン・合成観測の生成 |
| `twin.py` / `twin_fem.py` | 双子実験ドライバ(twin_fem=温度+変位版) |
| `node_locations.py` | 5ノードの3次元座標 |
| `plots.py` | 図の共通設定(日本語フォント・文字サイズ) |

### `daof/` … 【OpenFOAM版】実ソルバ連成
各メンバーが実際に chtMultiRegionFoam を回す、本物の field-space データ同化。
| ファイル | 中身 |
|----------|------|
| `of_io.py` | OpenFOAM フィールドの読み書き |
| `of_case.py` | メンバーケースの生成・実行・場の読み出し |
| `of_fem_twin.py` | 温度+変位観測の双子実験ドライバ |
| `of_plots.py` | 図の描画 |

### `fem/` … FrontISTR 変位の観測演算子
| ファイル | 中身 |
|----------|------|
| `fem_obs.py` | 温度場→FrontISTR熱膨張→上面変位(102_1のコードを再利用) |

### `run/` … 実行スクリプト(**ここを叩く**)
| ファイル | 何をする |
|----------|----------|
| `run_rom_fem.py` | 【ROM・数秒】0→600s の温度+変位データ同化 → `docs/img/rom_fem_*.png` |
| `run_openfoam_fem_enkf.py` | 【OpenFOAM・重い】実ソルバでデータ同化 → `docs/img/openfoam_fem_*.png` |
| `plot_da_temperature_field3d.py` | データ同化後の3D温度分布 → `da_temperature_field3d.png` |
| `plot_da_displacement_field.py` | データ同化された変位分布 → `da_displacement_field.png` |
| `plot_rom_field.py` | ROMの温度分布(点群) |
| `export_paraview.py` | ParaView用 VTU/PVD を書き出す → `paraview/` |

### `docs/` … 解説書と図 ★ここが読み物
| ファイル | 内容 |
|----------|------|
| `00_temp_displacement_da.md` | 温度+変位同化の仕組みと数式 |
| `01_method_walkthrough.md` | OpenFOAM版の1サイクルの流れ(実データで) |
| `02_beginner_guide.md` | **大学初学者向け超ていねい解説(まずこれ)** |
| `03_rom_derivation.md` | ROMの数学的導出(PDE→ODE、校正、なぜ5点で一致するか) |
| `04_paraview.md` | ParaViewでの見方(分布・時刻歴アニメ) |
| `05_folder_structure.md` | このファイル(フォルダ地図) |

#### `docs/img/` … 図(PNG、GitHubで表示される)
| 図 | 内容 |
|----|------|
| `rom_fem_temperature.png` | ROM 5点温度の時刻歴(0→600s, でたらめ/同化/真値) |
| `rom_fem_displacement.png` | ROM 上面変位の時刻歴 |
| `rom_fem_rmse.png` | ROM 誤差の収束(同化あり vs なし) |
| `da_temperature_field3d.png` | データ同化後の3D温度分布(実メッシュ断面) |
| `da_displacement_field.png` | データ同化された変位分布(変形形状) |
| `openfoam_fem_enkf_rmse.png` | OpenFOAM版 温度場RMSEの収束 |
| `openfoam_fem_enkf_Q.png` | OpenFOAM版 発熱量Qの推定 |
| `openfoam_fem_enkf_disp.png` | OpenFOAM版 変位観測の収束 |
| `openfoam_fem_enkf_field.png` | OpenFOAM版 温度場スナップショット |

### `results/` … 数値結果
| ファイル | 中身 |
|----------|------|
| `openfoam_fem_enkf_history.csv` | サイクルごとの RMSE・Q・変位 |
| `openfoam_fem_enkf_summary.yaml` | 最終結果の要約 |

### `openfoam/` … OpenFOAM ケース(**Git管理外・重い**)
`setup_base_case.sh` で 102_0 のメッシュを複製した `base_case/` と、
実行時に生成される `run_fem_enkf/`(真値・各メンバーのケース)が入る。
再実行で作り直せるのでリポジトリには含めない。

### `paraview/` … ParaView 用ファイル(**Git管理外・再生成可能**)
`run/export_paraview.py` で生成。
| ファイル | 中身 |
|----------|------|
| `temperature_fields.vtu` | データ同化後の3D温度分布(garbage/da/truth) |
| `displacement_fields.vtu` | データ同化された変位(Warp By Vectorで変形表示) |
| `temperature_timehistory.pvd` | **0→600s 時刻歴アニメーション**(これをParaViewで開く) |
| `timehistory/da_t*.vtu` | 上記アニメの各時刻の実体(31ステップ) |

---

## どこから始めるか(目的別)

- **とにかく理解したい** → `docs/02_beginner_guide.md` を読む
- **数式で納得したい** → `docs/03_rom_derivation.md` / `docs/00_temp_displacement_da.md`
- **すぐ動かしたい(数秒)** → `python3 run/run_rom_fem.py` → `docs/img/rom_fem_*.png`
- **ParaViewで見たい** → `python3 run/export_paraview.py` → `paraview/temperature_timehistory.pvd` を開く(手順 `docs/04_paraview.md`)
- **本物のOpenFOAMで回したい(重い)** → `docs/01_method_walkthrough.md` の手順
