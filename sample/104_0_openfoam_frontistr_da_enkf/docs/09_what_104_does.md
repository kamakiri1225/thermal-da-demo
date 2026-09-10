# 104 で行っていること(概要)

## 一言でいうと

> **中空円筒をヒータで温める実験を想定し、「数点の温度・変位の測定」だけから
> 円筒全体の温度分布と発熱量を推定する ― これをアンサンブルカルマンフィルタ
> (EnKF)で実現する。**

「初期状態も発熱量も分からない」でたらめな出発点からでも、数点の観測で
正解(真値)に補正できることを、実際に計算して示すのが104の目的。

---

## 何を推定するのか

| | 内容 |
|---|---|
| **推定したいもの(状態)** | 固体の温度分布 ＋ ヒータ発熱量 Q |
| **手に入る手がかり(観測)** | 表面2点の温度(hot/cold)＋ 上面2点の変位(FrontISTR) |
| **出発点** | でたらめな初期温度(10〜50℃バラバラ)＋ 間違った Q |
| **やること** | 観測が来るたびにEnKFで状態を補正し、真値へ寄せる |

温度は熱電対、変位は変位センサに対応。**測っていない場所(円筒内部・未観測点)まで、
観測どうしの相関を使って同時に補正される**のがポイント。

---

## どうやるのか(EnKFの1サイクル)

```
① 予報   … たくさんの「分身(アンサンブル)」を物理モデルで少し前進させる
② 観測   … 真値の数点にノイズを乗せた測定値が届く
③ 解析   … 分身のばらつき(標本共分散)を使い、観測へ引き寄せる(カルマン更新)
④ 書戻し … 補正した状態を次の予報の出発点にする
        ↑ これを観測のたびに繰り返す
```

分身のばらつきが「どの向きに直すべきか(相関)」を教えてくれるので、
数点の観測で全体を直せる。解析の中身は `dacore/enkf.py`。

---

## 2つのやり方で実施している

同じEnKFを、**前進モデル**を変えて2通りで回している。

> **前進モデル**とは「今の温度 → 少し後の温度」と**時間を前へ進める計算**のこと
> (予報モデル)。天気予報で「今の大気から明日を計算するシミュレーション」に相当。
> EnKFの①予報の工程がこれ。

| | ROM版 | OpenFOAM実機版 |
|---|---|---|
| 前進モデル(時間を進める計算) | 5点の集中定数モデル(NumPy) | chtMultiRegionFoam(実CFD, 20696セル) |
| 変位 | 温度→変位の線形写像 | FrontISTR(実FEM) |
| 規模 | 60メンバー × 30サイクル | 5メンバー × 2サイクル |
| 速さ | 数秒 | 数十分 |
| 用途 | 仕組みの理解・全時間(0-600s)の可視化 | 本物の精度・分布 |
| 解析(EnKF) | **`dacore/enkf.py`(共通)** | **`dacore/enkf.py`(共通)** |

ROMは本物(102_0/102_1)に校正済みなので挙動は本物そっくり。
「軽くて分かりやすいROM」と「重くて正確なOpenFOAM」を住み分けている。

### どのフォルダでやっているのか

両方とも同じ104フォルダの中。使うサブフォルダが違うだけ。

| 役割 | ROM版 | OpenFOAM実機版 |
|------|-------|----------------|
| 実行の入口 | `run/run_rom_fem.py` | `run/run_openfoam_fem_enkf.py` |
| 前進モデル | `dacore/cht_rom.py` | `daof/of_case.py` → `chtMultiRegionFoam` |
| 実ケースの置き場 | (不要, 純NumPy) | `openfoam/run_fem_enkf/`(Git管理外) |
| 変位 | `dacore/displacement.py` | `fem/fem_obs.py` → `fistr1` |
| 双子実験の司令塔 | `dacore/twin_fem.py` | `daof/of_fem_twin.py` |
| 解析(EnKF) | `dacore/enkf.py` | `dacore/enkf.py`(同じファイル) |

まとめると:
- **ROM版 = `dacore/` の中**(実ソルバを呼ばない)
- **OpenFOAM版 = `daof/` + `fem/` + `openfoam/`**
- **EnKF解析だけは両方とも `dacore/enkf.py` を共有**

```
ROM版      : run/run_rom_fem.py          → dacore/(cht_rom, displacement, twin_fem, enkf)
OpenFOAM版 : run/run_openfoam_fem_enkf.py → daof/ + fem/ + openfoam/  ＋  dacore/enkf.py(解析のみ共有)
```

---

## 何が示せたか(結果)

- **温度分布**: でたらめ初期(RMSE約10K)→ 数十秒で真値に一致(ROM 0.039K / OpenFOAM 0.020K)
- **発熱量 Q**: 間違った値 → 真値15Wへ回復
- **変位分布**: 温度が直ると、そこから計算される変形形状も真値に一致
- **未観測点も補正**: mid/top/core や円筒内部まで真値に一致
- **変位観測の効果**: 温度2点だけより、変位も足す方が温度場推定が鋭くなる
  (OpenFOAM: 0.085K → 0.020K)

代表的な図:
- `docs/img/rom_fem_temperature.png` … 温度の時刻歴(でたらめ/同化/真値)
- `docs/img/da_vs_truth_temperature.gif` … 温度分布の時刻歴アニメ(同化 vs 真値)
- `docs/img/da_temperature_field3d.png` … データ同化後の3D温度分布
- `docs/img/da_displacement_field.png` … データ同化された変位分布

---

## もっと知るには

- **理解優先** → [`02_beginner_guide.md`](02_beginner_guide.md)(初学者向け超ていねい)
- **作業の流れ・成果物の由来** → [`08_what_was_computed.md`](08_what_was_computed.md)
- **数式** → [`00_temp_displacement_da.md`](00_temp_displacement_da.md) /
  [`03_rom_derivation.md`](03_rom_derivation.md) /
  `../103_0_openfoam_da_enkf/docs/05_theory_to_code.md`
- **フォルダ地図** → [`05_folder_structure.md`](05_folder_structure.md)
