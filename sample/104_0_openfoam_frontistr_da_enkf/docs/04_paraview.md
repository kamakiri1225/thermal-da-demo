# ParaView で見る(温度分布・変位分布・時刻歴アニメーション)

`run/export_paraview.py` が、データ同化の結果を ParaView で開ける VTU / PVD に書き出す。

```
cd sample/104_0_openfoam_frontistr_da_enkf
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/export_paraview.py
```

出力(`paraview/`、再生成可能なので Git 管理外):

| ファイル | 中身 |
|----------|------|
| `temperature_fields.vtu` | 固体実メッシュ(20696セル)。セルデータ `T_garbage_degC / T_da_degC / T_truth_degC` |
| `displacement_fields.vtu` | 円筒メッシュ(5040節点)。点データ `U_garbage_m / U_da_m / U_truth_m`, `Uz_*_um` |
| `temperature_timehistory.pvd` | **0→600s の温度時刻歴アニメーション**(31ステップ)。`T_da / T_free / T_truth_degC` |

---

## 1. データ同化後の3D温度分布(静止画)

1. ParaView で `temperature_fields.vtu` を開く(Apply)
2. Coloring を **`T_da_degC`**(データ同化後)に
3. Filters → **Clip**(Plane, Y方向)で断面にすると内部分布が見える
4. `T_truth_degC` に切り替えると真値と一致することが確認できる。
   `T_garbage_degC` は一様(でたらめ)

## 2. データ同化された変位分布(変形形状)

1. `displacement_fields.vtu` を開く(Apply)
2. Filters → **Warp By Vector** → Vectors に **`U_da_m`**、Scale Factor を大きく(例 2000)
3. Coloring を **`Uz_da_um`** にするとヒータ側が持ち上がる傾きが見える
4. Vectors を `U_truth_m` に変えると真値の変形と一致。`U_garbage_m` は一様膨張で過大

## 3. 時刻歴アニメーション(0→600s)★

1. **`temperature_timehistory.pvd`** を開く(Apply)
2. Coloring を **`T_da_degC`** に、Clip で断面表示
3. 上部の **再生ボタン(▶)** で 0→600s を再生
   - ヒータ ON(0-300s)で温まり、OFF(300-600s)で冷える
   - `T_free_degC`(同化なし)に切り替えると、でたらめなまま真値からズレ続ける
   - `T_da_degC`(同化あり)は序盤で真値に吸い付き、以降ずっと一致

> 時刻歴の温度場は ROM(5ノード)を円筒メッシュへ逆距離補間して作っている。
> OpenFOAM 版のフル解像度時刻歴が要る場合は、`da_openfoam_config.yaml` の
> `t_end_s` を伸ばして `run/run_openfoam_fem_enkf.py` を回す(計算は重い)。

---

## 補足: 色スケールを揃える

3つの配列(garbage/da/truth)を見比べるときは、ParaView の
**Rescale to Custom Range** で全ケース同じ範囲にすると比較しやすい。
でたらめは桁違いに大きいので、da と truth を比較するときは真値の範囲に合わせるとよい。
