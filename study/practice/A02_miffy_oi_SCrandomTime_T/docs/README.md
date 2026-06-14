# A02_miffy_oi_T — OI デモ：ミッフィーが浮かび上がる

## ドキュメント一覧

| ファイル | 内容 |
|---------|------|
| `README.md`（このファイル） | フォルダ概要・問題設定・パラメータ |
| [`oi_method.md`](oi_method.md) | **最適内挿法（OI）の詳細説明** |

> EnKF は `../A02_miffy_enkf_T/`、カルマンフィルタ（熱伝導方程式）は `../A002_miffy_kf_T/` を参照してください。

---

## 概要

最適内挿法（OI: Optimal Interpolation）によるデータ同化デモ。
ランダムな初期推定値から、スパース観測センサを使って
ミッフィーの温度分布を再構築します。

```
初期状態（空間相関ノイズ、平均 50°C）
     ↓  センサ観測（60点/サイクル, ランダム配置）
  Cycle 1: かすかな輪郭が現れ始める
     ↓
  Cycle 10: ミッフィーの形が見え始める
     ↓
  Cycle 50: 正解分布に収束（RMSE ≈ 18.7°C）
```

---

## フォルダ構成

```
A02_miffy_oi_T/
├── python/
│   └── oi_miffy.py          ← OI 実装
├── img/
│   ├── anim_oi_miffy.gif       ← OI アニメーション
│   ├── fig03_oi_snapshots.png  ← 各サイクルのスナップショット
│   ├── fig04_oi_rmse.png       ← RMSE 収束カーブ
│   ├── fig05_enkf_vs_oi.png    ← EnKF との RMSE 比較（再生成時出力先）
│   └── fig_ex0*.png            ← 説明図
└── docs/
    ├── README.md    ← このファイル
    └── oi_method.md ← OI の理論・アルゴリズム詳細
```

---

## 問題設定

| 項目 | 値 |
|------|-----|
| 状態ベクトル次元 | $n = 3600$（60×60 グリッドの温度値） |
| 正解場 | ミッフィーの温度分布（T_BG=0°C, T_FUR=100°C, T_DRESS=55°C） |
| DA サイクル数 | 50 サイクル |

## OI パラメータ（oi_miffy.py）

| パラメータ | 値 | 意味 |
|-----------|-----|------|
| `M_OBS` | 60 | 観測点数/サイクル |
| `SIGMA_R` | 5.0 | 観測誤差 [°C] |
| `SIGMA_B0` | 30.0 | 初期背景誤差 [°C] |
| `SIGMA_B_MIN` | 3.0 | 最小背景誤差 [°C] |
| `DECAY` | 0.88 | σ_b 減衰率/サイクル |
| `L_CORR` | 5.0 | ガウス相関長 [cells] |
| `N_CYC` | 50 | DA サイクル数 |

---

## 実行方法

```bash
cd python
python oi_miffy.py    # OI 実行 → img/ に図を出力
                      # ※ EnKF との比較図は A02_miffy_enkf_T/python/enkf_miffy.py 実行後に生成
```
