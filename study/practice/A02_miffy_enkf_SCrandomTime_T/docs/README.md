# A02_miffy_enkf_T — EnKF デモ：ミッフィーが浮かび上がる

## ドキュメント一覧

| ファイル | 内容 |
|---------|------|
| `README.md`（このファイル） | フォルダ概要・問題設定・パラメータ |
| [`enkf_method.md`](enkf_method.md) | **アンサンブルカルマンフィルタ（EnKF）の詳細説明** |

> OI（最適内挿法）は `../A02_miffy_oi_T/` を参照してください。

---

## 概要

アンサンブルカルマンフィルタ（EnKF）によるデータ同化デモ。
ランダムな初期アンサンブルから、スパース観測センサを使って
ミッフィーの温度分布を再構築します。

```
初期アンサンブル（空間相関ノイズ、N_ENS=50 メンバー）
     ↓  センサ観測（30点/サイクル, ランダム配置）
  Cycle 1: かすかな輪郭が現れ始める
     ↓
  Cycle 10: ミッフィーの形が見え始める
     ↓
  Cycle 50: 正解分布に収束（RMSE ≈ 17.6°C）
```

---

## フォルダ構成

```
A02_miffy_enkf_T/
├── python/
│   ├── enkf_miffy.py           ← EnKF 実装
│   └── make_explanation_figs.py ← 説明図生成スクリプト
├── img/
│   ├── anim_enkf_miffy.gif        ← EnKF アニメーション
│   ├── fig01_enkf_snapshots.png   ← 各サイクルのスナップショット
│   ├── fig02_rmse_convergence.png ← RMSE 収束カーブ
│   ├── fig05_enkf_vs_oi.png       ← OI との RMSE 比較
│   ├── enkf_rmse.npy              ← RMSE データ（oi_miffy.py が参照）
│   └── fig_ex0*.png               ← 説明図（説明用）
└── docs/
    ├── README.md          ← このファイル
    └── enkf_method.md     ← EnKF の理論・アルゴリズム詳細
```

---

## 問題設定

| 項目 | 値 |
|------|-----|
| 状態ベクトル次元 | $n = 3600$（60×60 グリッドの温度値） |
| 正解場 | ミッフィーの温度分布（T_BG=0°C, T_FUR=100°C, T_DRESS=55°C） |
| DA サイクル数 | 50 サイクル |

## EnKF パラメータ（enkf_miffy.py）

| パラメータ | 値 | 意味 |
|-----------|-----|------|
| `N_ENS` | 50 | アンサンブルサイズ |
| `M_OBS` | 30 | 観測点数/サイクル |
| `SIGMA_R` | 5.0 | 観測誤差 [°C] |
| `R_LOC` | 8.0 | ローカライゼーション半径 [cells] |
| `INFL` | 1.05 | 乗算インフレーション係数 |
| `N_CYC` | 50 | DAサイクル数 |

---

## 実行方法

```bash
cd python
python enkf_miffy.py     # EnKF 実行 → img/ に図と enkf_rmse.npy を出力
```
