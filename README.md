# thermal-da-demo

熱データ同化を、Python の小さなモデルから OpenFOAM 連携へ段階的に試すためのデモプロジェクトです。

最終的な関心は、少数の温度センサと熱解析モデルを組み合わせて、測っていない場所の温度場を推定し、工作機械などの熱変位補償へつなげることです。

## 全体の目的

このリポジトリでは、次の流れを段階的に確認します。

```text
少数の温度センサ
  ↓
熱解析モデルによる予測
  ↓
データ同化による補正
  ↓
未測定点を含む温度場の推定
  ↓
熱変位・加工誤差補償への応用
```

実機で全点の温度を測ることは現実的ではありません。そこで、物理モデルの予測とセンサ観測を組み合わせ、測っていない場所の状態を推定する方法を検証します。

## 現在のサンプル

各サンプルの手法、目的、実行方法は [sample/README.md](sample/README.md) にまとめています。

| サンプル | 内容 | 主な手法 | スライド | 状態 |
|---|---|---|---|---|
| `sample/001_kalman_thermal_1d` | Python だけで動く 1 次元熱伝導デモ | Kalman filter | [▶ スライドを開く](https://kamakiri1225.github.io/thermal-da-demo/sample/001_kalman_thermal_1d/slides.html) | 実装済み |
| `sample/002-0_laplacian_da_1d` | OpenFOAM `laplacianFoam` と Python をつなぐデータ同化デモ | optimal interpolation | [▶ スライドを開く](https://kamakiri1225.github.io/thermal-da-demo/sample/002-0_laplacian_da_1d/slides.html) | 実装済み |
| `sample/002-1_laplacian_da_round_bar` | 丸棒の軸方向温度場を OpenFOAM / FrontISTR の2ソルバーで OI データ同化し比較するデモ | optimal interpolation | [▶ スライドを開く](https://kamakiri1225.github.io/thermal-da-demo/sample/002-1_laplacian_da_round_bar/slides.html) | 実装済み |

## リポジトリ構成

```text
thermal-da-demo/
├── README.md
├── sample/
│   ├── README.md
│   ├── 001_kalman_thermal_1d/
│   │   ├── main.py
│   │   ├── kalman_filter.py
│   │   ├── thermal_model.py
│   │   ├── slides.html
│   │   └── docs/
│   ├── 002-0_laplacian_da_1d/
│   │   ├── of_interface.py
│   │   ├── case_base/
│   │   ├── oi/
│   │   ├── kf/
│   │   ├── slides.html
│   │   └── docs/
│   ├── 002-1_laplacian_da_round_bar/     ← OpenFOAM 側
│   │   ├── README.md
│   │   ├── of_interface.py
│   │   ├── case_base/
│   │   ├── assets/                        ← スライド用画像・GIF
│   │   ├── slides.html
│   │   └── oi/
│   └── 002-1-FrontISTR_round_bar_da/     ← FrontISTR 側
│       ├── README.md
│       ├── fistr_interface.py
│       ├── case/
│       └── oi/
├── docs/
└── IoT/
```

## ロードマップ

| フェーズ | 内容 | 状態 |
|---|---|---|
| Phase 1 | Python のみで 1D 熱モデル + Kalman filter を確認 | 完了 |
| Phase 2a | OpenFOAM + OI データ同化 (1D 擬似モデル、2点センサ) | 完了 |
| Phase 2b | OpenFOAM + FrontISTR の 2 ソルバーで OI DA を実装・比較 (丸棒モデル) | 完了 |
| Phase 3 | アルミ丸棒片持ちはり実機実験 × EnKF パラメータ同定 DA (温度のみ vs 温度+変位) | 次に実施 |
| Phase 4 | 熱流体・熱構造連成解析 (OpenFOAM + FrontISTR/CalculiX) との DA 統合 | 計画中 |
| Phase 5 | DA 較正済みモデルによるサロゲートモデル構築 → リアルタイム熱変位予測 | 計画中 |

> 研究の核となる問い: **「温度センサを増やすより、変位センサを 1 点加える方が熱変位予測には効くか？」**  
> 詳細は [docs/research_roadmap.md](docs/research_roadmap.md) を参照。

## ドキュメント

- サンプル一覧と手法の説明: [sample/README.md](sample/README.md)
- 001 のスライド: [▶ スライドを開く (GitHub Pages)](https://kamakiri1225.github.io/thermal-da-demo/sample/001_kalman_thermal_1d/slides.html)
- 002 のスライド: [▶ スライドを開く (GitHub Pages)](https://kamakiri1225.github.io/thermal-da-demo/sample/002-0_laplacian_da_1d/slides.html)
- 002-1 のスライド: [▶ スライドを開く (GitHub Pages)](https://kamakiri1225.github.io/thermal-da-demo/sample/002-1_laplacian_da_round_bar/slides.html)
- 001 の詳細: [sample/001_kalman_thermal_1d/docs/00_index.md](sample/001_kalman_thermal_1d/docs/00_index.md)
- 002 の詳細: [sample/002-0_laplacian_da_1d/docs/00_index.md](sample/002-0_laplacian_da_1d/docs/00_index.md)
- 002-1 の詳細: [sample/002-1_laplacian_da_round_bar/README.md](sample/002-1_laplacian_da_round_bar/README.md)
- 先行研究メモ: [docs/](docs/)

## ライセンス

MIT License
