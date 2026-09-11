# データ同化・熱変位予測の研究資料

## 読む順番

1. [研究方針](overview/research_roadmap.md): 研究の問い、実験・解析計画、実装ロードマップ
2. [先行研究まとめ](overview/research_summary.md): 工作機械の熱誤差補償とセンサ配置の整理
3. [Kalman Filter論文レビュー](reviews/paper1_kalman_filter.md): 要点版
4. [Kalman Filter論文の詳細日本語資料](reviews/paper1_kalman_filter_ja.md): 図を含む詳細版
5. [センサ配置論文レビュー](reviews/paper2_sensor_placement.md): 要点版
6. [センサ配置論文の詳細日本語資料](reviews/paper2_sensor_placement_ja.md): 図を含む詳細版

## 構成

```text
research/
├── README.md
├── overview/               研究方針・先行研究全体
├── reviews/                論文の要点・日本語資料・紙面画像
├── papers/                 原著論文・参考PDF
│   ├── general/
│   ├── thermal_boundary/
│   ├── sensor_placement/
│   └── structural_contact/
└── topics/
    └── contact_boundary/   接触境界条件同定の研究案
```

## 参考論文

### データ同化・状態推定

- [Data Assimilation](<papers/general/Data Assimilation.pdf>)
- [Kalman Filter-Driven State Observer for Thermal Error Compensation](<papers/general/Kalman Filter-Driven State Observer for.pdf>)

### 熱パラメータ・熱境界条件の推定

- [Determination of Thermal Parameters using Thermography and Data Assimilation](<papers/thermal_boundary/Determination of Thermal Parameters using Thermography and Data Assimilation and its Application to the Convective Heat Transfer Coefficient.pdf>)
- [Estimation of Thermal Boundary Conditions in Rotating Bearing by Data Assimilation](<papers/thermal_boundary/Estimation of Thermal Boundary Conditions in Rotating Bearing by Data Assimilation.pdf>)

### センサ配置

- [Sensor Placement](papers/sensor_placement/SensorPlacement.pdf)

### 構造モデルの接触境界条件

- [切削時間時の時間応答に基づく接触部境界条件の推定](papers/structural_contact/切削時間時の時間応答に基づくデータ同化による構造モデルの接触部境界条件の推定.pdf)
- [時間応答に基づく接触面境界条件の逐次推定・更新](papers/structural_contact/時間応答に基づく構造モデルの接触面境界条件の逐次推定・更新システムの構築.pdf)

## 研究テーマ案

- [切削応答を使う接触境界条件同定](topics/contact_boundary/contact_boundary_da_cutting_response.md)
- [接触境界条件のオンライン更新システム](topics/contact_boundary/contact_boundary_da_online_update_system.md)

