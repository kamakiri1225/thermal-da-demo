# 索引 ― 「これが知りたい」から探す

**本編は `blog_001`〜`blog_005` の5本です。** 優位性・先行研究の確認範囲・追加検証は [研究レビュー](19_research_advantage_and_validation.md) を参照してください。 下の逆引き表から目的の節へ飛んでください。

Webで読む（数式・図・GIFつき）: <https://kamakiri1225.github.io/thermal-da-demo/>

---

## 逆引き ― 疑問からファイルへ

### データ同化の手続きが知りたい

| 知りたいこと | どこ |
|---|---|
| **この研究で実際に行った同化の手続き**（真値・観測・初期値・各サイクルの計算・Qが同定される仕組み） | **blog_004 §7-6 冒頭** |
| **発熱量 $Q$ をどうやって同定するのか**（手順・数式・実数値で1サイクル追う） | **blog_003 §4-0** |
| **変位計の値をどうやって同化に入れるのか**（状態に入れない理由・4ステップ・コード） | **blog_004 §7-0** |
| 観測データをどう作ったのか（双子実験・ノイズの載せ方） | blog_004 §7-0 の図 |
| \(D\) 行列とは何か・なぜFEMを呼ばずに済むのか | blog_004 §7-0・§7-2 |
| 放熱 $h$ はなぜ当たらないのか／どうすれば当たるのか | blog_003 §4-0-b |
| 温度しか測っていないのになぜ $Q$ が動くのか | blog_003 §4-0（ステップ4） |
| OI と EnKF は何が違うのか | blog_002 §3-6・blog_003 §1〜3 |
| EnKF の更新式を手計算で追いたい | blog_003 §3 |
| 粒子フィルタではだめなのか | blog_003 §5-2 |

### 低次元化（POD / ROM）が知りたい

| 知りたいこと | どこ |
|---|---|
| **Q-DEIM と gappy-POD は何が違うのか** | **blog_004 §4-0** |
| POD は何をしているのか（SVDの中身） | blog_004 §2 |
| 小さい行列で POD を手計算したい | blog_004 §3 |
| 代表点はどう選ばれるのか（貪欲法・枢軸付きQR） | blog_004 §4 |
| 5点から全温度場をどう復元するのか | blog_004 §5-4 |
| ROM の係数 $C,K,h$ はどう決めたのか | blog_004 §6-2 |
| 測っていないセルでも本当に合うのか | blog_004 §5-4（末尾） |

### センサをどこに置くかが知りたい

| 知りたいこと | どこ |
|---|---|
| **温度センサと変位センサで基準が違う**（2つの感度の使い分け） | **blog_005 §4-2b** |
| 観測の価値 $\Delta$ とは何か | blog_005 §1 |
| 推定対象が変われば最適観測も変わるのか | blog_005 §4-2c |
| 熱感度 $W=K_s^{-1}H_T$ はどこから来るのか | blog_005 §4-2b・blog_004 §7-2 |
| 観測を増やすとどれだけ良くなるのか | blog_004 §7-6 |
| **知りたい場所を測らずに当てられるのか**（循環の排除） | **blog_004 §7-4b** |
| **温度センサ0本、変位計だけで温度場を直せるのか** | **blog_004 §7-4c** |
| 傾き（反り）は工作機械で補正できるのか | blog_004 §7-4b 末尾 |

### 解析そのものが知りたい

| 知りたいこと | どこ |
|---|---|
| OpenFOAM の CHT で何を解いたのか | blog_001 §1〜2 |
| CFD の温度を FEM の節点へどう渡すのか（IDW補間） | blog_001 §2-3 |
| 熱伝達率はどう求めたのか | blog_001 §3 |
| 温度から変位はどう計算するのか（熱弾性FEM） | blog_004 §7-2 |
| $x,y,z$ の全成分はどうなっているのか | blog_004 §7-5 |
| なぜ工作機械で熱変位が問題なのか（文献つき） | blog_001 §0 |

---

## 本編5本

| # | ファイル | 内容 |
|---|---|---|
| 1 | [`blog_001_cht_thermal_expansion_and_htc.md`](blog_001_cht_thermal_expansion_and_htc.md) | OpenFOAM(CHT) → FrontISTR(熱膨張)、熱伝達率。**双子実験の真値の作り方** |
| 2 | [`blog_002_oi_data_assimilation.md`](blog_002_oi_data_assimilation.md) | OI（最適内挿）。更新式1本と3点モデルの手計算 |
| 3 | [`blog_003_ensemble_kalman_filter.md`](blog_003_ensemble_kalman_filter.md) | EnKF。**パラメータ同定の手順（§4-0）** |
| 4 | [`blog_004_pod_selected_rom.md`](blog_004_pod_selected_rom.md) | **POD／Q-DEIM／ROM の本編。変位の観測化（§7-0）** |
| 5 | [`blog_005_optimal_sensor_placement.md`](blog_005_optimal_sensor_placement.md) | 観測の価値とセンサ配置 |

## 補助資料（本編にない固有の内容だけ）

| ファイル | 固有の内容 |
|---|---|
| [`18_oi_sigma_t_validation.md`](18_oi_sigma_t_validation.md) | OI の $\sigma_T$ を振った検証・熱伝導の誤差伝播試算 |
| [`17_holdout_displacement_results.md`](17_holdout_displacement_results.md) | 未観測C/D変位の作業記録（2026-09-14）。**更新版は blog_004 §7-5** |
| [`16_unobserved_displacement_support.md`](16_unobserved_displacement_support.md) | 上の支援メモ |
| [`03_presentation_story.md`](03_presentation_story.md) | 発表ストーリーの下書き。**最新は `presentation/`** |
| [`research_handoff_prompt.md`](research_handoff_prompt.md) | 引き継ぎ指示書 |

## 発表資料は `../presentation/`

| ファイル | 内容 |
|---|---|
| `01_slides.html` | **発表スライド28枚**（reveal.js） |
| `02_positioning_novelty.md` | **研究の位置づけ・先行研究との差分・主張の設計** |
| `03_expected_questions.md` | 想定質問 Q1〜Q22 |
| `04_intro_literature.md` | 導入の文献調査・メーカー比較 |
| `06_conference_posters.html` | 学会別ポスター（4学会分） |

---

## 統合の記録（2026-09-30）

重複が大きかった3本を blog へ統合して削除しました（git 履歴には残っています）。

| 削除 | 行数 | 移行先 |
|---|---|---|
| `02_full_story.md` | 1,291 | `blog_004`。**§7-6 として RMSE の定義を吸収** |
| `00_data_assimilation_algorithm.md` | 1,123 | `blog_002`（OI）・`blog_003`（EnKF） |
| `01_pod_qdeim_algorithm.md` | 330 | `blog_004` §2〜4 |

復元手順:

```bash
git log --diff-filter=D -- sample/106_0_pod_selected_rom/docs/02_full_story.md
git show <commit>^:sample/106_0_pod_selected_rom/docs/02_full_story.md > 復元先.md
```
