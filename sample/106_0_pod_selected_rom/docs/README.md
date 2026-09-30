# 106 のドキュメント ― まず blog_001〜005 を読んでください

**本編は blog_001〜005 の5本です。** それ以外は補助資料で、内容は本編と重複します。

## 本編（ここから読む）

| # | ファイル | 内容 |
|---|---|---|
| 1 | [`blog_001_cht_thermal_expansion_and_htc.md`](blog_001_cht_thermal_expansion_and_htc.md) | OpenFOAM(CHT) → FrontISTR(熱膨張)、面ごとの熱伝達率。**双子実験の真値をどう作ったか** |
| 2 | [`blog_002_oi_data_assimilation.md`](blog_002_oi_data_assimilation.md) | OI（最適内挿）。更新式1本と3点モデルの手計算から、発熱量 $Q$ の逆推定まで |
| 3 | [`blog_003_ensemble_kalman_filter.md`](blog_003_ensemble_kalman_filter.md) | EnKF。分身の散らばりで共分散を作る。**パラメータ同定の手順（§4-0）と $h$ の可同定性（§4-0-b）** |
| 4 | [`blog_004_pod_selected_rom.md`](blog_004_pod_selected_rom.md) | **POD／Q-DEIM／ROM の本編**。$x,y,z$ 全成分の同化（§7-5）、観測構成の比較（§7-6） |
| 5 | [`blog_005_optimal_sensor_placement.md`](blog_005_optimal_sensor_placement.md) | 観測の価値 $\Delta$、センサ配置。**2つの感度の使い分け（§4-2b）** |

Webで読む（数式・図・GIF入り）: <https://kamakiri1225.github.io/thermal-da-demo/>

## 補助資料（本編にない固有の内容だけ）

| ファイル | 固有の内容 |
|---|---|
| [`18_oi_sigma_t_validation.md`](18_oi_sigma_t_validation.md) | OI の $\sigma_T=1.5$ K を振った検証。熱伝導による誤差伝播の試算 |
| [`17_holdout_displacement_results.md`](17_holdout_displacement_results.md) | 未観測C/D変位の検証（2026-09-14の作業記録）。**更新版は blog_004 §7-5** |
| [`16_unobserved_displacement_support.md`](16_unobserved_displacement_support.md) | 上の支援メモ |
| [`03_presentation_story.md`](03_presentation_story.md) | 発表ストーリーの下書き。**最新は `presentation/` 配下** |
| [`research_handoff_prompt.md`](research_handoff_prompt.md) | 別担当者・別セッションへの引き継ぎ指示書 |
| `slides_01_pod_qdeim.html` / `slides_02_full_story.html` | reveal.js スライド（社内向け） |

## 発表資料は `presentation/` に

| ファイル | 内容 |
|---|---|
| `opencae_2026_slides.html` | **発表スライド25枚**（reveal.js） |
| `opencae_2026_positioning.md` | 研究の位置づけ・先行研究との差分・主張の設計 |
| `opencae_2026_qa.md` | 想定質問 Q1〜Q22 |
| `opencae_2026_slide02_intro.md` | 導入の文献調査・メーカー比較 |
| `conference_posters.html` | 学会別ポスター（4学会分） |

---

## 統合の記録（2026-09-30）

重複が大きかった次の3本を **blog へ統合して削除**しました（git 履歴には残っています）。

| 削除したファイル | 行数 | 移行先 |
|---|---|---|
| `02_full_story.md` | 1,291 | `blog_004`（全体）。**§7-6「観測構成を変えて比べる」として RMSE の定義を吸収** |
| `00_data_assimilation_algorithm.md` | 1,123 | `blog_002`（OI）・`blog_003`（EnKF） |
| `01_pod_qdeim_algorithm.md` | 330 | `blog_004` §2〜4 |

吸収した固有内容：
- **平均RMSEの厳密な定義**（5 seed × 60メンバー × 全5点 × 10時刻の平均の取り方）→ blog_004 §7-6
- **OSSE**（Observing System Simulation Experiment＝双子実験の正式名）→ blog_001

復元するには `git log --diff-filter=D -- docs/02_full_story.md` で削除コミットを探し、
`git show <commit>^:sample/106_0_pod_selected_rom/docs/02_full_story.md` で取り出せます。
