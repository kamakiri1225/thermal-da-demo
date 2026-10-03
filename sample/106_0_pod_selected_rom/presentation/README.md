# 発表資料（オープンCAE学会シンポジウム B-16・2026年11月12日）

番号は**読む順**です。

| # | ファイル | 内容 | いつ使うか |
|---|---|---|---|
| **00** | [`00_summary.md`](00_summary.md) | **発表内容の整理（一言・7章の流れ・主な数値・主張できること/できないこと・未決定事項）** | **まずこれを読む** |
| **01** | [`01_slides.html`](01_slides.html) | 発表スライド（旧構成・35枚。背景→DA基礎→課題→ROM→結果→観測点設計の順） | 旧版（07 と内容が重なる） |
| **02** | [`02_positioning_novelty.md`](02_positioning_novelty.md) | **研究の位置づけ・先行研究との差分・何を主張すれば優位性が出るか・追加検証の結果** | 発表前の準備／質疑 |
| **03** | [`03_expected_questions.md`](03_expected_questions.md) | **想定質問 Q1〜Q22**（LETKF、ベイズ最適化、メーカー比較など） | 質疑 |
| **04** | [`04_intro_literature.md`](04_intro_literature.md) | 導入の文献調査（75%の出典）・メーカー4社比較・ISO 230-3・測定方法 | スライド2の根拠 |
| **05** | [`05_overview.md`](05_overview.md) | 発表内容6項目の概要 | 全体像の確認 |
| **06** | [`06_conference_posters.html`](06_conference_posters.html) / [`.pdf`](06_conference_posters.pdf) | 学会別ポスター4学会分（OpenCAE／計算工学／計算力学／実務応用） | 別学会への打診 |
| **07** | [`07_slides_displacement.html`](07_slides_displacement.html) / [`.pdf`](07_slides_displacement.pdf) | **発表スライド（新構成・25枚）**「測れない場所の熱変形を、測れる場所から当てる」。切り口A／B／C で組み立て直したもの。全スライドに**話す原稿**（reveal.js の **S** キーでノート表示） | **本番（新）** |

公開版スライド: <https://kamakiri1225.github.io/thermal-da-demo/opencae2026.html>

## 本編（技術内容）は `../docs/`

解説は blog_001〜005 にあります。**どこに何が書いてあるかは
[`../docs/README.md`](../docs/README.md) の逆引き表**から引いてください。

## 改名の記録（2026-09-30）

| 旧 | 新 |
|---|---|
| `opencae_2026_slides.html` | `01_slides.html` |
| `opencae_2026_positioning.md` | `02_positioning_novelty.md` |
| `opencae_2026_qa.md` | `03_expected_questions.md` |
| `opencae_2026_slide02_intro.md` | `04_intro_literature.md` |
| `opencae_2026_slide00_overview.md` | `05_overview.md` |
| `conference_posters.*` | `06_conference_posters.*` |

`git mv` なので履歴は `git log --follow <新ファイル名>` で追えます。

## 今回の研究レビュー（2026-10-01）

優位性・調査範囲・訂正点・追加試験は [研究レビュー](../docs/19_research_advantage_and_validation.md)。ブログ005と本番スライドにも要点を反映しました。学会別ポスターPDFは4枚のA2横向きです。

ポスターPDFはJavaScriptを使わず、TeXをMathJaxのSVGへ変換して組版できます。

```bash
python3 ../run/render_conference_posters_pdf.py --mathjax-root /path/to/node_modules/mathjax-full
```

Node.js、mathjax-full、beautifulsoup4、WeasyPrint、PyMuPDFが必要です。PDFではGIFは静止表示、HTMLでは動画として表示します。

## 07_slides_displacement.html について（2026-10-03 追加）

01 は「作ったもの」を順に説明する構成でした。07 は**主張**で組み直したものです。

| | 01（旧） | 07（新） |
|---|---|---|
| 軸 | 手法の説明順（DA 基礎 → ROM → 結果 → 観測設計） | **切り口A／B／C**（変位は情報源／測れない場所を当てる／条件が変わっても使える） |
| 目玉 | 温度場 0.16 K | **A・O を一度も測らずに反り 0.075 µm**（直接測るより正確。実ソルバ真値で3条件） |
| 先行研究 | 言及なし | **4つの流派と「空いている場所」**を1枚で |
| 配置設計 | 「悪い例」との比較のみ | **既存の配置法（A最適・D最適・Q-DEIM・素朴・ランダム）と同じ土俵で比較**（2026-10-03 追加） |
| 限界 | 1枚 | 1枚（**言えないこと**を7項目で明示） |
| 原稿 | なし | **全25枚に話す原稿**（S キー） |
| 1枚の密度 | 条件・表・カード・結論を全部載せる | **見出し（結論文）＋大きい図1つ＋結論1行**。細部は原稿とブログへ |

- **20分発表の想定で25枚**（1枚あたり約48秒）。検証ごとに「観測点 → 温度の時刻歴 → 変位の時刻歴 → 誤差」の順で見せます。会場のスクリーンで読めるよう、
  スライド用の図を別に作って使っています。
  - `docs/img/slide_*.png` … 誤差のまとめ（`run/make_slide_figs.py`。計算はせず results/*.json から描き直す）
  - `docs/img/ver{1..4}_detail.png` / `ver{1..4}_slide.png` … 検証ごとの**観測点・温度の時刻歴・変位の時刻歴**
    （`run/verification_timeseries.py`。EnKF を回し直して時刻歴を保存する）
- 図は `../docs/img/` を参照しています。フォルダごと移動する場合は相対パスに注意してください。
- reveal.js と MathJax は CDN から読み込みます。**会場でネットが使えない場合に備えて PDF 版**（`07_slides_displacement.pdf`、25ページ、16:9）を同梱しています。
- PDF は各スライドを 2560×1440 で描画して束ねたものです。アニメーション GIF は1コマ目で止まります。

公開版: <https://kamakiri1225.github.io/thermal-da-demo/opencae2026_b.html>
