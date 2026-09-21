# 元記事から作るTeX版PDF

`../blog_001_*.md`〜`../blog_005_*.md`を本文の原本とし、PandocでLaTeXへ変換、LuaLaTeXでPDFへ組版します。本文の要約・節の削除は行いません。

- 日本語本文：IPAex明朝、見出し・強調：IPAexゴシック
- 欧文：TeX Gyre Termes（Times系）
- 数式：TeX Gyre Termes Math
- **A4・二段組（論文形式）**、10pt、余白18mm、段間6.5mm。
- 表は `tables_twocol.lua` で booktabs の tabular に変換し、段幅を超える場合のみ縮小
  （二段組と相性の悪い longtable は使いません）。
- 別行数式は adjustbox で段幅以内に自動フィット。長いパス風リンクは `/` で改行可。
- GIFは中央フレームを静止画として掲載します。動画は元記事を参照してください。
- Mermaidはノードと矢印の接続を保ってGraphvizのベクトル図に変換します。

## 再生成

```bash
python3 build_pdf.py
# 1本だけ
python3 build_pdf.py blog_005
```

必要環境：Pandoc、LuaLaTeX、luatexja、IPAexフォント、TeX Gyreフォント、Graphviz、PythonのPillow。

`blog_001.tex`〜`blog_005.tex`が編集可能なTeXソース、`tex_assets/`が掲載画像とフロー図です。TeXのみから再組版する場合は、このディレクトリで `lualatex blog_001.tex` を2回実行します。`paper_header.tex`と`build_pdf.py`は再生成時の設定です。TeXの直接修正はMarkdownから再生成すると置き換わります。

`build_manifest.json`に原本のSHA-256と数式・図・表の変換件数を記録します。数式数はインラインと別行数式の合計です。
