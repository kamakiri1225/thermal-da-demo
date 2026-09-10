"""slides.json から reveal.js プレゼン(chapters=横スクロール, 章内=縦スクロール)を生成する。
数式は <eq>TeX</eq> を MathJax で表示。画像は presentation/ からの相対パス参照。
Usage: python3 presentation/build_reveal.py  ->  conference_reveal.html
"""
import html, json, re, subprocess
from pathlib import Path
HERE = Path(__file__).absolute().parent

def render_body(md):
    # <eq>TeX</eq> -> $$TeX$$（pandoc --mathjax が \[..\] を出し、MathJaxが描画）
    md = re.sub(r"<eq>(.*?)</eq>", lambda m: "\n\n$$"+m.group(1)+"$$\n\n", md, flags=re.S)
    return subprocess.check_output(
        ["pandoc","-f","markdown","-t","html5","--mathjax"], input=md, text=True)

def slide_html(s):
    eyebrow = f'<div class="eyebrow">{html.escape(s.get("part",""))}</div>' if s.get("part") else ""
    fig = ""
    if s.get("image"):
        fig = (f'<figure><img src="{html.escape(s["image"])}">'
               f'<figcaption>{html.escape(s.get("caption",""))}</figcaption></figure>')
    body = render_body(s.get("body",""))
    take = (f'<div class="takeaway">{html.escape(s["takeaway"])}</div>'
            if s.get("takeaway") else "")
    cls = "cover" if s.get("cover") else ("visual" if fig else "text")
    return (f'<section class="{cls}">{eyebrow}<h2>{html.escape(s["title"])}</h2>'
            f'<div class="content">{fig}<div class="body">{body}</div></div>{take}</section>')

def main():
    slides = json.loads((HERE/"slides.json").read_text())
    # 章にまとめる: 表紙 or 「第N部/補足資料」の cover を章の先頭にする
    chapters, cur = [], []
    for s in slides:
        is_div = s.get("cover") and (s["title"].startswith("第") or s["title"].startswith("補足"))
        is_title = s.get("cover") and not is_div
        if (is_div or is_title) and cur:
            chapters.append(cur); cur = []
        cur.append(s)
    if cur: chapters.append(cur)

    horiz = []
    for ch in chapters:
        inner = "".join(slide_html(s) for s in ch)
        horiz.append(f"<section>{inner}</section>")  # 外=横, 内=縦
    slides_html = "".join(horiz)

    css = """
    .reveal{font-family:'Noto Sans CJK JP',sans-serif;color:#172f43}
    .reveal h2{font-size:1.25em;color:#0c4a52;text-align:left;margin:0 0 .4em;text-transform:none}
    .reveal .eyebrow{font-size:.42em;font-weight:700;letter-spacing:1px;color:#007f86;text-align:left}
    .reveal .content{display:flex;gap:26px;align-items:flex-start;text-align:left}
    .reveal .body{font-size:.62em;line-height:1.5;flex:1;min-width:0}
    .reveal .visual .body{font-size:.55em}
    .reveal figure{margin:0;width:52%;flex:none}
    .reveal figure img{width:100%;max-height:60vh;object-fit:contain;background:#fff}
    .reveal figcaption{font-size:.38em;color:#53677a}
    .reveal .takeaway{margin-top:18px;padding:10px 16px;background:#e9f5f3;border-left:6px solid #007f86;
       font-weight:700;font-size:.6em;text-align:left;line-height:1.4}
    .reveal strong{color:#007f86}
    .reveal table{font-size:.55em;border-collapse:collapse}
    .reveal th,.reveal td{border-bottom:1px solid #c7d4df;padding:6px 10px}
    .reveal th{background:#edf4f6;color:#007f86}
    .reveal pre{font-size:.5em}
    .reveal .cover h2{font-size:1.5em;color:#0c4a52;text-align:center;margin-top:.5em}
    .reveal .cover .body{font-size:.8em;text-align:center}
    .reveal .cover .content{display:block}
    /* 数式と段落を詰めて見切れを防ぐ */
    .reveal .body p{margin:.35em 0}
    .reveal mjx-container{font-size:.82em !important;margin:.25em 0 !important}
    .reveal .body>*:first-child{margin-top:0}
    /* 万一あふれたら縦スクロールで全文読めるようにする */
    .reveal .slides section{max-height:92vh;overflow-y:auto}
    .reveal .slides section::-webkit-scrollbar{width:8px}
    .reveal .slides section::-webkit-scrollbar-thumb{background:#bcd; border-radius:4px}
    """
    doc = f"""<!doctype html><html lang="ja"><head><meta charset="utf-8">
<title>温度・変位観測を用いた熱状態推定 — reveal.js</title>
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/reveal.js/4.6.1/reveal.min.css">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/reveal.js/4.6.1/theme/white.min.css">
<style>{css}</style></head><body>
<div class="reveal"><div class="slides">{slides_html}</div></div>
<script src="https://cdnjs.cloudflare.com/ajax/libs/reveal.js/4.6.1/reveal.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/reveal.js/4.6.1/plugin/math/math.min.js"></script>
<script>
Reveal.initialize({{hash:true, controls:true, progress:true, slideNumber:'c/t',
  width:1280, height:720, margin:0.06,
  plugins:[RevealMath.MathJax3],
  math:{{mathjax:'https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-mml-chtml.js',
         tex:{{inlineMath:[['\\\\(','\\\\)']],displayMath:[['\\\\[','\\\\]']]}}}}}});
</script></body></html>"""
    (HERE/"conference_reveal.html").write_text(doc)
    print(f"reveal.js: {len(slides)}枚 / {len(chapters)}章 -> conference_reveal.html")
    print("章ごと=横(←→), 章内=縦(↑↓)。数式はMathJax(TeX)。")

if __name__=="__main__":
    main()
