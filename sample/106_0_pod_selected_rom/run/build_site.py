"""blog_001〜012（docs/blog_0NN_*.md）を、ブラウザで読む HTML にする（GitHub Pages の gh-pages 用）。

  出力：site/blog0NN.html と site/img/（記事が使う図だけを docs/img から写す）
  - 数式：MathJax（CDN）、目次つき、本文の幅 1120px。原稿の数式は GitHub 用の書き方（```math、$`…`$、run/gfm_math.py）なので、pandoc 用に戻してから変換する
  - ```mermaid のブロックは mermaid.js で図にする
  - 記事どうしのリンク blog_0NN_*.md → blog0NN.html、そのほかの .md やプログラムへのリンクは GitHub のページへ
  再現：python3 run/build_site.py
"""
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / 'docs'; SITE = ROOT / 'site'; (SITE / 'img').mkdir(parents=True, exist_ok=True)
GH = 'https://github.com/kamakiri1225/thermal-da-demo/blob/main/sample/106_0_pod_selected_rom/docs/'
HEAD = """<style>
   html body { max-width: 1120px !important; }
   body { padding: 0 24px; }
   img { max-width: 100%; height: auto; }
   table { font-size: 0.95em; }
   pre.mermaid { background: none; text-align: center; }
  </style>
<script type="module">
  import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10.9.1/dist/mermaid.esm.min.mjs';
  mermaid.initialize({ startOnLoad: true });
</script>
"""


def to_pandoc(text):
    """GitHub 用の数式の書き方（```math、$`…`$）を、pandoc が読める $$ … $$、$ … $ に戻す（run/gfm_math.py の逆）。"""
    text = re.sub(r'^([ \t>]*)```math[ \t]*\n(.*?)\n\1```[ \t]*$', lambda m: f'{m.group(1)}$$\n{m.group(2)}\n{m.group(1)}$$', text, flags=re.S | re.M)
    return re.sub(r'\$`(.+?)`\$', r'$\1$', text)


def fix_link(m):
    href = m.group(1)
    if href.startswith(('http', '#', 'mailto:')):
        return m.group(0)
    path, _, anchor = href.partition('#')
    b = re.match(r'blog_(\d{3})_.*\.md$', Path(path).name)
    if b and '/' not in path:
        new = f'blog{b.group(1)}.html'
    elif path.startswith('img/'):
        return m.group(0)
    else:
        new = GH + path
    return f'href="{new}' + (f'#{anchor}' if anchor else '') + '"'


def main():
    head = SITE / '_head.html'; head.write_text(HEAD, encoding='utf-8')
    used = set()
    for md in sorted(DOCS.glob('blog_0[0-9][0-9]_*.md')):
        num = md.name[5:8]
        text = md.read_text(encoding='utf-8')
        title = re.search(r'^# (.+)$', text, re.M).group(1).strip()
        used |= set(re.findall(r'\]\((img/[^)\s]+)\)', text)) | set(re.findall(r'src="(img/[^"]+)"', text))
        out = SITE / f'blog{num}.html'
        subprocess.run(['pandoc', '-', '-f', 'markdown+autolink_bare_uris', '--standalone', '--mathjax=https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js',
                        '--toc', '--metadata', f'pagetitle={title}', '--metadata', 'lang=ja',
                        f'--include-in-header={head}', '-o', str(out)], check=True, cwd=DOCS, input=to_pandoc(text), text=True)
        h = out.read_text(encoding='utf-8')
        h = re.sub(r'<pre class="mermaid"><code>(.*?)</code></pre>', r'<pre class="mermaid">\1</pre>', h, flags=re.S)
        h = re.sub(r'href="([^"]+)"', fix_link, h)
        out.write_text(h, encoding='utf-8')
        print('wrote', out.name, title)
    head.unlink()
    missing = []
    for rel in sorted(used):
        src = DOCS / rel
        if src.exists():
            shutil.copy2(src, SITE / rel)
        else:
            missing.append(rel)
    print('images', len(used), 'missing', missing)


if __name__ == '__main__':
    main()
