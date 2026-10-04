"""ブログ・ポスターをローカルの静的サイトとして組み立て、自分のPCだけで閲覧する.

GitHub Pages の「ログイン限定公開」は Enterprise Cloud 限定のため、
リポジトリを非公開にすると .io は使えなくなる。その代替。

使い方:
  python3 run/build_local_site.py            # site/ に生成するだけ
  python3 run/build_local_site.py --serve    # 生成して http://localhost:8765 で配信
"""
from __future__ import annotations
import os, re, sys, shutil, subprocess, http.server, socketserver, functools
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DOCS = ROOT / "docs"
OUT = ROOT / "site"
MJ = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"
PORT = 8765

TITLES = {
    "blog_001": "blog_001：OpenFOAM＋FrontISTR で熱流体固体連成→熱膨張",
    "blog_002": "blog_002：OI（最適内挿）でデータ同化",
    "blog_003": "blog_003：アンサンブルカルマンフィルタ（EnKF）",
    "blog_004": "blog_004：POD で代表点を選ぶ ROM",
    "blog_005": "blog_005：推定対象で最適センサは変わる",
    "blog_006": "blog_006：温度2点＋変位2点の同化を1ステップずつ追う",
    "blog_007": "blog_007：作った ROM はどこまで使えるか",
    "blog_008": "blog_008：温度だけでなく変位も同化に使う ― 本研究の主張",
}
LEAD = {
    "blog_001": "中空円筒のCHT解析からFrontISTRの熱膨張まで。熱伝達率の分布も出す。",
    "blog_002": "更新式1本と3点モデルの手計算から、熱感度Wによる観測点選びまで。",
    "blog_003": "分身の散らばりで共分散を作り、発熱量Qと放熱hを同時推定する。",
    "blog_004": "PODで場の型を出し、Q-DEIMで代表5点を選んで軽いROMを組む。",
    "blog_005": "観測の価値Δを1本の式で定義し、対象ごとに最適観測が変わることを示す。",
    "blog_006": "センサの読みから、状態・予報・観測演算子・共分散・ゲイン・更新まで、本番の数値で追う。",
    "blog_007": "発熱量・加熱のしかたが変わっても使えるか。使える条件と作り直しが要る条件を実計算で確かめる。",
    "blog_008": "先行研究を調べたうえで、変位を観測に入れることが研究の主張として成り立つかを整理する。",
}

CSS = """
<style>
 html body { max-width: 1120px !important; }
 body { padding: 0 24px; }
 img { max-width: 100%; height: auto; }
 table { font-size: 0.95em; }
</style>
"""

INDEX_CSS = """
body{font-family:"Noto Sans CJK JP",system-ui,sans-serif;max-width:980px;margin:0 auto;
     padding:32px 24px 64px;color:#1b2430;line-height:1.7}
h1{font-size:30px;border-bottom:3px solid #2E75D4;padding-bottom:10px}
.note{background:#f1f5fa;border-left:4px solid #2E75D4;padding:12px 16px;font-size:14.5px}
.card{display:block;border:1px solid #c9d4e0;border-radius:10px;padding:16px 18px;margin:12px 0;
      text-decoration:none;color:inherit;transition:.15s}
.card:hover{border-color:#2E75D4;box-shadow:0 2px 10px rgba(46,117,212,.18)}
.card b{font-size:17px;color:#14459c}
.card p{margin:6px 0 0;font-size:14.5px;color:#445}
.sub{font-size:13px;color:#667}
</style>
"""


REPO_BLOB = "https://github.com/kamakiri1225/thermal-da-demo/blob/main/"
DOCS_IN_REPO = "sample/106_0_pod_selected_rom/docs/"


def fix_links(html: Path) -> None:
    """Markdown 用の相対リンクを、公開ページで切れないリンクに書き換える.

    - blog_00N_*.md(#…)  → blog00N.html(#…)（サイト内の記事）
    - それ以外の相対パス（../run/*.py、21_*.md など）→ GitHub 上のファイル
    """
    import posixpath
    t = html.read_text(encoding="utf-8")

    def repl(m):
        href = m.group(1)
        if re.match(r"^(https?:|mailto:|#|img/|pdf/)", href):
            return m.group(0)
        path, _, frag = href.partition("#")
        mb = re.match(r"^blog_(00[1-8])_[^/]*\.md$", path)
        if mb:
            return f'href="blog{mb.group(1).replace("00", "00", 1)}.html' + (f"#{frag}" if frag else "") + '"'
        repo_path = posixpath.normpath(posixpath.join(DOCS_IN_REPO, path))
        return f'href="{REPO_BLOB}{repo_path}' + (f"#{frag}" if frag else "") + '"'

    t = re.sub(r'href="([^"]+)"', repl, t)
    html.write_text(t, encoding="utf-8")


def build():
    # 他の作業で追加されたHTMLや素材を消さず、生成対象だけ更新する。
    (OUT / "img").mkdir(parents=True, exist_ok=True)
    (OUT / "pdf").mkdir(parents=True, exist_ok=True)

    hdr = OUT / "_header.html"
    hdr.write_text(CSS, encoding="utf-8")

    mds = sorted(DOCS.glob("blog_00[1-8]_*.md"))
    used_imgs: set[str] = set()
    for md in mds:
        stem = md.name[:8]
        html = OUT / f"{stem.replace('_','')}.html"
        subprocess.run(["pandoc", str(md), "-f", "markdown+gfm_auto_identifiers", "--standalone", f"--mathjax={MJ}", "--toc",
                        "--include-in-header", str(hdr),
                        "--metadata", f"title={TITLES.get(stem, stem)}",
                        "-o", str(html)], check=True)
        used_imgs |= set(re.findall(r'img/[A-Za-z0-9_\-]+\.(?:png|gif|jpg)', md.read_text(encoding="utf-8")))
        fix_links(html)
        print("  built", html.name)

    # Markdown向けのブログリンクを、生成したHTML向けに読み替える。
    from bs4 import BeautifulSoup
    page_by_md = {p.name: f"{p.name[:8].replace('_', '')}.html" for p in mds}
    for page in sorted(OUT.glob("blog00[1-8].html")):
        soup = BeautifulSoup(page.read_text(encoding="utf-8"), "html.parser")
        for a in soup.find_all("a", href=True):
            href = a["href"]
            if href.startswith(("https:", "http:", "mailto:", "#")):
                continue
            dest, sep, fragment = href.partition("#")
            if Path(dest).name in page_by_md:
                target = page_by_md[Path(dest).name]
                if fragment:
                    target_soup = BeautifulSoup((OUT / target).read_text(encoding="utf-8"), "html.parser")
                    if not target_soup.find(id=fragment):
                        short = re.sub(r"^\d+(?:-\d+)*-", "", fragment)
                        if target_soup.find(id=short):
                            fragment = short
                a["href"] = target + ("#" + fragment if sep else "")
            elif dest.endswith(".md") and (DOCS / dest).exists():
                a["href"] = "../docs/" + dest + ("#" + fragment if sep else "")
        page.write_text(str(soup), encoding="utf-8")

    for rel in sorted(used_imgs):
        src = DOCS / rel
        if src.exists():
            shutil.copy2(src, OUT / rel)
    print(f"  copied {len(used_imgs)} images")

    for pdf in sorted((DOCS / "pdf").glob("blog_00[1-8].pdf")):
        shutil.copy2(pdf, OUT / "pdf" / pdf.name)
    poster = ROOT / "presentation" / "06_conference_posters.pdf"
    if poster.exists():
        shutil.copy2(poster, OUT / "pdf" / "posters.pdf")

    # ポスターHTML（画像パスを site 内へ書き換え）
    ph = ROOT / "presentation" / "06_conference_posters.html"
    if ph.exists():
        t = ph.read_text(encoding="utf-8").replace('src="../docs/img/', 'src="img/')
        (OUT / "posters.html").write_text(t, encoding="utf-8")
        for rel in set(re.findall(r'img/[A-Za-z0-9_\-]+\.(?:png|gif|jpg)', t)):
            src = DOCS / rel
            if src.exists():
                shutil.copy2(src, OUT / rel)

    stale_pdfs = [md for md in mds if (DOCS / "pdf" / (md.name[:8] + ".pdf")).exists()
                  and (DOCS / "pdf" / (md.name[:8] + ".pdf")).stat().st_mtime < md.stat().st_mtime]
    pdf_note = ('<p class="note">PDFには更新前の記事が含まれます。最新の説明は下のHTML記事を読んでください。</p>'
                if stale_pdfs else "")

    cards = "".join(
        f'<a class="card" href="{s.replace("_","")}.html"><b>{TITLES[s]}</b>'
        f'<p>{LEAD[s]}</p>'
        + (f'<p class="sub">PDF: <a href="pdf/{s}.pdf">{s}.pdf</a></p>' if (DOCS / "pdf" / f"{s}.pdf").exists() else '')
        + '</a>\n'
        for s in ["blog_001", "blog_002", "blog_003", "blog_004", "blog_005", "blog_006", "blog_007", "blog_008"])
    (OUT / "index.html").write_text(
        f'<!DOCTYPE html><html lang="ja"><head><meta charset="utf-8">'
        f'<title>熱データ同化 研究ノート（ローカル閲覧）</title>{INDEX_CSS}</head><body>'
        f'<h1>熱データ同化 研究ノート</h1>'
        f'<p class="note">このページは<b>あなたのPCの中だけ</b>で動いています'
        f'（<code>{OUT}</code>）。インターネットには公開されていません。</p>'
        f'{pdf_note}{cards}'
        f'<a class="card" href="posters.html"><b>学会発表ポスター（4件）</b>'
        f'<p>オープンCAE／計算工学／計算力学／精密工学。</p>'
        f'<p class="sub">PDF: <a href="pdf/posters.pdf">posters.pdf</a></p></a>'
        f'</body></html>', encoding="utf-8")
    hdr.unlink()
    print(f"\n✅ {OUT} に生成しました")


def serve():
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(OUT))
    socketserver.TCPServer.allow_reuse_address = True
    with socketserver.TCPServer(("127.0.0.1", PORT), handler) as httpd:
        print(f"\n🌐 http://localhost:{PORT}/ で閲覧できます（Ctrl+C で停止）")
        print("   ※ 127.0.0.1 のみ待ち受け。外部からは見えません。")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n停止しました")


if __name__ == "__main__":
    build()
    if "--serve" in sys.argv:
        serve()
