"""理論章から数式SVG埋め込みのHTML/PDFとMarkdownを生成する。
MathJaxに接続しなくても読み物版の数式を表示できる。
"""
from pathlib import Path
import base64,html,json,re,subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from weasyprint import HTML
HERE=Path(__file__).absolute().parent
ASSETS=HERE/'assets/theory'

def main():
 ASSETS.mkdir(exist_ok=True)
 s=json.loads((HERE/'theory_slides.json').read_text())
 articles=[];md=['# OI・KF・EnKF・PF：理論式の展開と104への接続\n\n共通の推定問題から導出し、具体例・実装へつなげる読み物版。過程ノイズ共分散は、ヒーター発熱量Qと区別してΩと表記する。\n']
 n=0
 for i,x in enumerate(s,1):
  def mathsvg(m):
   nonlocal n
   n+=1;raw=m.group(1)
   tex=raw.replace(r'\tfrac',r'\frac').replace(r'\mathbf1',r'\mathbf{1}')
   # MathJax/TeXは短縮分数を許すが、MathTextは波括弧が必要。
   tex=re.sub(r'\\frac1([A-Za-z])',r'\\frac{1}{\1}',tex)
   tex=tex.replace(r'\frac1{',r'\frac{1}{')
   tex=re.sub(r'\\le(?![A-Za-z])', r'\\leq', tex)
   tex=re.sub(r'\\ge(?![A-Za-z])', r'\\geq', tex)
   tex=re.sub(r'\\frac([0-9])([0-9A-Za-z])',r'\\frac{\1}{\2}',tex)
   f=ASSETS/f'equation_{n:03d}.svg'
   fig=plt.figure(figsize=(1,1));fig.text(0,0,'$'+tex+'$',fontsize=19)
   try:fig.savefig(f,bbox_inches='tight',pad_inches=.08,transparent=True)
   except Exception as e:raise RuntimeError(f'Formula {n}, slide {i}: {tex}') from e
   finally:plt.close(fig)
   data=base64.b64encode(f.read_bytes()).decode()
   return '\n<img class="equation" alt="'+html.escape(raw,quote=True)+'" src="data:image/svg+xml;base64,'+data+'">\n'
  body=re.sub(r'<eq>(.*?)</eq>',mathsvg,x['body'],flags=re.S)
  rendered=subprocess.check_output(['pandoc','-f','markdown','-t','html5'],input=body,text=True)
  article='<article id="theory-'+str(i)+'"><h2>'+html.escape(x['title'])+'</h2>'+rendered+'<p class="take"><strong>要点：</strong>'+html.escape(x['takeaway'])+'</p>'
  if x['notes']:article+='<p class="note">'+html.escape(x['notes'])+'</p>'
  article+='<p class="source">参考：<a href="'+x['source_url']+'">'+html.escape(x['source'])+'</a></p></article>'
  articles.append(article)
  md+=['## '+x['title']+'\n',x['body'].replace('<eq>','\n$$\n').replace('</eq>','\n$$\n')+'\n','**要点：'+x['takeaway']+'**\n',x['notes']+'\n','参考：['+x['source']+']('+x['source_url']+')\n']
 toc='<ol>'+''.join('<li><a href="#theory-'+str(i)+'">'+html.escape(x['title'])+'</a></li>' for i,x in enumerate(s,1))+'</ol>'
 css='''@page{size:A4;margin:18mm 17mm;@bottom-center{content:counter(page);font-size:9pt;color:#53677a}}body{font-family:"Noto Sans CJK JP",sans-serif;color:#172f43;font-size:11pt;line-height:1.7;max-width:960px;margin:auto}h1{font-size:24pt;color:#006c76}h2{font-size:16pt;color:#006c76;border-bottom:1px solid #bdd2d8;break-after:avoid;margin-top:28px}p{margin:8px 0}a{color:#006c76}.equation{display:block;max-width:100%;max-height:105px;margin:12px auto;break-inside:avoid}pre{font-family:"DejaVu Sans Mono","Noto Sans CJK JP",monospace;font-size:8pt;white-space:pre-wrap;overflow-wrap:anywhere;padding:12px;background:#f0f5f7;break-inside:avoid}table{border-collapse:collapse;width:100%;font-size:9pt}td,th{border-bottom:1px solid #c6d6dd;padding:7px;text-align:left}tr{break-inside:avoid}th{background:#eef6f5}.take{padding:9px;border-left:4px solid #007f86;background:#edf7f5}.note,.source{font-size:9pt;color:#4c6172}.toc{break-after:page;font-size:9pt}@media screen{body{padding:30px}.toc{columns:2}}'''
 doc='<!doctype html><html lang="ja"><meta charset="utf-8"><title>OI・KF・EnKF・PF 理論解説</title><style>'+css+'</style><body><h1>OI・KF・EnKF・PF<br>理論式の展開と104への接続</h1><p>記号 → 仮定 → 式の導出 → 数値例 → 実装の順で読む資料。数式は埋め込みSVGで、ネット接続なしでも表示できます。</p><div class="toc">'+toc+'</div>'+''.join(articles)+'</body></html>'
 (HERE/'theory_guide.md').write_text('\n'.join(md));(HERE/'theory_guide.html').write_text(doc)
 HTML(string=doc,base_url=str(HERE)).write_pdf(HERE/'theory_guide.pdf')
 print(f'Theory guide: {len(s)} topics, {n} rendered formulas, HTML/PDF/Markdown')
if __name__=='__main__':main()
