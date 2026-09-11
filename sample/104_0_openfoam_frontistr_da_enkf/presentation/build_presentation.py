"""Generate slides from slides.json without rerunning or changing solver cases.
Usage: python3 presentation/build_presentation.py
Dependencies: numpy, matplotlib, weasyprint, pandoc.
"""
from pathlib import Path
import pathlib
import base64, html, json, re, subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from weasyprint import HTML
from matplotlib import font_manager as _fm
import glob as _glob
for _fp in _glob.glob(str(pathlib.Path.home()/'.fonts'/'*CJK*.otf'))+_glob.glob('/usr/share/fonts/**/NotoSansCJK*.otf',recursive=True):
 try:_fm.fontManager.addfont(_fp)
 except Exception:pass
HERE=Path(__file__).absolute().parent
ROOT=HERE.parent
ASSETS=HERE/'assets'
plt.rcParams.update({'font.family':'Noto Sans CJK JP','font.size':15,'axes.spines.top':False,'axes.spines.right':False})
def charts():
 h=np.genfromtxt(ROOT/'results/openfoam_fem_enkf_history.csv',delimiter=',',names=True)
 for kind in ['rmse','q']:
  fig,ax=plt.subplots(figsize=(10.5,4.8),layout='constrained');t=h['time_s']
  if kind=='rmse':
   y=h['rmse_field_K'];ax.semilogy(t,y,'o-',lw=3,color='#007f86',ms=9)
   for x,v in zip(t,y):ax.annotate(f'{v:.4g} K',(x,v),xytext=(0,14),textcoords='offset points',ha='center')
   ax.set_ylim(.008,30);ax.set_ylabel('平均温度場と真値の RMSE [K]')
  else:
   y=h['q_mean_W'];s=h['q_std_W'];ax.fill_between(t,y-s,y+s,color='#007f86',alpha=.16,label='平均 ± 表示用標準偏差')
   ax.plot(t,y,'o-',lw=3,color='#007f86',ms=9,label='発熱量の推定平均');ax.axhline(15,color='#c15d35',ls='--',lw=2,label='真値 15 W')
   ax.set_ylim(3,23);ax.set_ylabel('発熱量 Q [W]');ax.legend(fontsize=11,loc='upper left')
  ax.set_xlabel('時刻 [s]（10 s・20 s は同化後）');ax.set_xticks(t);ax.set_xlim(-1,21);ax.grid(alpha=.2)
  fig.savefig(ASSETS/f'{kind}.png',dpi=170);plt.close(fig)
 u=np.array([[32.48,6.83,-10.70,-2.55,20.98],[.97,.71,.63,.81,.61]])
 fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
 for j,ax in enumerate(axs):
  ax.scatter(range(5),u[j],s=100,color='#007f86');ax.axhline([.409,.625][j],color='#c15d35',ls='--',label='合成観測（ノイズあり）')
  ax.set_title(f'{[10,20][j]}秒の同化前予報');ax.set_xticks(range(5),[f'{i:02d}' for i in range(5)]);ax.set_xlabel('メンバー');ax.set_ylabel('ヒータ側上面 Uz [µm]');ax.grid(alpha=.2);ax.legend(fontsize=10)
 fig.savefig(ASSETS/'forecast_displacement.png',dpi=170);plt.close(fig)
def uri(p):
 p=Path(p);mime={'.svg':'image/svg+xml','.gif':'image/gif','.jpg':'image/jpeg'}.get(p.suffix,'image/png')
 return 'data:'+mime+';base64,'+base64.b64encode(p.read_bytes()).decode()
def render(md,i):
 def equation(m):
  tex=m.group(1)
  tex=re.sub(r'\\le(?![a-zA-Z])','\\\\leq',tex);tex=re.sub(r'\\ge(?![a-zA-Z])','\\\\geq',tex)  # mathtext_fix
  path=ASSETS/f'eq_{i}_{m.start()}.svg'
  try:
   fig=plt.figure(figsize=(1,1));fig.text(0,0,'$'+tex+'$',fontsize=25)
   fig.savefig(path,bbox_inches='tight',pad_inches=.12,transparent=True);plt.close(fig)
  except Exception:
   plt.close('all')
   return '\n<pre class="equation-fallback">'+html.escape(tex)+'</pre>\n'  # mathtextで描けない式はTeX原文
  return '\n<img class="equation" src="'+uri(path)+'" alt="'+html.escape(tex,quote=True)+'">\n'
 md=re.sub(r'<eq>(.*?)</eq>',equation,md,flags=re.S)
 return subprocess.check_output(['pandoc','-f','markdown','-t','html5'],input=md,text=True)
def main():
 ASSETS.mkdir(exist_ok=True);charts();slides=json.loads((HERE/'slides.json').read_text());sections=[];notes=['# 発表者ノート\n\n所属・発表者名・学会名は未指定。表紙を発表前に設定する。\n'];mdout=['# 学会発表資料：温度・変位観測を用いた熱状態推定\n']
 for i,s in enumerate(slides,1):
  body=render(s.get('body',''),i);fig=''
  if s.get('image'):fig='<figure><img src="'+uri(HERE/s['image'])+'"><figcaption>'+html.escape(s.get('caption',''))+'</figcaption></figure>'
  layout=('visual' if fig else 'text')+(' wide' if s.get('wide') else '')+(' cover' if s.get('cover') else '')
  source=s.get('source','104ケースの設定・実装を基に作成')
  sections.append(f'<section class="slide {layout}" id="slide-{i}"><div class="eyebrow">{html.escape(s["part"])} · 104 / THERMAL DATA ASSIMILATION</div><h1>{html.escape(s["title"])}</h1><div class="content">{fig}<div class="body">{body}</div></div><div class="takeaway">{html.escape(s["takeaway"])}</div><footer><span>{html.escape(source)}</span><b>{i:02d} / {len(slides):02d}</b></footer></section>')
  notes.append(f'## {i:02d}. {s["title"]}\n\n{s["notes"]}\n\n出典：{source}\n')
  mdout.append(f'## {i:02d}. {s["title"]}\n\n{s.get("body", "")}\n\n'+(f'![{s.get("caption", "")} ]({s["image"]})\n\n' if fig else '')+f'**要点：{s["takeaway"]}**\n\n出典：{source}\n\n---\n')
 css=(HERE/'slides.css').read_text()
 doc='<!doctype html><html lang="ja"><head><meta charset="utf-8"><title>温度・変位観測を用いた熱状態推定 — 学会発表資料</title><style>'+css+'</style></head><body>'+''.join(sections)+'''<nav>← →：ページ移動 ／ F：全画面 ／ Ctrl+P：印刷</nav><script>let n=Number(location.hash.split('-')[1]||1);document.addEventListener('keydown',e=>{if(['ArrowRight','PageDown','ArrowLeft','PageUp'].includes(e.key)){e.preventDefault();n=Math.max(1,Math.min(document.querySelectorAll('.slide').length,n+(['ArrowRight','PageDown'].includes(e.key)?1:-1)));location.hash='slide-'+n;}if(e.key==='f'||e.key==='F'){document.fullscreenElement?document.exitFullscreen():document.documentElement.requestFullscreen();}});</script></body></html>'''
 (HERE/'conference_slides.html').write_text(doc);(HERE/'speaker_notes.md').write_text('\n'.join(notes));(HERE/'conference_slides.md').write_text('\n'.join(mdout))
 HTML(string=doc,base_url=str(HERE)).write_pdf(HERE/'conference_slides.pdf')
 print(f'Generated {len(slides)} slides: HTML / PDF / Markdown / notes')
if __name__=='__main__':main()
