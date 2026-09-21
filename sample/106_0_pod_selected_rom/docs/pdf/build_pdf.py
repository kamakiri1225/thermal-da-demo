#!/usr/bin/env python3
"""元MarkdownをPandocでTeXへ変換し、LuaLaTeXで論文調のPDFを生成する。"""
from pathlib import Path
import re, json, subprocess, hashlib, shutil, os, sys, html
from PIL import Image
HERE=Path(__file__).resolve().parent
DOC=HERE.parent
ASSETS=HERE/'tex_assets'
ASSETS.mkdir(exist_ok=True)
ENV=dict(os.environ,TEXMFVAR='/tmp/thermal-blog-texmf',TEXMFCACHE='/tmp/thermal-blog-texmf')
HEADER=r'''
\usepackage{luatexja-fontspec}
\setmainjfont[BoldFont=IPAexGothic,ItalicFont=IPAexMincho,BoldItalicFont=IPAexGothic]{IPAexMincho}
\setsansjfont[BoldFont=IPAexGothic,ItalicFont=IPAexMincho,BoldItalicFont=IPAexGothic]{IPAexGothic}
\setmonojfont[BoldFont=IPAexGothic,ItalicFont=IPAexGothic,BoldItalicFont=IPAexGothic]{IPAexGothic}
\setmathfont{TeX Gyre Termes Math}
\usepackage{fvextra}
\usepackage{seqsplit}
% 段幅フィット：実幅を測り、超過時のみresizebox。
% （adjustboxはpandocの\setkeys{Gin}既定と干渉して極小化するため使わない）
\makeatletter
\newcommand{\fitbox}[1]{\begingroup\setbox\z@\hbox{#1}%
 \ifdim\wd\z@>\linewidth \resizebox{\linewidth}{!}{\box\z@}\else\box\z@\fi\endgroup}
\makeatother
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{commandchars=\\\{\},breaklines,breakanywhere,fontsize=\small}
\usepackage{float}
\floatplacement{figure}{H}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf}
\usepackage{fancyhdr}
\pagestyle{fancy}\fancyhf{}\fancyfoot[C]{\thepage}
\renewcommand{\headrulewidth}{0pt}
\setlength{\emergencystretch}{3em}
\tolerance=2000
\setlength{\parindent}{1em}
\setlength{\parskip}{0.35em}
\usepackage{booktabs}
\setlength{\columnsep}{6.5mm}
\renewcommand{\arraystretch}{1.22}
\AtBeginDocument{\hypersetup{colorlinks=true,linkcolor=black,urlcolor=black}}
'''
(HERE/'paper_header.tex').write_text(HEADER)
def run(args, **kw):
 return subprocess.run(args,check=True,env=ENV,**kw)
def flow_pdf(code,stem):
 nodes={};edges=[]
 pattern=r'([A-Za-z][A-Za-z0-9_]*)(?:\[([^\]]*)\])?'
 for line in code.splitlines():
  if '-->' not in line:continue
  parts=line.split('-->')
  ids=[]
  for part in parts:
   m=re.search(pattern,part.strip())
   if not m:raise ValueError('Unsupported flow line '+line)
   key,label=m.groups();ids.append(key)
   if label is not None:nodes[key]=html.unescape(label.strip('"')).replace('<br/>','\n').replace('<br>','\n').replace('\\n','\n')
  edges.extend(zip(ids,ids[1:]))
 if not edges:raise ValueError('Empty flowchart')
 dot='digraph G {rankdir=LR; graph [bgcolor="white",pad="0.15",nodesep="0.3",ranksep="0.3"]; node [shape=box,style="rounded",fontname="IPAexMincho",fontsize=15,margin="0.15,0.12",color="#354858"]; edge [color="#354858"];\n'
 for k,v in nodes.items():dot+=f'{k} [label={json.dumps(v,ensure_ascii=False)}];\n'
 for a,b in edges:dot+=f'{a} -> {b};\n'
 dot+='}'
 target=ASSETS/(stem+'.pdf');(ASSETS/(stem+'.dot')).write_text(dot)
 run(['dot','-Tpdf',str(ASSETS/(stem+'.dot')),'-o',str(target)])
 return 'tex_assets/'+target.name

def make(source):
 stem=source.name[:8];raw=source.read_text();text=raw
 # HTML画像を本文から落とさずMarkdown画像へ。幅指定はPDFページに合わせる。
 def image_html(m):
  attrs=dict(re.findall(r'(\w+)="([^"]*)"',m.group()))
  return '\n\n!['+attrs.get('alt','図')+']('+attrs['src']+')\n\n'
 text=re.sub(r'<img\b[^>]*>',image_html,text)
 # Pandoc ASTを通して数式と図を構造的に処理。
 ast=json.loads(run(['pandoc','-f','markdown+tex_math_dollars+raw_html','-t','json'],input=text,text=True,capture_output=True).stdout)
 count={'math':0,'images':0,'flows':0,'gif':0,'tables':0}
 def visit(obj):
  if isinstance(obj,list):return [visit(v) for v in obj]
  if not isinstance(obj,dict):return obj
  typ=obj.get('t');c=obj.get('c')
  if typ=='Math':
   count['math']+=1
   math=c[1].replace(r'\gt', '>').replace(r'\lt','<')
   if c[0]['t']=='DisplayMath':
    return {'t':'RawInline','c':['latex',r'\[\fitbox{$\displaystyle '+math+r'$}\]']}
   obj['c'][1]=math
  elif typ=='Table':count['tables']+=1
  elif typ=='CodeBlock' and 'mermaid' in c[0][1]:
   count['flows']+=1;path=flow_pdf(c[1],f'{stem}_flow_{count["flows"]}')
   return {'t':'RawBlock','c':['latex',r'\begin{center}\includegraphics[width=\linewidth,height=0.30\textheight,keepaspectratio]{'+path+r'}\end{center}']}
  elif typ=='Image':
   count['images']+=1
   path=c[2][0];src=DOC/path
   if not src.exists():raise FileNotFoundError(src)
   if src.suffix.lower()=='.gif':
    count['gif']+=1;im=Image.open(src);im.seek(im.n_frames//2)
    out=ASSETS/(src.stem+'_frame.png');im.convert('RGB').save(out)
    c[1]+= [{'t':'Str','c':'（アニメーションの中央フレーム。動画は元記事参照）'}]
   else:
    out=ASSETS/src.name
    if src.resolve()!=out.resolve():shutil.copy2(src,out)
   c[2][0]='tex_assets/'+out.name
   c[0][2]=[['width','100%'],['height','0.36\\textheight']]
  elif typ=='Code':
   # 長いパス・CSVヘッダなどのインラインコードは任意位置で改行可にする。
   if len(c[1])>24:
    return {'t':'RawInline','c':['latex',r'\seqsplit{\texttt{\small '+c[1].replace('\\',r'\textbackslash{}').replace('_',r'\_').replace('%',r'\%').replace('&',r'\&').replace('#',r'\#').replace('$',r'\$')+'}}']}
  elif typ=='Link':
   # 長いパス風のリンク文字列は '/' '_' の後で改行を許可（二段組のはみ出し防止）。
   def breakable(inl):
    if isinstance(inl,dict) and inl.get('t')=='Str' and len(inl['c'])>20 and ('/' in inl['c'] or '_' in inl['c']):
     out=[];buf=''
     for ch in inl['c']:
      buf+=ch
      if ch in '/_':
       out.append({'t':'Str','c':buf});out.append({'t':'RawInline','c':['latex','\\allowbreak{}']});buf=''
     if buf:out.append({'t':'Str','c':buf})
     return out
    return [inl]
   c[1]=[x for inl in c[1] for x in breakable(inl)]
   target=c[2][0]
   if not target.startswith(('http','#','mailto:')):
    # PDFからも元リポジトリのコード・記事へ到達できるリンク。
    targetpath=(DOC/target).resolve()
    root=DOC.parents[2]
    try:rel=targetpath.relative_to(root)
    except ValueError:rel=None
    if rel is not None:c[2][0]='https://github.com/kamakiri1225/thermal-da-demo/blob/main/'+str(rel)
  return {k:visit(v) for k,v in obj.items()}
 ast=visit(ast)
 args=['pandoc','-f','json','-t','latex','-s','--top-level-division=section','--syntax-highlighting=tango',f'--lua-filter={HERE/"tables_twocol.lua"}','-V','documentclass=ltjsarticle','-V','classoption=10pt','-V','classoption=twocolumn','-V','geometry:margin=18mm','-V','mainfont=TeX Gyre Termes','-V','sansfont=TeX Gyre Heros','-V','monofont=IPAexGothic','-V','mathfont=TeX Gyre Termes Math','-V','linestretch=1.12','--include-in-header',str(HERE/'paper_header.tex')]
 # Pandoc 3.1 uses the old option name.
 args[args.index('--syntax-highlighting=tango')]='--highlight-style=tango'
 tex=run(args,input=json.dumps(ast),text=True,capture_output=True).stdout
 tex=re.sub(r'\\includegraphics\[[^\]]*\]', lambda m: r'\includegraphics[width=\linewidth,height=0.36\textheight,keepaspectratio]', tex)
 (HERE/(stem+'.tex')).write_text(tex)
 for passno in (1,2):
  with (HERE/(stem+f'.build{passno}.log')).open('w') as log:
   run(['lualatex','-interaction=nonstopmode','-halt-on-error',stem+'.tex'],cwd=HERE,stdout=log,stderr=subprocess.STDOUT)
 return {'source':source.name,'sha256':hashlib.sha256(raw.encode()).hexdigest(),**count}
if __name__=='__main__':
 selected=sys.argv[1:];records=[]
 for src in sorted(DOC.glob('blog_00[1-5]_*.md')):
  if selected and src.name[:8] not in selected:continue
  print('Building',src.name,flush=True);records.append(make(src));print('OK',records[-1],flush=True)
 (HERE/'build_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')
