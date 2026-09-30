"""Render the current four HTML posters with offline TeX SVGs and WeasyPrint.

Usage: python3 run/render_conference_posters_pdf.py --mathjax-root /path/to/mathjax-full
Dependencies: Node.js, mathjax-full, beautifulsoup4, weasyprint, PyMuPDF.
Animations become static frames in the PDF; use the HTML to see the animation.
"""
from pathlib import Path
import argparse
import json
import re
import subprocess
import tempfile
from bs4 import BeautifulSoup
from weasyprint import HTML
import fitz

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mathjax-root', required=True, type=Path)
    args = parser.parse_args()
    source = ROOT / 'presentation/06_conference_posters.html'
    soup = BeautifulSoup(source.read_text(), 'html.parser')
    # Remove scripts and navigation: PDF contains the four posters themselves.
    for element in soup.select('script,nav,.sources,dialog'):
        element.decompose()
    for element in list(soup.contents):
        if element.__class__.__name__ == 'Doctype':
            element.extract()
    expressions = []
    for node in list(soup.find_all(string=True)):
        if node.parent.name == 'style':
            continue
        def placeholder(match):
            expressions.append({'tex': match.group(1) or match.group(2),
                                'display': match.group(2) is not None})
            return f'__POSTER_MATH_{len(expressions)-1}__'
        node.replace_with(re.sub(r'\\\((.*?)\\\)|\\\[(.*?)\\\]', placeholder,
                                 str(node), flags=re.S))
    javascript = r'''
const fs=require('fs'),path=require('path');
const root=process.argv[2];
const {mathjax}=require(path.join(root,'js/mathjax.js'));
const {TeX}=require(path.join(root,'js/input/tex.js'));
const {SVG}=require(path.join(root,'js/output/svg.js'));
const {liteAdaptor}=require(path.join(root,'js/adaptors/liteAdaptor.js'));
const {RegisterHTMLHandler}=require(path.join(root,'js/handlers/html.js'));
const {AllPackages}=require(path.join(root,'js/input/tex/AllPackages.js'));
const adaptor=liteAdaptor(); RegisterHTMLHandler(adaptor);
const document=mathjax.document('',{InputJax:new TeX({packages:AllPackages}),OutputJax:new SVG({fontCache:'none'})});
const items=JSON.parse(fs.readFileSync(process.argv[3]));
fs.writeFileSync(process.argv[4],JSON.stringify(items.map(item=>adaptor.outerHTML(document.convert(item.tex,{display:item.display})))));
'''
    with tempfile.TemporaryDirectory(prefix='thermal-posters-') as directory:
        temp = Path(directory)
        (temp/'math.json').write_text(json.dumps(expressions))
        (temp/'math.js').write_text(javascript)
        subprocess.run(['node',str(temp/'math.js'),str(args.mathjax_root.resolve()),
                        str(temp/'math.json'),str(temp/'svg.json')],check=True)
        html = str(soup)
        for i,svg in enumerate(json.loads((temp/'svg.json').read_text())):
            if 'data-mjx-error' in svg:
                raise ValueError('Invalid TeX: '+expressions[i]['tex'])
            html = html.replace(f'__POSTER_MATH_{i}__',svg)
        html = html.replace('</head>', '''<style>
body{font-family:IPAexGothic,sans-serif}
mjx-container{display:inline-block;vertical-align:middle}
mjx-container[display="true"]{display:block;text-align:center}
.col{display:block}article{margin-bottom:9px}
@media print{.poster:last-child{break-after:auto}}
</style></head>''')
        pdf = temp/'posters.pdf'
        HTML(string=html,base_url=str(source.parent)).write_pdf(pdf)
        with fitz.open(pdf) as document:
            if len(document)!=4:
                raise ValueError(f'Expected 4 pages, got {len(document)}; existing PDF preserved')
        target = source.with_suffix('.pdf')
        target.write_bytes(pdf.read_bytes())
        print(f'Generated {target}: 4 pages, {len(expressions)} validated TeX expressions')

if __name__=='__main__':
    main()
