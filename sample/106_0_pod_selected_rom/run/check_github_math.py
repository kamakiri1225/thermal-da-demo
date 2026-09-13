r"""GitHubのMarkdown数式レンダリングを壊す書き方を検出する pre-push リンター.

GitHubの $...$ / $$...$$ 数式は、Markdown処理やMathJaxの制約で
特定の書き方をすると「認識されない」または「崩れて表示される」。
push前にこのスクリプトで全ルールを検査する。

使い方:
  python3 run/check_github_math.py <file1.md> <file2.md> ...
  （引数なしなら 106/docs と 104 の主要mdを検査）
終了コード: 問題があれば1、なければ0。

検出ルール（すべて実レンダリング検証で確立）:
 R1  $ の総数が奇数（閉じ忘れ／行跨ぎ）
 R2  引用ブロック(>)内の $$表示数式（GitHubで描画されない）
 R3  インライン$の外側が「非ASCII文字」or「ASCII英数字」に密着
     → GitHubが数式境界と認識せず失敗。空白かASCII句読点の境界が必要
     （開き側は開き括弧 （「【 も不可。閉じ側の全角句読点 、。 は可）
 R4  インライン$内にパイプ | / \|（Markdownが表罫線と誤認）→ \lVert \rVert 等へ
 R5  インライン$内に '' アポストロフィ（引用符と誤認）→ \prime へ
 R6  インライン$内に }_ ]_ )_（閉じ括弧直後の下線＝Markdown強調<em>を誘発）
     → 下線を文字直後へ並べ替え、または ```math フェンス化
 R7  数式内(inline/display/```math)にCJK文字（MathJaxに和文字体が無く□表示）
     → \text{...}は英語に、和文説明は数式の外へ
 R8  数式内に生% コメント（\% でなく %）→ 以降がコメント化して壊れる
 R9  数式の波括弧 {} が不一致
 R10 数式内に生Unicodeギリシャ文字（σ 等）→ \sigma 等へ
"""
from __future__ import annotations
import os, re, sys

INLINE=re.compile(r'(?<!\$)\$(?!\$)([^$\n]+?)\$(?!\$)')
DISP=re.compile(r'\$\$(.+?)\$\$', re.S)
MATHFENCE=re.compile(r'```math\n(.*?)```', re.S)
GREEK="αβγδεζηθικλμνξοπρστυφχψωΓΔΘΛΞΠΣΦΨΩ"
SAFE_OPEN=set("（「『【〔《〈｛［")

def is_cjk(c: str) -> bool:
    if not c: return False
    o=ord(c)
    return (0x3040<=o<=0x30FF or 0x4E00<=o<=0x9FFF or 0x3400<=o<=0x4DBF
            or 0x3000<=o<=0x303F or 0xFF00<=o<=0xFFEF)

# GitHubのMathJaxが禁止するマクロ（使うと "macros are not allowed" エラー表示になる）
DISALLOWED=["operatorname","newcommand","def","let","require","class",
            "cssId","style","href","renewcommand","providecommand"]

def check_expr_content(expr: str):
    """1つの数式内容の R7-R11 を返す（issue文字列のリスト）"""
    iss=[]
    if any(c in GREEK for c in expr): iss.append("R10 生ギリシャ文字")
    if any(is_cjk(c) for c in expr): iss.append("R7 数式内CJK")
    if re.search(r'(?<!\\)%', expr): iss.append("R8 生%コメント")
    if expr.count('{')!=expr.count('}'): iss.append("R9 波括弧不一致")
    for mac in DISALLOWED:
        if re.search(r'\\'+mac+r'[^A-Za-z]', expr+" "):
            iss.append(f"R11 GitHub禁止マクロ \\{mac}")
    return iss

def lint(path: str):
    problems=[]
    text=open(path, encoding="utf-8").read()
    lines=text.split("\n")
    in_fence=False; dollar_total=0
    for i,ln in enumerate(lines,1):
        s=ln.strip()
        if s.startswith("```"): in_fence=not in_fence; continue
        if in_fence: continue
        dollar_total+=len(re.findall(r'(?<!\\)\$', ln))
        if s.startswith(">") and "$$" in ln:
            problems.append((i,"R2 引用ブロック内の$$表示数式",s[:50]))
        for m in INLINE.finditer(ln):
            g=m.group(1)
            b=ln[m.start()-1] if m.start()>0 else ""
            a=ln[m.end()] if m.end()<len(ln) else ""
            # R3 境界
            if b and (ord(b)>=0x80 and b not in SAFE_OPEN) or (b and b.isascii() and b.isalnum()):
                problems.append((i,"R3 開き$の直前に境界なし(要空白)",f"…{b}${g[:18]}$"))
            if (a and is_cjk(a)) or (a and a.isascii() and a.isalnum()):
                problems.append((i,"R3 閉じ$の直後に境界なし(要空白)",f"${g[:18]}${a}…"))
            # R4-R6
            if "|" in g: problems.append((i,"R4 数式内パイプ|",f"${g[:24]}$"))
            if "''" in g: problems.append((i,"R5 アポストロフィ''",f"${g[:24]}$"))
            if re.search(r'[)\]}]_', g): problems.append((i,"R6 閉じ括弧直後の下線(強調誘発)",f"${g[:24]}$"))
            for it in check_expr_content(g): problems.append((i,it+"(inline)",f"${g[:24]}$"))
    # 表示数式・```math（R7-R10のみ）
    body=re.sub(r'```(?!math).*?```','',text,flags=re.S)
    for m in MATHFENCE.finditer(text):
        for it in check_expr_content(m.group(1)):
            problems.append(("?",it+"(```math)",m.group(1)[:30].replace("\n"," ")))
    t2=re.sub(r'```.*?```','',text,flags=re.S)
    for m in DISP.finditer(t2):
        for it in check_expr_content(m.group(1)):
            problems.append(("?",it+"(display)",m.group(1)[:30].replace("\n"," ")))
    if dollar_total%2==1:
        problems.append((0,"R1 $の総数が奇数",f"total={dollar_total}"))
    return problems

def main(argv):
    if argv:
        files=argv
    else:
        base=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        sample=os.path.dirname(base)
        files=[os.path.join(base,"docs",x) for x in
               ("00_data_assimilation_algorithm.md","01_pod_qdeim_algorithm.md",
                "02_full_story.md","03_presentation_story.md")]
        files.append(os.path.join(sample,"104_0_openfoam_frontistr_da_enkf","docs","16_presentation_story.md"))
    total=0
    for f in files:
        if not os.path.exists(f): print(f"[skip] {f}"); continue
        p=lint(f)
        name=os.path.basename(f)
        if not p:
            print(f"✓ {name}: 問題なし")
        else:
            print(f"✗ {name}: {len(p)}件")
            for ln,rule,ctx in p[:40]:
                print(f"    L{ln}: {rule}  {ctx}")
            total+=len(p)
    print(f"\n{'合格（全ルールOK）' if total==0 else f'不合格: 計{total}件'}")
    return 1 if total else 0

if __name__=="__main__":
    sys.exit(main(sys.argv[1:]))
