from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


DOC_DIR = Path(__file__).resolve().parent
OUT_DIR = DOC_DIR / "pdf"
BUILD_DIR = OUT_DIR / "_build"


def latex_escape(text: str) -> str:
    repl = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(repl.get(ch, ch) for ch in text)


def normalize_plain_text(text: str) -> str:
    """PDFの本文・コードで欠けやすいUnicode記号を安全なASCII表記に寄せる。"""
    sup_digits = {
        "⁰": "0",
        "¹": "1",
        "²": "2",
        "³": "3",
        "⁴": "4",
        "⁵": "5",
        "⁶": "6",
        "⁷": "7",
        "⁸": "8",
        "⁹": "9",
    }
    text = re.sub(r"⁻([⁰¹²³⁴⁵⁶⁷⁸⁹])", lambda m: "^-" + sup_digits[m.group(1)], text)

    repl = {
        "⁻": "^-",
        "⁰": "^0",
        "¹": "^1",
        "²": "^2",
        "³": "^3",
        "⁴": "^4",
        "⁵": "^5",
        "⁶": "^6",
        "⁷": "^7",
        "⁸": "^8",
        "⁹": "^9",
        "₀": "_0",
        "₁": "_1",
        "₂": "_2",
        "₃": "_3",
        "₄": "_4",
        "₅": "_5",
        "₆": "_6",
        "₇": "_7",
        "₈": "_8",
        "₉": "_9",
        "ᵀ": "^T",
        "↕": "|",
    }
    return "".join(repl.get(ch, ch) for ch in text)


def normalize_math_text(text: str) -> str:
    """LaTeX数式中でフォント欠けしやすい記号をLaTeX表記へ寄せる。"""
    text = text.replace("°\\text{C}", r"^\circ\text{C}")
    text = text.replace("°C", r"^\circ\mathrm{C}")
    return text


def inline_text_only(text: str) -> str:
    """数式ではない通常テキストだけをLaTeX用に変換する。"""
    text = normalize_plain_text(text)
    text = latex_escape(text)
    text = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", text)
    text = re.sub(r"`([^`]+)`", r"\\texttt{\1}", text)
    return text


def inline(text: str) -> str:
    """Markdownの文中数式をLaTeX数式として残し、それ以外をエスケープする。"""
    out: list[str] = []
    pos = 0

    for match in re.finditer(r"(?<!\\)\$(?!\$)(.+?)(?<!\\)\$", text):
        out.append(inline_text_only(text[pos : match.start()]))
        out.append(r"\(" + normalize_math_text(match.group(1)) + r"\)")
        pos = match.end()

    out.append(inline_text_only(text[pos:]))
    return "".join(out)


def display_math(content: str) -> str:
    """表示数式をLaTeXのdisplay math環境へ変換する。"""
    return "\\[\n" + normalize_math_text(content.strip()) + "\n\\]\n\n"


def image_latex(alt: str, rel: str, base: Path) -> str:
    img = (base / rel).resolve()
    if not img.exists():
        return inline(f"[image not found: {rel}]") + "\n\n"
    tex_path = img.as_posix()
    return (
        "\\begin{figure}[H]\n"
        "\\centering\n"
        f"\\includegraphics[width=0.95\\linewidth]{{{tex_path}}}\n"
        f"\\caption*{{{inline(alt)}}}\n"
        "\\end{figure}\n\n"
    )


def table_to_latex(rows: list[str]) -> str:
    parsed = []
    for row in rows:
        cells = [c.strip() for c in row.strip().strip("|").split("|")]
        parsed.append(cells)
    if len(parsed) < 2:
        return "\\begin{verbatim}\n" + "\n".join(rows) + "\n\\end{verbatim}\n\n"
    ncols = max(len(r) for r in parsed)
    colspec = "|" + "|".join(["p{%.3f\\linewidth}" % (0.92 / ncols)] * ncols) + "|"
    body = [f"\\begin{{longtable}}{{{colspec}}}", "\\hline"]
    for idx, cells in enumerate(parsed):
        if idx == 1 and all(re.match(r"^:?-{3,}:?$", c) for c in cells):
            body.append("\\hline")
            continue
        cells = cells + [""] * (ncols - len(cells))
        body.append(" & ".join(inline(c) for c in cells) + r" \\")
        body.append("\\hline")
    body.append("\\end{longtable}\n")
    return "\n".join(body) + "\n"


def markdown_to_latex(md: str, base: Path) -> str:
    lines = md.splitlines()
    out: list[str] = []
    i = 0
    in_code = False
    code_lang = ""
    code_buf: list[str] = []
    in_math = False
    math_buf: list[str] = []

    while i < len(lines):
        line = lines[i]

        if in_code:
            if line.strip().startswith("```"):
                if code_lang == "math":
                    out.append(display_math("\n".join(code_buf)))
                else:
                    out.append("\\begin{Verbatim}[fontsize=\\small]\n")
                    out.append("\n".join(normalize_plain_text(line) for line in code_buf))
                    out.append("\n\\end{Verbatim}\n\n")
                code_buf = []
                code_lang = ""
                in_code = False
            else:
                code_buf.append(line)
            i += 1
            continue

        if in_math:
            stripped = line.strip()
            if stripped.endswith("$$"):
                before_end = stripped[:-2].strip()
                if before_end:
                    math_buf.append(before_end)
                out.append(display_math("\n".join(math_buf)))
                math_buf = []
                in_math = False
            else:
                math_buf.append(line)
            i += 1
            continue

        stripped = line.strip()
        if stripped.startswith("```"):
            in_code = True
            code_lang = stripped[3:].strip().lower()
            code_buf = []
            i += 1
            continue
        single_line_math = re.match(r"^\$\$(.*)\$\$$", stripped)
        if single_line_math:
            out.append(display_math(single_line_math.group(1)))
            i += 1
            continue
        if stripped.startswith("$$"):
            in_math = True
            first = stripped[2:].strip()
            math_buf = [first] if first else []
            i += 1
            continue
        if stripped == "---":
            out.append("\\bigskip\\hrule\\bigskip\n")
            i += 1
            continue
        if not stripped:
            out.append("\n")
            i += 1
            continue

        img = re.match(r"!\[([^\]]*)\]\(([^)]+)\)", stripped)
        if img:
            out.append(image_latex(img.group(1), img.group(2), base))
            i += 1
            continue

        if stripped.startswith("|") and "|" in stripped[1:]:
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                rows.append(lines[i])
                i += 1
            out.append(table_to_latex(rows))
            continue

        m = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if m:
            level = len(m.group(1))
            title = inline(m.group(2))
            cmd = {
                1: "section*",
                2: "subsection*",
                3: "subsubsection*",
                4: "paragraph*",
                5: "paragraph*",
                6: "paragraph*",
            }[level]
            out.append(f"\\{cmd}{{{title}}}\n")
            i += 1
            continue

        if stripped.startswith("- "):
            out.append("\\begin{itemize}\n")
            while i < len(lines) and lines[i].strip().startswith("- "):
                out.append("\\item " + inline(lines[i].strip()[2:]) + "\n")
                i += 1
            out.append("\\end{itemize}\n\n")
            continue

        if re.match(r"^\d+\.\s+", stripped):
            out.append("\\begin{enumerate}\n")
            while i < len(lines) and re.match(r"^\d+\.\s+", lines[i].strip()):
                item = re.sub(r"^\d+\.\s+", "", lines[i].strip())
                out.append("\\item " + inline(item) + "\n")
                i += 1
            out.append("\\end{enumerate}\n\n")
            continue

        quote = stripped.startswith("> ")
        if quote:
            out.append("\\begin{quote}\n" + inline(stripped[2:]) + "\n\\end{quote}\n\n")
        else:
            out.append(inline(stripped) + "\n\n")
        i += 1

    return "".join(out)


def make_tex(md_path: Path) -> str:
    body = markdown_to_latex(md_path.read_text(encoding="utf-8"), md_path.parent)
    title = latex_escape(md_path.stem)
    return rf"""\documentclass[a4paper,11pt]{{ltjsarticle}}
\usepackage[margin=22mm]{{geometry}}
\usepackage{{amsmath,amssymb}}
\usepackage{{graphicx}}
\usepackage{{float}}
\usepackage{{caption}}
\usepackage{{longtable,array}}
\usepackage{{fancyvrb}}
\usepackage{{hyperref}}
\hypersetup{{colorlinks=true,linkcolor=blue,urlcolor=blue}}
\setlength{{\parindent}}{{0pt}}
\setlength{{\parskip}}{{0.6em}}
\title{{{title}}}
\date{{}}
\begin{{document}}
\maketitle
{body}
\end{{document}}
"""


def convert(md_path: Path) -> None:
    OUT_DIR.mkdir(exist_ok=True)
    BUILD_DIR.mkdir(exist_ok=True)
    tex_path = BUILD_DIR / f"{md_path.stem}.tex"
    tex_path.write_text(make_tex(md_path), encoding="utf-8")
    subprocess.run(
        [
            "lualatex",
            "-interaction=nonstopmode",
            "-halt-on-error",
            f"-output-directory={BUILD_DIR}",
            str(tex_path),
        ],
        cwd=DOC_DIR,
        check=True,
    )
    pdf_src = BUILD_DIR / f"{md_path.stem}.pdf"
    pdf_dst = OUT_DIR / f"{md_path.stem}.pdf"
    pdf_dst.write_bytes(pdf_src.read_bytes())


def main() -> int:
    md_files = sorted(p for p in DOC_DIR.glob("*.md") if p.name.lower() != "md_to_pdf.md")
    failed = []
    for md in md_files:
        print(f"converting {md.name} ...")
        try:
            convert(md)
        except Exception as exc:
            failed.append((md.name, str(exc)))
            print(f"FAILED {md.name}: {exc}", file=sys.stderr)
    if failed:
        print("\nfailed files:", file=sys.stderr)
        for name, err in failed:
            print(f"- {name}: {err}", file=sys.stderr)
        return 1
    print(f"\nPDF files written to: {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
