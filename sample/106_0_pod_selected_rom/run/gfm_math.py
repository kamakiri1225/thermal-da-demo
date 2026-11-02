"""Markdown の数式を GitHub で崩れない書き方に置き換える（中身は変えない）。

  ブロック：$$ … $$  →  ```math … ```
  文中　　：$ … $    →  $` … `$
  コード（``` の中、`…` の中）は触らない。pandoc の tex_math_dollars と同じ規則で $ を見つける。
  python3 gfm_math.py FILE...   （上書き）
"""
import re
import sys


def convert(text):
    out = []
    lines = text.split('\n')
    i = 0
    fence = None
    while i < len(lines):
        l = lines[i]
        s = l.strip()
        m = re.match(r'^(\s*)(`{3,}|~{3,})', l)
        if fence:
            out.append(l)
            if m and m.group(2)[0] == fence[0] and len(m.group(2)) >= len(fence) and s == m.group(2):
                fence = None
            i += 1; continue
        if m:
            fence = m.group(2); out.append(l); i += 1; continue
        ind = re.match(r'^(\s*(?:>\s*)?)', l).group(1)
        body = l[len(ind):]
        # 1 行の $$ … $$（行全体）
        mm = re.match(r'^\$\$(.+)\$\$\s*$', body)
        if mm and '$$' not in mm.group(1):
            out += [ind + '```math', ind + mm.group(1).strip(), ind + '```']; i += 1; continue
        # 複数行の $$ ブロック
        if body.startswith('$$') and body.count('$$') == 1:
            first = body[2:].strip()
            block = [first] if first else []
            j = i + 1
            while j < len(lines) and '$$' not in lines[j]:
                block.append(lines[j][len(ind):] if lines[j].startswith(ind) else lines[j]); j += 1
            if j >= len(lines):
                raise SystemExit(f'閉じていない $$ : 行 {i + 1}')
            last = lines[j][len(ind):] if lines[j].startswith(ind) else lines[j]
            pre, post = last.split('$$', 1)
            if pre.strip(): block.append(pre.rstrip())
            out += [ind + '```math'] + [ind + b if b else '' for b in block] + [ind + '```']
            if post.strip(): out.append(ind + post.strip())
            i = j + 1; continue
        out.append(convert_inline(l)); i += 1
    return '\n'.join(out)


def convert_inline(l):
    res = []; k = 0; n = len(l)
    while k < n:
        c = l[k]
        if c == '`':                                   # コードスパンは飛ばす
            run = re.match(r'`+', l[k:]).group(0); e = l.find(run, k + len(run))
            if e < 0: res.append(l[k:]); break
            res.append(l[k:e + len(run)]); k = e + len(run); continue
        if c == '\\' and k + 1 < n:
            res.append(l[k:k + 2]); k += 2; continue
        if c == '$' and l[k:k + 2] != '$$' and k + 1 < n and not l[k + 1].isspace() and l[k + 1] != '`':
            e = k + 1; found = -1
            while e < n:
                if l[e] == '\\': e += 2; continue
                if l[e] == '$':
                    if not l[e - 1].isspace() and not (e + 1 < n and l[e + 1].isdigit()):
                        found = e
                    break
                e += 1
            if found > 0:
                res.append('$`' + l[k + 1:found] + '`$'); k = found + 1; continue
        res.append(c); k += 1
    return ''.join(res)


if __name__ == '__main__':
    for p in sys.argv[1:]:
        t = open(p, encoding='utf-8').read()
        open(p, 'w', encoding='utf-8').write(convert(t))
        print('converted', p)
