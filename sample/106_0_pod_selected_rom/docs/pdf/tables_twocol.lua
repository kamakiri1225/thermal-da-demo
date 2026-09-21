-- 二段組対応：表を longtable ではなく booktabs の tabular にする。
-- 段幅を超える表だけ adjustbox(max width) で縮小（小さい表は拡大しない）。
function Table(t)
  local ncols = #t.colspecs
  local function cells_tex(row)
    local out = {}
    for _, cell in ipairs(row.cells) do
      local doc = pandoc.Pandoc(cell.contents)
      local tex = pandoc.write(doc, "latex")
      tex = tex:gsub("\n$", ""):gsub("\n", " ")
      table.insert(out, tex)
    end
    return table.concat(out, " & ")
  end
  local lines = {}
  table.insert(lines, "\\begin{center}")
  table.insert(lines, "\\fitbox{\\small\\begin{tabular}{" .. string.rep("l", ncols) .. "}")
  table.insert(lines, "\\toprule")
  for _, row in ipairs(t.head.rows) do
    table.insert(lines, cells_tex(row) .. " \\\\")
  end
  table.insert(lines, "\\midrule")
  for _, body in ipairs(t.bodies) do
    for _, row in ipairs(body.body) do
      table.insert(lines, cells_tex(row) .. " \\\\")
    end
  end
  table.insert(lines, "\\bottomrule")
  table.insert(lines, "\\end{tabular}}")
  table.insert(lines, "\\end{center}")
  return pandoc.RawBlock("latex", table.concat(lines, "\n"))
end
