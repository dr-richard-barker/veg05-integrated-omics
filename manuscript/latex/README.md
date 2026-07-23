# LaTeX manuscript (npj Microgravity / Springer Nature style)

Assembles the manuscript in the official **Springer Nature LaTeX template**
(`sn-jnl` class) — the format npj Microgravity accepts and typesets from.

```
latex/
├── main.tex          # assembled manuscript (sn-jnl class, sn-nature refs)
├── references.bib    # 22 references transcribed from the manuscript
├── figures/          # the 13 figures (PNG)
└── README.md         # this file
```

## How to compile

`sn-jnl.cls` / `sn-nature.bst` ship with Springer Nature's official template
(not vendored here).

- **Overleaf (recommended):** new project from the **"Springer Nature Article
  Template (sn-jnl)"** → replace its `main.tex` with this one, upload
  `references.bib` and `figures/`, compile with **pdfLaTeX**.
- **Local:** place `sn-jnl.cls` + `sn-nature.bst` here, then
  `pdflatex main` → `bibtex main` → `pdflatex main` → `pdflatex main`.

## Status / TODO before submission

- [ ] **Not yet compile-tested** — authored without a local TeX install; build
      once on Overleaf and fix any stragglers. Table 1 uses `longtable`
      (single-column); if the SN class is set to two-column, switch it to a
      standard `table`/`tabular` or `\onecolumn` around it.
- [ ] **Author block** — confirm author list, ORCID(s), and affiliation
      (currently a single-author placeholder; `% TODO` in `main.tex`).
- [ ] **Supplementary tables S1–S10** — referenced in-text; bundle the CSVs as
      supplementary data per journal instructions.
- [ ] **References** — 22 entries transcribed with journal/volume/pages; add
      DOIs where available.
- [ ] **Figures** — repo PNGs. npj prefers vector (PDF/EPS) or ≥300 dpi for
      final submission; swap files in `figures/` and paths resolve unchanged.

## Source

Ported from `../manuscript_with_figures.md`. Body text, Table 1, figures, and
references are the author's own content — nothing was invented.
