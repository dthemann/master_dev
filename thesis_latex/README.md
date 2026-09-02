# LaTeX thesis project

The short thesis is the working version. Build it in this directory with:

```bash
latexmk -f -pdf -interaction=nonstopmode -file-line-error Thesis_short.tex
```

`Thesis_short.tex` is the main document. Content is split into
`body_main_short.tex` and `body_appendix_short.tex`. Sources are in
`Literatur.bib`, the template is `twbook.cls`, and figures are under
`media/media/`. The current build is 142 pages.

`make_contact_sheets.py` tiles rendered pages from `tmp/pdf_final/pages` into
contact sheets for visual QA. Rasterise the PDF into that directory first.

The full-thesis build, its editing backups, the review reports and the stale
page renders were retired on 2026-09-02 into `obsolete/`, which is gitignored.
See `obsolete/MANIFEST.md` for the inventory and how to restore a file.
