# LaTeX thesis project

This project converts the consistency-reviewed Word thesis into the supplied
FH Technikum Wien `twbook` template.

Build in the project directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error Thesis.tex
```

The main document is `Thesis.tex`; the generated content is split into
`body_main.tex` and `body_appendix.tex`. The 91 reviewed source records are in
`Literatur.bib`. 89 are substantively used and printed in the thesis bibliography.
Extracted figures are under `media/media/`.

The earlier source-data review queries have been resolved in the current sources.
Mixed figure/table numbering from the Word document is handled through automatic LaTeX numbering and references.
The original Pandoc conversion is retained as `body_raw.tex` for auditability.
