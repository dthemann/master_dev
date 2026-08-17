# LaTeX thesis project

This project converts the consistency-reviewed Word thesis into the supplied
FH Technikum Wien `twbook` template.

Build in the project directory with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error Thesis.tex
```

The main document is `Thesis.tex`; the generated content is split into
`body_main.tex` and `body_appendix.tex`. The 72 reviewed references are in
`Literatur.bib`, and extracted figures are under `media/media/`.

Seven non-printing `% REVIEW QUERY:` comments preserve unresolved source-data
questions from the Word review. The eighth editorial query (mixed figure/table
numbering) is resolved here through automatic LaTeX numbering and references.
The original Pandoc conversion is retained as `body_raw.tex` for auditability.
