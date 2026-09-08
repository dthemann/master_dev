#!/usr/bin/env bash
# Phase 8: build the independent thesis copy with pdflatex via latexmk and rasterise the pages
# that hold the regenerated figures and changed tables. Writes ONLY inside thesis_latex_nearest/
# and Scripts/Analysis/nearest_copy_program/qa_raster/. The shipped thesis_latex/ is never touched.
set -euo pipefail
COPY=/home/manndo/master_dev/thesis_latex_nearest
QA=/home/manndo/master_dev/Scripts/Analysis/nearest_copy_program/qa_raster
cd "$COPY"
latexmk -pdf -f Thesis_short.tex > latexmk_build.log 2>&1 || true
echo "latexmk exit: ${PIPESTATUS[0]:-0}"; grep -c '^!' Thesis_short.log || true
grep -n '^!' Thesis_short.log | head -10 || true
grep -c 'Overfull' Thesis_short.log || true
grep -n 'undefined' Thesis_short.log | head -5 || true
mkdir -p "$QA"
/home/manndo/anaconda3/envs/vina/bin/python - <<'PY'
import fitz, os
doc = fitz.open("/home/manndo/master_dev/thesis_latex_nearest/Thesis_short.pdf")
print("pages:", doc.page_count)
qa = "/home/manndo/master_dev/Scripts/Analysis/nearest_copy_program/qa_raster"
# rasterise abstract/Kurzfassung pages and every page whose text mentions a regenerated figure or changed table
want = set()
for i, p in enumerate(doc):
    t = p.get_text()
    if any(k in t for k in ("Kurzfassung", "Abstract")) and i < 8: want.add(i)
    if any(k in t for k in ("nearest deposited copy", "nearest of all", "Reference convention", "single deposited instance")): want.add(i)
for i in sorted(want)[:40]:
    doc[i].get_pixmap(dpi=80).save(f"{qa}/page_{i+1:03d}.png")
print("rasterised pages:", [i+1 for i in sorted(want)[:40]])
PY
