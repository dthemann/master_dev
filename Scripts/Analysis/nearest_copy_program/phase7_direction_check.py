#!/usr/bin/env python
"""Direction rule for Phase 7: every value written into the yaml copy must be a value the thesis COPY prints.

Reads harness_signatures/phase7_patch_ledger.csv (key, line, col, old, new) and checks, for each patched token,
that the new value appears as a number token in thesis_latex_nearest/body_main_short.tex or body_appendix_short.tex
(or Thesis_short.tex for abstract numbers), and reports whether the old value has vanished from the same files.
Percentages may be printed with a following \\% and negative numbers with \\ensuremath{-}; both spellings are
searched. Values the copy does not print are listed for manual reconciliation (they are either derived quantities
the harness checks but the thesis never prints as such, or a genuine transcription miss). Read-only.
"""
from __future__ import annotations
import csv, re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
COPY = HERE.parents[2] / "thesis_latex_nearest"

def tokens(text: str) -> set[str]:
    text = text.replace("\\ensuremath{-}", "-").replace("\\%", "%").replace("--", " ")
    return set(re.findall(r"(?<![\w.])-?\d+(?:\.\d+)?(?:e-?\d+)?(?![\w.])", text))

def main():
    tex = "\n".join((COPY / f).read_text() for f in ("body_main_short.tex", "body_appendix_short.tex", "Thesis_short.tex"))
    toks = tokens(tex)
    rows = list(csv.DictReader(open(HERE / "harness_signatures" / "phase7_patch_ledger.csv")))
    miss, old_alive = [], []
    for r in rows:
        new, old = r["new"], r["old"]
        variants = {new, new.lstrip("+")}
        if "." in new:
            variants.add(new.rstrip("0").rstrip("."))
        if not variants & toks:
            miss.append(r)
        if old in toks and old not in variants:
            old_alive.append(r)
    print(f"{len(rows)} patched values; {len(rows) - len(miss)} printed in the copy; {len(miss)} NOT found; {len(old_alive)} old values still printed somewhere (may be legitimate reuse)")
    out = HERE / "harness_signatures" / "phase7_direction_misses.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) + ["status"]); w.writeheader()
        for r in miss: w.writerow({**r, "status": "new value not printed"})
        for r in old_alive: w.writerow({**r, "status": "old value still printed"})
    for r in miss[:60]:
        print(f"  NOT PRINTED  {r['key']:60s} {r['old']} -> {r['new']}")

if __name__ == "__main__":
    main()
