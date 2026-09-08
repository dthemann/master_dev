#!/usr/bin/env python
"""Phase 7 direction check: every value the program yaml asserts must be PRINTED in the thesis copy.

Reads a harness signature CSV (key, expected, actual, ok, source), takes the rows whose `ok` is False
(the rows that move under the switch), and for each numeric value inside `actual` checks that the
literal (with the thesis's formatting variants: 1,622 / 1622, 36.6 / 36.6\\%, a trailing zero dropped
or added, integers as counts) occurs in the copy's tex files named by the row's source. A value that
exists in the outputs but nowhere in the tex is a copy-from-outputs row and is reported. Read-only.
"""
from __future__ import annotations
import argparse, ast, re, sys
from pathlib import Path

REPO = Path("/home/manndo/master_dev")
COPY = REPO / "thesis_latex_nearest"

def variants(v) -> set[str]:
    out = set()
    if isinstance(v, bool):
        return out
    if isinstance(v, int):
        out.add(str(v)); out.add(f"{v:,}")
        return out
    if isinstance(v, float):
        for nd in (0, 1, 2, 3, 4):
            s = f"{v:.{nd}f}"
            out.add(s); out.add(s.replace(".", ","))
            if abs(v) >= 1000: out.add(f"{v:,.{nd}f}")
        return out
    return out

def flatten(x):
    if isinstance(x, (list, tuple)):
        for y in x: yield from flatten(y)
    elif isinstance(x, dict):
        for y in x.values(): yield from flatten(y)
    else:
        yield x

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--signature", type=Path, required=True)
    ap.add_argument("--tex", type=Path, nargs="*", default=[COPY / "body_main_short.tex", COPY / "body_appendix_short.tex", COPY / "Thesis_short.tex"])
    ap.add_argument("--only-failing", action="store_true", default=True)
    a = ap.parse_args()
    import csv
    text = "\n".join(p.read_text() for p in a.tex)
    text_c = text.replace("\\,", ",").replace("{,}", ",")
    rows = list(csv.DictReader(open(a.signature)))
    missing, checked = [], 0
    for r in rows:
        if a.only_failing and r["ok"] == "True":
            continue
        try:
            val = ast.literal_eval(r["actual"].replace("np.float64(", "(").replace("np.int64(", "("))
        except Exception:
            continue
        for v in flatten(val):
            vs = variants(v)
            if not vs: continue
            checked += 1
            if not any(s in text_c for s in vs):
                missing.append((r["key"], v))
    print(f"checked {checked} numeric values from {sum(r['ok']!='True' for r in rows)} moving rows; {len(missing)} not found in the copy's tex")
    for k, v in missing[:80]:
        print(f"  MISSING {k}: {v}")
    if len(missing) > 80: print(f"  ... {len(missing)-80} more")

if __name__ == "__main__":
    main()
