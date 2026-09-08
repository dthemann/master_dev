#!/usr/bin/env python
"""Change signature C: the re-transcribed (nearest) yaml values against the OLD canonical tree.

Builds a temporary spec from thesis_expected_values_nearest.yaml with every `_nearest` result path mapped back
to its canonical counterpart (the inverse of make_harness_copies.PATH_MAP) and runs the harness copy on it.
Reads the canonical tree only; writes harness_signatures/harness_signature_C_newvalues_oldtree.csv and the
temporary spec. Expected outcome: the same blocks that moved in signature A fail again, in the opposite
direction, and every block that did not move passes.
"""
from __future__ import annotations
import csv, importlib, sys
from pathlib import Path
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))  # repro_harness and friends live in Scripts/Analysis
ha = importlib.import_module("thesis_assertions_nearest")
from make_harness_copies import PATH_MAP

def main():
    src = (HERE / "thesis_expected_values_nearest.yaml").read_text()
    n = 0
    for old, new in PATH_MAP:
        if "thesis_latex" in old:          # keep the source pointers on the copy; only result paths are inverted
            continue
        c = src.count(new); n += c
        src = src.replace(new, old)
    tmp = HERE / "harness_signatures" / "thesis_expected_values_nearest.oldtree-spec.yaml"
    tmp.write_text(src)
    print(f"inverted {n} path occurrences -> {tmp.name}")
    ha.SPEC = tmp
    ha._spec = lambda: yaml.safe_load(tmp.read_text())
    rep = ha.run_all(verbose=False)
    out = HERE / "harness_signatures" / "harness_signature_C_newvalues_oldtree.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["key", "expected", "actual", "ok", "source"])
        for r in rep.results:
            w.writerow([r.key, r.expected, r.actual, r.ok, r.source])
    fails = rep.failures
    print(f"signature C: {len(rep.results)} checks, {len(fails)} differ against the old tree")
    import collections
    print(collections.Counter(f.key.split(".")[0].split("[")[0] for f in fails).most_common())

if __name__ == "__main__":
    main()
