#!/usr/bin/env python
"""Run the harness copy against a given spec yaml and write the full result table as CSV (signature runs).

Usage: run_harness_copy.py --spec <yaml> --out <csv>. Read-only on every input except --out.
"""
from __future__ import annotations
import argparse, csv, collections, importlib, sys
from pathlib import Path
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
ha = importlib.import_module("thesis_assertions_nearest")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", type=Path, default=HERE / "thesis_expected_values_nearest.yaml")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    ha.SPEC = a.spec
    ha._spec = lambda: yaml.safe_load(a.spec.read_text())
    rep = ha.run_all(verbose=False)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["key", "expected", "actual", "ok", "source"])
        for r in rep.results:
            w.writerow([r.key, r.expected, r.actual, r.ok, r.source])
    expected = tuple(f"{k}[" for k, v in yaml.safe_load(a.spec.read_text()).items() if isinstance(v, dict) and v.get("expected_to_fail"))
    print(rep.summary(expected))
    print(collections.Counter(f.key.split(".")[0].split("[")[0] for f in rep.failures).most_common())
    for f in rep.failures[:80]:
        print(f"  FAIL {f.key}: thesis {f.expected} recomputed {f.actual}")

if __name__ == "__main__":
    main()
