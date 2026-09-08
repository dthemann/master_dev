#!/usr/bin/env python
"""Phase 7 of the nearest-copy program: re-transcribe the moving values into the yaml COPY.

Step 1 (trace): run the harness copy with a position-tracking yaml loader so every Result knows the
(line, column) of the expected scalar it was compared against. Step 2 (patch): for every failing check with a
traced position, replace that scalar token in thesis_expected_values_nearest.yaml by the recomputed value,
formatted at the precision of the token it replaces (same decimals, or the same significant figures for values
below 0.01). Failures without a traceable token are written to a ledger for manual transcription. Step 3
(signature B): re-run the harness copy on the patched yaml. Only files inside the program folder are written;
the canonical yaml and harness are never read for writing. A backup of the pre-Phase-7 yaml copy is kept.
"""
from __future__ import annotations
import argparse, csv, importlib, re, shutil, sys
from pathlib import Path
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))  # repro_harness and friends live in Scripts/Analysis
ha = importlib.import_module("thesis_assertions_nearest")

class PInt(int):
    _pos = None
class PFloat(float):
    _pos = None
class PStr(str):
    _pos = None

class PosLoader(yaml.SafeLoader):
    pass

def _wrap(cls, base):
    def construct(loader, node):
        v = base(loader, node)
        try:
            w = cls(v)
        except Exception:
            return v
        w._pos = (node.start_mark.line, node.start_mark.column)
        return w
    return construct

PosLoader.add_constructor("tag:yaml.org,2002:int", _wrap(PInt, yaml.SafeLoader.construct_yaml_int))
PosLoader.add_constructor("tag:yaml.org,2002:float", _wrap(PFloat, yaml.SafeLoader.construct_yaml_float))
PosLoader.add_constructor("tag:yaml.org,2002:str", _wrap(PStr, yaml.SafeLoader.construct_yaml_str))

TOKEN = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")

def run_traced(spec_path: Path):
    ha.SPEC = spec_path
    ha._spec = lambda: yaml.load(spec_path.read_text(), Loader=PosLoader)
    positions = []
    orig_add = ha.Report.add
    def add(self, *a, **k):
        r = orig_add(self, *a, **k)
        exp = a[1] if len(a) > 1 else k.get("expected")
        positions.append(getattr(exp, "_pos", None))
        return r
    ha.Report.add = add
    try:
        rep = ha.run_all(verbose=False)
    finally:
        ha.Report.add = orig_add
    return rep, positions

def fmt_like(token: str, actual) -> str | None:
    """Format `actual` at the precision of `token`; None when the value cannot be expressed."""
    s = str(actual)
    m = re.fullmatch(r"np\.float64\((.*)\)|np\.int64\((.*)\)", s)
    if m:
        s = m.group(1) or m.group(2)
    try:
        val = float(s)
    except ValueError:
        return None
    if "e" in token.lower():
        mant = token.lower().split("e")[0]
        d = len(mant.split(".")[1]) if "." in mant else 0
        return f"{val:.{d}e}"
    if "." not in token:
        return str(int(round(val))) if abs(val - round(val)) < 1e-9 else None
    d = len(token.split(".")[1])
    out = f"{val:.{d}f}"
    if abs(val) < 0.01 and val != 0 and float(out) == 0.0 or (abs(val) < 0.01 and val != 0 and abs(float(out) - val) / abs(val) > 0.02):
        sig = len(token.split(".")[1].lstrip("0")) or 1
        out = f"{val:.{sig}g}"
    return out

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--yaml", type=Path, default=HERE / "thesis_expected_values_nearest.yaml")
    ap.add_argument("--out-dir", type=Path, default=HERE / "harness_signatures")
    ap.add_argument("--apply", action="store_true", help="write the patched yaml (default: dry run, report only)")
    a = ap.parse_args()
    a.out_dir.mkdir(exist_ok=True)
    rep, pos = run_traced(a.yaml)
    rows = []
    for r, p in zip(rep.results, pos):
        rows.append(dict(key=r.key, expected=r.expected, actual=r.actual, ok=r.ok, source=r.source,
                         line=None if p is None else p[0], col=None if p is None else p[1]))
    with open(a.out_dir / "phase7_trace.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    fails = [r for r in rows if not r["ok"]]
    print(f"trace: {len(rows)} checks, {len(fails)} failing, {sum(r['line'] is not None for r in fails)} with a yaml position")
    lines = a.yaml.read_text().split("\n")
    patches, ledger = {}, []
    for r in fails:
        if r["line"] is None:
            ledger.append({**r, "reason": "no traced position (expected value is derived, not a yaml scalar)"}); continue
        ln, col = r["line"], r["col"]
        m = TOKEN.match(lines[ln], col)
        if not m:
            ledger.append({**r, "reason": f"no numeric token at {ln+1}:{col+1}: {lines[ln][col:col+20]!r}"}); continue
        tok = m.group(0)
        try:
            if abs(float(tok) - float(str(r["expected"]))) > 1e-12:
                ledger.append({**r, "reason": f"token {tok} != expected {r['expected']}"}); continue
        except ValueError:
            ledger.append({**r, "reason": f"expected {r['expected']!r} is not numeric"}); continue
        new = fmt_like(tok, r["actual"])
        if new is None:
            ledger.append({**r, "reason": f"cannot format actual {r['actual']!r} like token {tok}"}); continue
        patches.setdefault(ln, []).append((col, tok, new, r["key"]))
    # one yaml scalar can back several checks (the same cell compared under two keys): patch it once, and refuse
    # to proceed when two checks want different values for the same token
    for ln, items in list(patches.items()):
        by_col = {}
        for col, tok, new, key in items:
            if col in by_col and by_col[col][1] != new:
                sys.exit(f"CONFLICT at {ln+1}:{col+1}: {by_col[col][2]} wants {by_col[col][1]}, {key} wants {new}")
            by_col.setdefault(col, (tok, new, key))
        patches[ln] = [(col, tok, new, key) for col, (tok, new, key) in by_col.items()]
    n_patch = sum(len(v) for v in patches.values())
    print(f"patchable: {n_patch}; manual ledger: {len(ledger)}")
    with open(a.out_dir / "phase7_patch_ledger.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["key", "line", "col", "old", "new"])
        for ln in sorted(patches):
            for col, tok, new, key in sorted(patches[ln]):
                w.writerow([key, ln + 1, col + 1, tok, new])
    with open(a.out_dir / "phase7_manual_ledger.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ledger[0]) if ledger else ["key"]); w.writeheader(); w.writerows(ledger)
    if not a.apply:
        print("dry run; pass --apply to write the yaml"); return
    backup = a.yaml.with_suffix(".yaml.pre-phase7")
    if not backup.exists():
        shutil.copy2(a.yaml, backup)
    for ln, items in patches.items():
        line = lines[ln]
        for col, tok, new, _ in sorted(items, reverse=True):
            assert line[col:col + len(tok)] == tok
            line = line[:col] + new + line[col + len(tok):]
        lines[ln] = line
    a.yaml.write_text("\n".join(lines))
    print(f"patched {n_patch} tokens on {len(patches)} lines -> {a.yaml.name}; backup {backup.name}")
    # signature B: plain harness copy run on the patched yaml
    ha.SPEC = a.yaml
    ha._spec = lambda: yaml.safe_load(a.yaml.read_text())
    rep2 = ha.run_all(verbose=False)
    with open(a.out_dir / "harness_signature_B_newvalues_newtree.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["key", "expected", "actual", "ok", "source"])
        for r in rep2.results:
            w.writerow([r.key, r.expected, r.actual, r.ok, r.source])
    expected = tuple(f"{k}[" for k, v in yaml.safe_load(a.yaml.read_text()).items() if isinstance(v, dict) and v.get("expected_to_fail"))
    print("signature B:", rep2.summary(expected))
    for f in rep2.failures:
        print(f"  STILL FAILING {f.key}: thesis {f.expected} recomputed {f.actual}")

if __name__ == "__main__":
    main()
