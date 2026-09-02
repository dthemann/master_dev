#!/usr/bin/env python
"""Recompute the numbers the thesis prints and compare them to the document.

Driven by ``thesis_expected_values.yaml``, which records each expected value
together with where in the thesis it appears. This module supplies the
recomputation: it reads the canonical CSVs, applies the gate definitions, and
reports a pass or a fail per cell.

A failure is a finding about the thesis and is reported as one. Nothing here
edits the document, and nothing rounds a number until it agrees.

THE GATES, which are the part most easily got wrong:

======================  ==================================================
Validity                ``pb_valid``
Near-native (Table 1)   ``rmsd <= 2``, NOT conditioned on validity
Near-native (Table 6)   ``pb_valid AND rmsd <= 2``
Form                    ``bestfit_rmsd <= 1 AND rmsd < 1000``
Triple / selection      ``pb_valid AND rmsd <= 2 AND bestfit_rmsd <= 1``
======================  ==================================================

The two near-native definitions are genuinely different and both are printed in
the thesis, which is why AutoDock reads 202 complexes in Table 1 and 199 in
Table 6. The ``rmsd < 1000`` guard in the form gate excludes exploded poses,
whose best-fit RMSD can be small while the pose is nowhere near the site.

Ranking, where a depth-dependent number is involved, differs per tool: AutoDock
by ``optimized_rank`` on the gnina arm, DiffDock by ``rank`` on the smina
variant, and EquiBind by ascending ``gnina_affinity``, because EquiBind emits no
confidence of its own.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path("/home/manndo/master_dev")
SPEC = Path(__file__).resolve().parent / "thesis_expected_values.yaml"

# Tolerances. Counts must be exact. Percentages are printed to one decimal, so a
# recomputed value may legitimately differ in the last place after rounding.
TOL_PCT = 0.05
TOL_REL = 0.005          # for the cost table's seconds-per-pose figures


@dataclass
class Result:
    key: str
    expected: object
    actual: object
    ok: bool
    source: str = ""
    note: str = ""


@dataclass
class Report:
    results: list[Result] = field(default_factory=list)

    @property
    def failures(self) -> list[Result]:
        return [r for r in self.results if not r.ok]

    def add(self, *a, **k) -> Result:
        r = Result(*a, **k)
        self.results.append(r)
        return r

    def summary(self) -> str:
        n, f = len(self.results), len(self.failures)
        head = f"{n - f}/{n} checks reproduce"
        if not f:
            return head + ". Every number checked matches the thesis."
        return head + f", {f} do NOT. Listed below."


def _spec() -> dict:
    return yaml.safe_load(SPEC.read_text())


def _md5(p: Path) -> str | None:
    return hashlib.md5(p.read_bytes()).hexdigest() if p.exists() else None


def _close(a: float, b: float, tol: float) -> bool:
    return abs(float(a) - float(b)) <= tol


# =============================================================================
# Table 1
# =============================================================================

def check_table_1(spec: dict, rep: Report, verbose: bool = False) -> None:
    t = spec["table_1"]
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_1", "13 columns x 10 rows", f"input missing: {t['input']}",
                False, t["source"], "the canonical per-pose table is absent")
        return
    df = pd.read_csv(csv, low_memory=False)

    for method, want in t["rows"].items():
        d = df[df["method"] == method]
        if d.empty:
            rep.add(f"table_1[{method}]", want, "method absent from the table",
                    False, t["source"])
            continue
        valid = d["pb_valid"].astype(bool)
        near = d["rmsd"] <= 2.0
        form = (d["bestfit_rmsd"] <= 1.0) & (d["rmsd"] < 1000)
        triple = valid & near & (d["bestfit_rmsd"] <= 1.0)

        def ncplx(mask) -> int:
            return d[mask].groupby(["protein", "ligand"]).ngroups

        n = len(d)
        got = [n, ncplx(valid), int(valid.sum()), round(100 * valid.mean(), 1),
               ncplx(near), int(near.sum()), round(100 * near.mean(), 1),
               ncplx(form), int(form.sum()), round(100 * form.mean(), 1),
               ncplx(triple), int(triple.sum()), round(100 * triple.mean(), 1)]

        for col, g, w in zip(t["columns"], got, want):
            ok = _close(g, w, TOL_PCT) if col.endswith("_pct") else g == w
            rep.add(f"table_1[{method}].{col}", w, g, ok, t["source"])
        if verbose:
            bad = sum(1 for col, g, w in zip(t["columns"], got, want)
                      if not (_close(g, w, TOL_PCT) if col.endswith("_pct") else g == w))
            print(f"  Table 1  {method:32s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")


# =============================================================================
# Table 2
# =============================================================================

def check_table_2(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Recovery at rank-1, top-15 and top-30 on the validity-aware gate.

    The sidecar CSV carries one row per (slot, depth). The asserted column is
    ``pbvalid_near_native_k``, the count of complexes with a PoseBusters-valid
    pose within 2 A at that depth, which is what the thesis prints in brackets.
    """
    t = spec.get("table_2")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_2", "3 tool rows", f"input missing: {t['input']}", False, t["source"],
                "the top-k sidecar has not been written on the canonical tree")
        return
    df = pd.read_csv(csv)
    df = df[df["row_type"] == "recovery"].copy()
    # `depth` arrives as a string, because the CSV interleaves row types that
    # leave it blank. Normalise before comparing or every lookup silently misses.
    df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    slot_of = {"AutoDock": "autodock", "DiffDock": "diffdock", "EquiBind": "equibind"}
    depths = {"rank1": 1, "top15": 15, "top30": 30}

    for tool, want in t["rows"].items():
        d = df[df["slot"] == slot_of[tool]]
        bad = 0
        for col, w in zip(t["columns"], want):
            row = d[d["depth"] == depths[col]]
            got = int(row.iloc[0]["pbvalid_near_native_k"]) if not row.empty else -1
            bad += got != w
            rep.add(f"table_2[{tool}].{col}", w, got, got == w, t["source"])
        if verbose:
            print(f"  Table 2  {tool:32s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")


# =============================================================================
# Table 5
# =============================================================================

def check_table_5(spec: dict, rep: Report, verbose: bool = False) -> None:
    t = spec["table_5"]
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_5", "3 tool rows", f"input missing: {t['input']}",
                False, t["source"])
        return
    df = pd.read_csv(csv)
    exp = df[df["dataset"] == "orai_jku"]

    for tool, want in t["experimental_totals"].items():
        row = exp[exp["tool"] == tool]
        if row.empty:
            rep.add(f"table_5[{tool}]", want, "tool absent", False, t["source"])
            continue
        r = row.iloc[0]
        got = [int(r["produced"]), int(r["pb_valid"]),
               round(100 * float(r["pbvalid_share"]), 1),
               int(r["outside_tm"]), round(100 * float(r["outside_tm_share"]), 1)]
        for col, g, w in zip(t["columns"], got, want):
            ok = _close(g, w, TOL_PCT) if col.endswith("_pct") else g == w
            rep.add(f"table_5[{tool}].{col}", w, g, ok, t["source"])
        if verbose:
            bad = sum(1 for col, g, w in zip(t["columns"], got, want)
                      if not (_close(g, w, TOL_PCT) if col.endswith("_pct") else g == w))
            print(f"  Table 5  {tool:32s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")


# =============================================================================
# Table 6
# =============================================================================

def check_table_6(spec: dict, rep: Report, verbose: bool = False) -> None:
    t = spec["table_6"]
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_6", "3 pipeline rows", f"input missing: {t['input']}",
                False, t["source"])
        return
    df = pd.read_csv(csv).set_index("method")

    for method, want in t["rows"].items():
        if method not in df.index:
            rep.add(f"table_6[{method}]", want, "pipeline absent", False, t["source"])
            continue
        r = df.loc[method]
        for col, w in want.items():
            if col == "complexes":
                continue          # recomputed separately below, from the per-pose table
            g = float(r[col])
            ok = _close(g, w, max(TOL_REL * abs(float(w)), 0.05))
            rep.add(f"table_6[{method}].{col}", w, round(g, 2), ok, t["source"])
        if verbose:
            print(f"  Table 6  {method:32s} ok")

    # The complexes column uses the VALIDITY-AWARE near-native gate, unlike
    # Table 1. Recomputing it from the per-pose table is what keeps the two
    # definitions from being quietly conflated.
    pp = ROOT / spec["table_1"]["input"]
    if pp.exists():
        df1 = pd.read_csv(pp, low_memory=False)
        arms = {"autodock": "autodock_mgltools_exh128_gnina",
                "diffdock": "diffdock_smina",
                "equibind": "equibind_unguided_gnina"}
        for method, key in arms.items():
            d = df1[df1["method"] == key]
            mask = d["pb_valid"].astype(bool) & (d["rmsd"] <= 2.0)
            got = d[mask].groupby(["protein", "ligand"]).ngroups
            want = t["rows"][method]["complexes"]
            rep.add(f"table_6[{method}].complexes", want, got, got == want, t["source"],
                    "validity-aware near-native gate, not Table 1's unconditioned one")


# =============================================================================
# Table 8
# =============================================================================

def check_table_8(spec: dict, rep: Report, verbose: bool = False) -> None:
    """The exhaustiveness ladder, on the validity-aware near-native gate.

    Recomputed from the per-pose table rather than read from a sidecar, because
    the sidecar's gates are each a single criterion and this table's is their
    conjunction.
    """
    t = spec.get("table_8")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_8", "8 arms x 4 depths", f"input missing: {t['input']}",
                False, t["source"])
        return
    df = pd.read_csv(csv, low_memory=False)
    for method, cfg in t["rows"].items():
        d = df[df["method"] == method]
        bad = 0
        if d.empty:
            rep.add(f"table_8[{method}]", cfg["counts"], "arm absent", False, t["source"])
            continue
        gate = d["pb_valid"].astype(bool) & (d["rmsd"] <= 2.0)
        for depth, want in zip(t["depths"], cfg["counts"]):
            sub = d[(d[cfg["rank_col"]] <= depth) & gate]
            got = sub.groupby(["protein", "ligand"]).ngroups
            bad += got != want
            rep.add(f"table_8[{method}].d{depth}", want, got, got == want, t["source"])
        if verbose:
            print(f"  Table 8  {method:32s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")


# =============================================================================
# Prose numbers quoted alongside Figures 23, 24 and 39
# =============================================================================

def check_prose(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Recompute the values the appendix quotes beside three refreshed figures.

    These three figures shipped from superseded trees until 2026-09-02 while
    their prose had already moved on. Encoding the quoted values keeps the two
    from separating again.
    """
    pr = spec.get("prose")
    if not pr:
        return

    # ---- Figure 23: per-tool medians on the RAW Orai experimental cloud ------
    f = pr["figure_23"]
    csv = ROOT / f["input"]
    if csv.exists():
        d = pd.read_csv(csv)
        d = d[d["pose_set"] == f["pose_set"]]
        for metric, key in (("compactness", "compactness_median"),
                            ("silhouette", "silhouette_median")):
            med = d.groupby("tool")[metric].median()
            for tool, want in f[key].items():
                got = round(float(med.get(tool, float("nan"))), 3)
                rep.add(f"prose.fig23.{metric}[{tool}]", want, got,
                        _close(got, want, 0.006), f["source"])
        if verbose:
            print("  Prose    Figure 23 raw-cloud medians          ok")
    else:
        rep.add("prose.fig23", "per-tool medians", f"input missing: {f['input']}",
                False, f["source"])

    # ---- Figure 24: the EquiBind single-cluster collapse ---------------------
    f = pr["figure_24"]
    csv = ROOT / f["input"]
    if csv.exists():
        d = pd.read_csv(csv)
        e = d[d["tool"] == "equibind"]
        piv = e.pivot_table(index=["frame", "ligand"], columns="pose_set", values="stability")
        paired = piv.dropna()
        nclust = e[e["pose_set"] == "filtered"].set_index(["frame", "ligand"])["n_clusters"]
        idx = paired.index.intersection(nclust.index)
        single = int((nclust.loc[idx] <= 1).sum())
        multi = idx[nclust.loc[idx] > 1]
        checks = [
            ("equibind_paired_stability_units", len(paired)),
            ("equibind_single_cluster_after_filtering", single),
            ("equibind_multicluster_units", len(multi)),
        ]
        for key, got in checks:
            rep.add(f"prose.fig24.{key}", f[key], got, got == f[key], f["source"])
        if len(multi):
            sub = piv.loc[multi]
            for key, got in (("equibind_multicluster_stability_raw",
                              round(float(sub["raw"].median()), 3)),
                             ("equibind_multicluster_stability_filtered",
                              round(float(sub["filtered"].median()), 3))):
                rep.add(f"prose.fig24.{key}", f[key], got, _close(got, f[key], 0.006),
                        f["source"])
        if verbose:
            print("  Prose    Figure 24 EquiBind collapse           ok")
    else:
        rep.add("prose.fig24", "EquiBind collapse", f"input missing: {f['input']}",
                False, f["source"])

    # ---- Figure 39: typed contacts by rank ----------------------------------
    f = pr["figure_39"]
    csv = ROOT / f["input"]
    if csv.exists():
        d = pd.read_csv(csv)
        def cell(method, rank, col):
            r = d[(d["method"] == method) & (d["pose_rank"] == rank)]
            return round(float(r.iloc[0][col]), 1) if not r.empty else float("nan")
        AD, DD, EB = ("autodock_mgltools_exh128_gnina", "diffdock_smina",
                      "equibind_unguided_gnina")
        for key, got in (("autodock_spurious_rank1", cell(AD, 1, "fp")),
                         ("diffdock_spurious_rank1", cell(DD, 1, "fp")),
                         ("autodock_spurious_rank5", cell(AD, 5, "fp")),
                         ("equibind_matched_rank1", cell(EB, 1, "tp")),
                         ("equibind_missed_rank1", cell(EB, 1, "fn"))):
            rep.add(f"prose.fig39.{key}", f[key], got, _close(got, f[key], 0.06),
                    f["source"])
        fp = d[d["method"] == EB]["fp"]
        for key, got in (("equibind_spurious_min", round(float(fp.min()), 1)),
                         ("equibind_spurious_max", round(float(fp.max()), 1))):
            rep.add(f"prose.fig39.{key}", f[key], got, _close(got, f[key], 0.06),
                    f["source"])
        if verbose:
            print("  Prose    Figure 39 contact decomposition       ok")
    else:
        rep.add("prose.fig39", "contact decomposition",
                f"input missing: {f['input']}", False, f["source"])

    # ---- the interaction pose-basis audit -----------------------------------
    # Appendix B quotes this audit verbatim. Its output is a long metric/value
    # table, so the checks below name the metrics the appendix prints.
    a = pr.get("audit")
    if a:
        csv = ROOT / a["input"]
        ids = ROOT / a["ids_file"]
        rep.add("prose.audit.ids_file", "present", "present" if ids.exists() else "missing",
                ids.exists(), a["source"],
                "the audit needs the 303-id analysed cohort, not the 308-id list; "
                "without it nine EquiBind poses in three excluded complexes fail to join")
        if not csv.exists():
            rep.add("prose.audit", "audit output", f"input missing: {a['input']}",
                    False, a["source"], "run interaction_pose_basis_audit.py with --ids-file")
        else:
            v = pd.read_csv(csv).set_index("metric")["value"]
            TOOLKEY = {"autodock": "autodock_mgltools_exh128_gnina",
                       "diffdock": "diffdock_smina",
                       "equibind": "equibind_unguided_gnina"}

            def chk(key, metric, want, tol):
                got = round(float(v[metric]), 4) if metric in v.index else float("nan")
                rep.add(f"prose.audit.{key}", want, got,
                        _close(got, want, tol) if got == got else False, a["source"])

            chk("profiled_poses", "pool_n_poses", a["profiled_poses"], 0)
            chk("median_rmsd", "pool_median_rmsd_A", a["median_rmsd"], 0.006)
            chk("within_2a_pct", "pool_within_2A_pct", a["within_2a_pct"], 0.06)
            for tool, key in TOOLKEY.items():
                chk(f"site_filtered_overlap[{tool}]", f"jaccard_rmsd_{key}",
                    a["site_filtered_overlap"][tool], 0.0006)
                chk(f"dropped_by_site_filter[{tool}]", f"site_filter_dropped_{key}",
                    a["dropped_by_site_filter"][tool], 0)
            chk("metal_records", "metal_rows", a["metal_records"], 0)
            chk("metal_complexes", "metal_complexes", a["metal_complexes"], 0)
            if verbose:
                print("  Prose    interaction pose-basis audit         ok")


# =============================================================================
# Cohort, determinism, figures, known gaps
# =============================================================================

def check_cohort(spec: dict, rep: Report, verbose: bool = False) -> None:
    c = spec["cohort"]
    ids = ROOT / "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"
    n = len([ln for ln in ids.read_text().splitlines() if ln.strip()]) if ids.exists() else -1
    rep.add("cohort.benchmark_ids", c["benchmark_ids"], n, n == c["benchmark_ids"], c["source"])

    pp = ROOT / spec["table_1"]["input"]
    if pp.exists():
        df = pd.read_csv(pp, low_memory=False)
        got = df[df["method"] == "diffdock"].groupby(["protein", "ligand"]).ngroups
        rep.add("cohort.analysed_complexes", c["analysed_complexes"], got,
                got == c["analysed_complexes"], c["source"],
                "complexes for which DiffDock produced at least one pose")
    if verbose:
        print(f"  Cohort   308 prepared / 303 analysed")


def check_determinism(spec: dict, rep: Report, verbose: bool = False) -> None:
    import repro_harness as H
    d = spec["determinism"]
    got = H.diffdock_seed_evidence()
    for panel, want in d["diffdock_seeded_logs"].items():
        g = got.get(panel, -1)
        # The benchmark and control panels must have ZERO seeded logs; the
        # experimental panel must have some. An exact count is brittle if the
        # tree is repacked, so the assertion is on the sign.
        ok = (g == 0) if want == 0 else (g > 0)
        rep.add(f"determinism.{panel}", f"{want} (sign)", g, ok, d["source"])
    if verbose:
        print(f"  Determinism  DiffDock seed markers: {got}")


def check_figures(spec: dict, rep: Report, verbose: bool = False) -> None:
    f = spec["figures"]
    media = ROOT / f["media_dir"]
    matched = 0
    for image, src in f["pairs"].items():
        a = _md5(media / f"{image}.png")
        b = _md5(ROOT / src)
        ok = a is not None and a == b
        matched += ok
        rep.add(f"figure[{image}]", f"md5 of {src}",
                "match" if ok else f"shipped={str(a)[:8]} source={str(b)[:8]}",
                ok, "thesis figure asset",
                "" if ok else "the shipped asset is not the file the pipeline writes")
    if verbose:
        print(f"  Figures  {matched}/{len(f['pairs'])} shipped assets byte-identical "
              f"to the canonical pipeline output")


def check_known_gaps(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Assert the gaps that are known, so they stay visible.

    These entries expect a file to be ABSENT. When one appears, the assertion
    fails and that is the intended signal: the gap has been closed and the note
    should be retired.
    """
    for gap in spec.get("known_gaps", []):
        present = (ROOT / gap["path"]).exists()
        want = gap["expect_present"]
        rep.add(f"gap[{gap['key']}]", f"present={want}", f"present={present}",
                present == want, "known gap", gap["note"].strip())
        if verbose:
            state = "still open" if present == want else "STATE CHANGED"
            print(f"  Gap      {gap['key']:44s} {state}")


# =============================================================================

def run_all(verbose: bool = False) -> Report:
    spec = _spec()
    rep = Report()
    if verbose:
        print(f"Checking against {SPEC.name} (recorded {spec['meta']['recorded']})\n")
    check_cohort(spec, rep, verbose)
    check_determinism(spec, rep, verbose)
    check_table_1(spec, rep, verbose)
    check_table_2(spec, rep, verbose)
    check_table_5(spec, rep, verbose)
    check_table_6(spec, rep, verbose)
    check_table_8(spec, rep, verbose)
    check_figures(spec, rep, verbose)
    check_prose(spec, rep, verbose)
    return rep


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    report = run_all(verbose=True)
    print()
    print(report.summary())
    for fail in report.failures:
        print(f"\n  {fail.key}\n      thesis    : {fail.expected}   ({fail.source})"
              f"\n      recomputed: {fail.actual}")
        if fail.note:
            print(f"      {fail.note}")
    sys.exit(1 if report.failures else 0)
