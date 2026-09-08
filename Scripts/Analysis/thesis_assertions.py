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
import sys
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

    def summary(self, expected_fail_prefixes: tuple[str, ...] = ()) -> str:
        n, fails = len(self.results), self.failures
        expected = [r for r in fails
                    if any(r.key.startswith(p) for p in expected_fail_prefixes)]
        unexpected = [r for r in fails if r not in expected]
        head = f"{n - len(fails)}/{n} checks reproduce"
        if not fails:
            return head + ". Every number checked matches the thesis."
        parts = [head]
        if expected:
            parts.append(f"{len(expected)} are RECORDED DEFECTS the spec expects to fail")
        if unexpected:
            parts.append(f"{len(unexpected)} are NOT expected and need attention")
        return ", ".join(parts) + ". Listed below."


def _spec() -> dict:
    return yaml.safe_load(SPEC.read_text())


def _md5(p: Path) -> str | None:
    return hashlib.md5(p.read_bytes()).hexdigest() if p.exists() else None


def _close(a: float, b: float, tol: float) -> bool:
    # 1e-9 hair: a value printed from x.x5 sits exactly on the window edge in binary
    return abs(float(a) - float(b)) <= tol + 1e-9


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
# Tables 18 and 19
# =============================================================================

def check_tables_18_19(spec: dict, rep: Report, verbose: bool = False) -> None:
    """The top-k recovery table and the within-tool optimisation contrasts.

    Both were extended on 2026-09-02: Table 18 gained AutoDock Vina (raw) and
    DiffDock + smina, Table 19 gained the AutoDock rescoring and DiffDock smina
    contrasts. Table 19's Holm family was split from the between-tool one at the
    same time, so `family_size` is asserted as well; a silent change there would
    move every star in the table.
    """
    t18 = spec.get("table_18")
    if t18:
        csv = ROOT / t18["input"]
        if not csv.exists():
            rep.add("table_18", "7 variants", f"input missing: {t18['input']}",
                    False, t18["source"])
        else:
            import math
            d = pd.read_csv(csv)
            for variant, cells in t18["rows"].items():
                bad = 0
                for k, (want_near, want_valid) in zip(t18["depths"], cells):
                    row = d[(d["variant"] == variant) & (d["k"] == k)]
                    if row.empty:
                        rep.add(f"table_18[{variant}].k{k}", [want_near, want_valid],
                                "variant absent", False, t18["source"])
                        bad += 1
                        continue
                    r = row.iloc[0]
                    n = int(r["n_complexes"])
                    got = [round(100 * r["near_k"] / n, 1), round(100 * r["valid_k"] / n, 1)]
                    for lbl, g, w in (("near", got[0], want_near), ("valid", got[1], want_valid)):
                        ok = _close(g, w, TOL_PCT)
                        bad += not ok
                        rep.add(f"table_18[{variant}].k{k}.{lbl}", w, g, ok, t18["source"])
                if verbose:
                    print(f"  Table 18 {variant:24s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")

    t19 = spec.get("table_19")
    if t19:
        csv = ROOT / t19["input"]
        if not csv.exists():
            rep.add("table_19", "4 contrasts", f"input missing: {t19['input']}",
                    False, t19["source"])
            return
        d = pd.read_csv(csv)
        if "family" not in d.columns:
            rep.add("table_19.family", "family column present", "absent", False,
                    t19["source"], "the split Holm families are not recorded in the CSV")
            return
        d = d[d["family"] == t19["family"]]
        rep.add("table_19.family_size", t19["family_size"], int(d["family_size"].iloc[0]),
                int(d["family_size"].iloc[0]) == t19["family_size"], t19["source"],
                "a changed family size moves every star in the table")
        for label, cells in t19["rows"].items():
            a, b = [x.strip() for x in label.split("->")]
            bad = 0
            for k, want in zip(t19["depths"], cells):
                row = d[(d["baseline"] == a) & (d["comparison"] == b) & (d["k"] == k)]
                if row.empty:
                    rep.add(f"table_19[{label}].k{k}", want, "contrast absent", False,
                            t19["source"])
                    bad += 1
                    continue
                r = row.iloc[0]
                got = [round(r["baseline_rate_%"], 1), round(r["comparison_rate_%"], 1),
                       int(r["baseline_only_wins"]), int(r["comparison_only_wins"])]
                for i, lbl in enumerate(("raw", "opt", "raw_wins", "opt_wins")):
                    ok = (_close(got[i], want[i], TOL_PCT) if i < 2 else got[i] == want[i])
                    bad += not ok
                    rep.add(f"table_19[{label}].k{k}.{lbl}", want[i], got[i], ok, t19["source"])
            if verbose:
                print(f"  Table 19 {label:46s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")


# =============================================================================
# Table 21
# =============================================================================

def check_table_21(spec: dict, rep: Report, verbose: bool = False) -> None:
    """The cross-tool recovery contrasts and the Holm p-values quoted in the note.

    This table was the one appendix float with no assertion behind it, and that is
    exactly where the pipeline had drifted from the thesis. The appendix declares
    the family is computed on the gnina-rescored AutoDock arm and prints that arm's
    rates; the script tested raw Vina instead. The rates, the discordant counts and
    the corrected p-values are all pinned here, because the note quotes all three.
    """
    t21 = spec.get("table_21")
    if not t21:
        return
    csv = ROOT / t21["input"]
    if not csv.exists():
        rep.add("table_21", "4 contrasts", f"input missing: {t21['input']}",
                False, t21["source"])
        return
    d = pd.read_csv(csv)
    if "family" not in d.columns:
        rep.add("table_21.family", "family column present", "absent", False,
                t21["source"], "the split Holm families are not recorded in the CSV")
        return
    d = d[d["family"] == t21["family"]]
    if d.empty:
        rep.add("table_21.family", t21["family"], "no rows in that family", False,
                t21["source"])
        return
    rep.add("table_21.family_size", t21["family_size"], int(d["family_size"].iloc[0]),
            int(d["family_size"].iloc[0]) == t21["family_size"], t21["source"],
            "a changed family size moves every corrected p in the note")
    for label, cells in t21["rows"].items():
        a, b = [x.strip() for x in label.split("->")]
        bad = 0
        for k, want in zip(t21["depths"], cells):
            row = d[(d["baseline"] == a) & (d["comparison"] == b) & (d["k"] == k)]
            if row.empty:
                rep.add(f"table_21[{label}].k{k}", want, "contrast absent", False,
                        t21["source"],
                        "the cross-tool baseline no longer matches the declared arm")
                bad += 1
                continue
            r = row.iloc[0]
            got = [round(r["baseline_rate_%"], 1), round(r["comparison_rate_%"], 1),
                   int(r["baseline_only_wins"]), int(r["comparison_only_wins"]),
                   float(f"{r['p_holm']:.4g}")]
            for i, lbl in enumerate(("base", "comp", "base_wins", "comp_wins", "p_holm")):
                if i < 2:
                    ok = _close(got[i], want[i], TOL_PCT)
                elif i < 4:
                    ok = got[i] == want[i]
                else:
                    ok = _close(got[i], want[i], abs(want[i]) * 0.02 + 1e-12)
                bad += not ok
                rep.add(f"table_21[{label}].k{k}.{lbl}", want[i], got[i], ok, t21["source"])
        if verbose:
            print(f"  Table 21 {label:46s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")


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

def check_prose_refiners(spec: dict, rep: Report, verbose: bool = False) -> None:
    """The Discussion's within-family refiner contrasts.

    Nothing in the pipeline computes these, so before this check the whole
    passage was unverifiable. It also carried a false blanket claim, that
    neither within-family difference was separable, which holds for DiffDock
    but for EquiBind only on the selection endpoint.

    RANKING TRAP: both EquiBind arms have rank == 999 on every row and each
    populates only its own refiner's affinity column, so each must be ranked by
    its own column. Ranking the smina arm by gnina_affinity gives an empty frame
    and a silently meaningless comparison.
    """
    t = spec.get("prose", {}).get("refiner_contrasts")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("prose[refiner_contrasts]", "5 contrasts",
                f"input missing: {t['input']}", False, t["source"])
        return
    import numpy as np
    from stats_utils import mcnemar_exact
    pp = pd.read_csv(csv, low_memory=False)
    rank_by = {"equibind_unguided_smina": "smina_affinity",
               "equibind_unguided_gnina": "gnina_affinity"}
    cx = sorted(set(map(tuple, pp[pp["method"] == "diffdock"][["protein", "ligand"]]
                        .drop_duplicates().values)))

    def vec(method, k, gate):
        d = pp[pp["method"] == method].copy()
        col = rank_by.get(method)
        if col:
            d = d[d[col].notna()]
            d["_r"] = d.groupby(["protein", "ligand"])[col].rank(method="first")
        else:
            d["_r"] = d["rank"]
        m = (d["_r"] <= k) & (d["rmsd"] <= 2.0)
        if gate in ("double", "triple"):
            m &= d["pb_valid"].astype(bool)
        if gate == "triple":
            m &= (d["bestfit_rmsd"] <= 1.0) & (d["rmsd"] < 1000)
        hit = set(map(tuple, d[m][["protein", "ligand"]].drop_duplicates().values))
        return np.array([1 if c in hit else 0 for c in cx])

    cases = [
        ("equibind_near_top15", "equibind_unguided_gnina", "equibind_unguided_smina",
         "near", "gnina", "smina"),
        ("equibind_selection_top15", "equibind_unguided_gnina", "equibind_unguided_smina",
         "triple", "gnina", "smina"),
        ("diffdock_triple_top15", "diffdock_smina", "diffdock_gnina",
         "triple", "smina", "gnina"),
        ("diffdock_double_top15", "diffdock_smina", "diffdock_gnina",
         "double", "smina", "gnina"),
        ("diffdock_near_top15", "diffdock_smina", "diffdock_gnina",
         "near", "smina", "gnina"),
    ]
    for key, ma, mb, gate, ka, kb in cases:
        want = t.get(key)
        if not want:
            continue
        a, b = vec(ma, 15, gate), vec(mb, 15, gate)
        _, _, pval = mcnemar_exact(a, b)
        got = {ka: int(a.sum()), kb: int(b.sum()), "mcnemar_p": float(f"{pval:.4g}")}
        bad = 0
        for fld in (ka, kb):
            ok = got[fld] == want[fld]
            bad += not ok
            rep.add(f"prose.refiners.{key}.{fld}", want[fld], got[fld], ok, t["source"])
        ok = _close(got["mcnemar_p"], want["mcnemar_p"],
                    max(abs(want["mcnemar_p"]) * 0.02, 1e-9))
        bad += not ok
        rep.add(f"prose.refiners.{key}.mcnemar_p", want["mcnemar_p"], got["mcnemar_p"],
                ok, t["source"])
        if verbose:
            print(f"  Prose    refiners {key:26s} {'ok' if not bad else f'{bad} DIFFER'}")


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

    # The two assets with no generator. Recorded, not hidden. Each caption
    # disclaims measurement in the document itself. The assertion is that they
    # are still absent from the tree: if one acquires a source, this fails and
    # the entry should be promoted into `pairs` rather than left as a note.
    for image, meta in f.get("unreproducible", {}).items():
        shipped = _md5(media / f"{image}.png")
        found = _find_by_md5(shipped, skip=("thesis_latex/media", "thesis_latex_nearest/media", "obsolete", ".backup"))
        rep.add(f"figure[{image}]", "no generator in the repository",
                "still none" if found is None else f"now produced by {found}",
                found is None, f"Figure {meta['figure']}, {meta['what']}",
                meta["why"].strip())
        if verbose:
            state = "no generator (as declared)" if found is None else f"NOW HAS ONE: {found}"
            print(f"  Figures  {image:8s} Figure {meta['figure']:<3d} {state}")

    # Completeness. Every image the document actually includes must be either
    # paired or declared. This is the assertion that makes the other two mean
    # something: without it, adding a figure to the thesis silently adds an
    # unchecked asset and the count above still reads as full coverage.
    used = _images_used_by_thesis(spec)
    accounted = set(f["pairs"]) | set(f.get("unreproducible", {}))
    unaccounted = sorted(used - accounted, key=_image_sort_key)
    rep.add("figures[coverage]", f"all {len(used)} included images accounted for",
            "complete" if not unaccounted else f"{len(unaccounted)} unaccounted: "
            + ", ".join(unaccounted),
            not unaccounted, "\\includegraphics in the two body files",
            "an image the thesis prints with neither a provenance pair nor a "
            "declared reason is an unchecked asset")
    if verbose:
        print(f"  Figures  {len(used)} included, {len(accounted & used)} accounted for"
              + ("" if not unaccounted else f", MISSING {', '.join(unaccounted)}"))


def _image_sort_key(name: str) -> int:
    digits = "".join(c for c in name if c.isdigit())
    return int(digits) if digits else 0


def _images_used_by_thesis(spec: dict) -> set[str]:
    """The image stems the two body files actually \\includegraphics.

    Read from the document rather than listed here, so the coverage assertion
    cannot drift away from what is printed.
    """
    import re
    used: set[str] = set()
    for key in ("thesis_main", "thesis_appendix"):
        p = ROOT / spec["meta"][key]
        if p.exists():
            used |= set(re.findall(r"media/media/(image\d+)\.png", p.read_text()))
    return used


def _find_by_md5(digest: str | None, skip: tuple[str, ...] = ()) -> str | None:
    """Repo-relative path of any PNG whose md5 equals `digest`, else None.

    Size is checked before hashing, so this walks the tree once cheaply. The
    `skip` fragments keep the shipped copy and the parked backups from counting
    as a generator for themselves.
    """
    import os
    if digest is None:
        return None
    # The shipped copy is the one file whose digest we already know, so its size
    # is the size any match must have. Comparing sizes first means the walk
    # hashes a handful of candidates instead of every PNG in the repository.
    target_size = next(
        (c.stat().st_size for c in (ROOT / "thesis_latex/media/media").glob("*.png")
         if _md5(c) == digest), None)
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in (".git", "__pycache__", "node_modules")]
        rel_dir = os.path.relpath(dirpath, ROOT)
        if any(s in rel_dir for s in skip):
            continue
        for fn in filenames:
            if not fn.lower().endswith(".png"):
                continue
            p = Path(dirpath) / fn
            try:
                if target_size is not None and p.stat().st_size != target_size:
                    continue
                if _md5(p) == digest:
                    return os.path.relpath(p, ROOT)
            except OSError:
                continue
    return None


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



def check_table_23(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Table 23: both reference conventions per arm and depth, from the sensitivity CSV."""
    t = spec.get("table_23")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_23", "3 arms x 4 depths x 8 cells", f"input missing: {t['input']}", False, t["source"])
        return
    df = pd.read_csv(csv)
    cols = [("recovered_nearest", "recovered_nearest"), ("recovered_instance", "recovered_instance"), ("gain", "gain"),
            ("form_nearest", "form_complexes_nearest"), ("form_instance", "form_complexes_instance"),
            ("combined_nearest", "combined_complexes_nearest"), ("combined_instance", "combined_complexes_instance")]
    for tool, rows in t["rows"].items():
        key = t["method_keys"][tool]
        bad = 0
        for depth, want in zip(t["depths"], rows):
            allr = df[(df.arm == key) & (df.depth == depth) & (df.stratum == "all")]
            sing = df[(df.arm == key) & (df.depth == depth) & (df.stratum == "single_copy")]
            if allr.empty or sing.empty:
                rep.add(f"table_23[{tool}].d{depth}", want, "row absent", False, t["source"])
                bad += 1
                continue
            got = [int(allr.iloc[0][c]) for _, c in cols] + [int(sing.iloc[0]["recovered_nearest"])]
            ok = got == [int(x) for x in want]
            bad += not ok
            rep.add(f"table_23[{tool}].d{depth}", want, got, ok, t["source"])
        if verbose:
            print(f"  Table 23 {tool:24s} {'ok' if not bad else f'{bad} ROWS DIFFER'}")


def check_convention(spec: dict, rep: Report, verbose: bool = False) -> None:
    """The per-pose table's reference convention, asserted by name."""
    c = spec.get("convention")
    if not c:
        return
    import json as _json
    man = ROOT / c["manifest"]
    if man.exists():
        m = _json.loads(man.read_text())
        sig = m.get("signature", {})
        rep.add("convention.manifest", c["reference_convention"], m.get("reference_convention"),
                m.get("reference_convention") == c["reference_convention"], c["source"])
        rep.add("convention.schema", c["schema"], sig.get("schema"), sig.get("schema") == c["schema"], c["source"])
    else:
        rep.add("convention.manifest", c["reference_convention"], f"missing: {c['manifest']}", False, c["source"])
    csv = ROOT / c["per_pose_csv"]
    if csv.exists():
        col = pd.read_csv(csv, usecols=["reference_convention"], low_memory=False)["reference_convention"]
        vals = sorted(col.dropna().astype(str).unique().tolist())
        rep.add("convention.per_pose_column", [c["reference_convention"]], vals, vals == [c["reference_convention"]], c["source"])
    summ = ROOT / c["summary"]
    if summ.exists():
        d = _json.loads(summ.read_text())
        rep.add("convention.multi_copy_303", c["multi_copy_303"], d.get("multi_copy_303"), d.get("multi_copy_303") == c["multi_copy_303"], c["source"])
        rep.add("convention.single_copy_303", c["single_copy_303"], d.get("single_copy_303"), d.get("single_copy_303") == c["single_copy_303"], c["source"])
        c308 = d.get("count_308") or {}
        got308 = c308.get("multi_copy") if isinstance(c308, dict) else None
        rep.add("convention.multi_copy_308", c["multi_copy_308"], got308, got308 == c["multi_copy_308"], c["source"])
    if verbose:
        print("  Convention  nearest-copy manifest, column and copy counts  checked")


def check_nearest_program(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Sidecars of the nearest-copy program that the appendix and the footnote :144 print."""
    p = (spec.get("prose") or {}).get("nearest_copy_program")
    if not p:
        return
    import json as _json
    src = p["source"]
    sb = ROOT / p["selection_bias"]
    if sb.exists():
        d = pd.read_csv(sb)
        top = d[(d.gate == "triple") & (d.depth == 15)].set_index("family")["bias_pp"]
        for fam, want in p["winners_curse_top15_triple_pp"].items():
            got = round(float(top[fam]), 1) if fam in top.index else None
            rep.add(f"prose.nearest_copy_program.winners_curse[{fam}]", want, got, got == want, src)
        lad = d[d.family.str.startswith("AutoDock ladder")]
        for key, (gate, depth) in {"rank1_triple": ("triple", 1), "top15_double": ("double", 15), "rank1_double": ("double", 1)}.items():
            row = lad[(lad.gate == gate) & (lad.depth == depth)]
            got = round(float(row.iloc[0]["bias_pp"]), 1) if len(row) else None
            want = p["winners_curse_ladder_pp"][key]
            rep.add(f"prose.nearest_copy_program.winners_curse_ladder.{key}", want, got, got == want, src)
    else:
        rep.add("prose.nearest_copy_program.winners_curse", "sidecar", f"missing: {p['selection_bias']}", False, src)
    ic = ROOT / p["identity_cost"]
    if ic.exists():
        cost = _json.loads(ic.read_text())["cost"]
        spec_ic = p["identity_cost_diffdock_smina"]
        for conv in ("nearest", "instance"):
            got = [cost.get(f"diffdock_smina|{conv}|{dpt}") for dpt in spec_ic["depths"]]
            rep.add(f"prose.nearest_copy_program.identity_cost.diffdock_smina.{conv}", spec_ic[conv], got, got == spec_ic[conv], src)
        for arm in p["identity_cost_zero_arms"]:
            got = sorted({v for k, v in cost.items() if k.startswith(arm + "|")})
            rep.add(f"prose.nearest_copy_program.identity_cost.{arm}", [0], got, got == [0], src)
    else:
        rep.add("prose.nearest_copy_program.identity_cost", "sidecar", f"missing: {p['identity_cost']}", False, src)
    it = ROOT / p["itt"]
    if it.exists():
        d = _json.loads(it.read_text())
        got = (d.get("conventions", {}).get("nearest", {}) or {}).get("dropped_recovered")
        rep.add("prose.nearest_copy_program.itt_dropped_recovered", p["itt_dropped_recovered"], got, got == p["itt_dropped_recovered"], src)
    else:
        rep.add("prose.nearest_copy_program.itt", "sidecar", f"missing: {p['itt']}", False, src)
    ac = ROOT / p["alternate_copies"]
    if ac.exists():
        d = _json.loads(ac.read_text())
        for k, want in p["alternate_copies_in_cube"].items():
            rep.add(f"prose.nearest_copy_program.alternate_copies.{k}", want, d.get(k), d.get(k) == want, src)
    else:
        rep.add("prose.nearest_copy_program.alternate_copies", "sidecar", f"missing: {p['alternate_copies']}", False, src)
    if verbose:
        print("  Prose    nearest-copy program sidecars           checked")


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
    check_table_3(spec, rep, verbose)
    check_table_4(spec, rep, verbose)
    check_table_5(spec, rep, verbose)
    check_table_6(spec, rep, verbose)
    check_table_8(spec, rep, verbose)
    check_table_9(spec, rep, verbose)
    check_table_10(spec, rep, verbose)
    check_table_11(spec, rep, verbose)
    check_table_12(spec, rep, verbose)
    check_table_14(spec, rep, verbose)
    check_table_16(spec, rep, verbose)
    check_table_17(spec, rep, verbose)
    check_tables_18_19(spec, rep, verbose)
    check_table_20(spec, rep, verbose)
    check_table_21(spec, rep, verbose)
    check_table_22(spec, rep, verbose)
    check_table_23(spec, rep, verbose)
    check_convention(spec, rep, verbose)
    check_nearest_program(spec, rep, verbose)
    check_table_24(spec, rep, verbose)
    check_table_25(spec, rep, verbose)
    check_table_26(spec, rep, verbose)
    check_table_27(spec, rep, verbose)
    check_figures(spec, rep, verbose)
    check_prose(spec, rep, verbose)
    check_prose_refiners(spec, rep, verbose)
    check_prose_appendix_h9(spec, rep, verbose)
    return rep



# =============================================================================
# Tables added 2026-09-03
#
# Six floats that the harness did not read. Table 14 was the sharpest case: no
# script anywhere in the repository produced it, so its twenty numbers lived
# only in the document. The others each had an artifact on disk that nothing
# compared against, which is the condition under which Table 21 drifted.
# =============================================================================

def check_table_3(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Crystal-cluster reach and co-reach, parsed from the figure's stats sidecar.

    The sidecar is the plotting run's own record of what it computed. Reading it
    is a stronger check than re-deriving the reach rates here would be, because a
    re-derivation could agree with the thesis while disagreeing with the figure
    printed beside it.
    """
    t = spec.get("table_3")
    if not t:
        return
    src = ROOT / t["input"]
    if not src.exists():
        rep.add("table_3", "7 statistics x 4 depths", f"input missing: {t['input']}",
                False, t["source"])
        return
    import re
    text = src.read_text()
    bad = 0

    for tool, want_rows in t["reach_pct"].items():
        m = re.search(rf"^{tool}\*?\s+.*?reach rate(.*)$", text, re.M)
        got = re.findall(r"(\d+)%\s+\[(\d+)%.(\d+)%\]", m.group(1)) if m else []
        for depth, want, cell in zip(t["depths"], want_rows, got):
            trip = [int(x) for x in cell]
            ok = trip == list(want)
            bad += not ok
            rep.add(f"table_3[{tool} reach].top{depth}", want, trip, ok, t["source"])
        if len(got) != len(t["depths"]):
            bad += 1
            rep.add(f"table_3[{tool} reach]", f"{len(t['depths'])} depths",
                    f"parsed {len(got)}", False, t["source"])

    for pair, want_rows in t["co_reach_phi"].items():
        a, b = pair.split("+")
        m = re.search(rf"^{a}\*?\s*\+\s*{b}\*?\s+.*?co-reach(.*)$", text, re.M)
        got = [float(x) for x in re.findall(r"[φf]=([+-]\d+\.\d+)", m.group(1))] if m else []
        for depth, want, val in zip(t["depths"], want_rows, got):
            ok = _close(val, want, 0.005)
            bad += not ok
            rep.add(f"table_3[{pair} phi].top{depth}", want, val, ok, t["source"])

    m = re.search(r"^All three reach.*$", text, re.M)
    obs_exp = re.findall(r"(\d+) obs / (\d+) exp", m.group(0)) if m else []
    for depth, wo, we, cell in zip(t["depths"], t["all_three_observed"],
                                   t["all_three_expected"], obs_exp):
        got = (int(cell[0]), int(cell[1]))
        ok = got == (wo, we)
        bad += not ok
        rep.add(f"table_3[all three].top{depth}", f"{wo}/{we}", f"{got[0]}/{got[1]}",
                ok, t["source"])

    if verbose:
        print(f"  Table 3  reach, co-reach phi and all-three counts   "
              f"{'ok' if not bad else f'{bad} CELLS DIFFER'}")


def check_table_4(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Native-interaction recovery by pose rank, all sixty cells.

    The complex count n is asserted per tool and per rank, not just the three
    rates. It is the evidence for the table's own caveat that the tools are
    compared on overlapping rather than identical cohorts, so a silent change in
    n would move every rate beneath it without any rate looking wrong.
    """
    t = spec.get("table_4")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_4", "3 tools x 5 ranks x 4 rows", f"input missing: {t['input']}",
                False, t["source"])
        return
    df = pd.read_csv(csv)
    for tool, want in t["rows"].items():
        key = t["method_keys"][tool]
        d = df[df["method"] == key].set_index("rank")
        bad = 0
        # Half of the last printed place, plus a hair for binary representation.
        HALF_ULP = 0.00051
        for field, col, tol in (("n", "n", 0), ("precision", "precision", HALF_ULP),
                                ("recall", "recall", HALF_ULP), ("f1", "f1", HALF_ULP)):
            for rank, exp in zip(t["ranks"], want[field]):
                if rank not in d.index:
                    bad += 1
                    rep.add(f"table_4[{tool}.{field}].k{rank}", exp, "rank absent",
                            False, t["source"])
                    continue
                got = d.loc[rank, col]
                # Compared unrounded. Rounding the stored value first and then
                # demanding equality rejects a cell that is correct: the CSV
                # carries 0.5035 and 0.4635, which the thesis prints as 0.504 and
                # 0.463, and Python's round() takes both the other way.
                ok = (int(got) == exp) if field == "n" else _close(float(got), exp, tol)
                bad += not ok
                rep.add(f"table_4[{tool}.{field}].k{rank}", exp,
                        int(got) if field == "n" else round(float(got), 4), ok, t["source"])
        if verbose:
            print(f"  Table 4  {tool:10s} n, precision, recall, F1   "
                  f"{'ok' if not bad else f'{bad} CELLS DIFFER'}")


def check_table_10(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Pooled cost per generated pose, plus the derived quantities in its note.

    The charged column is recomputed from the two hardware currencies on the
    single basis the Results chapter defines, rather than read from a stored
    field. That is the point of the check: the appendix table and the Results
    convention have to stay the same convention.
    """
    t = spec.get("table_10")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_10", "3 pipelines x 3 currencies", f"input missing: {t['input']}",
                False, t["source"])
        return
    df = pd.read_csv(csv).set_index("method")
    threads = t["threads"]
    for method, (w_charged, w_gpu, w_cpu) in t["rows"].items():
        r = df.loc[method]
        cpu, gpu = float(r["cpu_core_s_per_generated"]), float(r["gpu_s_per_generated"])
        charged = cpu / threads + gpu
        bad = 0
        for label, got, want in (("charged", charged, w_charged), ("gpu", gpu, w_gpu),
                                 ("cpu_core", cpu, w_cpu)):
            ok = _close(round(got, 2), want, 0.005)
            bad += not ok
            rep.add(f"table_10[{method}.{label}]", want, round(got, 3), ok, t["source"])
        if verbose:
            print(f"  Table 10 {method:10s} charged, GPU and CPU-core per pose   "
                  f"{'ok' if not bad else f'{bad} CELLS DIFFER'}")

    fn, ad = t["footnote"], df.loc["autodock"]
    search_h = float(ad["total_cpu_core_h"]) / threads
    charged_h = search_h + float(ad["total_gpu_h"])
    derived = {
        "autodock_poses_generated": int(ad["poses_generated"]),
        "autodock_cpu_core_h": round(float(ad["total_cpu_core_h"]), 2),
        "autodock_gpu_h": round(float(ad["total_gpu_h"]), 2),
        "autodock_search_wall_h": round(search_h, 2),
        "autodock_charged_h": round(charged_h, 2),
        "autodock_rescoring_share_pct": round(100 * float(ad["total_gpu_h"]) / charged_h, 1),
        "equibind_full_run_wall_h": round(float(df.loc["equibind", "total_wall_full_h"]), 2),
    }
    nbad = 0
    for k, want in fn.items():
        got = derived[k]
        ok = (got == want) if isinstance(want, int) else _close(got, want, 0.005)
        nbad += not ok
        rep.add(f"table_10[note.{k}]", want, got, ok, t["source"])
    if verbose:
        print(f"  Table 10 footnote derivations                       "
              f"{'ok' if not nbad else f'{nbad} DIFFER'}")


def check_table_14(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Receptor geometry across the four frames, measured from the PDBs.

    This float had no generator. The measurement is reconstructed here from the
    copy-only receptor inputs: the pore axis is the first principal component of
    the CA cloud, the origin is that cloud's centroid, and each diagnostic atom's
    perpendicular distance is averaged over the six subunits. The subunits are
    found by residue number and atom name because the chain column of these files
    was flattened to a single value and the six copies live in SEGID.
    """
    t = spec.get("table_14")
    if not t:
        return
    import numpy as np

    def measure(pdb: Path) -> tuple[dict, int]:
        names, seqs, xyz = [], [], []
        for line in pdb.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                names.append((line[12:16].strip(), line[17:20].strip()))
                seqs.append(int(line[22:26]))
                xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        xyz = np.asarray(xyz)
        ca = xyz[[i for i, (n, _) in enumerate(names) if n == "CA"]]
        centre = ca.mean(axis=0)
        axis = np.linalg.svd(ca - centre, full_matrices=False)[2][0]
        axis = axis / np.linalg.norm(axis)
        out, copies = {}, {}
        for key in t["rows"]:
            res, seq, atom = key.split()
            idx = [i for i, (n, r) in enumerate(names)
                   if n == atom and r == res and seqs[i] == int(seq)]
            v = xyz[idx] - centre
            perp = v - np.outer(v @ axis, axis)
            out[key] = float(np.linalg.norm(perp, axis=1).mean()) if idx else float("nan")
            copies[key] = len(idx)
        return out, copies

    measured, bad = {}, 0
    for frame, rel in t["receptors"].items():
        p = ROOT / rel
        if not p.exists():
            rep.add(f"table_14[{frame}]", "receptor present", f"missing: {rel}",
                    False, t["source"])
            return
        measured[frame], copies = measure(p)
        n_sub = min(copies.values())
        ok = n_sub == t["subunits"]
        bad += not ok
        rep.add(f"table_14[{frame}.subunits]", t["subunits"], n_sub, ok, t["source"],
                "the six subunits are what the table averages over")

    frames = list(t["receptors"])
    for key, want_row in t["rows"].items():
        for frame, want in zip(frames, want_row):
            got = round(measured[frame][key], 2)
            ok = _close(got, want, 0.005)
            bad += not ok
            rep.add(f"table_14[{key}].{frame}", want, got, ok, t["source"])
    if verbose:
        print(f"  Table 14 five atoms x four frames, measured from the PDBs   "
              f"{'ok' if not bad else f'{bad} CELLS DIFFER'}")

    pr = t["prose"]
    later = frames[1:]
    spreads = {
        "later_frame_spread_arg91": ("ARG 91 CZ",),
        "later_frame_spread_glu106": ("GLU 106 CD",),
        "later_frame_spread_asp114": ("ASP 114 CG",),
    }
    pbad = 0
    for field, (key,) in spreads.items():
        vals = [measured[f][key] for f in later]
        got = round(max(vals) - min(vals), 2)
        ok = _close(got, pr[field], 0.005)
        pbad += not ok
        rep.add(f"table_14[prose.{field}]", pr[field], got, ok, pr["source"])
    for field, keys in (("fr0_offset_aspartates",
                         ["ASP 110 CG", "ASP 112 CG", "ASP 114 CG"]),
                        ("fr0_offset_gate_filter", ["ARG 91 CZ", "GLU 106 CD"])):
        offs = [sum(measured[f][k] for f in later) / len(later) - measured["Fr0"][k]
                for k in keys]
        got = [round(min(offs), 1), round(max(offs), 1)]
        ok = all(_close(g, w, 0.05) for g, w in zip(got, pr[field]))
        pbad += not ok
        rep.add(f"table_14[prose.{field}]", pr[field], got, ok, pr["source"])
    if verbose:
        print(f"  Table 14 prose spreads and Fr0 offsets                      "
              f"{'ok' if not pbad else f'{pbad} DIFFER'}")


def check_table_20(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Fixed-rank optimisation gain by pose-rank band.

    Every pose stays at its own rank, so this is a pooled rate difference within
    each band and not a re-ranking. The validity rows are printed as whole points
    and two of them land near a half, so the comparison is against the unrounded
    gain: rounding first would accept either neighbouring integer.
    """
    t = spec.get("table_20")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_20", "6 rows x 6 bands", f"input missing: {t['input']}",
                False, t["source"])
        return
    df = pd.read_csv(csv)
    col_for = {"near-native": "near_%", "PoseBusters-valid": "pbv_%",
               "near-native and PB-valid": "both_%"}
    for row_key, want in t["rows"].items():
        tool, criterion = row_key.split("/")
        col = col_for[criterion]
        bad, i = 0, 0
        for band, (lo, hi) in t["bands"].items():
            for opt in ("smina", "gnina"):
                raw = df[(df.tool == tool) & (df.optimizer == "raw")
                         & df["rank"].between(lo, hi)]
                cur = df[(df.tool == tool) & (df.optimizer == opt)
                         & df["rank"].between(lo, hi)]
                if raw.empty or cur.empty:
                    rep.add(f"table_20[{row_key}].{band}/{opt}", want[i], "band absent",
                            False, t["source"])
                    i += 1
                    continue
                base = (raw[col] * raw.n).sum() / raw.n.sum()
                got = (cur[col] * cur.n).sum() / cur.n.sum() - base
                # Whole-point rows carry a 0.5 window, one-decimal rows 0.05.
                tol = 0.5 if float(want[i]).is_integer() and abs(want[i]) >= 10 else 0.05
                ok = abs(got - want[i]) <= tol + 1e-9
                bad += not ok
                rep.add(f"table_20[{row_key}].{band}/{opt}", want[i], round(got, 2),
                        ok, t["source"])
                i += 1
        if verbose:
            print(f"  Table 20 {row_key:36s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")


def check_table_22(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Accurate-but-invalid share by ranking depth, from counts over the 303.

    Derived from the same sidecar Table 18 reads, which carries near_k and
    valid_k as counts. That is the right source twice over. It is the definition
    the table states, the complex-level gap between Table 18's two columns, and it
    avoids re-deriving a per-tool ranking here: the raw EquiBind arm has no usable
    rank of its own and no gnina affinity to stand in for one, so a re-derivation
    silently returns an empty pool and every cell reads zero.

    Deliberately NOT computed by subtracting the printed cells of Table 18. Those
    are already rounded to one decimal, and differencing them moves three of the
    twenty entries by a tenth, which would then read as a defect in this table
    rather than as the rounding artefact it is.
    """
    t = spec.get("table_22")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_22", "5 variants x 4 depths", f"input missing: {t['input']}",
                False, t["source"])
        return
    try:
        from stats_utils import wilson_ci
    except ImportError as exc:                       # pragma: no cover
        rep.add("table_22", "5 variants x 4 depths", f"cannot import: {exc}",
                False, t["source"])
        return

    df = pd.read_csv(csv)
    rows = dict(t["rows"])
    rows[t["footnote_variant"]] = t["footnote_autodock_raw"]

    labels = t.get("variant_labels", {})
    for variant, want_rows in rows.items():
        d = df[df["variant"] == labels.get(variant, variant)].set_index("k")
        bad = 0
        for depth, want in zip(t["depths"], want_rows):
            if depth not in d.index:
                bad += 1
                rep.add(f"table_22[{variant}].k{depth}", list(want), "depth absent",
                        False, t["source"])
                continue
            r = d.loc[depth]
            n = int(r["n_complexes"])
            k = int(r["near_k"]) - int(r["valid_k"])
            lo, hi = (100 * x for x in wilson_ci(k, n))
            got = [round(100.0 * k / n, 1), round(lo, 1), round(hi, 1)]
            ok = n == t["n"] and all(_close(g, w, 0.05) for g, w in zip(got, want))
            bad += not ok
            rep.add(f"table_22[{variant}].k{depth}", list(want), got, ok, t["source"])
        if verbose:
            label = variant + (", footnote" if variant == t["footnote_variant"] else "")
            print(f"  Table 22 {label:34s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")


def _arm_sets(df, method, rank_col, depth):
    """(near-native, validity-aware) complex sets for one arm at one depth."""
    import thesis_endpoint_diagnostics as TED
    near = TED.recovered(df, method, depth, rank_col, False)
    valid = TED.recovered(df, method, depth, rank_col, True)
    return set(near[near].index), set(valid[valid].index)


def check_table_25(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Decomposition of the rank-1 to best-of-top-15 pass-all gain.

    The four rows are an identity, and asserting it is half the point: net has to
    equal gained minus stay-invalid plus rescued, or the decomposition is not a
    decomposition. Both the parts and the identity are checked, so a change that
    keeps the arithmetic while moving the parts still fails.
    """
    t = spec.get("table_25")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_25", "3 arms x 4 rows", f"input missing: {t['input']}",
                False, t["source"])
        return
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import thesis_endpoint_diagnostics as TED
    df = TED.load(str(csv))
    for arm, want in t["rows"].items():
        method, rank_col = t["arms"][arm]
        n1, v1 = _arm_sets(df, method, rank_col, 1)
        n15, v15 = _arm_sets(df, method, rank_col, 15)
        gained = n15 - n1
        got = [len(gained), len(gained - v15), len((n1 & v15) - v1), len(v15) - len(v1)]
        bad = 0
        for label, g, w in zip(("gained", "stay_invalid", "rescued", "net"), got, want):
            bad += g != w
            rep.add(f"table_25[{arm}].{label}", w, g, g == w, t["source"])
        # The arithmetic identity gained - stay_invalid + rescued == net is a
        # TAUTOLOGY given these four definitions and can never fail, so asserting
        # it would be decoration. What CAN fail is the nesting the identity rests
        # on: the pools are cumulative, so a complex recovered at rank-1 must
        # still be recovered at top-15, and validity-aware recovery is a strict
        # subset of near-native recovery. Either would break if a ranking column
        # or a gate changed, and the decomposition would stop being one.
        for label, sub, sup in (("V1 within V15", v1, v15),
                                ("N1 within N15", n1, n15),
                                ("V15 within N15", v15, n15),
                                ("V1 within N1", v1, n1)):
            leaked = sub - sup
            bad += bool(leaked)
            rep.add(f"table_25[{arm}].{label}", "subset holds",
                    "holds" if not leaked else f"{len(leaked)} complexes outside",
                    not leaked, t["source"],
                    "the pools are nested and validity-aware recovery is a subset "
                    "of near-native recovery; the decomposition assumes both")
        if verbose:
            print(f"  Table 25 {arm:24s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")


def check_table_26(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Near-nativeness and form recovery across ten thresholds and three depths.

    Both halves are validity-aware. The gate is pb_valid AND distance <= t, not a
    bare distance gate, which is what makes the 2 A in-place column reproduce
    Table 2. That equality is asserted at the end rather than assumed, because it
    is the one join between these two tables and nothing else would notice if it
    broke.
    """
    t = spec.get("table_26")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_26", "180 cells", f"input missing: {t['input']}", False, t["source"])
        return
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import thesis_endpoint_diagnostics as TED
    df = TED.load(str(csv))
    n = t["n"]
    got_2A = {}
    for metric, col in t["metrics"].items():
        for arm, per_depth in t["rows"][metric].items():
            method, rank_col = t["arms"][arm]
            x = df[df.method == method].copy()
            x["_rk"] = TED._effective_rank(x, rank_col)
            bad = 0
            for depth, want_row in per_depth.items():
                s = x[x._rk <= int(depth)]
                valid = s.pb_valid.astype("boolean").fillna(False)
                dist = pd.to_numeric(s[col], errors="coerce")
                for thr, want in zip(t["thresholds"], want_row):
                    hit = s.assign(ok=valid & (dist <= thr)).groupby("cid").ok.any()
                    pct = round(100.0 * int(hit.sum()) / n, 1)
                    ok = _close(pct, want, 0.05)
                    bad += not ok
                    rep.add(f"table_26[{metric}/{arm}].d{depth}.t{thr}", want, pct,
                            ok, t["source"])
                    if metric == "RMSD" and thr == 2:
                        got_2A[(arm, int(depth))] = int(hit.sum())
            if verbose:
                print(f"  Table 26 {metric:12s} {arm:24s} "
                      f"{'ok' if not bad else f'{bad} CELLS DIFFER'}")

    # The join to Table 2. Same gate, same depths, so the counts must agree.
    t2 = spec.get("table_2")
    if t2 and "cross_check" in t:
        short = {"AutoDock Vina + gnina": "AutoDock", "DiffDock + smina": "DiffDock",
                 "EquiBind + gnina": "EquiBind"}
        cols = {d: i for i, d in enumerate((1, 15, 30))}
        nbad = 0
        for arm, tool in short.items():
            for depth in t["cross_check"]["table_2_depths"]:
                want = t2["rows"][tool][cols[depth]]
                got = got_2A.get((arm, depth))
                ok = got == want
                nbad += not ok
                rep.add(f"table_26[join Table 2].{tool}.d{depth}", want, got, ok,
                        t["cross_check"]["note"].strip())
        if verbose:
            print(f"  Table 26 2 A in-place column joins Table 2      "
                  f"{'ok' if not nbad else f'{nbad} DIFFER'}")


def check_table_27(spec: dict, rep: Report, verbose: bool = False) -> None:
    """In-place-RMSD-versus-form distributions, from the three-depth sidecar.

    The input is pinned to the sidecar built at rank-1, top-5 and top-15. Its
    sibling filmstrip_stats__per_tool_depth.csv carries the same columns at depths
    1, 3 and 5, so reading that one would supply no top-15 row and would line its
    depth-3 values up against this table's top-5 without erroring.
    """
    t = spec.get("table_27")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_27", "9 rows x 13 columns", f"input missing: {t['input']}",
                False, t["source"])
        return
    df = pd.read_csv(csv)
    for key, want in t["rows"].items():
        tool, depth = key.split("/")
        d = df[(df.tool == tool) & (df.depth == int(depth))]
        if d.empty:
            rep.add(f"table_27[{key}]", want, "row absent", False, t["source"])
            continue
        r = d.iloc[0]
        bad = 0
        for col, w in zip(t["columns"], want):
            got = float(r[col])
            # Counts exact; medians and quartiles to three places; percentages to one.
            if isinstance(w, int) and col in ("n_poses", "n_valid_complexes"):
                ok = int(got) == w
                shown = int(got)
            else:
                places = 3 if "median" in col or "_q" in col or col == "median_r" else 1
                shown = round(got, places + 1)
                ok = _close(got, w, 0.5 * 10 ** -places + 1e-9)
            bad += not ok
            rep.add(f"table_27[{key}].{col}", w, shown, ok, t["source"])
        if verbose:
            print(f"  Table 27 {key:14s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")


def check_table_16(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Principal ligand descriptors of the benchmark, over all 308 prepared ligands.

    The denominator is 308 and not the 303 analysed. This table describes the
    library that was docked rather than the cohort that survived to the
    three-tool comparison, and asserting the count keeps the two apart.
    """
    t = spec.get("table_16")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("table_16", "10 descriptors x 5 statistics",
                f"input missing: {t['input']}", False, t["source"])
        return
    df = pd.read_csv(csv, index_col=0)
    cols = {"median": "50%", "q1": "25%", "q3": "75%", "min": "min", "max": "max"}
    for key, want in t["rows"].items():
        if key not in df.index:
            rep.add(f"table_16[{key}]", want, "descriptor absent", False, t["source"])
            continue
        r = df.loc[key]
        places = t["decimals"][key]
        tol = 0.5 * 10 ** -places + 1e-9
        bad = 0
        n = int(r["count"])
        ok_n = n == t["n"]
        bad += not ok_n
        rep.add(f"table_16[{key}].n", t["n"], n, ok_n, t["source"])
        for (label, col), w in zip(cols.items(), want):
            got = float(r[col])
            ok = _close(got, w, tol)
            bad += not ok
            rep.add(f"table_16[{key}].{label}", w, round(got, places + 2), ok, t["source"])
        if verbose:
            print(f"  Table 16 {key:14s} median, IQR and range   "
                  f"{'ok' if not bad else f'{bad} CELLS DIFFER'}")


def check_table_24(spec: dict, rep: Report, verbose: bool = False) -> None:
    """PoseBusters failure decomposition by check group, all eight variants.

    Joined on pose_file. The PoseBusters table's docking_method names the ENGINE,
    so autodock_mgltools_exh128 carries the raw and the gnina rows together and a
    join on it would pool them and double the pose count without erroring.
    pose_file is unique in both tables and identifies the variant implicitly.
    """
    t = spec.get("table_24")
    if not t:
        return
    import numpy as np
    pb_csv, met_csv = ROOT / t["pb_input"], ROOT / t["metrics_input"]
    if not pb_csv.exists() or not met_csv.exists():
        rep.add("table_24", "8 variants x 12 columns", "input missing",
                False, t["source"])
        return
    checks = [c for grp in t["groups"].values() for c in grp]
    rep.add("table_24[applied checks]", t["n_applied_checks"], len(checks),
            len(checks) == t["n_applied_checks"], t["source"],
            "the three groups must partition the twenty-two applied checks")
    key = t["join_on"]
    pb = pd.read_csv(pb_csv, low_memory=False, usecols=checks + [key])
    met = pd.read_csv(met_csv, usecols=["method", key, "rmsd"], low_memory=False)

    def tb(col) -> "np.ndarray":
        return col.astype("boolean").fillna(False).to_numpy()

    pb["_valid"] = np.logical_and.reduce([tb(pb[c]) for c in checks])
    j = met.merge(pb, on=key, how="left")
    unjoined = int(j["_valid"].isna().sum())
    rep.add("table_24[join]", "every pose joined", f"{unjoined} unjoined",
            unjoined == 0, t["source"],
            "a partial join would silently drop poses from a variant's denominator")

    def block(x) -> list:
        share = lambda cols: round(
            100 * (~np.logical_and.reduce([tb(x[c]) for c in cols])).mean(), 1)
        return [len(x), round(100 * tb(x["_valid"]).mean(), 1),
                share(t["groups"]["chemical"]), share(t["groups"]["intramolecular"]),
                share(t["groups"]["intermolecular"]),
                round(100 * (~tb(x[t["single_check"]])).mean(), 1)]

    labels = ["poses", "valid_pct", "chem_pct", "intra_pct", "inter_pct", "mindist_pct"]
    for half, wanted in (("all", t["all_poses"]), ("near_native", t["near_native"])):
        for arm, want in wanted.items():
            d = j[j.method == t["arms"][arm]]
            if half == "near_native":
                d = d[d.rmsd <= t["near_native_rmsd"]]
            if d.empty:
                rep.add(f"table_24[{half}/{arm}]", want, "variant absent",
                        False, t["source"])
                continue
            got, bad = block(d), 0
            for label, g, w in zip(labels, got, want):
                ok = (g == w) if label == "poses" else _close(g, w, 0.05)
                bad += not ok
                rep.add(f"table_24[{half}/{arm}].{label}", w, g, ok, t["source"])
            if verbose:
                print(f"  Table 24 {half:11s} {arm:18s} "
                      f"{'ok' if not bad else f'{bad} CELLS DIFFER'}")

    # The joins back to Table 1. Same poses, same validity gate, so a divergence
    # means the decomposition has stopped describing the summary table.
    t1 = spec.get("table_1")
    if t1:
        idx = {c: i for i, c in enumerate(t1["columns"])}
        nbad = 0
        for arm, method in t["arms"].items():
            row = t1["rows"].get(method)
            if not row:
                continue
            for label, half, col in (("poses", "all_poses", "poses"),
                                     ("valid_pct", "all_poses", "valid_pct"),
                                     ("near_poses", "near_native", "near_poses")):
                w1 = row[idx["poses" if col == "poses" else
                             ("valid_pct" if col == "valid_pct" else "near_poses")]]
                got = t[half][arm][0 if col == "near_poses" else
                                   (0 if col == "poses" else 1)]
                ok = _close(got, w1, 0.05)
                nbad += not ok
                rep.add(f"table_24[join Table 1].{arm}.{label}", w1, got, ok,
                        t["note"].strip())
        if verbose:
            print(f"  Table 24 joins Table 1 on poses, validity and near-native   "
                  f"{'ok' if not nbad else f'{nbad} DIFFER'}")


def check_prose_appendix_h9(spec: dict, rep: Report, verbose: bool = False) -> None:
    """The two-pipeline contrast, a whole prose section with no float behind it.

    It is where the thesis answers its own headline question, and none of it was
    asserted. Two traps are encoded rather than left to be rediscovered.

    Orientation: bounded_claim() computes p_b - p_a with AutoDock as a, so its
    difference and interval are the negation of the printed ones.

    Equivalence: TOST at alpha 0.05 reads a 90 per cent interval, not the 95 per
    cent one quoted beside it, and the smallest passing margin is the larger
    absolute bound rounded UP. Rounding it down names a margin that fails.
    """
    t = spec.get("prose", {}).get("appendix_h9")
    if not t:
        return
    csv = ROOT / t["input"]
    if not csv.exists():
        rep.add("prose[appendix_h9]", "12 statistics x 4 depths",
                f"input missing: {t['input']}", False, t["source"])
        return
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        import thesis_endpoint_diagnostics as TED
        from stats_utils import (newcombe_paired_diff_ci, mcnemar_exact,
                                 mcnemar_power, tost_paired_proportions)
        from scipy.stats import norm
    except ImportError as exc:                       # pragma: no cover
        rep.add("prose[appendix_h9]", "12 statistics", f"cannot import: {exc}",
                False, t["source"])
        return
    import math

    df = TED.load(str(csv))
    a_meth, a_rank = t["arm_a"]
    b_meth, b_rank = t["arm_b"]
    eq = t["equivalence"]
    z90 = norm.ppf(0.95)
    bad = 0
    for i, depth in enumerate(t["depths"]):
        A = TED.recovered(df, a_meth, depth, a_rank, True)
        B = TED.recovered(df, b_meth, depth, b_rank, True)
        idx = A.index.union(B.index)
        a = A.reindex(idx, fill_value=False).astype(int).values
        b = B.reindex(idx, fill_value=False).astype(int).values

        ci = newcombe_paired_diff_ci(a, b)                 # p_b - p_a
        got = {"diff_pp": -100 * ci["diff"],
               "ci95_lo": -100 * ci["hi"], "ci95_hi": -100 * ci["lo"]}
        for field, value in got.items():
            want = t[field][i]
            ok = _close(round(value, 1), want, 0.05)
            bad += not ok
            rep.add(f"prose[h9].{field}.k{depth}", want, round(value, 2), ok, t["source"])

        _, _, p = mcnemar_exact(a, b)
        want_p = t["mcnemar_p"][i]
        ok = _close(p, want_p, max(5e-5, 0.02 * want_p))
        bad += not ok
        rep.add(f"prose[h9].mcnemar_p.k{depth}", want_p, round(p, 6), ok, t["source"])

        c90 = newcombe_paired_diff_ci(a, b, z=z90)
        lo90, hi90 = -100 * c90["hi"], -100 * c90["lo"]
        for field, value in (("ci90_lo", lo90), ("ci90_hi", hi90)):
            want = eq[field][i]
            ok = _close(round(value, 1), want, 0.05)
            bad += not ok
            rep.add(f"prose[h9].{field}.k{depth}", want, round(value, 2), ok, t["source"])
        margin = math.ceil(max(abs(lo90), abs(hi90)) * 10) / 10
        want_m = eq["smallest_passing_margin"][i]
        ok = _close(margin, want_m, 1e-9)
        bad += not ok
        rep.add(f"prose[h9].smallest_margin.k{depth}", want_m, margin, ok, t["source"],
                "the larger absolute 90% bound, rounded UP to one decimal")

        if depth == 1:
            r = t["rank1"]
            pw = mcnemar_power(a, b)
            n = ci["n"]
            checks = {
                "both_recover": ci["n11"], "neither_recovers": ci["n00"],
                "discordant": ci["n_discordant"],
                "discordant_pct": round(100 * ci["n_discordant"] / n, 1),
                "mde_pp": round(100 * pw["mde"], 1),
                "observed_power": round(pw["observed_power"], 2),
            }
            # observed_power is printed to two places and is only 0.12, so the
            # generic one-decimal window would accept anything from 0.07 to 0.17.
            tol = {"observed_power": 0.005}
            for field, value in checks.items():
                want = r[field]
                ok = (value == want) if isinstance(want, int) \
                    else _close(value, want, tol.get(field, 0.05))
                bad += not ok
                rep.add(f"prose[h9].{field}", want, value, ok, t["source"])
            got_n = int(round(pw["n_needed"], -2))
            ok = got_n == r["n_needed"]
            bad += not ok
            rep.add("prose[h9].n_needed", r["n_needed"], got_n, ok, t["source"],
                    "printed as 'roughly 4,200'; compared to the nearest hundred")

            # "Even the 428 entries of the earlier preprint release would leave
            # the MDE near 8.2 points." Same discordant proportion, larger n.
            scaled = math.sqrt(ci["n"] / 428.0)
            got_428 = round(100 * pw["mde"] * scaled, 1)
            ok = _close(got_428, r["mde_at_428"], 0.05)
            bad += not ok
            rep.add("prose[h9].mde_at_428", r["mde_at_428"], got_428, ok, t["source"],
                    "the MDE scales as 1/sqrt(n) at a fixed discordant proportion")

        # Equivalence verdicts. These were encoded and unread until the mutation
        # pass caught it: four values apiece backing "only rank-1 is equivalent,
        # within 10 points" and "no depth is equivalent within 5 points".
        for margin_pp, field in ((10, "equivalent_within_10pp"),
                                 (5, "equivalent_within_5pp")):
            want = eq[field][i]
            got = bool(tost_paired_proportions(a, b, margin=margin_pp / 100)["equivalent"])
            ok = got == want
            bad += not ok
            rep.add(f"prose[h9].{field}.k{depth}", want, got, ok, t["source"],
                    f"TOST at a {margin_pp}-point margin")
    if verbose:
        print(f"  Prose    Appendix H.9 two-pipeline contrast   "
              f"{'ok' if not bad else f'{bad} VALUES DIFFER'}")


def check_table_12(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Interaction-type mean counts on Orai1, plus the table's footnote.

    Parsed from the compare run's own stats file, whose "mean/pose" line is the
    mean of the per-complex means, which is the basis the caption states.

    History worth keeping: until 2026-09-03 the printed AutoDock column was the
    PRE_FR0 generation and this check was pinned as an expected failure. The table
    has since been regenerated onto the canonical run and two rows added, so the
    check is now expected to PASS. If it starts failing again, the first thing to
    test is whether the table has been rebuilt from a superseded tree; the giveaway
    is that every failing cell is an AutoDock one, because the DiffDock block is
    bit-identical across all four generations on disk.

    The footnote block is asserted too. It previously sat in the spec unread, which
    is how its stale "eleven" and 0.142 survived every run while the row cells were
    failing loudly. A spec value no checker reads is worse than no check at all.
    """
    t = spec.get("table_12")
    if not t:
        return
    src = ROOT / t["input"]
    if not src.exists():
        rep.add("table_12", "9 rows x 5 columns", f"input missing: {t['input']}",
                False, t["source"])
        return
    import re
    text = src.read_text()

    def section(tool: str) -> str:
        m = re.search(rf"── {re.escape(tool)} ──\n(.*?)(?=\n\s*── |\n=====)", text, re.S)
        return m.group(1) if m else ""

    parsed = {}
    for short, tool in t["tools"].items():
        for row, label in t["stats_labels"].items():
            m = re.search(rf"^\s*{re.escape(label)}\S*\s+mean/pose Bench\s+([\d.]+)\s+vs\s+Exp\s+([\d.]+).*?"
                          rf"Cliff [^=]*=([+-][\d.]+)", section(tool), re.M)
            if m:
                parsed[(short, row)] = (float(m.group(1)), float(m.group(2)), float(m.group(3)))

    nfail = 0
    for row, want in t["rows"].items():
        got_ad = parsed.get(("AD", row))
        got_dd = parsed.get(("DD", row))
        cells = [("AD bench", 0, got_ad, 0), ("AD exp", 1, got_ad, 1),
                 ("AD delta", 2, got_ad, 2), ("DD bench", 3, got_dd, 0),
                 ("DD exp", 4, got_dd, 1)]
        bad = 0
        for label, wi, tup, gi in cells:
            w = want[wi]
            got = None if tup is None else tup[gi]
            ok = got is not None and _close(got, w, 0.005)
            bad += not ok
            rep.add(f"table_12[{row}].{label}", w,
                    "row not parsed" if got is None else got, ok, t["source"])
        nfail += bad
        if verbose:
            print(f"  Table 12 {row:18s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")
    if verbose and nfail:
        print(f"  Table 12: {nfail} cells differ. If every one is an AutoDock cell, "
              f"suspect a rebuild from a superseded tree rather than a pipeline change.")

    _check_table_12_footnote(t, rep, verbose)


def _check_table_12_footnote(t: dict, rep: Report, verbose: bool = False) -> None:
    """The footnote's BH count and its ligand-unit minima.

    Read from the ligand-level contrasts sidecar rather than the stats file. The
    stats file writes "<1e-4" for the strongest cells, which cannot be ranked
    numerically, and it carries no ligand-unit collapse at all.

    The count is over the twenty-six distinct tool-by-type cells, meaning thirteen
    interaction types for each of the two mature tools. EquiBind is excluded from
    the table and from this count because it retains a single usable experimental
    pose.
    """
    fn = t.get("footnote")
    src = ROOT / t["footnote_input"] if t.get("footnote_input") else None
    if not fn or src is None:
        return
    if not src.exists():
        rep.add("table_12.footnote", "BH count and ligand-unit minima",
                f"input missing: {t['footnote_input']}", False, t["source"])
        return

    import csv as _csv

    def _f(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return None

    rows = [r for r in _csv.DictReader(src.open())
            if "interaction-type" in (r.get("family") or "")
            and r.get("tool") in ("autodock", "diffdock")
            and _f(r.get("p_adj")) is not None]

    n_clear = sum(1 for r in rows
                  if r["unit"] != "ligand" and _f(r["p_adj"]) < 0.05)
    rep.add("table_12.footnote.cells_clearing_bh", fn["cells_clearing_bh"], n_clear,
            n_clear == fn["cells_clearing_bh"], t["source"])

    fx_rows = [r for r in rows if r["unit"] != "ligand"]
    # The denominator the footnote names. Guards a regenerated sidecar quietly
    # shrinking the pool, which would leave the numerator looking plausible.
    if "total_cells" in fn:
        rep.add("table_12.footnote.total_cells", fn["total_cells"], len(fx_rows),
                len(fx_rows) == fn["total_cells"], t["source"])

    # Rows where BOTH tools clear, which the footnote names explicitly.
    per_type = {}
    for r in fx_rows:
        per_type.setdefault(r["metric"], {})[r["tool"]] = _f(r["p_adj"])
    both = sum(1 for v in per_type.values()
               if len(v) == 2 and all(p < 0.05 for p in v.values()))
    if "both_tools_clear" in fn:
        rep.add("table_12.footnote.both_tools_clear", fn["both_tools_clear"], both,
                both == fn["both_tools_clear"], t["source"])

    # COVERAGE, not just cell values. The footnote asserts that every cell clearing
    # correction appears as a printed row. Checking the cells one by one cannot catch
    # a row being DELETED from the spec, because the loop is over the spec itself: drop
    # a row and its five checks vanish silently while the harness still reports success.
    # That is the exact shape of the defect this table carried, so tie the printed row
    # count to the data instead. The fourteen clearing cells span nine distinct types.
    clearing_types = {m for m, v in per_type.items() if any(p < 0.05 for p in v.values())}
    rep.add("table_12.footnote.clearing_types_all_printed",
            f"{len(clearing_types)} types, one row each",
            f"{len(t['rows'])} rows printed",
            len(clearing_types) == len(t["rows"]), t["source"])

    for tool, key in (("autodock", "smallest_adjusted_autodock"),
                      ("diffdock", "smallest_adjusted_diffdock")):
        vals = [_f(r["p_adj"]) for r in rows
                if r["unit"] == "ligand" and r["tool"] == tool]
        got = round(min(vals), 3) if vals else None
        rep.add(f"table_12.footnote.{key}", fn[key],
                "no ligand-unit rows" if got is None else got,
                got is not None and _close(got, fn[key], 0.0005), t["source"])
        # The footnote's qualitative claim: nothing clears at the ligand unit.
        # Read the expectation from the spec so it can be perturbed like any other.
        if vals and "none_clears_at_ligand_unit" in fn:
            want = fn["none_clears_at_ligand_unit"]
            got = min(vals) >= 0.05
            rep.add(f"table_12.footnote.{tool}_none_clears_at_ligand_unit",
                    want, got, got == want, t["source"])

    if verbose:
        print(f"  Table 12 footnote: {n_clear} of 26 clear BH, {both} rows clear for both tools")


def check_table_9(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Docking-protocol parameters, read back from the committed configs.

    Only the settings are checkable; the rest of the table is prose. These are
    worth asserting anyway, because a config edited for one arm while the appendix
    keeps describing the old one is invisible until someone reads both files side
    by side, and nothing else in the harness reads a config at all.
    """
    t = spec.get("table_9")
    if not t:
        return
    cfg_dir = ROOT / t["config_dir"]

    def load(name: str) -> dict | None:
        p_cfg = cfg_dir / name
        return yaml.safe_load(p_cfg.read_text()) if p_cfg.exists() else None

    bad = 0
    for tool in ("autodock", "diffdock", "equibind"):
        block = t[tool]
        cfg = load(block["config"])
        if cfg is None:
            rep.add(f"table_9[{tool}]", "config present",
                    f"missing: {block['config']}", False, t["source"])
            bad += 1
            continue
        for field, want in block.items():
            if field == "config":
                continue
            if field == "n_rdkit_seeds":
                got = len(cfg.get("rdkit_seeds") or [])
            else:
                got = cfg.get(field)
            ok = got == want
            bad += not ok
            rep.add(f"table_9[{tool}].{field}", want, got, ok, t["source"],
                    f"read from {block['config']}")
        if verbose:
            print(f"  Table 9  {tool:10s} protocol settings   "
                  f"{'ok' if not bad else 'SOME DIFFER'}")

    lbad = 0
    for name, want in t["ladder"].items():
        cfg = load(name)
        got = None if cfg is None else cfg.get("exhaustiveness")
        ok = got == want
        lbad += not ok
        rep.add(f"table_9[ladder].{want}", want, got, ok, t["source"],
                f"the ladder rung declared by {name}")
    if verbose:
        print(f"  Table 9  exhaustiveness ladder 18/32/64/92/128   "
              f"{'ok' if not lbad else f'{lbad} DIFFER'}")


def check_table_11(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Orai1 protocol parameters and the control panel's pose and unit accounting.

    The accounting is the substantive half. The two panels were sampled at
    different budgets and the appendix spends a page bounding what that does to
    the comparison; every number in that argument was unasserted, including the
    12,315 control poses the top-ten-cap analysis rests on.
    """
    t = spec.get("table_11")
    if not t:
        return
    cfg_dir = ROOT / t["config_dir"]
    bad = 0
    for name, block in t["configs"].items():
        p_cfg = cfg_dir / block["config"]
        if not p_cfg.exists():
            rep.add(f"table_11[{name}]", "config present",
                    f"missing: {block['config']}", False, t["source"])
            bad += 1
            continue
        cfg = yaml.safe_load(p_cfg.read_text())
        for field, want in block.items():
            if field == "config":
                continue
            got = cfg.get(field)
            ok = got == want
            bad += not ok
            rep.add(f"table_11[{name}].{field}", want, got, ok, t["source"],
                    f"read from {block['config']}")
    if verbose:
        print(f"  Table 11 protocol settings, both panels   "
              f"{'ok' if not bad else f'{bad} DIFFER'}")

    csv = ROOT / t["control_pb"]
    if not csv.exists():
        rep.add("table_11[accounting]", "control PoseBusters table",
                f"missing: {t['control_pb']}", False, t["source"])
        return
    df = pd.read_csv(csv, low_memory=False,
                     usecols=["docking_method", "optimizer", "pocket_source",
                              "protein", "ligand", "pose_name"])
    sel = {}
    for tool, (meth, opt, pocket) in t["arms"].items():
        x = df[df.docking_method == meth]
        if opt is not None:
            x = x[x.optimizer == opt]
        if pocket is not None:
            x = x[x.pocket_source == pocket]
        sel[tool] = x

    abad = 0
    for tool, want in t["control_poses"].items():
        got = len(sel[tool])
        ok = got == want
        abad += not ok
        rep.add(f"table_11[poses].{tool}", want, got, ok, t["source"])
    units = {}
    for tool, want in t["control_units"].items():
        units[tool] = set(map(tuple, sel[tool][["protein", "ligand"]]
                              .drop_duplicates().values))
        got = len(units[tool])
        ok = got == want
        abad += not ok
        rep.add(f"table_11[units].{tool}", want, got, ok, t["source"])

    g = t["diffdock_gap"]
    missing = units["autodock"] - units["diffdock"]
    ok = len(missing) == g["units_missing"]
    abad += not ok
    rep.add("table_11[gap].units_missing", g["units_missing"], len(missing), ok,
            t["source"], "units AutoDock covers and DiffDock does not")
    import collections
    by_frame = collections.Counter(protein for protein, _ in missing)
    for frame, want in g["by_frame"].items():
        got = by_frame.get(frame, 0)
        ok = got == want
        abad += not ok
        rep.add(f"table_11[gap].{frame}", want, got, ok, t["source"])
    per = sel["diffdock"].groupby(["protein", "ligand"]).size()
    counts = {"units_full_ten": int((per == 10).sum()),
              "units_under_ten": int((per < 10).sum()),
              "units_with_nine": int((per == 9).sum()),
              "units_with_one": int((per == 1).sum())}
    for field, got in counts.items():
        want = g[field]
        ok = got == want
        abad += not ok
        rep.add(f"table_11[gap].{field}", want, got, ok, t["source"])
    if verbose:
        print(f"  Table 11 control pose and unit accounting   "
              f"{'ok' if not abad else f'{abad} DIFFER'}")


def check_table_17(spec: dict, rep: Report, verbose: bool = False) -> None:
    """Physicochemical descriptors of the three experimental Orai1 ligands.

    Recomputed with RDKit from the optimised SDFs. Two traps are enforced rather
    than trusted: the GSK row must come from the neutral file and not from the
    phenolate beside it, and rotatable bonds must be counted hydrogen-suppressed.
    Both would otherwise change printed cells while still looking reasonable.
    """
    t = spec.get("table_17")
    if not t:
        return
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
    except ImportError as exc:                       # pragma: no cover
        rep.add("table_17", "3 ligands x 10 descriptors", f"cannot import RDKit: {exc}",
                False, t["source"])
        return
    RDLogger.DisableLog("rdApp.*")

    def load(rel: str):
        f = ROOT / rel
        if not f.exists():
            return None
        m = Chem.MolFromMolFile(str(f), removeHs=t["remove_hs"])
        return None if m is None else Chem.RemoveHs(m)

    for name, want in t["ligands"].items():
        m = load(want["input"])
        if m is None:
            rep.add(f"table_17[{name}]", "molecule loads",
                    f"could not read {want['input']}", False, t["source"])
            continue
        mw = Descriptors.MolWt(m)
        hbd, hba = rdMolDescriptors.CalcNumHBD(m), rdMolDescriptors.CalcNumHBA(m)
        clogp = Crippen.MolLogP(m)
        got = {
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "mw": mw,
            "heavy_atoms": m.GetNumHeavyAtoms(),
            "rot_bonds": rdMolDescriptors.CalcNumRotatableBonds(m),
            "arom_rings": rdMolDescriptors.CalcNumAromaticRings(m),
            "hbd": hbd,
            "hba": hba,
            "tpsa": rdMolDescriptors.CalcTPSA(m),
            "clogp": clogp,
            "ro5_violations": sum([hbd > 5, hba > 10, mw > 500, clogp > 5]),
        }
        # Half of the last printed place. Compared unrounded, so an exact half
        # such as 2abp-NH2's TPSA of 35.25 matches the printed 35.3.
        tol = {"mw": 0.05, "tpsa": 0.05, "clogp": 0.005}
        bad = 0
        for field, w in want.items():
            if field == "input":
                continue
            g = got[field]
            ok = (g == w) if isinstance(w, (str, int)) and field not in tol \
                else _close(g, w, tol.get(field, 0.05))
            bad += not ok
            rep.add(f"table_17[{name}].{field}", w,
                    g if isinstance(g, (str, int)) else round(g, 4), ok, t["source"])
        if verbose:
            print(f"  Table 17 {name:22s} {'ok' if not bad else f'{bad} CELLS DIFFER'}")

    # The phenolate must stay a different species from the neutral row above.
    w = t["wrong_gsk_file"]
    m = load(w["path"])
    got = "file absent" if m is None else rdMolDescriptors.CalcMolFormula(m)
    ok = got == w["formula"]
    rep.add("table_17[wrong GSK file stays distinct]", w["formula"], got, ok,
            t["source"], w["note"].strip())
    if verbose:
        print(f"  Table 17 phenolate file is a distinct species   "
              f"{'ok' if ok else 'CHANGED'}")


# =============================================================================

if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    report = run_all(verbose=True)
    print()
    expected = tuple(f"{k}[" for k, v in _spec().items()
                     if isinstance(v, dict) and v.get("expected_to_fail"))
    print(report.summary(expected))
    for fail in report.failures:
        print(f"\n  {fail.key}\n      thesis    : {fail.expected}   ({fail.source})"
              f"\n      recomputed: {fail.actual}")
        if fail.note:
            print(f"      {fail.note}")
    sys.exit(1 if report.failures else 0)
