#!/usr/bin/env python
"""Build ``XRay_PoseBusters_Control.ipynb``.

The notebook is generated rather than hand-edited for the same reason
``_build_reproduction_notebook.py`` generates the thesis reproduction notebook:
a notebook edited in place accumulates execution-order damage that its page
order hides. Edit this file, re-run it, then execute the notebook.

    /home/manndo/anaconda3/envs/vina/bin/python \
        Scripts/Analysis/_build_xray_control_notebook.py
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path("/home/manndo/master_dev")
OUT = ROOT / "XRay_PoseBusters_Control.ipynb"

cells: list[dict] = []


def _cid() -> str:
    return f"cell-{len(cells):02d}"


def md(text: str) -> None:
    cells.append({"cell_type": "markdown", "id": _cid(), "metadata": {},
                  "source": text.strip("\n").splitlines(keepends=True)})


def code(text: str) -> None:
    cells.append({"cell_type": "code", "id": _cid(), "execution_count": None,
                  "metadata": {}, "outputs": [],
                  "source": text.strip("\n").splitlines(keepends=True)})


# =============================================================================
md(r"""
# PoseBusters on the crystal ligands of the PoseBusters benchmark set

This notebook runs the PoseBusters validity battery on the **experimental
X-ray ligand geometries** of the benchmark set, each paired with **its own
benchmark protein**. Nothing here is docked. The ligand handed to the checker
is the deposited crystallographic conformer, so the result is the ceiling that
the docking arms of the thesis are measured against, measured with the same
instrument rather than assumed.

It is a self-contained control experiment. It reads only the benchmark set and
the receptor tree of the dominant AutoDock arm, writes only into its own output
directory, and changes nothing the thesis reads.

### Why it is worth running

The thesis reports a PoseBusters pass rate for every docking arm and treats a
failure as a defect of the pose. That reading holds only if the battery passes
the crystal geometry in the first place. Three things can break it, and only a
measurement separates them.

- A **crystal ligand can fail its own battery**. Deposited geometries carry
  refinement restraints, partial occupancy and, in low-resolution structures,
  real strain. A check that the deposited conformer fails is not a check a
  docking program can be asked to pass.
- The **receptor file decides the intermolecular verdict**. Nine of the
  twenty-two applied checks measure the ligand against the protein, its
  cofactors and its waters. The distributed crystal protein still carries
  cofactors; the prepared docking receptor has them stripped. Whether that
  difference matters is measurable, and it is measured below.
- The **battery is not free of its own artefacts**. Six entries of the
  benchmark set still carry stray atoms of the cognate ligand inside the
  protein file, which is a data defect rather than a geometry defect.

### The two arms

Both arms bust the identical crystal ligand file. Only the protein differs.

| Arm | Protein handed to the checker |
| --- | --- |
| `crystal` | `Data/PoseBuster Benchmark Set/<ID>/<ID>_protein.pdb`, as distributed |
| `prepared` | the receptor the dominant AutoDock arm was docked and validated against |

The second arm is what makes the comparison to the thesis legitimate: every
docked pose in the thesis was busted against that receptor, so the ceiling has
to be measured against it too.

### What this notebook is not

It does not compute a root-mean-square deviation, and it supplies no reference
ligand to the checker. That matches the thesis, whose validity battery is
reference-free in every dataset. A crystal ligand scored against itself would
return zero by construction and would measure nothing.
""")

# =============================================================================
code(r"""
# ── Setup ───────────────────────────────────────────────────────────────────
import hashlib, json, os, pickle, re, select, signal, struct, sys, time
from pathlib import Path

ROOT = Path("/home/manndo/master_dev")
os.chdir(ROOT)

BENCH   = ROOT / "Data/PoseBuster Benchmark Set"
IDS     = BENCH / "posebusters_pdb_ccd_ids.txt"
PREPTREE= ROOT / "Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128"
THESIS_PB = ROOT / "posebusters_results/benchmark_matched_equibind/dock/posebusters_filtered_results.csv"
EXPECTED  = ROOT / "Scripts/Analysis/thesis_expected_values.yaml"

OUTDIR  = ROOT / "posebusters_results/xray_crystal_control"
FIGDIR  = OUTDIR / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

# ── run knobs ───────────────────────────────────────────────────────────────
# A full pass is 616 calls, one per complex per arm. Cost is dominated by the
# internal-energy check, which embeds and minimises a fifty-conformer ensemble:
# a small ligand takes under a second and the largest took ninety-three, for
# about seventy-three minutes of processor time in total. On thirty-two workers
# that is roughly ten minutes of wall clock, most of it the tail of a handful of
# large flexible ligands. Results are cached in the raw CSV below, so a re-run
# of this notebook costs nothing unless OVERWRITE is set.
N_WORKERS    = min(32, os.cpu_count() or 8)
POSE_TIMEOUT = 600      # seconds; a single bust() that wedges is killed, not waited on
OVERWRITE    = False    # True re-busts everything, ignoring the cached raw CSV

RAW_CSV = OUTDIR / "xray_posebusters_raw.csv"

sys.path.insert(0, str(ROOT / "Scripts/Analysis"))
sys.path.insert(0, str(ROOT / "Scripts/Docking/Posebusters"))

print(f"project root   : {ROOT}")
print(f"benchmark set  : {BENCH.relative_to(ROOT)}")
print(f"output         : {OUTDIR.relative_to(ROOT)}")
print(f"workers        : {N_WORKERS}   per-complex timeout: {POSE_TIMEOUT}s")
""")

# =============================================================================
md(r"""
## 1. Environment and patch guard

The numbers below are tied to the versions in this table, and to one patch that
lives outside the repository and is therefore invisible to `git`.

The **PoseBusters patch** repairs the handler guarding the check for missing
force-field parameters. As shipped, that handler reads the second argument of
an assertion carrying only one, raises `IndexError`, and aborts the whole
`bust()` call rather than the one check. Every ligand the force field cannot
parameterise, which is to say anything with hypervalent sulfur, a phosphate or
a coordinated metal, would be dropped from the results with no error reported
and no row written. Silent loss of exactly the hardest ligands would bias this
control upwards, so the guard is a hard stop rather than a warning.

The canonical list of pass or fail checks is imported from the thesis runner
rather than retyped, then compared against a literal written here. The result
table carries about a hundred columns under `full_report`, several of them
boolean diagnostics that a heuristic would mistake for verdicts. One of them,
the most extreme protein clash flag, is `False` on every row and would fail
every complex on its own. Restricting the verdict to the allowlist is what
keeps it correct.
""")

code(r"""
import importlib.metadata as _md
import numpy as np, pandas as pd

# The twenty-two stock PoseBusters verdicts under config="dock", v0.6.x.
# Imported from the thesis runner so it cannot drift from the rest of the
# pipeline, and compared against the literal below so a drift is loud.
_LITERAL_CHECKS = (
    "mol_pred_loaded", "mol_cond_loaded", "sanitization", "inchi_convertible",
    "all_atoms_connected", "no_radicals",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "internal_energy",
    "protein-ligand_maximum_distance",
    "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors",
    "minimum_distance_to_waters",
    "volume_overlap_with_protein",
    "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors",
    "volume_overlap_with_waters",
)
try:
    from run_posebusters import CANONICAL_TEST_COLUMNS as PB_CHECKS
    _check_source = "imported from Scripts/Docking/Posebusters/run_posebusters.py"
except Exception as exc:                       # pragma: no cover
    PB_CHECKS, _check_source = _LITERAL_CHECKS, f"literal fallback ({exc})"
PB_CHECKS = tuple(PB_CHECKS)
assert PB_CHECKS == _LITERAL_CHECKS, (
    "the thesis runner's verdict allowlist has drifted from this notebook's literal:\n"
    f"  runner   : {PB_CHECKS}\n  notebook : {_LITERAL_CHECKS}")

# Group taxonomy, identical to the one Table 24 of the thesis partitions on.
CHECK_GROUPS = {
    "chemical": ["mol_pred_loaded", "mol_cond_loaded", "sanitization",
                 "inchi_convertible", "all_atoms_connected", "no_radicals"],
    "intramolecular": ["bond_lengths", "bond_angles", "aromatic_ring_flatness",
                       "non-aromatic_ring_non-flatness", "double_bond_flatness",
                       "internal_steric_clash", "internal_energy"],
    "intermolecular": ["protein-ligand_maximum_distance", "minimum_distance_to_protein",
                       "minimum_distance_to_organic_cofactors",
                       "minimum_distance_to_inorganic_cofactors",
                       "minimum_distance_to_waters", "volume_overlap_with_protein",
                       "volume_overlap_with_organic_cofactors",
                       "volume_overlap_with_inorganic_cofactors",
                       "volume_overlap_with_waters"],
}
assert sorted(c for g in CHECK_GROUPS.values() for c in g) == sorted(PB_CHECKS), \
    "the three groups must partition the twenty-two applied checks"

# ── version table ───────────────────────────────────────────────────────────
WANT = {"python": "3.12.12", "posebusters": "0.6.3", "rdkit": "2025.09.1",
        "pandas": "2.3.3"}
rows = []
found_py = ".".join(map(str, sys.version_info[:3]))
rows.append(("python", found_py, WANT["python"], found_py == WANT["python"]))
for pkg in ("posebusters", "rdkit", "pandas", "numpy", "scipy", "matplotlib"):
    try:
        got = _md.version(pkg if pkg != "rdkit" else "rdkit")
    except Exception as exc:
        got = f"<{exc}>"
    want = WANT.get(pkg, "")
    ok = (not want) or got == want or got.replace(".0", ".") == want.replace(".0", ".")
    rows.append((pkg, got, want or "(unpinned)", ok))

# ── the energy-ratio patch ──────────────────────────────────────────────────
import posebusters as _pb_pkg
PATCH = Path(_pb_pkg.__file__).parent / "modules/energy_ratio.py"
patched, note = False, "energy_ratio.py not found"
if PATCH.exists():
    block = re.search(r"UFFHasAllMoleculeParams.*?except Exception as e:\s*\n(.*?)return _empty_results",
                      PATCH.read_text(), re.S)
    patched = bool(block) and "e.args[1]" not in block.group(1)
    note = "handler logs the exception" if patched else "handler still indexes e.args[1]"
rows.append(("posebusters energy-ratio patch", note, "patched", patched))

w = max(len(r[0]) for r in rows)
for name, got, want, ok in rows:
    print(f"  {'ok   ' if ok else 'DRIFT'} {name:{w}s}  {str(got)[:34]:36s} want: {want}")
print()
print(f"verdict allowlist: {len(PB_CHECKS)} checks, {_check_source}")

bad = [r[0] for r in rows if not r[3]]
assert not bad, f"environment drift, the numbers below would not be comparable: {bad}"
print("environment matches the thesis pipeline")
""")

# =============================================================================
md(r"""
## 2. Cohort and inputs

The cohort is the 308 identifiers the thesis prepared and docked, read from the
same list file the docking pipeline read. Three files are resolved per entry
and every one of them is hashed, so a later re-run that disagrees can be traced
to an input rather than to the checker.

The prepared receptor is the `.pdb` twin the docking stage kept beside the
`.pdbqt` it actually docked into. It is heavy-atom identical to that `.pdbqt`,
which is what allows a pose validated against it to be compared with a crystal
ligand validated against it.
""")

code(r"""
def sha256(path: Path, _cache={}) -> str:
    key = (str(path), path.stat().st_mtime_ns, path.stat().st_size)
    if key not in _cache:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        _cache[key] = h.hexdigest()
    return _cache[key]


ids = [ln.strip() for ln in IDS.read_text().splitlines() if ln.strip()]
assert len(ids) == 308, f"expected the 308-identifier cohort, found {len(ids)}"

manifest = []
for i in ids:
    lig  = BENCH / i / f"{i}_ligand.sdf"
    cry  = BENCH / i / f"{i}_protein.pdb"
    prep = PREPTREE / i / "_staging/receptors/pdbqt" / f"{i}_protein_mgl_tools.pdb"
    manifest.append({
        "entry": i, "pdb": i.split("_")[0], "ccd": i.split("_", 1)[1],
        "ligand_file": str(lig), "ligand_sha256": sha256(lig) if lig.exists() else "",
        "crystal_receptor": str(cry), "crystal_sha256": sha256(cry) if cry.exists() else "",
        "prepared_receptor": str(prep), "prepared_sha256": sha256(prep) if prep.exists() else "",
    })
mf = pd.DataFrame(manifest)

gaps = {c: int((mf[c] == "").sum()) for c in ("ligand_sha256", "crystal_sha256", "prepared_sha256")}
print(f"cohort                      : {len(mf)} complexes")
for c, n in gaps.items():
    print(f"missing {c.replace('_sha256',''):<20s}: {n}")
assert sum(gaps.values()) == 0, f"inputs missing: {gaps}"

mf.to_csv(OUTDIR / "xray_input_manifest.csv", index=False)
print(f"\nmanifest written to {(OUTDIR / 'xray_input_manifest.csv').relative_to(ROOT)}")
print(f"distinct ligand files       : {mf.ligand_sha256.nunique()}")
print(f"distinct crystal receptors  : {mf.crystal_sha256.nunique()}")
print(f"distinct prepared receptors : {mf.prepared_sha256.nunique()}")
""")

# =============================================================================
md(r"""
## 3. What is actually in the two protein files

The intermolecular half of the battery measures the ligand against whatever the
protein file contains. Reading that content before running is the difference
between interpreting a cofactor clash and inventing one.

Two things are worth knowing before the verdict arrives. The distributed
crystal proteins carry no water at all, so the two water checks cannot fail for
a reason connected to water. And a handful of entries still carry stray atoms of
the cognate ligand inside the protein file, which the checker will classify as
an organic cofactor sitting on top of the ligand it is validating.
""")

code(r"""
def survey(pdb: Path) -> dict:
    atoms = het = wat = 0
    resnames: dict[str, int] = {}
    elements: dict[str, int] = {}
    for line in pdb.read_text().splitlines():
        if line.startswith("ATOM"):
            atoms += 1
        elif line.startswith("HETATM"):
            het += 1
            rn = line[17:20].strip()
            resnames[rn] = resnames.get(rn, 0) + 1
            if rn == "HOH":
                wat += 1
        else:
            continue
        el = line[76:78].strip() or "?"
        elements[el] = elements.get(el, 0) + 1
    return {"atom": atoms, "hetatm": het, "water": wat,
            "resnames": resnames, "elements": elements}


aud = []
for r in mf.itertuples():
    c = survey(Path(r.crystal_receptor))
    p = survey(Path(r.prepared_receptor))
    aud.append({
        "entry": r.entry,
        "crystal_polymer_atoms": c["atom"], "crystal_hetatm": c["hetatm"],
        "crystal_waters": c["water"],
        "crystal_cognate_atoms": c["resnames"].get(r.ccd, 0),
        "crystal_hydrogens": c["elements"].get("H", 0),
        "prepared_polymer_atoms": p["atom"], "prepared_hetatm": p["hetatm"],
        "prepared_hydrogens": p["elements"].get("H", 0),
        "cofactor_residues": ";".join(sorted(k for k in c["resnames"] if k not in ("HOH",))),
    })
ad = pd.DataFrame(aud)
ad.to_csv(OUTDIR / "xray_receptor_audit.csv", index=False)

print("crystal protein files, as distributed")
print(f"  entries carrying any heteroatom record : {(ad.crystal_hetatm > 0).sum()} of {len(ad)}")
print(f"  entries carrying crystallographic water: {(ad.crystal_waters > 0).sum()} of {len(ad)}")
print(f"  entries carrying explicit hydrogen     : {(ad.crystal_hydrogens > 0).sum()} of {len(ad)}")
print()
print("prepared docking receptors")
print(f"  entries carrying any heteroatom record : {(ad.prepared_hetatm > 0).sum()} of {len(ad)}")
print(f"  entries carrying explicit hydrogen     : {(ad.prepared_hydrogens > 0).sum()} of {len(ad)}")
print(f"  polymer atom count identical to crystal: "
      f"{(ad.crystal_polymer_atoms == ad.prepared_polymer_atoms).sum()} of {len(ad)}")
print()

leftover = ad[ad.crystal_cognate_atoms > 0][["entry", "crystal_cognate_atoms"]]
print(f"entries whose protein file still contains atoms of the cognate ligand: {len(leftover)}")
if len(leftover):
    for r in leftover.itertuples():
        print(f"    {r.entry:12s} {r.crystal_cognate_atoms} atom(s)")
    print("  these are a defect of the distributed file, not of the ligand geometry;")
    print("  they are flagged on every table below rather than silently removed.")

top_cof = (ad.cofactor_residues.str.split(";").explode().replace("", np.nan).dropna()
             .value_counts().head(12))
print("\nmost frequent heteroatom residues in the crystal protein files")
for name, n in top_cof.items():
    print(f"    {name:6s} {n:4d} entries")
""")

# =============================================================================
md(r"""
## 4. The checker

One `bust()` call per complex per arm, under `config="dock"`, with
`full_report=True` and no reference ligand. That is the identical call the
thesis pipeline makes for every docked pose.

Two guards are carried over from that pipeline because they are not optional in
a parallel run. Each worker is pinned to a single thread, since the internal
energy check embeds and minimises a conformer ensemble with the force field's
default request for every core, and one such pool per worker is a thread count
that collapses throughput on a loaded machine. And each call runs in a hand
forked grandchild so that a wedge inside the conformer generator can be killed
by the kernel. A timer inside the interpreter cannot break a tight loop in C++,
and a pool worker is not allowed to spawn a child by the ordinary route.
""")

code(r"""
import multiprocessing as mp

_buster = None


def _worker_init():
    global _buster
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ.setdefault(v, "1")
    os.environ.setdefault("PYTHONWARNINGS", "ignore")
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        pass
    # Force the internal-energy ensemble onto one thread. get_energies() reads
    # new_conformation as a module global, so replacing the attribute is enough.
    try:
        from posebusters.modules import energy_ratio as _er
        _orig = _er.new_conformation
        if getattr(_orig, "__name__", "") != "_single_thread_new_conformation":
            def _single_thread_new_conformation(mol, n_confs=1, num_threads=0,
                                                energy_minimization=True):
                return _orig(mol, n_confs, 1, energy_minimization)
            _er.new_conformation = _single_thread_new_conformation
    except Exception:
        pass
    from posebusters import PoseBusters
    _buster = PoseBusters(config="dock")


class _Timeout(Exception):
    pass


_PR_SET_PDEATHSIG = 1
try:
    import ctypes as _ctypes
    _LIBC = _ctypes.CDLL("libc.so.6", use_errno=True)
except Exception:
    _LIBC = None


def _write_all(fd, data):
    while data:
        data = data[os.write(fd, data):]


def _read_exact(fd, n, deadline):
    buf = b""
    while len(buf) < n:
        left = deadline - time.monotonic()
        if left <= 0 or not select.select([fd], [], [], left)[0]:
            return None
        chunk = os.read(fd, n - len(buf))
        if not chunk:
            return None
        buf += chunk
    return buf


def _reap(pid, killed):
    try:
        if not killed:
            os.waitpid(pid, 0)
            return
    except (ChildProcessError, OSError):
        return
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.kill(pid, sig)
        except OSError:
            return
        for _ in range(30):
            time.sleep(0.05)
            try:
                if os.waitpid(pid, os.WNOHANG)[0]:
                    return
            except (ChildProcessError, OSError):
                return


def _bust_guarded(lig, prot, timeout):
    # bust() in a forked grandchild, hard-killed if it outruns *timeout*.
    r, w = os.pipe()
    pid = os.fork()
    if pid == 0:                                   # child
        try:
            os.close(r)
            if _LIBC is not None:
                try:
                    _LIBC.prctl(_PR_SET_PDEATHSIG, signal.SIGKILL, 0, 0, 0)
                except Exception:
                    pass
            try:
                df = _buster.bust(lig, None, prot, full_report=True)
                blob = pickle.dumps(("ok", df), protocol=pickle.HIGHEST_PROTOCOL)
            except Exception as e:
                blob = pickle.dumps(("err", f"{type(e).__name__}: {e}"))
            _write_all(w, struct.pack(">Q", len(blob)))
            _write_all(w, blob)
        except Exception:
            pass
        finally:
            os._exit(0)
    os.close(w)                                    # parent
    deadline = time.monotonic() + timeout
    got = False
    try:
        head = _read_exact(r, 8, deadline)
        if head is None:
            raise _Timeout()
        (n,) = struct.unpack(">Q", head)
        body = _read_exact(r, n, deadline)
        if body is None:
            raise _Timeout()
        status, payload = pickle.loads(body)
        got = True
    finally:
        os.close(r)
        _reap(pid, killed=not got)
    if status == "err":
        raise RuntimeError(payload)
    return payload


def _one(task):
    entry, arm, lig, prot, timeout = task
    t0 = time.monotonic()
    try:
        df = _bust_guarded(lig, prot, timeout)
    except _Timeout:
        return entry, arm, None, f"TIMEOUT after {timeout}s"
    except Exception as e:
        return entry, arm, None, f"{type(e).__name__}: {e}"
    df = df.reset_index(drop=True)
    df["entry"] = entry
    df["arm"] = arm
    df["ligand_file"] = lig
    df["receptor_file"] = prot
    df["elapsed_s"] = round(time.monotonic() - t0, 3)
    return entry, arm, df, None


def run_arms(manifest_df, arms, timeout=POSE_TIMEOUT, workers=N_WORKERS):
    tasks = []
    for r in manifest_df.itertuples():
        for arm, col in arms.items():
            tasks.append((r.entry, arm, r.ligand_file, getattr(r, col), timeout))
    frames, errors = [], []
    t0 = time.monotonic()
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=workers, initializer=_worker_init) as pool:
        for k, (entry, arm, df, err) in enumerate(pool.imap_unordered(_one, tasks, chunksize=1), 1):
            if df is not None:
                frames.append(df)
            else:
                errors.append({"entry": entry, "arm": arm, "error": err})
            if k % 100 == 0 or k == len(tasks):
                print(f"    {k:4d}/{len(tasks)}  {time.monotonic() - t0:6.1f}s  "
                      f"{len(errors)} error(s)", flush=True)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return out, pd.DataFrame(errors)


print("checker defined:  PoseBusters(config='dock').bust(ligand, None, protein, full_report=True)")
print(f"                  no reference ligand, {len(PB_CHECKS)} verdict columns")
""")

# =============================================================================
code(r"""
ARMS = {"crystal": "crystal_receptor", "prepared": "prepared_receptor"}

if RAW_CSV.exists() and not OVERWRITE:
    raw = pd.read_csv(RAW_CSV, low_memory=False)
    errs = pd.DataFrame()
    print(f"loaded cached results from {RAW_CSV.relative_to(ROOT)}  ({len(raw)} rows)")
    print("set OVERWRITE = True in the setup cell to re-bust from scratch")
else:
    print(f"busting {len(mf)} complexes x {len(ARMS)} arms on {N_WORKERS} workers")
    raw, errs = run_arms(mf, ARMS)
    raw.to_csv(RAW_CSV, index=False)
    print(f"\nwritten to {RAW_CSV.relative_to(ROOT)}")

print(f"\nrows            : {len(raw)}")
print(f"complexes       : {raw.entry.nunique()}")
print(f"rows per arm    : {dict(raw.arm.value_counts())}")
if len(errs):
    print(f"\nERRORS ({len(errs)}):")
    print(errs.to_string(index=False))
assert len(raw) == len(mf) * len(ARMS), (
    f"expected {len(mf) * len(ARMS)} rows, got {len(raw)}; an entry was dropped, "
    "which is the failure mode the energy-ratio patch exists to prevent")
print("\nevery complex returned a row in both arms")
""")

# =============================================================================
md(r"""
## 5. The verdict

A complex is valid when it passes all twenty-two checks. A missing value is a
failure, not a pass, which is the same strict reading the thesis pipeline uses:
the checks that return a missing value are the ones whose module bailed out,
and a module that could not run has not cleared anything.
""")

code(r"""
_BOOL_STR = {"True": True, "true": True, "False": False, "false": False}


def strict_bool(v):
    # Missing or unparseable is a failure, never a pass.
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, str):
        return _BOOL_STR.get(v.strip(), False)
    if v is None:
        return False
    try:
        if pd.isna(v):
            return False
    except (TypeError, ValueError):
        pass
    if isinstance(v, (int, float, np.integer, np.floating)):
        return bool(v)
    return False


res = raw.copy()
missing = [c for c in PB_CHECKS if c not in res.columns]
assert not missing, f"result table is missing stock verdict columns: {missing}"
for c in PB_CHECKS:
    res[c] = res[c].map(strict_bool).astype(bool)
res["pb_valid"] = res[list(PB_CHECKS)].all(axis=1)
for g, cols in CHECK_GROUPS.items():
    res[f"pass_{g}"] = res[cols].all(axis=1)

# a complex is flagged if its protein file still holds atoms of its own ligand
flagged = set(ad.loc[ad.crystal_cognate_atoms > 0, "entry"])
res["cognate_contamination"] = res.entry.isin(flagged)


def wilson_ci(k, n, z=1.96):
    if not n:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = (z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / d
    return (max(0.0, c - h), min(1.0, c + h))


summary = []
for arm in ARMS:
    for label, sub in (("all 308", res[res.arm == arm]),
                       ("clean 302", res[(res.arm == arm) & ~res.cognate_contamination])):
        k, n = int(sub.pb_valid.sum()), len(sub)
        lo, hi = wilson_ci(k, n)
        summary.append({
            "arm": arm, "cohort": label, "complexes": n, "valid": k,
            "valid_pct": 100 * k / n,
            "ci95_low_pct": 100 * lo, "ci95_high_pct": 100 * hi,
            "chemical_pct": 100 * sub.pass_chemical.mean(),
            "intramolecular_pct": 100 * sub.pass_intramolecular.mean(),
            "intermolecular_pct": 100 * sub.pass_intermolecular.mean(),
        })
sm = pd.DataFrame(summary)
sm.to_csv(OUTDIR / "xray_validity_summary.csv", index=False)

print("PoseBusters validity of the deposited crystal ligands\n")
hdr = f"{'arm':10s} {'cohort':10s} {'n':>4s} {'valid':>6s} {'rate':>8s} {'95% CI':>16s}   " \
      f"{'chem':>6s} {'intra':>6s} {'inter':>6s}"
print(hdr)
print("-" * len(hdr))
for r in sm.itertuples():
    print(f"{r.arm:10s} {r.cohort:10s} {r.complexes:4d} {r.valid:6d} {r.valid_pct:7.2f}% "
          f" [{r.ci95_low_pct:5.2f}, {r.ci95_high_pct:6.2f}]   "
          f"{r.chemical_pct:5.1f}% {r.intramolecular_pct:5.1f}% {r.intermolecular_pct:5.1f}%")
print("\nchem / intra / inter are the share of complexes passing every check in that group.")
""")

# =============================================================================
md(r"""
## 6. Which checks the crystal ligands fail

A pass rate on its own says nothing about whether a failure is a property of
the geometry or of the file it was compared with. The per-check table separates
them, and the named entries make each failure auditable one structure at a
time.
""")

code(r"""
GROUP_OF = {c: g for g, cols in CHECK_GROUPS.items() for c in cols}

per_check = []
for arm in ARMS:
    sub = res[res.arm == arm]
    for c in PB_CHECKS:
        fail = sub.loc[~sub[c], "entry"].tolist()
        per_check.append({
            "arm": arm, "check": c, "group": GROUP_OF[c],
            "n": len(sub), "passed": int(sub[c].sum()), "failed": len(fail),
            "fail_pct": 100 * len(fail) / len(sub),
            "failing_entries": ";".join(sorted(fail)),
        })
pc = pd.DataFrame(per_check)
pc.to_csv(OUTDIR / "xray_check_failure_rates.csv", index=False)

for arm in ARMS:
    sub = pc[(pc.arm == arm) & (pc.failed > 0)].sort_values("failed", ascending=False)
    print(f"\n{arm} receptor arm — checks failed by at least one crystal ligand")
    if sub.empty:
        print("    none; every crystal ligand passes every check in this arm")
        continue
    for r in sub.itertuples():
        ent = r.failing_entries.split(";")
        shown = ", ".join(ent[:8]) + (f", +{len(ent) - 8} more" if len(ent) > 8 else "")
        print(f"    {r.check:42s} {r.group:15s} {r.failed:3d} ({r.fail_pct:5.2f}%)  {shown}")

never = sorted(set(PB_CHECKS) - set(pc.loc[pc.failed > 0, 'check']))
print(f"\nchecks failed by no crystal ligand in either arm: {len(never)} of {len(PB_CHECKS)}")
""")

code(r"""
# Every complex that fails somewhere, with the reason and the diagnostics behind it.
DIAG = ["energy_ratio", "shortest_bond_relative_length", "longest_bond_relative_length",
        "most_extreme_relative_angle", "shortest_noncovalent_relative_distance",
        "smallest_distance_protein", "volume_overlap_protein",
        "smallest_distance_organic_cofactors", "volume_overlap_organic_cofactors"]
DIAG = [c for c in DIAG if c in res.columns]

bad = res[~res.pb_valid].copy()
bad["failed_checks"] = [";".join(c for c in PB_CHECKS if not row[c])
                        for _, row in bad.iterrows()]
cols = ["entry", "arm", "cognate_contamination", "failed_checks"] + DIAG
bad[cols].sort_values(["entry", "arm"]).to_csv(OUTDIR / "xray_failing_entries.csv", index=False)

print(f"complex-arm pairs failing at least one check: {len(bad)} of {len(res)}\n")
if len(bad):
    show = bad[["entry", "arm", "cognate_contamination", "failed_checks"]]
    show = show.sort_values(["entry", "arm"])
    print(show.to_string(index=False, max_colwidth=64))
else:
    print("none")
print(f"\nfull table with diagnostics: {(OUTDIR / 'xray_failing_entries.csv').relative_to(ROOT)}")
""")

# =============================================================================
md(r"""
## 7. Resolving every failure to the atom pair that caused it

A failed check names a test, not a cause. The four intermolecular distance
checks each record the atom pair that came closest, so each failure can be
resolved to a named residue in the protein file and read for what it is.

There is a trap in reading those four checks as four independent findings. The
atom-type filter that is supposed to restrict each check to its own group does
not remove organic cofactors when the protein file carries residue information:
the branch that would drop them fires only when residue information is absent.
So a single close contact with an organic cofactor is counted in the protein
bucket and in the water bucket as well as its own, and one contact can fail
three checks at once. The cell below does not take that on trust. It compares
the reported atom pairs across the four groups and says outright when they are
the same pair.
""")

code(r"""
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

DIST_CHECKS = {
    "minimum_distance_to_protein": "protein",
    "minimum_distance_to_organic_cofactors": "organic_cofactors",
    "minimum_distance_to_inorganic_cofactors": "inorganic_cofactors",
    "minimum_distance_to_waters": "waters",
}
CLASH_CUTOFF = 0.75          # dock.yml, identical for all four distance checks

_prot_cache = {}


def residue_of(receptor_path, idx):
    # Load exactly as PoseBusters does, so atom indices are the same ones it reported.
    if receptor_path not in _prot_cache:
        _prot_cache[receptor_path] = Chem.MolFromPDBFile(
            receptor_path, sanitize=False, removeHs=False, proximityBonding=False)
    m = _prot_cache[receptor_path]
    if m is None or idx is None or not np.isfinite(idx) or int(idx) >= m.GetNumAtoms():
        return None
    a = m.GetAtomWithIdx(int(idx))
    info = a.GetPDBResidueInfo()
    if info is None:
        return {"element": a.GetSymbol(), "atom_name": "?", "residue": "?",
                "residue_number": -1, "chain": "?", "hetero": None}
    return {"element": a.GetSymbol(), "atom_name": info.GetName().strip(),
            "residue": info.GetResidueName().strip(),
            "residue_number": info.GetResidueNumber(),
            "chain": info.GetChainId().strip(), "hetero": info.GetIsHeteroAtom()}


contacts = []
for _, row in res[~res.pb_valid].iterrows():
    for check, suf in DIST_CHECKS.items():
        if row[check]:
            continue
        idx = row.get(f"most_extreme_protein_atom_id_{suf}")
        part = residue_of(row["receptor_file"], idx)
        contacts.append({
            "entry": row["entry"], "arm": row["arm"], "check": check,
            "ligand_atom_id": row.get(f"most_extreme_ligand_atom_id_{suf}"),
            "ligand_element": row.get(f"most_extreme_ligand_element_{suf}"),
            "partner_atom_id": idx,
            "partner_element": (part or {}).get("element"),
            "partner_atom_name": (part or {}).get("atom_name"),
            "partner_residue": (part or {}).get("residue"),
            "partner_residue_number": (part or {}).get("residue_number"),
            "partner_chain": (part or {}).get("chain"),
            "partner_is_heteroatom": (part or {}).get("hetero"),
            "distance_A": row.get(f"most_extreme_distance_{suf}"),
            "sum_of_radii_A": row.get(f"most_extreme_sum_radii_{suf}"),
            "relative_distance": row.get(f"most_extreme_relative_distance_{suf}"),
            "clash_cutoff": CLASH_CUTOFF,
        })
ct = pd.DataFrame(contacts)
ct.to_csv(OUTDIR / "xray_failing_contacts.csv", index=False)

if ct.empty:
    print("no intermolecular distance check failed in either arm")
else:
    print("every failed intermolecular distance check, resolved to its atom pair\n")
    for r in ct.itertuples():
        print(f"  {r.entry}  [{r.arm}]  {r.check}")
        print(f"      ligand {r.ligand_element} atom {int(r.ligand_atom_id)}"
              f"  to  {r.partner_element} {r.partner_atom_name}"
              f" of {r.partner_residue} {int(r.partner_residue_number)} chain {r.partner_chain}"
              f"   ({'heteroatom record' if r.partner_is_heteroatom else 'polymer record'})")
        print(f"      {r.distance_A:.3f} A against a radius sum of {r.sum_of_radii_A:.2f} A"
              f"  =  {r.relative_distance:.3f} of it, cutoff {r.clash_cutoff}")

    # Is one contact being counted more than once?
    dup = (ct.groupby(["entry", "arm", "ligand_atom_id", "partner_atom_id"])
             .agg(checks=("check", lambda s: ", ".join(sorted(s))), n=("check", "size"))
             .reset_index())
    multi = dup[dup.n > 1]
    print()
    if multi.empty:
        print("each failed check reports a different atom pair; no double counting")
    else:
        for r in multi.itertuples():
            print(f"  {r.entry} [{r.arm}]: ONE contact, atom pair "
                  f"({int(r.ligand_atom_id)}, {int(r.partner_atom_id)}), fails {r.n} checks")
            print(f"      {r.checks}")
        print()
        print("  The same pair in more than one bucket is the classification leak described")
        print("  above, not three independent defects. Counting it once, the crystal cohort")
        print(f"  loses {multi.entry.nunique()} complex(es) to intermolecular contact, not "
              f"{len(ct)} check failures.")
""")

# =============================================================================
md(r"""
## 8. Does the receptor file change the verdict

The two arms are the same 308 ligands measured against two protein files, so
the comparison is paired and the exact McNemar test is the right one. The
question it answers is narrow and load-bearing. If stripping cofactors changes
the crystal ligand's verdict, then a docked pose validated against the stripped
receptor is being held to a different standard than the crystal, and the
ceiling has to be quoted from the prepared arm.
""")

code(r"""
from scipy.stats import binomtest

piv = res.pivot_table(index="entry", columns="arm", values="pb_valid", aggfunc="first")
a, b = piv["crystal"].astype(bool), piv["prepared"].astype(bool)
n01 = int((~a & b).sum())     # crystal fails, prepared passes
n10 = int((a & ~b).sum())     # crystal passes, prepared fails
n11 = int((a & b).sum())
n00 = int((~a & ~b).sum())

print(f"paired on {len(piv)} complexes\n")
print(f"    both arms valid                    : {n11}")
print(f"    valid on the crystal protein only  : {n10}")
print(f"    valid on the prepared receptor only: {n01}")
print(f"    invalid in both arms               : {n00}")

if n01 + n10 == 0:
    print("\nno complex changes verdict between the two receptor files;")
    print("the exact McNemar test is undefined with zero discordant pairs and is not reported.")
else:
    p = binomtest(n10, n01 + n10, 0.5).pvalue
    print(f"\nexact McNemar on {n01 + n10} discordant pairs: p = {p:.4g}")

print("\nper-check difference between the arms (prepared minus crystal, in failures)")
d = (pc.pivot_table(index="check", columns="arm", values="failed")
       .fillna(0).astype(int))
d["delta"] = d["prepared"] - d["crystal"]
moved = d[d.delta != 0].sort_values("delta")
if moved.empty:
    print("    no check changes its failure count between the arms")
else:
    print(moved.to_string())
""")

# =============================================================================
md(r"""
## 9. The ceiling next to the docking arms

The docked totals come from the thesis's own pinned expected-values file, so
this comparison cannot disagree with the appendix table it sits beside. The
per-check decomposition is not in that file and has to be recomputed from the
result table, which is where the arms can go wrong: the engine column is not
the arm. One engine key holds both the raw and the rescored poses, with the
optimiser in a separate column, and the three EquiBind arms are one engine key
split by pocket source and refinement variant. Joining on the engine alone
pools arms and doubles their pose counts.

The reconstruction is therefore gated. Each arm's pose count and overall pass
rate must reproduce the pinned values before any per-check number is reported,
and the cell stops if they do not.

The two bases differ on purpose. The docking rows pool every produced pose of
an arm; the crystal rows are one deposited geometry per complex, because a
crystal conformer has no distribution of ranks to pool over.
""")

code(r"""
import yaml

t24 = yaml.safe_load(EXPECTED.read_text())["table_24"]
assert sorted(c for g in t24["groups"].values() for c in g) == sorted(PB_CHECKS), \
    "the pinned group taxonomy has drifted from this notebook's"

rows = []
for arm in ARMS:
    s_ = sm[(sm.arm == arm) & (sm.cohort == "all 308")].iloc[0]
    rows.append({
        "arm": f"crystal ligand ({arm} receptor)", "basis": "308 complexes",
        "n": int(s_.complexes), "valid_pct": s_.valid_pct,
        "chem_fail_pct": 100 - s_.chemical_pct,
        "intra_fail_pct": 100 - s_.intramolecular_pct,
        "inter_fail_pct": 100 - s_.intermolecular_pct,
    })
for label, (poses, valid, chem, intra, inter, _mind) in t24["all_poses"].items():
    rows.append({"arm": label, "basis": "all produced poses", "n": poses,
                 "valid_pct": valid, "chem_fail_pct": chem,
                 "intra_fail_pct": intra, "inter_fail_pct": inter})
cmp_df = pd.DataFrame(rows)
cmp_df.to_csv(OUTDIR / "xray_vs_docked_comparison.csv", index=False)

hdr = f"{'arm':36s} {'basis':20s} {'n':>6s} {'valid':>8s} {'chem':>7s} {'intra':>7s} {'inter':>7s}"
print(hdr)
print("-" * len(hdr))
for r in cmp_df.itertuples():
    print(f"{r.arm:36s} {r.basis:20s} {r.n:6d} {r.valid_pct:7.1f}% "
          f"{r.chem_fail_pct:6.1f}% {r.intra_fail_pct:6.1f}% {r.inter_fail_pct:6.1f}%")
print()
print("valid is the share passing all twenty-two checks; the three group columns are")
print("the share FAILING somewhere in that group, so they can sum past the invalid share.")
print(f"docked rows read from {EXPECTED.relative_to(ROOT)}, table_24.all_poses")
""")

code(r"""
# Per-check decomposition of the docking arms, recomputed and then gated on the
# pinned totals before any of it is used.
print(f"reading {THESIS_PB.relative_to(ROOT)} ...", flush=True)
td = pd.read_csv(THESIS_PB, low_memory=False)
for c in PB_CHECKS:
    td[c] = td[c].map(strict_bool).astype(bool)
_opt = td["optimizer"].fillna("raw")

ARM_SELECTOR = {
    "AutoDock (raw)":   (td.docking_method == "autodock_mgltools_exh128") & (_opt == "original"),
    "AutoDock + gnina": (td.docking_method == "autodock_mgltools_exh128") & (_opt == "gnina"),
    "DiffDock (raw)":   (td.docking_method == "diffdock") & (_opt == "original"),
    "DiffDock + smina": (td.docking_method == "diffdock") & (_opt == "smina"),
    "DiffDock + gnina": (td.docking_method == "diffdock") & (_opt == "gnina"),
    "EquiBind (raw)":   (td.docking_method == "equibind_guided") & (td.pocket_source == "unguided") & (td.refine_variant == "raw"),
    "EquiBind + smina": (td.docking_method == "equibind_guided") & (td.pocket_source == "unguided") & (td.refine_variant == "smina"),
    "EquiBind + gnina": (td.docking_method == "equibind_guided") & (td.pocket_source == "unguided") & (td.refine_variant == "gnina"),
}
assert set(ARM_SELECTOR) == set(t24["arms"]), "arm labels differ from the pinned set"

recon, per_check_docked = [], []
for label, mask in ARM_SELECTOR.items():
    sub = td[mask]
    n = len(sub)
    valid_pct = 100 * sub[list(PB_CHECKS)].all(axis=1).mean()
    want_n, want_pct = t24["all_poses"][label][0], t24["all_poses"][label][1]
    recon.append({"arm": label, "poses": n, "pinned_poses": want_n,
                  "valid_pct": round(valid_pct, 1), "pinned_valid_pct": want_pct,
                  "reconciles": n == want_n and abs(round(valid_pct, 1) - want_pct) < 0.051})
    for c in PB_CHECKS:
        per_check_docked.append({"arm": label, "check": c, "group": GROUP_OF[c],
                                 "n": n, "failed": int((~sub[c]).sum()),
                                 "fail_pct": 100 * (~sub[c]).mean()})
rc = pd.DataFrame(recon)
print()
print(rc.to_string(index=False))
assert rc.reconciles.all(), (
    "the arm reconstruction does not reproduce the pinned pose counts and pass rates; "
    "the per-check decomposition below would be measuring different arms and is not reported")
print("\nall eight arms reproduce their pinned pose count and pass rate exactly")

pcd = pd.DataFrame(per_check_docked)
crystal_rows = pc[pc.arm == "crystal"].assign(arm="Crystal ligand")[
    ["arm", "check", "group", "n", "failed", "fail_pct"]]
combined = pd.concat([crystal_rows, pcd], ignore_index=True)
combined.to_csv(OUTDIR / "xray_vs_docked_per_check.csv", index=False)

wide = combined.pivot(index="arm", columns="check", values="fail_pct")
ARM_ORDER = ["Crystal ligand"] + list(ARM_SELECTOR)
wide = wide.loc[ARM_ORDER, [c for g in ("chemical", "intramolecular", "intermolecular")
                            for c in CHECK_GROUPS[g]]]
print("\nchecks no arm ever fails, crystal or docked:")
never_any = [c for c in wide.columns if wide[c].max() == 0]
for c in never_any:
    print(f"    {c}")
print(f"  {len(never_any)} of {len(PB_CHECKS)}")
""")


# =============================================================================
md(r"""
## 10. How much headroom the crystal geometries leave

A pass is a threshold crossing, and a threshold says nothing about how close
the crystal geometry sits to it. The bounds are read out of the installed
`dock.yml` rather than typed here, because two of them are not what a reader
would guess: the internal energy ratio is judged against one hundred, not
against the function's own default of seven, and the ensemble it is compared
with holds fifty conformers, not one hundred.

Four of the applied checks compare a single reported number against a fixed
bound, and those four can be drawn with the bound on the axis. The bond length,
bond angle and internal clash checks cannot. They compare each measurement
against a distance-geometry bound computed for that molecule and then widened
by a factor, so there is no one line on a shared axis that separates a pass
from a failure. Those are reported as a table of the widening factors instead
of being drawn with a bound that would not mean anything.
""")

code(r"""
import yaml as _yaml
import posebusters as _pbp

_dock_cfg = _yaml.safe_load((Path(_pbp.__file__).parent / "config/dock.yml").read_text())
_params = {}
for _m in _dock_cfg["modules"]:
    _params.setdefault(_m["name"], _m.get("parameters") or {})

# The four checks whose verdict is one reported number against one fixed bound.
BOUNDED = {
    "energy_ratio": (
        _params["Energy ratio"]["threshold_energy_ratio"], "above",
        "Internal Energy, as a Ratio to the Conformer Ensemble Average"),
    "most_extreme_relative_distance_protein": (
        _params["Distance to protein"]["clash_cutoff"], "below",
        "Closest Approach to the Protein, Relative to the Sum of van der Waals Radii"),
    "volume_overlap_protein": (
        _params["Volume overlap with protein"]["clash_cutoff"], "above",
        "Volume Shared With the Protein, as a Fraction of the Ligand Volume"),
    "aromatic_ring_maximum_distance_from_plane": (
        _params["Ring flatness"]["threshold_flatness"], "above",
        "Aromatic Ring Deviation, From Its Best-Fit Plane (Angstrom)"),
}
BOUNDED = {k: v for k, v in BOUNDED.items() if k in res.columns}

# Checks judged against a per-molecule distance-geometry bound, which no single
# line on a shared axis can represent.
UNBOUNDED_NOTE = {
    "Bond Lengths Within Bounds": ("shortest_bond_relative_length / longest_bond_relative_length",
                                   _params["Geometry"]["threshold_bad_bond_length"]),
    "Bond Angles Within Bounds": ("most_extreme_relative_angle",
                                  _params["Geometry"]["threshold_bad_angle"]),
    "No Internal Steric Clash": ("shortest_noncovalent_relative_distance",
                                 _params["Geometry"]["threshold_clash"]),
}

print("thresholds read from the installed PoseBusters dock configuration")
print(f"  {Path(_pbp.__file__).parent / 'config/dock.yml'}")
print(f"  internal energy ratio judged against "
      f"{_params['Energy ratio']['threshold_energy_ratio']}, over an ensemble of "
      f"{_params['Energy ratio']['ensemble_number_conformations']} conformers")
print()

cry = res[res.arm == "crystal"]
stat_rows = []
for c, (lim, side, _nice) in BOUNDED.items():
    v = pd.to_numeric(cry[c], errors="coerce").dropna()
    beyond = int((v > lim).sum()) if side == "above" else int((v < lim).sum())
    stat_rows.append({
        "diagnostic": c, "n": len(v), "min": v.min(), "q05": v.quantile(0.05),
        "median": v.median(), "q95": v.quantile(0.95), "max": v.max(),
        "bound": lim, "fails_when": side, "beyond_the_bound": beyond,
    })
st = pd.DataFrame(stat_rows)
st.to_csv(OUTDIR / "xray_diagnostic_headroom.csv", index=False)
print("checks with a fixed bound, measured over the 308 crystal ligands\n")
print(st.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

# The count beyond each bound must reconcile with the verdicts already computed.
_link = {"energy_ratio": "internal_energy",
         "most_extreme_relative_distance_protein": "minimum_distance_to_protein",
         "volume_overlap_protein": "volume_overlap_with_protein",
         "aromatic_ring_maximum_distance_from_plane": "aromatic_ring_flatness"}
print()
for r in st.itertuples():
    check = _link[r.diagnostic]
    failed = int((~cry[check]).sum())
    flag = "ok" if failed == r.beyond_the_bound else "MISMATCH"
    print(f"  {flag:8s} {r.diagnostic:44s} beyond the bound {r.beyond_the_bound:3d}"
          f"   {check} failures {failed:3d}")

print()
print("checks judged against a per-molecule distance-geometry bound, no fixed line")
for name, (col, factor) in UNBOUNDED_NOTE.items():
    print(f"  {name:32s} widens the computed bound by {factor}   ({col})")
""")

# =============================================================================
md(r"""
## 11. Figures

Colour follows the thesis house style: each docking tool keeps its own hue, and
the crystal reference is drawn in a neutral dark grey so it never reads as a
fourth tool. Every axis names the quantity it measures and nothing is
abbreviated.
""")

code(r"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 200, "font.size": 9,
    "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
})

CRYSTAL_GREY = "#3f3f3f"
ARM_COLOUR = {"crystal": CRYSTAL_GREY, "prepared": "#9E9E9E"}
TOOL_COLOUR = {"AutoDock": "#1f77b4", "DiffDock": "#ff7f0e", "EquiBind": "#2ca02c"}

NICE_CHECK = {
    "mol_pred_loaded": "Ligand File Loaded",
    "mol_cond_loaded": "Protein File Loaded",
    "sanitization": "Chemical Sanitization",
    "inchi_convertible": "Convertible to an International Chemical Identifier",
    "all_atoms_connected": "All Atoms Connected",
    "no_radicals": "No Radical Electrons",
    "bond_lengths": "Bond Lengths Within Bounds",
    "bond_angles": "Bond Angles Within Bounds",
    "internal_steric_clash": "No Internal Steric Clash",
    "aromatic_ring_flatness": "Aromatic Rings Flat",
    "non-aromatic_ring_non-flatness": "Non-Aromatic Rings Not Flat",
    "double_bond_flatness": "Double Bonds Flat",
    "internal_energy": "Internal Energy Within Bounds",
    "protein-ligand_maximum_distance": "Ligand Not Too Far From the Protein",
    "minimum_distance_to_protein": "Minimum Distance to the Protein",
    "minimum_distance_to_organic_cofactors": "Minimum Distance to Organic Cofactors",
    "minimum_distance_to_inorganic_cofactors": "Minimum Distance to Inorganic Cofactors",
    "minimum_distance_to_waters": "Minimum Distance to Waters",
    "volume_overlap_with_protein": "Volume Overlap With the Protein",
    "volume_overlap_with_organic_cofactors": "Volume Overlap With Organic Cofactors",
    "volume_overlap_with_inorganic_cofactors": "Volume Overlap With Inorganic Cofactors",
    "volume_overlap_with_waters": "Volume Overlap With Waters",
}
NICE_GROUP = {"chemical": "Chemical Validity",
              "intramolecular": "Ligand Internal Geometry",
              "intermolecular": "Ligand Against the Protein"}


def label_panels(axes):
    # (A), (B), (C) ... in the margin outside each panel, kept out of the layout.
    import string
    axes = [ax for ax in np.ravel(np.asarray(axes, dtype=object)) if ax is not None]
    labels = []
    for ax, letter in zip(axes, string.ascii_uppercase):
        t = ax.annotate(f"({letter})", xy=(0, 1), xycoords="axes fraction",
                        xytext=(-6, 6), textcoords="offset points",
                        ha="right", va="bottom", fontweight="bold", fontsize=11,
                        gid="panel_label")
        t.set_in_layout(False)
        labels.append(t)
    return labels


def save(fig, name, extra=()):
    path = FIGDIR / name
    art = list(fig.get_default_bbox_extra_artists()) + list(extra)
    fig.savefig(path, bbox_inches="tight", bbox_extra_artists=art)
    plt.close(fig)
    print(f"    {path.relative_to(ROOT)}")
    return path


print("figure helpers ready")
""")

code(r"""
# Figure 1 - every check, both receptor arms, on an axis wide enough to see one
# complex in three hundred.
order = [c for g in ("chemical", "intramolecular", "intermolecular") for c in CHECK_GROUPS[g]]
y = np.arange(len(order))

worst = 100 - pc.fail_pct.max()
xlo = min(99.0, np.floor(worst * 10) / 10 - 0.1)

fig, ax = plt.subplots(figsize=(8.6, 7.6))
for yy in y:
    ax.plot([xlo, 100], [yy, yy], color="#e6e6e6", lw=0.8, zorder=0)
# The two arms agree on twenty-one of the twenty-two checks, so they are drawn on
# slightly offset rows; otherwise the second marker hides the first entirely.
for arm, marker, size, dy in (("crystal", "o", 42, -0.17), ("prepared", "D", 30, 0.17)):
    d = pc[pc.arm == arm].set_index("check")
    vals = [100 - d.loc[c, "fail_pct"] for c in order]
    ax.scatter(vals, y + dy, s=size, marker=marker, color=ARM_COLOUR[arm], zorder=3,
               edgecolor="white", linewidth=0.8,
               label=("Distributed crystal protein" if arm == "crystal"
                      else "Prepared docking receptor"))
    for yy, v in zip(y, vals):
        if v < 100:
            ax.annotate(f"{v:.2f}%", (v, yy + dy), textcoords="offset points",
                        xytext=(-8, 0), ha="right", va="center", fontsize=7.5,
                        color=ARM_COLOUR[arm])

ax.set_yticks(y)
ax.set_yticklabels([NICE_CHECK[c] for c in order])
ax.invert_yaxis()
ax.set_xlim(xlo, 100 + (100 - xlo) * 0.04)
ax.set_ylim(len(order) - 0.4, -0.6)
ax.set_xlabel("Share of Crystal Ligands Passing the Check (percent)")
ax.set_ylabel("PoseBusters Check")
ax.axvline(100, color="#bbbbbb", lw=0.9, zorder=1)

start = 0
for g in ("chemical", "intramolecular", "intermolecular"):
    n = len(CHECK_GROUPS[g])
    if g == "intramolecular":
        ax.axhspan(start - 0.5, start + n - 0.5, color="#000000", alpha=0.045, zorder=0)
    ax.text(xlo + (100 - xlo) * 0.015, start - 0.30, NICE_GROUP[g],
            fontsize=8, style="italic", color="#555555")
    start += n

leg = ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=2, frameon=False)
ttl = fig.suptitle("Every PoseBusters Check Applied to the Deposited Crystal Ligands\n"
                   "of the 308 Benchmark Complexes", y=1.09)
save(fig, "01_per_check_pass_rate.png", extra=[leg, ttl])
""")

code(r"""
# Figure 2 - the same twenty-two checks, crystal geometry against every docking arm.
M = wide.values.astype(float)
rows_ = list(wide.index)
cols_ = list(wide.columns)

fig, ax = plt.subplots(figsize=(11.6, 4.6))
im = ax.imshow(M, aspect="auto", cmap="magma_r", vmin=0, vmax=100)

for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        if v == 0:
            continue
        ax.text(j, i, ("<0.1" if v < 0.1 else f"{v:.0f}" if v >= 10 else f"{v:.1f}"),
                ha="center", va="center", fontsize=6.6,
                color="white" if v > 55 else "#222222")

ax.set_xticks(np.arange(len(cols_)))
ax.set_xticklabels([NICE_CHECK[c] for c in cols_], rotation=45, ha="right",
                   rotation_mode="anchor", fontsize=7.5)
ax.set_yticks(np.arange(len(rows_)))
ax.set_yticklabels(rows_, fontsize=8.5)
for tick, name in zip(ax.get_yticklabels(), rows_):
    tick.set_color(CRYSTAL_GREY if name == "Crystal ligand"
                   else TOOL_COLOUR[name.split()[0]])
    if name == "Crystal ligand":
        tick.set_fontweight("bold")

# group separators along the check axis
b = 0
for g in ("chemical", "intramolecular"):
    b += len(CHECK_GROUPS[g])
    ax.axvline(b - 0.5, color="white", lw=2.4)
ax.axhline(0.5, color="white", lw=2.4)

cb = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.012)
cb.set_label("Share of Poses Failing the Check (percent)", fontsize=8.5)
ax.set_xlabel("PoseBusters Check")
ax.set_ylabel("Arm")
ttl = fig.suptitle("Which Checks Separate a Docked Pose From the Deposited Crystal Geometry\n"
                   "Unlabelled cells are checks that arm never fails. The crystal row is one geometry\n"
                   "per complex; "
                   "the docking rows pool all produced poses",
                   y=1.10, fontsize=10)
save(fig, "02_crystal_vs_docked_per_check.png", extra=[ttl])
""")

code(r"""
# Figure 3 - headroom between the crystal geometries and the bounds they are judged by.
show = list(BOUNDED)
fig, axes = plt.subplots(1, len(show), figsize=(3.3 * len(show), 4.1))
axes = np.atleast_1d(axes)
for ax, c in zip(axes, show):
    v = pd.to_numeric(res.loc[res.arm == "crystal", c], errors="coerce").dropna()
    lim, side, nice = BOUNDED[c]

    # A few complexes carry an extreme value that would squash the bulk into one
    # bin. Draw the central range and say how many were left outside it.
    lo, hi = float(v.quantile(0.005)), float(v.quantile(0.995))
    span = hi - lo
    if span <= 0:
        lo, hi, span = float(v.min()) - 1.0, float(v.max()) + 1.0, 2.0
    outside = int(((v < lo) | (v > hi)).sum())

    # Bring the bound onto the axis when it is anywhere near the data.
    near = (lo - 0.5 * span) <= lim <= (hi + 0.5 * span)
    xlim = ((min(lo, lim) - 0.06 * span, max(hi, lim) + 0.06 * span) if near
            else (lo - 0.06 * span, hi + 0.06 * span))

    ax.hist(v.clip(*xlim), bins=40, range=xlim, color=CRYSTAL_GREY, alpha=0.85,
            edgecolor="white", linewidth=0.4)
    ax.set_xlim(*xlim)
    # Reserve a clear band at the top so the two annotations never sit on a bar.
    top = ax.get_ylim()[1]
    ax.set_ylim(0, top * 1.38)
    top = ax.get_ylim()[1]
    if near:
        ax.axvline(lim, color="#c0392b", lw=1.7, ls="--")
        # Put the label on whichever side of the line has room inside the axes.
        room_right = xlim[1] - lim
        ha, dx = (("left", 0.015 * span) if room_right > 0.32 * span
                  else ("right", -0.015 * span))
        ax.text(lim + dx, top * 0.98, f"bound {lim:g}", color="#c0392b",
                fontsize=8, va="top", ha=ha)
    else:
        ax.text(0.5, 0.98, f"bound {lim:g}, far outside this range",
                transform=ax.transAxes, ha="center", va="top", fontsize=8, color="#c0392b")
    note = f"a value {side} it fails"
    if outside:
        note += f"\n{outside} complex(es) fall outside the drawn range"
    ax.text(0.5, 0.905, note, transform=ax.transAxes, ha="center", va="top",
            fontsize=7, color="#666666")

    ax.set_xlabel("\n".join(nice.split(", ")), fontsize=8)
    ax.set_ylabel("Number of Complexes")
labels = label_panels(axes)
ttl = fig.suptitle("How Far the Deposited Crystal Geometries Sit From the Bounds They Are Judged By",
                   y=1.05)
fig.tight_layout()
save(fig, "03_diagnostic_headroom.png", extra=labels + [ttl])
""")

code(r"""
# Figure 4 — the crystal ceiling next to the docking arms of the thesis.
lab_order = list(t24["all_poses"].keys())
vals = [t24["all_poses"][k][1] for k in lab_order]


def tool_of(label):
    for t in TOOL_COLOUR:
        if label.startswith(t):
            return t
    return "AutoDock"


def tint(hexcol, f):
    hexcol = hexcol.lstrip("#")
    r, g, b = (int(hexcol[i:i + 2], 16) for i in (0, 2, 4))
    r, g, b = (int(c + (255 - c) * f) for c in (r, g, b))
    return f"#{r:02x}{g:02x}{b:02x}"


colours, seen = [], {}
for lab in lab_order:
    t = tool_of(lab)
    k = seen.get(t, 0)
    seen[t] = k + 1
    colours.append(tint(TOOL_COLOUR[t], 0.30 * k))

fig, ax = plt.subplots(figsize=(9.0, 4.4))
x = np.arange(len(lab_order))
bars = ax.bar(x, vals, color=colours, edgecolor="white", linewidth=0.8)
ax.bar_label(bars, fmt="%.1f", fontsize=8, padding=2)

for arm, ls in (("crystal", "-"), ("prepared", "--")):
    v = sm[(sm.arm == arm) & (sm.cohort == "all 308")].valid_pct.iloc[0]
    ax.axhline(v, color=CRYSTAL_GREY, lw=1.7, ls=ls,
               label=(f"Crystal ligand, distributed protein  ({v:.1f}%)" if arm == "crystal"
                      else f"Crystal ligand, prepared receptor  ({v:.1f}%)"))

ax.set_xticks(x)
ax.set_xticklabels(lab_order, rotation=30, ha="right", rotation_mode="anchor", fontsize=8)
ax.set_ylim(0, 105)
ax.set_xlabel("Docking Arm Reported in the Thesis")
ax.set_ylabel("Share Passing Every PoseBusters Check (percent)")

handles = [Patch(facecolor=TOOL_COLOUR[t], label=t) for t in TOOL_COLOUR]
handles += ax.get_legend_handles_labels()[0]
leg = ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.02),
                ncol=3, frameon=False, fontsize=8)
ttl = fig.suptitle("The Crystal Ceiling Against Every Docking Arm\n"
                   "Docking bars pool all produced poses; the crystal lines are one deposited geometry per complex",
                   y=1.20, fontsize=10)
save(fig, "04_ceiling_vs_docking_arms.png", extra=[leg, ttl])
""")

# =============================================================================
md(r"""
## 12. What this run produced
""")

code(r"""
headline_c = sm[(sm.arm == "crystal") & (sm.cohort == "all 308")].iloc[0]
headline_p = sm[(sm.arm == "prepared") & (sm.cohort == "all 308")].iloc[0]

print("HEADLINE")
print(f"  {headline_c.valid:d} of {headline_c.complexes:d} deposited crystal ligands "
      f"({headline_c.valid_pct:.2f}%) pass all twenty-two PoseBusters checks")
print(f"  against their own distributed protein file.")
print(f"  Against the prepared docking receptor the figure is "
      f"{headline_p.valid:d} of {headline_p.complexes:d} ({headline_p.valid_pct:.2f}%).")
print()
print("FILES")
for p in sorted(OUTDIR.rglob("*")):
    if p.is_file():
        print(f"  {p.relative_to(ROOT)}  ({p.stat().st_size / 1024:.0f} KB)")
""")

# =============================================================================
nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "vina", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12.12"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}
OUT.write_text(json.dumps(nb, indent=1) + "\n")
print(f"wrote {OUT}  ({len(cells)} cells)")
