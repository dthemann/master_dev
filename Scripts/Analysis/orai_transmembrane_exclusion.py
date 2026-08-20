#!/usr/bin/env python3
"""Define the Orai1 transmembrane (TM1 pore) exclusion layer and filter out
docked poses that land inside it.

Biology
-------
Orai1 is a hexameric Ca2+ channel; its pore is lined by the TM1 helix of each of
the 6 subunits (chains A-F here). Going along the pore axis, TM1 runs from the
cytosolic basic region (Arg91, "R91") to the extracellular selectivity filter
(Glu106, "E106").  NOTE: residue 106 in these structures is GLUTAMATE (E106), not
glutamine -- glutamine at that end of TM1 is Q108. The membrane-spanning band
between the R91 ring and the E106 ring is the transmembrane region: it is where
the channel conducts ions, not a druggable surface pocket. Poses that dock there
are therefore not biologically meaningful binding sites and are excluded.

Geometry ("the layer")
----------------------
Per receptor frame we take the six R91 Calpha atoms and the six E106 Calpha atoms
(one per subunit). Their centroids define the two ends of the pore axis:

    p_R91  = mean of the 6 Arg91  CA          (t = 0        , axial_frac = 0)
    p_E106 = mean of the 6 Glu106 CA          (t = L        , axial_frac = 1)
    axis   = unit(p_E106 - p_R91)             (~ the pore / membrane normal)

For any point x, t(x) = (x - p_R91) . axis  is its position along the axis and
radial(x) = |(x - p_R91) - t*axis| its distance from the axis. The exclusion
LAYER is the axial slab  -pad <= t <= L + pad  (a slab perpendicular to the pore
axis, i.e. the membrane-spanning band). r_pore, the mean radial distance of the
12 ring CA atoms, describes the pore/backbone radius and is used only for the
optional "in_pore_lumen" sub-flag; the primary transmembrane test is slab
membership at any radius (a ligand on the lipid-facing TM surface is still in the
membrane).

Because the slab is radially unbounded, drawing it at r_pore would misleadingly
suggest only the conducting channel is excluded. So the layer also records r_outer:
the 99th-percentile radial distance of all protein heavy atoms whose axial position
falls inside the slab -- i.e. how far the hexamer's OUTER SHELL reaches from the axis
in this membrane band (~35 A vs ~8 A for the pore). The visualization (PDB cage /
PyMOL CGO) draws the exclusion boundary out to r_outer, with the pore cylinder kept
as an inner reference and explicit R91 / E106 boundary planes, so the picture matches
the actual test: the whole cross-section of the membrane-spanning band is excluded.

A pose is flagged in_transmembrane when its ligand heavy-atom CENTROID falls in
the slab. frac_atoms_in_layer and majority_in_layer are also reported so a
stricter atom-level rule can be applied without recomputing.

Because the 4 Orai frames (START-Fr0 + MD snapshots Fr300/400/499) are different
conformations in different coordinate systems, the layer is rebuilt per frame from
that frame's own receptor file (the exact `protein_file_used` the poses were
docked into), so poses and layer always share a coordinate system.

Inputs / outputs
----------------
Pose source: the Orai PoseBusters per-pose CSV (`posebusters_filtered_results.csv`),
which lists every generated pose (AutoDock Vina, DiffDock, EquiBind) with its
`pose_file` (per-pose SDF in receptor coords) and `protein_file_used` receptor.

Per dataset out-dir (default posebusters_results/<ds>/transmembrane_filter/):
  tm_exclusion_zones.json     -- per-frame layer parameters (axis, ends, L, r_pore, r_outer)
  <frame>_tm_zone.pdb         -- pseudoatom cage of the band: outer shell + pore walls,
                                 centreline, R91/E106 end-cap rings (load w/ receptor)
  <frame>_tm_zone.pml         -- PyMOL: receptor + translucent outer-shell exclusion
                                 cylinder (r_outer), inner pore reference, R91/E106 planes
  tm_pose_classification.csv  -- EVERY pose + geometry + flags
  tm_excluded_poses.csv       -- poses inside the TM layer (the filtered-out list)
  tm_kept_poses.csv           -- poses that survive the filter
  tm_invalid_poses.csv        -- poses dropped as not placed on the protein
  summary.json                -- counts by method / frame

  fig_tool_tm_loss.png        -- per tool: of all produced poses, how many PB-valid
                                 ones are lost to the TM vs survive (+ tool_tm_loss_summary.csv).
                                 Descriptive only; all inferential tests live in
                                 fig_tool_tm_loss_stats.txt (and tm_stats.json).
  fig_tool_tm_loss_stats.txt  -- the statistics that used to be printed on the figure
                                 panel (per-tool loss + Wilson CI, toolchain x fate
                                 G-test, per-(frame,ligand) paired loss test).
  fig_pbvalid_lost_matrix.png -- Orai frame x ligand matrix of PB-valid poses lost to
                                 the TM (+ pbvalid_lost_matrix.csv)
  fig_rank_survival.png       -- pose rank vs number of PB-valid poses OUTSIDE the TM
                                 (the potentially-valuable poses; + rank_survival.csv)

Cross-dataset combined figures (<out-root>/transmembrane_filter_combined/, built from
both datasets' persisted tm_pose_classification.csv unless --no-combined):
  fig_tool_tm_loss.png         -- two panels, Exp. Ligands (Orai x JKU) and
                                  Benchmark x Orai, of the per-tool PB-valid fate
  fig_tool_tm_survival_box.png -- box+whisker (jittered per-unit points) of TM-filter
                                  survival per (frame x ligand): survival rate + count,
                                  per tool x category
  fig_tool_tm_loss_stats.txt   -- statistics for both combined figures
  fig_rank_recovery_comparison.png   -- how well each tool's NATIVE pose ranking surfaces
                                  valuable poses (PB-valid & outside TM): cumulative
                                  recovery@N per complex (fraction with >=1 valuable pose
                                  in the top-N), Exp. vs Benchmark, Wilson 95% CI +
                                  availability ceiling. Per-complex normalised so the two
                                  datasets are comparable despite the ~100x count gap.
  fig_rank_enrichment_comparison.png -- marginal per-rank hit-rate (is the k-th ranked pose
                                  valuable?) + base-rate line; isolates ranking SKILL from
                                  overall yield (declining = score concentrates valid poses
                                  at the top; flat = validity-blind score)
  fig_rank_comparison_stats.txt      -- recoverable counts/CIs + small-n / pseudoreplication
                                  caveats for the two rank-comparison figures

Toolchains / optimizer variants
-------------------------------
DiffDock and EquiBind each carry raw + optimizer variants (DiffDock raw/smina/gnina;
EquiBind smina/gnina). The per-pose CSVs keep EVERY variant, each tagged with its
`variant` and a display `toolchain` (an asterisk marks a gnina-optimised toolchain:
DiffDock*, EquiBind*). The FIGURES / TABLE / summary report ONE variant per tool,
chosen by --diffdock-variant / --equibind-variant (default gnina), so DiffDock is
not double-counted. EquiBind has no native score, so its poses are ranked by gnina
energy (`gnina_affinity`, most-negative = rank 1) for the rank-survival figure.

Run in the analysis conda env (`vina`):
    conda run -n vina python Scripts/Analysis/orai_transmembrane_exclusion.py
    # gnina toolchains are the default; restrict ligands / pick variants:
    conda run -n vina python Scripts/Analysis/orai_transmembrane_exclusion.py \
        --dataset orai_jku --ligands gsk7975a Synta-66 \
        --diffdock-variant gnina --equibind-variant gnina
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# Shared, unit-tested statistical helpers (Holm/BH, Wilson CI, G-test of
# independence, Cochran-Armitage trend, ...). Robust path insert so the module
# imports whether the script is run from the repo root or its own directory.
try:
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
    import stats_utils as su
    _HAVE_SU = True
except Exception:                                   # pragma: no cover
    su = None
    _HAVE_SU = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    _HAVE_MPL = True
except Exception:                                   # pragma: no cover
    _HAVE_MPL = False

# Shared (A)/(B)/... panel labeller for the multi-panel combined figures.
try:
    from pocket_comparison_report import _label_panels  # noqa: E402
except Exception:                                   # pragma: no cover
    def _label_panels(axes, fontsize: int = 13) -> None:
        import string
        flat = list(axes.ravel()) if hasattr(axes, "ravel") else (
            list(axes) if isinstance(axes, (list, tuple)) else [axes])
        i = 0
        for ax in flat:
            if ax is None:
                continue
            ax.set_title(f"({string.ascii_uppercase[i % 26]})", loc="left",
                         fontweight="bold", fontsize=fontsize)
            i += 1

# Per-tool identity colours (fixed order). Matches the tool palette of the
# pose_comparison 09f_pbvalid_yield_boxplot figure (posebusters_pose_comparison.py
# TOOL_COLORS): AutoDock=tab:blue, DiffDock=tab:orange, EquiBind=tab:green (the
# equibind_unguided green the 09f box uses).
TOOL_STYLE = {
    "autodock": ("AutoDock Vina", "#1f77b4", "o"),
    "diffdock": ("DiffDock", "#ff7f0e", "^"),
    "equibind": ("EquiBind", "#2ca02c", "s"),
}
# Fate of a produced pose (semantic outcome palette; hatch on the "lost" segment
# so the survive/lost distinction is not carried by colour alone — CVD-safe).
_FATE = [
    ("survive",  "PB-valid · outside TM (kept)",       "#55A868", ""),
    ("lost_tm",  "PB-valid · inside TM (lost)",        "#C44E52", "///"),
    ("unplaced", "PB-valid · off-protein/unevaluable", "#DD8452", ""),
    ("nonvalid", "not PB-valid",                       "#CFCFCF", ""),
]


def _frame_sort_key(protein: str) -> int:
    """MD-time ordinal for an Orai frame label (START-Fr0 < MDSnap-Fr300 < Fr400
    < Fr499). Used to order the frames for the Cochran-Armitage trend test."""
    m = re.search(r"Fr(\d+)", str(protein))
    return int(m.group(1)) if m else 10 ** 9


def _tool_key(method: str) -> str:
    m = str(method).lower()
    for k in ("autodock", "diffdock", "equibind"):
        if m.startswith(k):
            return k
    return m


def _pose_variant(method, optimizer, refine_variant, pose_name) -> str:
    """The optimizer/refinement variant of a pose: 'gnina' | 'smina' | 'raw' (and
    'none' for an unoptimised AutoDock pose). DiffDock optimisation is post-hoc
    (`optimizer` col / optimized_<tool>/ path); EquiBind refinement is inline
    (`refine_variant` col); AutoDock optimisation is post-hoc like DiffDock's and
    is carried by the same `optimizer` column, so a run that holds both raw and
    gnina-optimised AutoDock poses keeps them apart instead of merging them into
    one doubled pose cloud."""
    base = _tool_key(method)
    if base == "autodock":
        v = str(optimizer).strip().lower()
        if v in ("gnina", "smina", "gnina_refinement"):
            return v
        n = str(pose_name)
        for tok in ("gnina_refinement", "gnina", "smina"):
            if f"optimized_{tok}" in n or n.endswith(f"_{tok}.sdf"):
                return tok
        return "none"
    if base == "diffdock":
        v = str(optimizer).strip().lower()
        if v in ("gnina", "smina"):
            return v
        n = str(pose_name)
        if "optimized_gnina" in n or n.endswith("_gnina.sdf"):
            return "gnina"
        if "optimized_smina" in n or n.endswith("_smina.sdf"):
            return "smina"
        return "raw"                          # 'original' / NaN
    if base == "equibind":
        v = str(refine_variant).strip().lower()
        return v if v in ("gnina", "smina") else "raw"
    return "none"


def _toolchain_label(base: str, variant: str) -> str:
    """Display label; an asterisk marks a gnina-optimised toolchain (DiffDock*,
    EquiBind*). Raw = the plain tool name; smina named explicitly."""
    name = TOOL_STYLE.get(base, (str(base).title(),))[0]
    if base == "autodock" or variant in ("none", ""):
        return name
    if variant == "gnina":
        return f"{name}*"
    if variant == "raw":
        return name
    return f"{name} ({variant})"


# Reported AutoDock optimizer variant, set once in main() from --autodock-variant.
# Module-level for the same reason _TOPN_ALLOW is in the sister script: it is a
# cross-cutting reporting filter read by several figure builders, and threading it
# through five signatures would add churn without adding clarity. "all" keeps every
# AutoDock variant, which is the historical behaviour for runs that have only one.
_AD_VARIANT: str = "all"


def _variant_matches(variant: str, selector: str) -> bool:
    """Does a pose variant satisfy a --<tool>-variant selector ('all' = any)."""
    sel = str(selector).lower()
    if sel == "all":
        return True
    if sel in ("raw", "original", "none"):
        # DiffDock/EquiBind label an unrefined pose "raw"; an unoptimised AutoDock
        # pose is labelled "none". Both mean "the tool's own untouched output".
        return variant in ("raw", "none")
    return variant == sel

# Reuse the pipeline's single source of truth for the PoseBusters validity verdict
# (the canonical 20-test allowlist + its bool coercion) so pb_valid here agrees
# EXACTLY with run_posebusters.py / the validity report. The "any boolean-like
# column" heuristic is deliberately NOT used: it sweeps in diagnostic columns such
# as `most_extreme_clash_protein` (always False) and mis-coerces object-dtype
# bool/NaN columns, silently flipping real failures to passes.
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_PB_DIR = _PROJECT_ROOT / "Scripts" / "Docking" / "Posebusters"
for _p in (str(_PROJECT_ROOT), str(_PB_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
try:
    from run_posebusters import CANONICAL_TEST_COLUMNS, coerce_test_cols_to_bool
    _HAVE_CANON = True
except Exception:                                   # pragma: no cover
    CANONICAL_TEST_COLUMNS, coerce_test_cols_to_bool = (), None
    _HAVE_CANON = False

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
R91_RESNUM, R91_RESNAME = 91, "ARG"
E106_RESNUM, E106_RESNAME = 106, "GLU"
RING_ATOM = "CA"

# metadata columns carried into the output CSVs when present.
# ``optimized_rank`` is load-bearing, not decorative: _pose_rank() prefers it for the
# optimised AutoDock arm, and it reads that column off `out` (the carried frame), so
# dropping it here silently demotes every gnina pose back to its pre-optimisation Vina
# order — the exact failure _pose_rank's docstring warns about.
_CARRY = [
    "docking_method", "protein", "ligand", "pose_name", "pose_file",
    "autodock_rank", "optimized_rank", "autodock_affinity", "diffdock_confidence",
    "smina_affinity", "gnina_affinity",
    "refine_variant", "optimizer", "pocket_source", "protein_file_used",
    "protein-ligand_maximum_distance",
]

# Physical-placement gate. A pose only counts as a real binding pose if the ligand
# actually sits on the protein. DiffDock occasionally emits rigid-body-translated
# poses hundreds-to-thousands of Å out of the pocket (internally intact but not
# placed — e.g. the JKU 2abp/boron worst-rank poses). geom_ok (|coord|<=1e4) is too
# loose to catch these, so they would otherwise fall into tm_kept_poses.csv as
# "outside the TM" with absurd axial_frac. We drop them via PoseBusters' own
# protein-ligand_maximum_distance verdict (authoritative), backed by a geometric
# envelope in case that column is absent. Legit poses sit within radial ~58 Å and
# axial_frac ~[-2.2, 2.9]; the corrupt ones run to radial 12000 Å / axial_frac 400.
PLACEMENT_FLAG = "protein-ligand_maximum_distance"
RADIAL_MAX = 80.0                    # Å from the pore axis (receptor envelope + margin)
AXIAL_FRAC_LO, AXIAL_FRAC_HI = -4.0, 5.0   # centroid position along R91→E106 axis

DATASETS = {
    "orai_benchmark": "posebusters_results/orai_benchmark/dock/posebusters_filtered_results.csv",
    "orai_jku": "posebusters_results/orai_jku/dock/posebusters_filtered_results.csv",
}

# Display categories for the cross-dataset (combined) figures. Orai × JKU are the
# experimental ("Exp.") ligands; orai_benchmark are the PoseBuster-benchmark ligands
# docked into Orai. Ordered Exp. first, Benchmark second.
CATEGORY_LABELS = {
    "orai_jku": "Exp. Ligands",
    "orai_benchmark": "Benchmark ligands × Orai",
}
CATEGORY_ORDER = ["orai_jku", "orai_benchmark"]
# category identity colours (distinct from the fate palette used inside each bar)
CATEGORY_COLOUR = {
    "Exp. Ligands":              "#8172B3",   # muted purple
    "Benchmark ligands × Orai":  "#64B5CD",   # muted cyan
}


# ---------------------------------------------------------------------------
# Receptor parsing / zone construction
# ---------------------------------------------------------------------------
def _ring_ca(pdb_path: Path, resnum: int, resname: str, atom: str = RING_ATOM) -> np.ndarray:
    """All `atom` coords for residue (resname, resnum) across every chain."""
    pts = []
    for ln in pdb_path.read_text(errors="ignore").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        if ln[12:16].strip() != atom:
            continue
        if ln[17:20].strip() != resname:
            continue
        try:
            if int(ln[22:26]) != resnum:
                continue
            pts.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        except ValueError:
            continue
    return np.asarray(pts, dtype=float)


def _pca_axis(pdb_path: Path) -> Optional[np.ndarray]:
    """Principal (largest-variance) axis of all CA atoms -- a sanity reference for
    the R91->E106 axis (they should be near-parallel for a channel)."""
    xs = []
    for ln in pdb_path.read_text(errors="ignore").splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "CA":
            try:
                xs.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            except ValueError:
                continue
    if len(xs) < 10:
        return None
    a = np.asarray(xs)
    _, vecs = np.linalg.eigh(np.cov((a - a.mean(0)).T))
    return vecs[:, -1]


def _protein_heavy_xyz(pdb_path: Path) -> np.ndarray:
    """All protein heavy-atom coordinates (ATOM records, hydrogens excluded)."""
    pts = []
    for ln in pdb_path.read_text(errors="ignore").splitlines():
        if not ln.startswith("ATOM"):
            continue
        el = ln[76:78].strip() or ln[12:16].strip().lstrip("0123456789")[:1]
        if el.upper() == "H":
            continue
        try:
            pts.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        except ValueError:
            continue
    return np.asarray(pts, dtype=float)


def build_zone(pdb_path: Path, pad: float) -> Dict:
    """Construct the TM exclusion layer for one receptor frame."""
    r91 = _ring_ca(pdb_path, R91_RESNUM, R91_RESNAME)
    e106 = _ring_ca(pdb_path, E106_RESNUM, E106_RESNAME)
    if len(r91) == 0 or len(e106) == 0:
        raise ValueError(f"{pdb_path}: found {len(r91)} R91 and {len(e106)} E106 CA atoms")
    p_lo = r91.mean(0)                      # R91 end  (axial_frac 0)
    p_hi = e106.mean(0)                     # E106 end (axial_frac 1)
    vec = p_hi - p_lo
    length = float(np.linalg.norm(vec))
    axis = vec / length
    # radial pore radius from the 12 ring CA atoms
    ring = np.vstack([r91, e106])
    proj = (ring - p_lo) @ axis
    radial = np.linalg.norm((ring - p_lo) - np.outer(proj, axis), axis=1)
    r_pore = float(radial.mean())
    # OUTER SHELL radius of the transmembrane band. The exclusion test is slab
    # membership at ANY radius (a ligand on the lipid-facing TM surface is still in
    # the membrane), so the pore radius r_pore describes only the conducting channel,
    # NOT the true radial extent of the excluded band. For an honest visualization we
    # measure how far the hexamer reaches out from the axis WITHIN this axial slab:
    # the radial distance of every protein heavy atom whose axial position falls in
    # [t_lo, t_hi]. r_outer uses the 99th percentile so a few flexible loop tips that
    # happen to sit axially in-band (but are not part of the compact TM bundle) do not
    # balloon the drawn boundary; r_outer_max keeps the absolute extent for reference.
    r_outer = r_outer_max = None
    n_heavy_in_slab = 0
    heavy = _protein_heavy_xyz(pdb_path)
    if len(heavy):
        dh = heavy - p_lo
        th = dh @ axis
        rh = np.linalg.norm(dh - np.outer(th, axis), axis=1)
        in_slab = (th >= -pad) & (th <= length + pad)
        if in_slab.any():
            rs = rh[in_slab]
            n_heavy_in_slab = int(in_slab.sum())
            r_outer = round(float(np.percentile(rs, 99)), 3)
            r_outer_max = round(float(rs.max()), 3)
    # sanity: angle between the residue-defined axis and the CA-PCA axis
    pca = _pca_axis(pdb_path)
    axis_pca_deg = None
    if pca is not None:
        c = abs(float(np.clip(np.dot(axis, pca / np.linalg.norm(pca)), -1, 1)))
        axis_pca_deg = round(math.degrees(math.acos(c)), 2)
    return {
        "receptor_file": str(pdb_path),
        "p_R91": p_lo.tolist(), "p_E106": p_hi.tolist(),
        "axis": axis.tolist(), "length": length, "pad": pad,
        "t_lo": -pad, "t_hi": length + pad,
        "r_pore_mean": r_pore, "r_pore_max": float(radial.max()),
        "r_outer": r_outer, "r_outer_max": r_outer_max,
        "n_heavy_in_slab": n_heavy_in_slab,
        "n_R91_subunits": int(len(r91)), "n_E106_subunits": int(len(e106)),
        "axis_vs_pca_deg": axis_pca_deg,
        "residues": {"lower": f"{R91_RESNAME}{R91_RESNUM}", "upper": f"{E106_RESNAME}{E106_RESNUM}"},
    }


def classify_points(zone: Dict, xyz: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return (t, radial) along/from the pore axis for an (n,3) array of points."""
    p_lo = np.asarray(zone["p_R91"]); axis = np.asarray(zone["axis"])
    d = xyz - p_lo
    t = d @ axis
    radial = np.linalg.norm(d - np.outer(t, axis), axis=1)
    return t, radial


# ---------------------------------------------------------------------------
# Ligand pose parsing (V2000 SDF; heavy atoms only, no RDKit needed)
# ---------------------------------------------------------------------------
def _read_v3000(lines: List[str]) -> List[list]:
    out, in_atom = [], False
    for ln in lines:
        s = ln.strip()
        if s.startswith("M  V30 BEGIN ATOM"):
            in_atom = True; continue
        if s.startswith("M  V30 END ATOM"):
            break
        if in_atom and s.startswith("M  V30"):
            parts = s.split("V30", 1)[1].split()          # idx element x y z aamap ...
            if len(parts) >= 5:
                el = parts[1]
                if el not in ("H", "D"):
                    out.append([float(parts[2]), float(parts[3]), float(parts[4])])
    return out


def _read_heavy_xyz(path: str) -> Optional[np.ndarray]:
    """Heavy-atom coordinates from a single-molecule SDF (V2000 or V3000)."""
    try:
        lines = Path(path).read_text(errors="ignore").splitlines()
        if len(lines) < 4:
            return None
        if "V3000" in lines[3]:
            out = _read_v3000(lines)
        else:
            natoms = int(lines[3][:3])
            out = []
            for ln in lines[4:4 + natoms]:
                el = ln[31:34].strip()
                if el in ("H", "D"):
                    continue
                out.append([float(ln[0:10]), float(ln[10:20]), float(ln[20:30])])
        if not out:
            return None
        return np.asarray(out, dtype=float)
    except (ValueError, IndexError, OSError):
        return None


# worker globals (avoid re-sending zones per task)
_ZONES: Dict[str, Dict] = {}


def _init_worker(zones: Dict[str, Dict]):
    global _ZONES
    _ZONES = zones


def _process_pose(task: Tuple[int, str, str]) -> Dict:
    idx, frame, pose_file = task
    zone = _ZONES.get(frame)
    res = {"_idx": idx, "n_heavy_atoms": 0, "cx": None, "cy": None, "cz": None,
           "axial_t": None, "axial_frac": None, "radial_dist": None,
           "frac_atoms_in_layer": None, "in_transmembrane": None,
           "in_pore_lumen": None, "majority_in_layer": None, "geom_ok": False}
    if zone is None:
        res["geom_error"] = "no_zone_for_frame"
        return res
    xyz = _read_heavy_xyz(pose_file)
    if xyz is None or not np.all(np.isfinite(xyz)) or np.max(np.abs(xyz)) > 1e4:
        res["geom_error"] = "unreadable_or_corrupt"
        return res
    cent = xyz.mean(0)
    t_c, rad_c = classify_points(zone, cent[None, :])
    t_c, rad_c = float(t_c[0]), float(rad_c[0])
    t_atoms, _ = classify_points(zone, xyz)
    in_layer_atoms = (t_atoms >= zone["t_lo"]) & (t_atoms <= zone["t_hi"])
    frac = float(in_layer_atoms.mean())
    in_tm = bool(zone["t_lo"] <= t_c <= zone["t_hi"])
    res.update({
        "n_heavy_atoms": int(len(xyz)),
        "cx": round(float(cent[0]), 3), "cy": round(float(cent[1]), 3),
        "cz": round(float(cent[2]), 3),
        "axial_t": round(t_c, 3),
        "axial_frac": round(t_c / zone["length"], 4),
        "radial_dist": round(rad_c, 3),
        "frac_atoms_in_layer": round(frac, 4),
        "in_transmembrane": in_tm,
        "in_pore_lumen": bool(in_tm and rad_c <= zone["r_pore_mean"]),
        "majority_in_layer": bool(frac >= 0.5),
        "geom_ok": True,
    })
    return res


# ---------------------------------------------------------------------------
# pb_valid derivation
# ---------------------------------------------------------------------------
def _derive_pb_valid(df: pd.DataFrame) -> pd.Series:
    """Canonical PoseBusters validity: passes every canonical test present.
    Uses the exact allowlist + coercion from run_posebusters.py."""
    if not _HAVE_CANON:
        return pd.Series(np.nan, index=df.index)
    checks = [c for c in CANONICAL_TEST_COLUMNS if c in df.columns]
    if not checks:
        return pd.Series(np.nan, index=df.index)
    sub = df[checks].copy()
    coerce_test_cols_to_bool(sub, checks)
    return sub.all(axis=1)


# ---------------------------------------------------------------------------
# Visualization writers
# ---------------------------------------------------------------------------
def _perp_basis(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    tmp = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(axis, tmp); u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    return u, v


def write_zone_pdb(zone: Dict, out: Path, n_ring: int = 48, n_layers: int = 9):
    """Pseudoatom cage tracing the TM exclusion band. The band is an axial slab
    between the R91 and E106 rings, unbounded in radius; here we draw it out to the
    hexamer's OUTER SHELL so the visualization reflects that the exclusion covers the
    whole membrane-spanning cross-section (the lipid-facing surface), not just the
    conducting pore. Elements/chains separate the parts for easy colouring in PyMOL:
      chain O (resn OUT, carbon)   -- outer shell wall  (r_outer)
      chain P (resn POR, oxygen)   -- inner pore wall    (r_pore, reference)
      chain Y (resn AXS, nitrogen) -- pore centreline
      chain B (resn CAP, sulphur)  -- end-cap rings at the R91 and E106 boundary planes
    """
    p_lo = np.asarray(zone["p_R91"]); axis = np.asarray(zone["axis"])
    L = zone["length"]
    r_pore = float(zone["r_pore_mean"])
    r_out = float(zone.get("r_outer") or r_pore)
    u, v = _perp_basis(axis)
    recs = []
    serial = 0

    def _wall(radius: float, el: str, resn: str, ch: str):
        nonlocal serial
        for k in range(n_layers):
            frac = k / (n_layers - 1)
            center = p_lo + frac * L * axis
            for j in range(n_ring):
                th = 2 * math.pi * j / n_ring
                p = center + radius * (math.cos(th) * u + math.sin(th) * v)
                serial += 1
                recs.append((serial, el, resn, ch, serial, p, frac * 100))

    _wall(r_out, "C", "OUT", "O")                 # outer shell wall (the exclusion boundary)
    _wall(r_pore, "O", "POR", "P")                # inner pore wall (reference)
    # pore centreline
    for k in range(n_layers):
        frac = k / (n_layers - 1)
        p = p_lo + frac * L * axis
        serial += 1
        recs.append((serial, "N", "AXS", "Y", serial, p, frac * 100))
    # end-cap rings: concentric rings from the pore out to the shell at the R91 (frac 0)
    # and E106 (frac 1) planes, so the two axial boundaries of the band read as disks
    n_cap_rings = 4
    for frac in (0.0, 1.0):
        center = p_lo + frac * L * axis
        for m in range(n_cap_rings + 1):
            radius = r_pore + (r_out - r_pore) * m / n_cap_rings
            for j in range(n_ring):
                th = 2 * math.pi * j / n_ring
                p = center + radius * (math.cos(th) * u + math.sin(th) * v)
                serial += 1
                recs.append((serial, "S", "CAP", "B", serial, p, frac * 100))
    lines = ["REMARK  Orai1 transmembrane exclusion band (R91 ring <-> E106 ring)",
             f"REMARK  axis {np.round(axis,4).tolist()}  length {L:.2f}  "
             f"r_pore {r_pore:.2f}  r_outer {r_out:.2f}"]
    for serial, el, resn, ch, resi, p, b in recs:
        lines.append(
            f"HETATM{serial % 100000:>5d} {el:>2s}   {resn:>3s} {ch}{resi % 10000:>4d}    "
            f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00{b:6.2f}          {el:>2s}"
        )
    lines.append("END")
    out.write_text("\n".join(lines) + "\n")


def write_zone_pml(zone: Dict, frame: str, out: Path):
    p_lo = np.asarray(zone["p_R91"]); p_hi = np.asarray(zone["p_E106"])
    axis = np.asarray(zone["axis"])
    r_pore = float(zone["r_pore_mean"])
    r_out = float(zone.get("r_outer") or r_pore)
    receptor = zone["receptor_file"]
    zone_pdb = out.with_suffix(".pdb").name
    # thin end-cap disks (short, fat CGO cylinders) marking the R91 and E106 planes,
    # each spanning the full band cross-section (radius r_outer) so the two axial
    # boundaries of the excluded band are explicit, colour-matched to the residues.
    cap = 0.5
    r91_a = p_lo - 0.5 * cap * axis;  r91_b = p_lo + 0.5 * cap * axis
    e106_a = p_hi - 0.5 * cap * axis; e106_b = p_hi + 0.5 * cap * axis
    # The whole scene is wrapped in a `python ... python end` block so PyMOL runs it as
    # genuine Python: the multi-line CGO list literals with injected coordinates parse
    # reliably (the .pml command interpreter mishandles multi-line list continuation).
    txt = f"""# PyMOL: Orai1 transmembrane exclusion band for {frame}
# usage:  pymol {out.name}
# The excluded band is the axial slab between the R91 and E106 CA rings, unbounded in
# radius. It is drawn out to the hexamer OUTER SHELL (r_outer={r_out:.1f} A) so it is
# clear the exclusion covers the whole membrane-spanning cross-section, not just the
# pore (r_pore={r_pore:.1f} A, drawn as an inner reference).
python
from pymol import cmd
from pymol.cgo import CYLINDER

cmd.load(r"{receptor}", "receptor")
cmd.hide("everything", "receptor")
cmd.show("cartoon", "receptor")
cmd.color("grey80", "receptor")
cmd.set("cartoon_transparency", 0.15, "receptor")
cmd.select("R91", "receptor and resi 91 and resn ARG")
cmd.select("E106", "receptor and resi 106 and resn GLU")
cmd.show("sticks", "R91 or E106")
cmd.color("blue", "R91")
cmd.color("red", "E106")

# OUTER shell exclusion boundary -- spans the whole hexamer cross-section, R91 -> E106
tm_outer = [CYLINDER, {p_lo[0]:.3f}, {p_lo[1]:.3f}, {p_lo[2]:.3f}, {p_hi[0]:.3f}, {p_hi[1]:.3f}, {p_hi[2]:.3f}, {r_out:.3f}, 1.0, 0.55, 0.0, 1.0, 0.55, 0.0]
cmd.load_cgo(tm_outer, "TM_outer_shell")
cmd.set("cgo_transparency", 0.70, "TM_outer_shell")

# INNER pore / conducting channel (reference only)
tm_pore = [CYLINDER, {p_lo[0]:.3f}, {p_lo[1]:.3f}, {p_lo[2]:.3f}, {p_hi[0]:.3f}, {p_hi[1]:.3f}, {p_hi[2]:.3f}, {r_pore:.3f}, 0.60, 0.0, 0.0, 0.60, 0.0, 0.0]
cmd.load_cgo(tm_pore, "TM_pore_channel")
cmd.set("cgo_transparency", 0.35, "TM_pore_channel")

# R91 (blue) and E106 (red) boundary planes -- the axial limits of the excluded band
r91_plane = [CYLINDER, {r91_a[0]:.3f}, {r91_a[1]:.3f}, {r91_a[2]:.3f}, {r91_b[0]:.3f}, {r91_b[1]:.3f}, {r91_b[2]:.3f}, {r_out:.3f}, 0.10, 0.30, 1.0, 0.10, 0.30, 1.0]
cmd.load_cgo(r91_plane, "R91_boundary_plane")
cmd.set("cgo_transparency", 0.80, "R91_boundary_plane")
e106_plane = [CYLINDER, {e106_a[0]:.3f}, {e106_a[1]:.3f}, {e106_a[2]:.3f}, {e106_b[0]:.3f}, {e106_b[1]:.3f}, {e106_b[2]:.3f}, {r_out:.3f}, 1.0, 0.10, 0.10, 1.0, 0.10, 0.10]
cmd.load_cgo(e106_plane, "E106_boundary_plane")
cmd.set("cgo_transparency", 0.80, "E106_boundary_plane")

cmd.zoom("receptor")
cmd.set("two_sided_lighting", 1)
# For a pseudoatom cage instead of solid CGO, load {zone_pdb}:
#   cmd.load(r"{zone_pdb}", "tm_cage")   # chains O=outer P=pore Y=axis B=end-caps
python end
"""
    out.write_text(txt)


# ---------------------------------------------------------------------------
# Pose rank + figures
# ---------------------------------------------------------------------------
def _pose_rank(method, pose_name, pose_file, autodock_rank=None,
               optimized_rank=None) -> Optional[int]:
    """Native rank of a pose (1 = top-scored). AutoDock = Vina mode / affinity rank,
    DiffDock = confidence rank, EquiBind = unguided sample index. None if unknown.

    For an optimised AutoDock pose the RE-RANKED order (``optimized_rank``) is the
    tool's actual output order, so it wins over the pre-optimisation Vina rank.
    Using ``autodock_rank`` there would rank gnina poses by the order gnina was
    brought in to replace."""
    m = _tool_key(method)
    pf = str(pose_file)
    name = Path(pf).name if pf and pf != "nan" else str(pose_name)
    if m == "autodock":
        if optimized_rank is not None and pd.notna(optimized_rank):
            try:
                return int(float(optimized_rank))
            except (TypeError, ValueError):
                pass
        if autodock_rank is not None and pd.notna(autodock_rank):
            try:
                return int(float(autodock_rank))
            except (TypeError, ValueError):
                pass
        mt = re.search(r"_(?:model|pose)(\d+)", name)
        return int(mt.group(1)) if mt else None
    if m == "diffdock":
        mt = re.search(r"rank(\d+)", name)
        return int(mt.group(1)) if mt else None
    mt = re.search(r"(?:unguided|guided|pose|rank|model)_?(\d+)", name)
    return int(mt.group(1)) if mt else None


def _apply_equibind_gnina_rank(out: pd.DataFrame) -> int:
    """Override EquiBind pose_rank with a gnina-energy rank (most-negative
    gnina_affinity = rank 1) within each (frame, ligand). EquiBind emits no native
    confidence score, so its unguided sample index is not a quality rank; the gnina
    minimisation energy is. Poses with no gnina_affinity get pose_rank = None (not
    gnina-rankable). Returns the number of poses ranked. Non-EquiBind untouched."""
    if "gnina_affinity" not in out.columns:
        return 0
    eb = out["docking_method"].astype(str).str.startswith("equibind")
    if not eb.any():
        return 0
    gaff = pd.to_numeric(out["gnina_affinity"], errors="coerce")
    n = 0
    for _, idx in out[eb].groupby(["protein", "ligand"]).groups.items():
        sub = gaff.loc[idx]
        ranked = sub.dropna().sort_values(kind="stable")     # ascending: most-negative first
        for r, i in enumerate(ranked.index, 1):
            out.at[i, "pose_rank"] = r
            n += 1
        for i in sub[sub.isna()].index:                      # no gnina energy → unrankable
            out.at[i, "pose_rank"] = None
    return n


def _coerce_bool(s: pd.Series) -> pd.Series:
    """Coerce a possibly-string pb_valid column to real bools (CSV round-trips can
    yield 'True'/'False' object dtype, which breaks `== True` masks). NaN/unknown
    → False (matches the old `== True` semantics, and never invents validity)."""
    if s.dtype == bool:
        return s
    m = {"True": True, "true": True, "1": True, "1.0": True, True: True, 1: True, 1.0: True,
         "False": False, "false": False, "0": False, "0.0": False, False: False, 0: False, 0.0: False}

    def _one(v):
        if isinstance(v, float) and v != v:          # NaN
            return False
        return m.get(v, bool(v) if isinstance(v, (int, float, bool)) else False)
    return s.map(_one).astype(bool)


def _prep_fate(out: pd.DataFrame, dd_variant: str, eb_variant: str):
    """Select the reported optimizer variant per tool, tag every pose with its fate
    (survive / lost_tm / unplaced / nonvalid) and build the plotting scaffolding
    shared by the per-dataset and combined figures.

    Returns (df, pv, st, order, tc_style, base_of, counts, totals), or None when the
    selection is empty. `df` is a copy carrying _base/variant/toolchain/_fate.
    """
    df = out.copy()
    df["_base"] = df["docking_method"].map(_tool_key)
    if "variant" not in df.columns:
        df["variant"] = [_pose_variant(m, o, r, pn) for m, o, r, pn in zip(
            df["docking_method"], df.get("optimizer", ""), df.get("refine_variant", ""),
            df.get("pose_name", ""))]
    if "toolchain" not in df.columns:
        df["toolchain"] = [_toolchain_label(b, v) for b, v in zip(df["_base"], df["variant"])]
    # restrict the REPORTING to the selected optimizer variant per tool (default gnina
    # for the learned tools, --autodock-variant for AutoDock); per-pose CSVs still
    # hold every variant.
    keep = ((df["_base"] == "autodock") & df["variant"].map(lambda v: _variant_matches(v, _AD_VARIANT))
            | ((df["_base"] == "diffdock") & df["variant"].map(lambda v: _variant_matches(v, dd_variant)))
            | ((df["_base"] == "equibind") & df["variant"].map(lambda v: _variant_matches(v, eb_variant))))
    df = df[keep].copy()
    if df.empty:
        return None
    if "pb_valid" in df.columns:
        df["pb_valid"] = _coerce_bool(df["pb_valid"])
    else:
        df["pb_valid"] = False
    pv = df["pb_valid"] == True                                    # noqa: E712
    st = df["analysis_status"]
    df["_fate"] = np.where(~pv, "nonvalid",
                  np.where(st == "kept", "survive",
                  np.where(st == "transmembrane", "lost_tm", "unplaced")))
    # ordered toolchains present + a (colour, marker) per toolchain (colour = base tool)
    base_of = dict(zip(df["toolchain"], df["_base"]))
    order = [tc for base in ("autodock", "diffdock", "equibind")
             for tc in sorted(df.loc[df["_base"] == base, "toolchain"].unique())]
    _EXTRA_MK = ["D", "v", "P", "X", "*"]
    tc_style = {}
    for base in ("autodock", "diffdock", "equibind"):
        for i, tc in enumerate([t for t in order if base_of.get(t) == base]):
            mk = TOOL_STYLE[base][2] if i == 0 else _EXTRA_MK[(i - 1) % len(_EXTRA_MK)]
            tc_style[tc] = (TOOL_STYLE[base][1], mk)
    if not order:
        return None
    counts = {tc: df[df["toolchain"] == tc]["_fate"].value_counts() for tc in order}
    totals = {tc: int((df["toolchain"] == tc).sum()) for tc in order}
    return df, pv, st, order, tc_style, base_of, counts, totals


def _compute_fig1_stats(df: pd.DataFrame, pv: pd.Series, order: List[str],
                        counts: Dict, name: str):
    """Toolchain × fate statistics for the PB-valid-loss figure. Returns
    (f1_dict, small_n). Numbers are identical to the previous in-figure test; they
    now live in the sidecar txt/json instead of being drawn on the panel. Never
    raises — a failure degrades to a minimal dict."""
    st = df["analysis_status"]
    f1: Dict = {
        "dataset": name,
        "unit_of_analysis_caveat": (
            "Counts below are individual POSES — pseudoreplicated: each tool emits "
            "many correlated poses per (frame, ligand) unit. The G-test / "
            "Cochran-Armitage p-values are computed on pose counts and therefore "
            "OVERSTATE significance; the honest independent unit is the "
            "(frame, ligand) pair. Per-(frame,ligand) loss-fraction aggregates and "
            "n_units are reported alongside so the pooled result can be read against "
            "the real unit."),
    }
    try:
        _u = df.loc[pv, ["protein", "ligand"]].drop_duplicates()
        n_units = int(len(_u))
        n_ligands = int(df.loc[pv, "ligand"].nunique())
    except Exception:
        n_units = n_ligands = 0
    small_n = (n_units < 12) or (n_ligands < 10)
    f1["n_frame_ligand_units_with_pbvalid"] = n_units
    f1["n_distinct_ligands"] = n_ligands
    f1["small_n_exploratory"] = bool(small_n)
    if not _HAVE_SU:
        return f1, small_n
    _fate_keys = [k for k, _lab, _c, _h in _FATE]     # survive, lost_tm, unplaced, nonvalid
    try:
        f1["unit"] = "poses (pseudoreplicated — see unit_of_analysis_caveat)"
        f1["fate_categories"] = _fate_keys
        f1["toolchains"] = list(order)
        # Wilson CI on the per-tool TM-loss fraction (lost_tm / PB-valid)
        per_tool = {}
        for tc in order:
            nvalid = int(((df["toolchain"] == tc) & pv).sum())
            nlost = int(counts[tc].get("lost_tm", 0))
            lo, hi = su.wilson_ci(nlost, nvalid) if nvalid else (float("nan"), float("nan"))
            per_tool[tc] = {"pb_valid": nvalid, "lost_tm": nlost,
                            "loss_fraction": (nlost / nvalid) if nvalid else None,
                            "wilson_ci": [lo, hi]}
        f1["per_toolchain_loss_fraction"] = per_tool
        # G-test on the toolchain × fate table (drop all-zero rows/cols)
        tbl = np.array([[int(counts[tc].get(k, 0)) for k in _fate_keys]
                        for tc in order], float)
        keep_r = tbl.sum(1) > 0
        keep_c = tbl.sum(0) > 0
        rows_k = [tc for tc, kr in zip(order, keep_r) if kr]
        cols_k = [k for k, kc in zip(_fate_keys, keep_c) if kc]
        sub = tbl[np.ix_(keep_r, keep_c)]
        if sub.shape[0] >= 2 and sub.shape[1] >= 2:
            g = su.gtest_independence(sub)
            resid = np.asarray(g["residuals"])
            f1["gtest"] = {"G": g["G"], "p": g["p"], "df": g["df"],
                           "cramers_v": g["cramers_v"],
                           "n_low_expected": g["n_low_expected"],
                           "rows": rows_k, "cols": cols_k,
                           "adjusted_residuals": resid.tolist()}
            if "lost_tm" in cols_k:
                jc = cols_k.index("lost_tm")
                f1["lost_tm_residual_flagged"] = [rows_k[ir] for ir in range(len(rows_k))
                                                  if abs(float(resid[ir, jc])) > 2.0]
        else:
            f1["gtest"] = {"skipped": "table has <2 non-empty rows/cols"}
        # pseudoreplication-free view: per-(frame,ligand) loss fraction per tool
        agg = {}
        for tc in order:
            s = df[(df["toolchain"] == tc) & pv]
            if s.empty:
                continue
            fr = (s.assign(_l=(s["analysis_status"] == "transmembrane").astype(float))
                    .groupby(["protein", "ligand"])["_l"].mean())
            agg[tc] = {"n_units": int(fr.size),
                       "mean_per_unit_loss_fraction": float(fr.mean()) if fr.size else None,
                       "per_unit_loss_fraction": {f"{a}|{b}": float(v)
                                                  for (a, b), v in fr.items()}}
        f1["per_frame_ligand_loss_fraction"] = agg
        # pseudoreplication-FREE paired test across tools on the per-unit loss fraction
        try:
            unit_mat = pd.DataFrame({tc: pd.Series(agg[tc]["per_unit_loss_fraction"])
                                     for tc in order if tc in agg})
            put: Dict = {"test": "friedman + wilcoxon (paired across tools)",
                         "quantity": "per-(frame,ligand) TM-loss fraction",
                         "unit": "(frame,ligand) pair"}
            boot = {}
            for tc in order:
                vals = np.array(list(agg.get(tc, {}).get("per_unit_loss_fraction", {}).values()), float)
                if vals.size:
                    est, blo, bhi = su.bootstrap_ci(vals, statistic=np.mean)
                    boot[tc] = {"mean": float(est), "ci": [float(blo), float(bhi)],
                                "n_units": int(vals.size)}
            put["per_tool_mean_ci"] = boot
            complete = unit_mat.dropna()
            put["n_units_complete"] = int(len(complete))
            if len(complete) >= 5 and complete.shape[1] >= 3:
                put.update(su.paired_continuous(complete, labels=list(complete.columns)))
            elif len(complete) >= 5 and complete.shape[1] == 2:
                a, b = list(complete.columns)
                rb, p, npair = su.wilcoxon_rankbiserial(complete[a].to_numpy(float),
                                                        complete[b].to_numpy(float))
                put["pairwise"] = [{"a": a, "b": b, "rank_biserial": float(rb),
                                    "p_raw": float(p), "star": su.p_stars(p), "n": int(npair)}]
            else:
                put["note"] = "n too small — exploratory"
            f1["per_frame_ligand_loss_test"] = put
        except Exception as e:                        # pragma: no cover
            f1["per_frame_ligand_loss_test"] = {"error": str(e)}
        if small_n:
            f1["note"] = "exploratory — few independent (frame,ligand) units"
    except Exception as e:                            # pragma: no cover
        f1["error"] = str(e)
    return f1, small_n


def _fig1_stats_text(f1: Dict) -> str:
    """Render the toolchain×fate statistics (previously drawn on the panel) as a
    human-readable sidecar. Consumes the dict from _compute_fig1_stats."""
    L: List[str] = []
    name = f1.get("dataset", "")
    L.append(f"PB-valid poses lost to the transmembrane vs surviving — statistics")
    L.append(f"dataset / category: {name}")
    L.append("=" * 72)
    L.append("")
    L.append(f1.get("unit_of_analysis_caveat", ""))
    L.append("")
    L.append(f"independent (frame×ligand) units with ≥1 PB-valid pose : "
             f"{f1.get('n_frame_ligand_units_with_pbvalid', '?')}")
    L.append(f"distinct ligands                                       : "
             f"{f1.get('n_distinct_ligands', '?')}")
    if f1.get("small_n_exploratory"):
        L.append("→ EXPLORATORY: too few independent units/ligands for formal inference.")
    L.append("")

    ptl = f1.get("per_toolchain_loss_fraction", {})
    if ptl:
        L.append("Per-toolchain TM loss of PB-valid poses (pose-level, pseudoreplicated):")
        L.append(f"  {'toolchain':16s} {'PB-valid':>9s} {'lost→TM':>8s} {'loss%':>7s}  Wilson-95%-CI")
        for tc, d in ptl.items():
            nv = d.get("pb_valid", 0); nl = d.get("lost_tm", 0)
            lf = d.get("loss_fraction"); ci = d.get("wilson_ci", [None, None])
            lfp = f"{100 * lf:.0f}%" if lf is not None else "n/a"
            cip = (f"[{100 * ci[0]:.0f}–{100 * ci[1]:.0f}]%"
                   if ci and ci[0] is not None and ci[0] == ci[0] else "")
            L.append(f"  {tc:16s} {nv:9d} {nl:8d} {lfp:>7s}  {cip}")
        L.append("")

    g = f1.get("gtest", {})
    if g and "G" in g:
        L.append("Toolchain × fate contingency — G-test of independence (pose-level):")
        line = (f"  G={g['G']:.2f}, df={g['df']}, {su.fmt_p(g['p']) if _HAVE_SU else 'p='+str(g['p'])} "
                f"{su.p_stars(g['p']) if _HAVE_SU else ''}, Cramér's V={g['cramers_v']:.3f}")
        if g.get("n_low_expected"):
            line += f", {g['n_low_expected']} low-expected cell(s)"
        L.append(line)
        flagged = f1.get("lost_tm_residual_flagged", [])
        if flagged:
            L.append(f"  |adjusted residual|>2 on the lost-to-TM cell: {', '.join(flagged)}")
        L.append("")
    elif g.get("skipped"):
        L.append(f"Toolchain × fate G-test skipped: {g['skipped']}")
        L.append("")

    put = f1.get("per_frame_ligand_loss_test", {})
    if put and "error" not in put:
        L.append("Pseudoreplication-free view — per-(frame×ligand) TM-loss fraction, paired across tools:")
        for tc, b in put.get("per_tool_mean_ci", {}).items():
            L.append(f"  {tc:16s} mean loss {100 * b['mean']:.0f}%  "
                     f"CI[{100 * b['ci'][0]:.0f}–{100 * b['ci'][1]:.0f}]%  (n_units={b['n_units']})")
        om = put.get("omnibus")
        if om and om.get("p") == om.get("p"):
            L.append(f"  omnibus Friedman: χ²={om['chi2']:.2f}, "
                     f"{su.fmt_p(om['p']) if _HAVE_SU else 'p='+str(om['p'])} "
                     f"{su.p_stars(om['p']) if _HAVE_SU else ''}, Kendall W={om['kendall_w']:.2f}, n={om['n']}")
        for pr in put.get("pairwise", []) or []:
            star = pr.get("star", su.p_stars(pr.get("p_holm", pr.get("p_raw"))) if _HAVE_SU else "")
            pval = pr.get("p_holm", pr.get("p_raw"))
            L.append(f"  {pr['a']} vs {pr['b']}: rank-biserial={pr['rank_biserial']:.2f}, "
                     f"p={pval:.3g} {star} (n={pr['n']})")
        if put.get("note"):
            L.append(f"  note: {put['note']}")
        L.append("")
    return "\n".join(L).rstrip() + "\n"


def _draw_fate_bars(ax, df: pd.DataFrame, pv: pd.Series, order: List[str],
                    counts: Dict, totals: Dict, title: str, ymax=None,
                    annotate: bool = True) -> None:
    """Descriptive stacked fate bars for one dataset on a supplied axis. Draws the
    survive/lost/unplaced/nonvalid segments + a produced/PB-valid/lost annotation.
    Carries NO inferential statistics (those live in the sidecar txt)."""
    if ymax is None:
        ymax = max(totals.values()) if totals else 1
    x = np.arange(len(order)); bottoms = np.zeros(len(order))
    for key, label, colour, hatch in _FATE:
        vals = np.array([int(counts[tc].get(key, 0)) for tc in order], float)
        ax.bar(x, vals, bottom=bottoms, color=colour, width=0.6, label=label,
               hatch=hatch, edgecolor="white", linewidth=0.7)
        for xi, (v, b) in enumerate(zip(vals, bottoms)):
            if v >= 0.03 * ymax:
                ax.text(xi, b + v / 2, f"{int(v)}", ha="center", va="center",
                        fontsize=8, color="#222")
        bottoms += vals
    if annotate:
        for xi, tc in enumerate(order):
            nvalid = int(((df["toolchain"] == tc) & pv).sum())
            nlost = int(counts[tc].get("lost_tm", 0))
            pct = 100 * nlost / nvalid if nvalid else 0
            ax.text(xi, totals[tc] + 0.015 * ymax,
                    f"produced {totals[tc]}\nPB-valid {nvalid}\n"
                    f"lost to TM {nlost} ({pct:.0f}%)",
                    ha="center", va="bottom", fontsize=8, color="#222")
    ax.set_xticks(x); ax.set_xticklabels(order)
    ax.set_ylabel("Number of poses produced")
    ax.margins(y=0.18)
    ax.set_title(title, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)


def make_figures(out: pd.DataFrame, out_dir: Path, name: str,
                 dd_variant: str = "gnina", eb_variant: str = "gnina") -> None:
    """Three PB-valid-focused figures + their data CSVs (see module docstring)."""
    if not _HAVE_MPL:
        print("  matplotlib unavailable — skipping figures")
        return
    prep = _prep_fate(out, dd_variant, eb_variant)
    if prep is None:
        print("  no poses for the selected variant — skipping figures")
        return
    df, pv, st, order, tc_style, base_of, counts, totals = prep
    ymax = max(totals.values()) if totals else 1

    # ── Statistics: computed once, then written to sidecars (fig_tool_tm_loss_stats.txt
    #    + tm_stats.json). They are no longer drawn on the panel — the figure stays
    #    purely descriptive (per user request). ───────────────────────────────────
    tm_stats: Dict = {"dataset": name}
    f1, small_n = _compute_fig1_stats(df, pv, order, counts, name)
    tm_stats["unit_of_analysis_caveat"] = f1.get("unit_of_analysis_caveat")
    tm_stats["n_frame_ligand_units_with_pbvalid"] = f1.get("n_frame_ligand_units_with_pbvalid")
    tm_stats["n_distinct_ligands"] = f1.get("n_distinct_ligands")
    tm_stats["small_n_exploratory"] = f1.get("small_n_exploratory")
    tm_stats["fig_tool_tm_loss"] = f1
    try:
        (out_dir / "fig_tool_tm_loss_stats.txt").write_text(_fig1_stats_text(f1))
        print(f"  stats → {out_dir/'fig_tool_tm_loss_stats.txt'}")
    except Exception as e:                            # pragma: no cover
        print(f"  [warn] writing fig_tool_tm_loss_stats.txt failed: {e}")

    # ── Figure 1: per-toolchain PB-valid fate (produced → survive vs lost-to-TM).
    #    Descriptive only; no inferential tests, no footnote. ──────────────────────
    fig, ax = plt.subplots(figsize=(1.9 * len(order) + 3.5, 6.2))
    _draw_fate_bars(ax, df, pv, order, counts, totals,
                    f"PB-valid poses lost to the transmembrane vs surviving, per toolchain  "
                    f"({CATEGORY_LABELS.get(name, name)})",
                    ymax=ymax)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=2,
              frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_tool_tm_loss.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ── Figure 2: matrix of PB-valid poses lost to TM, per Orai frame × ligand ──
    lost = df[pv & (st == "transmembrane")]
    mat = lost.groupby(["protein", "ligand"]).size().unstack("ligand", fill_value=0)
    if mat.size == 0:
        mat = pd.DataFrame(0, index=sorted(df["protein"].unique()),
                           columns=sorted(df["ligand"].unique()))
    validm = (df[pv].groupby(["protein", "ligand"]).size().unstack("ligand", fill_value=0)
              .reindex(index=mat.index, columns=mat.columns, fill_value=0))
    MAX_COLS = 24
    col_note = ""
    if mat.shape[1] > MAX_COLS:
        top = mat.sum(0).sort_values(ascending=False).head(MAX_COLS).index
        mat, validm = mat[top], validm[top]
        col_note = f"  (top {MAX_COLS} ligands by loss of {df['ligand'].nunique()})"
    frames, ligs = list(mat.index), list(mat.columns)
    cmap = LinearSegmentedColormap.from_list("loss", ["#f7f7f7", "#C44E52"])
    fig, ax = plt.subplots(figsize=(max(6, 0.62 * len(ligs) + 3),
                                    max(3, 0.6 * len(frames) + 2)))
    im = ax.imshow(mat.values, cmap=cmap, aspect="auto",
                   vmin=0, vmax=max(1, int(mat.values.max())))
    ax.set_xticks(range(len(ligs)))
    ax.set_xticklabels(ligs, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(frames)))
    ax.set_yticklabels([f.replace("Orai1WT-", "") for f in frames], fontsize=9)
    vmax = mat.values.max() if mat.values.size else 1
    for i in range(len(frames)):
        for j in range(len(ligs)):
            nl, nv = int(mat.values[i, j]), int(validm.values[i, j])
            if nv:
                cell = f"{nl}/{nv}\n{100 * nl / nv:.0f}%"
            else:
                cell = f"{nl}"
            ax.text(j, i, cell, ha="center", va="center", linespacing=1.15,
                    fontsize=7.5, color="white" if nl > 0.6 * vmax else "#333")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03).set_label(
        "PB-valid poses lost to transmembrane", fontsize=9)

    # Cochran-Armitage trend: does the fraction of PB-valid poses lost to the TM
    # rise/fall across the ORDERED MD frames (START-Fr0 → Fr300 → Fr400 → Fr499)?
    # Pooled across toolchains + per toolchain. Scores = the MD-time frame numbers.
    ca_title = ""
    if _HAVE_SU:
        try:
            f2: Dict = {"unit": "poses (pseudoreplicated — see unit_of_analysis_caveat)",
                        "frame_order_scores": "Fr number (MD time)"}
            frames_ord = sorted(df.loc[pv, "protein"].unique(), key=_frame_sort_key)
            scores = [_frame_sort_key(f) for f in frames_ord]

            def _ca(mask):
                succ = [int((mask & pv & (df["protein"] == f) & (st == "transmembrane")).sum())
                        for f in frames_ord]
                tot = [int((mask & pv & (df["protein"] == f)).sum()) for f in frames_ord]
                z, p, sgn = su.cochran_armitage(succ, tot, scores=scores)
                return {"frames": [str(f) for f in frames_ord], "scores": scores,
                        "lost": succ, "pb_valid": tot,
                        "z": None if z != z else float(z),
                        "p": None if p != p else float(p), "sign": int(sgn)}

            true_all = pd.Series(True, index=df.index)
            if len(frames_ord) >= 3:
                pooled = _ca(true_all)
                f2["pooled"] = pooled
                per_tc = {}
                for tc in order:
                    per_tc[tc] = _ca(df["toolchain"] == tc)
                f2["per_toolchain"] = per_tc
                if pooled["p"] is not None:
                    trend = ("↑" if pooled["sign"] > 0 else "↓" if pooled["sign"] < 0 else "flat")
                    ca_title = (f"Cochran–Armitage frame trend (pooled): z={pooled['z']:.2f}, "
                                f"{su.fmt_p(pooled['p'])} {su.p_stars(pooled['p'])} {trend}")
            else:
                f2["skipped"] = f"only {len(frames_ord)} ordered frame(s) — no trend test"
            if small_n:
                f2["note"] = "exploratory — few independent (frame,ligand) units"
            tm_stats["fig_pbvalid_lost_matrix"] = f2
        except Exception as e:                        # pragma: no cover
            print(f"  [warn] fig_pbvalid_lost_matrix stats failed: {e}")
            tm_stats["fig_pbvalid_lost_matrix"] = {"error": str(e)}

    _t2 = (f"PB-valid poses lost to the transmembrane, per Orai frame × ligand{col_note}\n"
           f"cell = lost / total PB-valid (% of PB-valid lost)   "
           f"({CATEGORY_LABELS.get(name, name)})")
    if ca_title:
        _t2 += ("\n" + ca_title
                + ("  · exploratory" if small_n else "")
                + "  · poses (pseudorepl.)")
    ax.set_title(_t2, fontsize=10)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_pbvalid_lost_matrix.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ── Figure 3: rank survival — PB-valid poses outside the TM, by rank ──
    # AutoDock = Vina affinity rank, DiffDock = confidence rank, EquiBind = gnina-
    # energy rank (assigned upstream, since EquiBind has no native score).
    rk = df[pv & df["pose_rank"].notna()].copy()
    fig, ax = plt.subplots(figsize=(8, 5))
    rank_rows = []
    for tc in order:
        sub = rk[rk["toolchain"] == tc]
        if sub.empty:
            continue
        sub = sub.assign(pose_rank=sub["pose_rank"].astype(int))
        total = sub.groupby("pose_rank").size()
        valuable = (sub[sub["analysis_status"] == "kept"].groupby("pose_rank").size()
                    .reindex(total.index, fill_value=0))
        colour, mk = tc_style[tc]
        ax.plot(total.index, valuable.values, "-", marker=mk, color=colour, lw=2, ms=6,
                label=tc)
        ax.plot(total.index, total.values, "--", color=colour, lw=1.2, alpha=0.4)
        for rr in total.index:
            rank_rows.append({"toolchain": tc, "rank": int(rr),
                              "pb_valid": int(total[rr]),
                              "pb_valid_outside_tm": int(valuable[rr])})
    ax.set_xlabel("Pose rank (1 = best-scored: AutoDock affinity / DiffDock confidence / "
                  "EquiBind gnina energy)")
    ax.set_ylabel("Number of PB-valid poses outside the transmembrane")
    ax.set_title(f"Potentially valuable poses by rank\n"
                 f"solid = PB-valid & outside TM   ·   dashed = all PB-valid   "
                 f"({CATEGORY_LABELS.get(name, name)})",
                 fontsize=11)
    # cap the x-axis at the AutoDock/DiffDock scoring depth so those tools stay
    # readable; EquiBind's deeper gnina-energy-ranked tail is clipped.
    if not rk.empty:
        ranked = rk[rk["_base"].isin(["autodock", "diffdock"])]["pose_rank"]
        xmax = int((ranked if not ranked.empty else rk["pose_rank"]).astype(int).max())
        ax.set_xlim(0.5, xmax + 0.5)
        ax.set_xticks(range(1, xmax + 1))
    ax.set_ylim(bottom=0)
    ax.grid(True, axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(title="Toolchain", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_rank_survival.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ── supporting data CSVs ──
    rows = []
    for tc in order:
        sub = df[df["toolchain"] == tc]; c = sub["_fate"].value_counts()
        nvalid = int((sub["pb_valid"] == True).sum())                # noqa: E712
        rows.append({"toolchain": tc, "produced": int(len(sub)), "pb_valid": nvalid,
                     "pbvalid_survive_outside_tm": int(c.get("survive", 0)),
                     "pbvalid_lost_to_tm": int(c.get("lost_tm", 0)),
                     "pbvalid_offprotein_or_uneval": int(c.get("unplaced", 0)),
                     "pct_pbvalid_lost_to_tm": round(100 * c.get("lost_tm", 0) / nvalid, 1)
                     if nvalid else 0.0})
    pd.DataFrame(rows).to_csv(out_dir / "tool_tm_loss_summary.csv", index=False)
    mat.to_csv(out_dir / "pbvalid_lost_matrix.csv")
    pd.DataFrame(rank_rows).to_csv(out_dir / "rank_survival.csv", index=False)

    # fig_rank_survival: descriptive only (no inferential test drawn on the figure,
    # per validation plan §8 — a tool-vs-tool AUC permutation would need many more
    # independent units than Orai x JKU offers). Record the per-rank survival counts
    # so the descriptive claim is recoverable from the sidecar.
    if _HAVE_SU:
        try:
            tm_stats["fig_rank_survival"] = {
                "test": "none — descriptive",
                "note": ("PB-valid poses outside the TM by native pose rank; solid = "
                         "outside TM, dashed = all PB-valid. See rank_survival.csv."),
                "rows": rank_rows}
        except Exception as e:                        # pragma: no cover
            print(f"  [warn] fig_rank_survival stats note failed: {e}")

    # sidecar: full numeric statistical results next to the figures
    if _HAVE_SU:
        try:
            def _san(o):
                if isinstance(o, dict):
                    return {k: _san(v) for k, v in o.items()}
                if isinstance(o, (list, tuple)):
                    return [_san(v) for v in o]
                if isinstance(o, (np.floating,)):
                    return float(o)
                if isinstance(o, (np.integer,)):
                    return int(o)
                if isinstance(o, np.ndarray):
                    return o.tolist()
                return o
            (out_dir / "tm_stats.json").write_text(json.dumps(_san(tm_stats), indent=2))
            print(f"  stats → {out_dir/'tm_stats.json'}")
        except Exception as e:                        # pragma: no cover
            print(f"  [warn] writing tm_stats.json failed: {e}")

    print("  figures → fig_tool_tm_loss.png, fig_pbvalid_lost_matrix.png, fig_rank_survival.png")


# ---------------------------------------------------------------------------
# Cross-dataset (combined) figures: Exp. Ligands vs Benchmark × Orai
# ---------------------------------------------------------------------------
def _survival_units(df: pd.DataFrame, pv: pd.Series) -> Dict[str, pd.DataFrame]:
    """Per-(frame×ligand) TM-survival of PB-valid poses, per toolchain. The TM
    filter only acts on PLACED poses (kept + transmembrane); off-protein poses are
    excluded from the denominator. Returns {toolchain: DataFrame[protein, ligand,
    n_placed, n_survive, rate]}; units with no placed PB-valid pose are dropped."""
    st = df["analysis_status"]
    sub = df[pv & st.isin(["kept", "transmembrane"])].copy()
    if sub.empty:
        return {}
    sub["_keep"] = (sub["analysis_status"] == "kept").astype(float)
    out: Dict[str, pd.DataFrame] = {}
    for tc, g in sub.groupby("toolchain"):
        agg = (g.groupby(["protein", "ligand"])["_keep"]
                 .agg(n_placed="size", n_survive="sum").reset_index())
        agg["rate"] = agg["n_survive"] / agg["n_placed"].clip(lower=1)
        out[str(tc)] = agg
    return out


def _draw_survival_boxes(ax, data_by_cat: Dict[str, Dict[str, np.ndarray]],
                         tools: List[str], cat_labels: List[str], ylabel: str,
                         ymax: Optional[float] = None) -> None:
    """Grouped box+whisker with jittered points. For each tool, one box per category
    side by side (coloured by category). Empty groups are skipped but keep their slot."""
    rng = np.random.default_rng(0)
    n_cat = max(len(cat_labels), 1)
    width = 0.8 / n_cat
    box_data, positions, colours = [], [], []
    for ti, tool in enumerate(tools):
        for ci, cat in enumerate(cat_labels):
            vals = np.asarray(data_by_cat.get(cat, {}).get(tool, []), float)
            vals = vals[~np.isnan(vals)]
            pos = ti + (ci - (n_cat - 1) / 2.0) * width
            colour = CATEGORY_COLOUR.get(cat, "#888888")
            if vals.size:
                box_data.append(vals); positions.append(pos); colours.append(colour)
                # jittered per-unit points. Small-n groups (the Exp. category) get
                # prominent markers — the points, not the quartiles, are the signal
                # there; large-n groups (Benchmark) get small low-alpha dots so they
                # read as a density cloud instead of an overplotted slab.
                jit = (rng.random(vals.size) - 0.5) * width * 0.6
                if vals.size <= 40:
                    ax.scatter(pos + jit, vals, s=13, color=colour, edgecolor="white",
                               linewidth=0.3, alpha=0.9, zorder=3)
                else:
                    ax.scatter(pos + jit, vals, s=4, color=colour, edgecolor="none",
                               alpha=0.22, zorder=2)
    if box_data:
        bp = ax.boxplot(box_data, positions=positions, widths=width * 0.85,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color="#222", linewidth=1.4),
                        whiskerprops=dict(color="#555"), capprops=dict(color="#555"))
        for patch, c in zip(bp["boxes"], colours):
            patch.set_facecolor(c); patch.set_alpha(0.35); patch.set_edgecolor(c)
    ax.set_xticks(range(len(tools)))
    ax.set_xticklabels(tools)
    ax.set_xlim(-0.6, len(tools) - 0.4)
    ax.set_ylabel(ylabel)
    if ymax is not None:
        ax.set_ylim(0, ymax)
    ax.grid(True, axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)


def make_category_figures(class_paths: Dict[str, Path], out_dir: Path,
                          dd_variant: str = "gnina", eb_variant: str = "gnina") -> None:
    """Cross-dataset figures splitting the per-tool TM analysis into the two ligand
    CATEGORIES — Exp. Ligands (Orai × JKU) and Benchmark × Orai — reading each
    dataset's persisted tm_pose_classification.csv:

      fig_tool_tm_loss.png          -- two category panels of the per-tool PB-valid
                                       fate (produced → survive vs lost-to-TM)
      fig_tool_tm_survival_box.png  -- box+whisker (jittered points) of TM-survival
                                       per (frame×ligand) unit, tool × category
      fig_tool_tm_loss_stats.txt    -- all statistics for both figures (moved off the
                                       panels)
    """
    if not _HAVE_MPL:
        print("  matplotlib unavailable — skipping combined figures")
        return
    usecols = ["dataset", "docking_method", "protein", "ligand", "pose_file",
               "gnina_affinity", "refine_variant", "pb_valid", "variant",
               "toolchain", "pose_rank", "analysis_status"]
    cats = []
    for ds in CATEGORY_ORDER:
        p = class_paths.get(ds)
        if not p or not Path(p).exists():
            print(f"  combined: {CATEGORY_LABELS.get(ds, ds)} classification not found ({p}) — skipping it")
            continue
        try:
            raw = pd.read_csv(p, usecols=lambda c: c in usecols, low_memory=False)
        except Exception as e:                        # pragma: no cover
            print(f"  [warn] combined: cannot read {p}: {e}")
            continue
        prep = _prep_fate(raw, dd_variant, eb_variant)
        if prep is None:
            print(f"  combined: no poses for the selected variant in {ds} — skipping it")
            continue
        df, pv, st, order, tc_style, base_of, counts, totals = prep
        f1, small_n = _compute_fig1_stats(df, pv, order, counts, CATEGORY_LABELS[ds])
        cats.append({"ds": ds, "label": CATEGORY_LABELS[ds], "df": df, "pv": pv,
                     "order": order, "counts": counts, "totals": totals,
                     "f1": f1, "small_n": small_n,
                     "survival": _survival_units(df, pv)})
    if not cats:
        print("  combined: no dataset classification CSVs available — skipping combined figures")
        return
    out_dir.mkdir(parents=True, exist_ok=True)

    # canonical tool order across categories (variant already selected)
    tools: List[str] = []
    for c in cats:
        for tc in c["order"]:
            if tc not in tools:
                tools.append(tc)
    _rank = {"AutoDock Vina": 0, "DiffDock": 1, "DiffDock*": 1, "EquiBind": 2, "EquiBind*": 2}
    tools.sort(key=lambda t: (_rank.get(t, 9), t))
    cat_labels = [c["label"] for c in cats]

    # ── Combined figure 1: two category panels of the per-tool PB-valid fate ──
    ncol = len(cats)
    fig, axes = plt.subplots(1, ncol, figsize=(6.2 * ncol, 6.6), squeeze=False)
    axes = axes.ravel()
    for ax, c in zip(axes, cats):
        _draw_fate_bars(ax, c["df"], c["pv"], c["order"], c["counts"], c["totals"],
                        c["label"])
    _label_panels(axes[:len(cats)])
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(_FATE),
               frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("PB-valid poses lost to the transmembrane vs surviving, per toolchain",
                 fontsize=13, y=0.99)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(out_dir / "fig_tool_tm_loss.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ── Combined figure 2: TM-survival per (frame×ligand) unit, tool × category ──
    #    Panel A = survival RATE (normalised — recommended); Panel B = survival COUNT
    #    (the literal "how many survive"). Box+whisker + jittered per-unit points; for
    #    the small-n Exp. category the individual points carry the real information.
    rate_by_cat = {c["label"]: {tc: sv["rate"].to_numpy(float) for tc, sv in c["survival"].items()}
                   for c in cats}
    count_by_cat = {c["label"]: {tc: sv["n_survive"].to_numpy(float) for tc, sv in c["survival"].items()}
                    for c in cats}
    fig, (axR, axC) = plt.subplots(1, 2, figsize=(6.6 * 2, 5.6))
    _draw_survival_boxes(axR, rate_by_cat, tools, cat_labels,
                         "TM-survival rate of PB-valid poses\n(kept / [kept + inside-TM], per frame×ligand)",
                         ymax=1.02)
    axR.set_title("Survival rate (normalised)", fontsize=11)
    cmax = 0.0
    for cat in cat_labels:
        for tc in tools:
            v = count_by_cat.get(cat, {}).get(tc, np.array([]))
            if len(v):
                cmax = max(cmax, float(np.nanmax(v)))
    _draw_survival_boxes(axC, count_by_cat, tools, cat_labels,
                         "PB-valid poses surviving the TM filter\n(count per frame×ligand unit)",
                         ymax=cmax * 1.08 + 1)
    axC.set_title("Survival count", fontsize=11)
    _label_panels([axR, axC])
    cat_handles = [plt.Line2D([0], [0], marker="s", linestyle="none",
                              markerfacecolor=CATEGORY_COLOUR.get(cat, "#888"),
                              markeredgecolor="white", markersize=10, label=cat)
                   for cat in cat_labels]
    fig.legend(handles=cat_handles, loc="lower center", ncol=len(cat_labels),
               frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("TM-filter survival of PB-valid poses, per toolchain "
                 "(each point = one frame × ligand)", fontsize=13, y=1.0)
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    fig.savefig(out_dir / "fig_tool_tm_survival_box.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ── stats sidecar (both figures) ──
    lines: List[str] = []
    lines.append("PB-valid TM loss / survival — statistics (combined Exp. vs Benchmark)")
    lines.append("#" * 72)
    lines.append("")
    for c in cats:
        lines.append(f"########## {c['label']}  ({c['ds']}) ##########")
        lines.append("")
        lines.append(_fig1_stats_text(c["f1"]))
        lines.append("")
    lines.append("=" * 72)
    lines.append("TM-filter survival per (frame × ligand) unit  (box-plot figure)")
    lines.append("rate = PB-valid poses kept / (kept + inside-TM); count = PB-valid poses kept")
    lines.append("=" * 72)
    lines.append("")
    lines.append(f"  {'category':18s} {'tool':16s} {'n_units':>7s} "
                 f"{'median_rate':>11s} {'IQR_rate':>16s} {'median_count':>12s}")
    for c in cats:
        for tc in tools:
            sv = c["survival"].get(tc)
            if sv is None or sv.empty:
                continue
            r = sv["rate"].to_numpy(float); n = sv["n_survive"].to_numpy(float)
            q1, q3 = np.percentile(r, [25, 75])
            lines.append(f"  {c['label']:18s} {tc:16s} {len(sv):7d} "
                         f"{np.median(r):11.2f} [{q1:.2f}–{q3:.2f}]{'':6s} {np.median(n):12.1f}")
    lines.append("")
    # cross-category test per tool (Exp vs Benchmark) on the survival RATE
    if _HAVE_SU and len(cats) == 2:
        a_lab, b_lab = cats[0]["label"], cats[1]["label"]
        lines.append(f"Cross-category Mann–Whitney (survival rate: {a_lab} vs {b_lab}), per tool:")
        for tc in tools:
            sa = cats[0]["survival"].get(tc)
            sb = cats[1]["survival"].get(tc)
            if sa is None or sb is None or sa.empty or sb.empty:
                continue
            res = su.mannwhitney_cliffs(sa["rate"].to_numpy(float), sb["rate"].to_numpy(float))
            if not res:
                continue
            lines.append(f"  {tc:16s} U={res['U']:.0f}, {su.fmt_p(res['p'])} "
                         f"{su.p_stars(res['p'])}, Cliff's δ={res['cliffs_delta']:+.2f} "
                         f"(n_{a_lab.split()[0]}={res['n_a']}, n_{b_lab.split()[0]}={res['n_b']})")
        lines.append("")
        lines.append("NOTE: Exp. Ligands has only ~3 ligands × 4 frames — its boxes summarise very")
        lines.append("few units, so the jittered points (not the quartiles) carry the real signal.")
    try:
        (out_dir / "fig_tool_tm_loss_stats.txt").write_text("\n".join(lines).rstrip() + "\n")
        print(f"  combined stats → {out_dir/'fig_tool_tm_loss_stats.txt'}")
    except Exception as e:                            # pragma: no cover
        print(f"  [warn] writing combined stats txt failed: {e}")
    print(f"  combined figures → {out_dir/'fig_tool_tm_loss.png'}, "
          f"{out_dir/'fig_tool_tm_survival_box.png'}")


# ---------------------------------------------------------------------------
# Cross-dataset RANK-quality comparison: how well does each tool's native pose
# ranking surface valuable poses (PB-valid & outside TM), Exp. Ligands vs
# Benchmark ligands × Orai?  Normalised per-complex views, so the two datasets
# are comparable despite the ~100× complex-count gap (1202 vs 12).
# ---------------------------------------------------------------------------
# dataset marker/linestyle (colour comes from CATEGORY_COLOUR; the redundant
# marker+dash keeps the two series distinguishable without relying on colour).
CATEGORY_MARK = {"orai_jku": ("^", "--"), "orai_benchmark": ("o", "-")}
# panel order for the rank-comparison figures: strongest → weakest overall yield,
# so the (near-)empty EquiBind panel lands last and reads as the deliberate
# "catastrophe" panel rather than as missing data.
_RANK_PANEL_ORDER = {"diffdock": 0, "autodock": 1, "equibind": 2}


def _first_valuable_rank(tool_df: pd.DataFrame, universe: pd.MultiIndex) -> pd.Series:
    """Smallest native rank carrying a VALUABLE pose (PB-valid & outside TM), per
    complex, reindexed onto the full shared `universe` of (protein, ligand) pairs.
    A complex the tool never produced / never got a valuable pose for → np.inf, so
    it counts as a MISS at every rank depth (the full-universe denominator refuses
    to reward a tool for silently dropping hard complexes)."""
    val = tool_df[tool_df["valuable"]]
    if val.empty:
        mr = pd.Series(dtype=float)
    else:
        mr = val.groupby(["protein", "ligand"])["pose_rank"].min()
    return mr.reindex(universe).astype(float).fillna(np.inf)


def _recovery_stats(tool_df: pd.DataFrame, universe: pd.MultiIndex, max_rank: int) -> Dict:
    """Cumulative recovery@N (fraction of the shared universe with ≥1 valuable pose
    among the tool's top-N ranked poses) + Wilson 95% CI, plus the availability
    ceiling (≥1 valuable pose ANYWHERE in the full ranked list, incl. ranks >N)."""
    first = _first_valuable_rank(tool_df, universe)
    n = int(len(first))
    ns = list(range(1, max_rank + 1))
    rec_k = [int((first <= N).sum()) for N in ns]
    lo, hi = [], []
    for k in rec_k:
        a, b = su.wilson_ci(k, n) if (n and _HAVE_SU) else (float("nan"), float("nan"))
        lo.append(a); hi.append(b)
    ceil_k = int(np.isfinite(first).sum())
    ca, cb = su.wilson_ci(ceil_k, n) if (n and _HAVE_SU) else (float("nan"), float("nan"))
    produced = int(tool_df.drop_duplicates(["protein", "ligand"]).shape[0])
    return {"N": ns, "n": n, "k": rec_k,
            "frac": [k / n if n else float("nan") for k in rec_k],
            "lo": lo, "hi": hi, "ceiling_k": ceil_k,
            "ceiling": ceil_k / n if n else float("nan"),
            "ceiling_ci": (ca, cb), "produced_complexes": produced}


def _marginal_stats(tool_df: pd.DataFrame, max_rank: int) -> Dict:
    """Marginal per-rank hit-rate: P(the pose at EXACTLY rank k is valuable | the
    tool produced a pose at rank k) + Wilson CI, and the pooled base rate over
    ranks 1..max_rank. A declining profile = the native score concentrates valuable
    poses at the top (real ranking skill); a flat profile = a validity-blind score."""
    ns, frac, lo, hi, kk, dd = [], [], [], [], [], []
    for k in range(1, max_rank + 1):
        at = tool_df[tool_df["pose_rank"] == k]
        d = int(at.drop_duplicates(["protein", "ligand"]).shape[0])
        h = int(at[at["valuable"]].drop_duplicates(["protein", "ligand"]).shape[0])
        a, b = su.wilson_ci(h, d) if (d and _HAVE_SU) else (float("nan"), float("nan"))
        ns.append(k); kk.append(h); dd.append(d)
        frac.append(h / d if d else float("nan")); lo.append(a); hi.append(b)
    top = tool_df[tool_df["pose_rank"] <= max_rank]
    base_d = int(len(top)); base_h = int(top["valuable"].sum())
    return {"rank": ns, "k": kk, "d": dd, "frac": frac, "lo": lo, "hi": hi,
            "base_rate": base_h / base_d if base_d else float("nan")}


def make_rank_comparison_figures(class_paths: Dict[str, Path], out_dir: Path,
                                 dd_variant: str = "gnina", eb_variant: str = "gnina",
                                 max_rank: int = 10) -> None:
    """Two normalised, per-complex RANK-quality comparison figures + a stats sidecar,
    contrasting Exp. Ligands vs Benchmark ligands × Orai on how well each tool's
    NATIVE pose ranking surfaces valuable poses (PB-valid & outside the TM):

      fig_rank_recovery_comparison.png   -- cumulative recovery@N (≥1 valuable pose in
                                            the top-N ranked poses), per tool, both
                                            datasets, Wilson 95% CI + availability ceiling
      fig_rank_enrichment_comparison.png -- marginal per-rank hit-rate (is the k-th ranked
                                            pose valuable?) + per-series base-rate line;
                                            isolates ranking SKILL from overall yield
      fig_rank_comparison_stats.txt      -- recoverable counts/CIs + honesty caveats
                                            (inferential dataset contrasts kept off-panel)
    """
    if not _HAVE_MPL:
        print("  matplotlib unavailable — skipping rank-comparison figures")
        return
    usecols = ["docking_method", "protein", "ligand", "pb_valid", "variant",
               "toolchain", "pose_rank", "analysis_status"]
    cats: List[Dict] = []
    for ds in CATEGORY_ORDER:
        p = class_paths.get(ds)
        if not p or not Path(p).exists():
            print(f"  rank-cmp: {CATEGORY_LABELS.get(ds, ds)} classification not found ({p}) — skipping it")
            continue
        try:
            raw = pd.read_csv(p, usecols=lambda c: c in usecols, low_memory=False)
        except Exception as e:                            # pragma: no cover
            print(f"  [warn] rank-cmp: cannot read {p}: {e}")
            continue
        prep = _prep_fate(raw, dd_variant, eb_variant)
        if prep is None:
            print(f"  rank-cmp: no poses for the selected variant in {ds} — skipping it")
            continue
        df = prep[0].copy()
        df = df[df["pose_rank"].notna()].copy()
        if df.empty:
            continue
        df["pose_rank"] = df["pose_rank"].astype(int)
        df["valuable"] = (df["pb_valid"] == True) & (df["analysis_status"] == "kept")  # noqa: E712
        universe = (df.drop_duplicates(["protein", "ligand"])
                      .set_index(["protein", "ligand"]).index)
        cats.append({"ds": ds, "label": CATEGORY_LABELS[ds], "df": df,
                     "order": prep[3], "universe": universe, "n": int(len(universe))})
    if not cats:
        print("  rank-cmp: no dataset classification CSVs available — skipping rank-comparison figures")
        return
    out_dir.mkdir(parents=True, exist_ok=True)

    # tools present across categories, ordered strongest → weakest
    tools: List[str] = []
    base_of: Dict[str, str] = {}
    for c in cats:
        cdf = c["df"]
        for tc in c["order"]:
            # `order` (from _prep_fate) is computed before the pose_rank.notna()
            # filter, so a toolchain with only NaN-ranked poses is listed here yet
            # absent from cdf — guard the empty selection so it is skipped rather
            # than crashing (which the blanket try/except would swallow, silently
            # dropping the whole comparison including the valid AutoDock/DiffDock).
            sel = cdf.loc[cdf["toolchain"] == tc, "_base"]
            if sel.empty:
                continue
            if tc not in tools:
                tools.append(tc)
            base_of.setdefault(tc, str(sel.iloc[0]))
    tools.sort(key=lambda t: (_RANK_PANEL_ORDER.get(base_of.get(t, ""), 9), t))

    # per (tool, category) recovery + marginal stats
    rec: Dict[str, Dict[str, Dict]] = {}
    mrg: Dict[str, Dict[str, Dict]] = {}
    for tc in tools:
        rec[tc], mrg[tc] = {}, {}
        for c in cats:
            sub = c["df"][c["df"]["toolchain"] == tc]
            if sub.empty:
                continue
            rec[tc][c["ds"]] = _recovery_stats(sub, c["universe"], max_rank)
            mrg[tc][c["ds"]] = _marginal_stats(sub, max_rank)

    def _legend(fig, extra: str):
        h = [plt.Line2D([0], [0], color=CATEGORY_COLOUR.get(c["label"], "#888"),
                        marker=CATEGORY_MARK[c["ds"]][0], linestyle=CATEGORY_MARK[c["ds"]][1],
                        lw=2, ms=7, label=f"{c['label']} (n={c['n']})") for c in cats]
        fig.legend(handles=h, loc="lower center", ncol=len(cats), frameon=False,
                   fontsize=9.5, bbox_to_anchor=(0.5, 0.055))
        fig.text(0.5, 0.012, extra, ha="center", va="bottom", fontsize=7.6, color="#666")

    xs = list(range(1, max_rank + 1))
    ncol = len(tools)

    # ── Figure A: cumulative recovery@N ──────────────────────────────────────
    fig, axes = plt.subplots(1, ncol, figsize=(5.0 * ncol, 5.2), sharey=True, squeeze=False)
    axes = axes.ravel()
    for ax, tc in zip(axes, tools):
        ax.axvline(1, color="#999", lw=0.8, ls=":", alpha=0.5, zorder=0)  # recovery@1 guide
        for c in cats:
            r = rec[tc].get(c["ds"])
            if not r:
                continue
            colour = CATEGORY_COLOUR.get(c["label"], "#888")
            mk, dash = CATEGORY_MARK[c["ds"]]
            ax.fill_between(r["N"], r["lo"], r["hi"], color=colour, alpha=0.16,
                            hatch=("///" if c["ds"] == "orai_jku" else None),
                            edgecolor=colour, linewidth=0)
            ax.plot(r["N"], r["frac"], color=colour, marker=mk, ls=dash, lw=1.8, ms=6,
                    zorder=4)
            # availability ceiling (≥1 valuable pose anywhere in the full ranked list).
            # Numeric label only in the EquiBind panel, where the curve sits well below
            # the ceiling (the informative gap); on the near-saturated AutoDock/DiffDock
            # panels the dotted line alone avoids label pile-up (exact counts → sidecar).
            ax.axhline(r["ceiling"], color=colour, lw=1.0, ls=(0, (1, 2)), alpha=0.7)
            if base_of.get(tc) == "equibind":
                _va = "top" if c["ds"] == "orai_jku" else "bottom"
                ax.text(max_rank + 0.15, r["ceiling"], f"ceiling {r['ceiling_k']}/{r['n']}",
                        color=colour, fontsize=7.5, va=_va, ha="left")
            # raw k/n for the tiny Exp. set at recovery@1 (the headline value), so the
            # rate is never disguised as a smooth curve
            if c["ds"] == "orai_jku":
                ax.annotate(f"{r['k'][0]}/{r['n']}", (1, r["frac"][0]),
                            textcoords="offset points", xytext=(-3, 7),
                            ha="right", fontsize=7.5, color=colour, zorder=5)
        base = base_of.get(tc, "")
        # make EquiBind's benchmark zero visceral, not "empty"
        if base == "equibind":
            rb = rec[tc].get("orai_benchmark")
            if rb and rb["ceiling_k"] == 0:
                ax.text(0.5, 0.55,
                        f"Benchmark: 0/{rb['n']} complexes ever valuable\n"
                        "EquiBind geometry catastrophically distorted\n"
                        "(bond angles/lengths fail ~80%)",
                        transform=ax.transAxes, ha="center", va="center", fontsize=8,
                        color=CATEGORY_COLOUR.get("Benchmark ligands × Orai", "#888"),
                        style="italic")
        ax.set_title(tc, fontsize=11)
        ax.set_xlabel("Rank depth N (top-N native-ranked poses)")
        ax.set_xlim(0.5, max_rank + 0.5)
        ax.set_xticks(xs)
        ax.set_ylim(-0.02, 1.03)
        ax.grid(True, axis="y", alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Complexes with ≥1 valuable pose\nin the top-N (fraction)")
    _label_panels(axes[:len(tools)])
    fig.suptitle("Ranking quality: does the native pose rank surface a valuable pose near the top?\n"
                 "cumulative recovery@N   ·   valuable = PoseBusters-valid & outside the "
                 "transmembrane   ·   dotted = availability ceiling (any rank)",
                 fontsize=12.5, y=1.0)
    _legend(fig, "bands = Wilson 95% CI · Exp. Ligands is exploratory (n=12; hatched band, k/n labels) · "
                 "N is a fixed triage budget capped at 10 (EquiBind emits ~30) · overlapping bands ⇒ "
                 "datasets statistically indistinguishable")
    fig.tight_layout(rect=(0, 0.13, 1, 0.94))
    fig.savefig(out_dir / "fig_rank_recovery_comparison.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ── Figure B: marginal per-rank hit-rate (ranking-skill decomposition) ────
    fig, axes = plt.subplots(1, ncol, figsize=(5.0 * ncol, 5.2), sharey=True, squeeze=False)
    axes = axes.ravel()
    for ax, tc in zip(axes, tools):
        for c in cats:
            m = mrg[tc].get(c["ds"])
            if not m:
                continue
            colour = CATEGORY_COLOUR.get(c["label"], "#888")
            mk, dash = CATEGORY_MARK[c["ds"]]
            ax.fill_between(m["rank"], m["lo"], m["hi"], color=colour, alpha=0.16,
                            hatch=("///" if c["ds"] == "orai_jku" else None),
                            edgecolor=colour, linewidth=0)
            ax.plot(m["rank"], m["frac"], color=colour, marker=mk, ls=dash, lw=1.8, ms=6,
                    zorder=4)
            ax.axhline(m["base_rate"], color=colour, lw=1.0, ls=(0, (1, 2)), alpha=0.7)
        base = base_of.get(tc, "")
        if base == "equibind":
            rb = rec[tc].get("orai_benchmark")
            if rb and rb["ceiling_k"] == 0:
                ax.text(0.5, 0.55, f"Benchmark: 0/{rb['n']} valuable at any rank",
                        transform=ax.transAxes, ha="center", va="center", fontsize=8.5,
                        color=CATEGORY_COLOUR.get("Benchmark ligands × Orai", "#888"),
                        style="italic")
        ax.set_title(tc, fontsize=11)
        ax.set_xlabel("Native pose rank (1 = best-scored)")
        ax.set_xlim(0.5, max_rank + 0.5)
        ax.set_xticks(xs)
        ax.set_ylim(-0.02, 1.03)
        ax.grid(True, axis="y", alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Pose at exactly this rank is valuable\n(fraction of complexes)")
    _label_panels(axes[:len(tools)])
    fig.suptitle("Ranking SKILL: is the pose at each rank position valuable?\n"
                 "marginal per-rank hit-rate   ·   dotted = each series' base rate over ranks 1–10   ·   "
                 "declining ⇒ score concentrates valid poses at the top; flat ⇒ validity-blind",
                 fontsize=12.5, y=1.0)
    _legend(fig, "bands = Wilson 95% CI · Exp. Ligands is exploratory (n=12; hatched band) · "
                 "points above the dotted base-rate line = enrichment at that rank")
    fig.tight_layout(rect=(0, 0.13, 1, 0.94))
    fig.savefig(out_dir / "fig_rank_enrichment_comparison.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # ── stats sidecar (tests kept off the panels) ────────────────────────────
    L: List[str] = []
    L.append("Rank-quality comparison — Exp. Ligands vs Benchmark ligands × Orai")
    L.append("valuable pose = PoseBusters-valid AND outside the transmembrane conduction pore")
    L.append("denominator = full shared universe of attempted (frame × ligand) complexes per")
    L.append("dataset (a complex the tool produced no pose for counts as a MISS).")
    L.append("#" * 74)
    L.append("")
    L.append("RECOVERY@N  (fraction of complexes with ≥1 valuable pose in the top-N ranked)")
    L.append("-" * 74)
    L.append(f"  {'category':26s} {'tool':16s} {'n':>5s} {'produced':>8s} "
             f"{'rec@1':>14s} {'rec@10':>14s} {'ceiling':>14s}")
    for tc in tools:
        for c in cats:
            r = rec[tc].get(c["ds"])
            if not r:
                continue
            def _ci(k, lo, hi):
                return f"{k}/{r['n']} [{lo*100:.0f}-{hi*100:.0f}%]"
            L.append(f"  {c['label']:26s} {tc:16s} {r['n']:5d} {r['produced_complexes']:8d} "
                     f"{_ci(r['k'][0], r['lo'][0], r['hi'][0]):>14s} "
                     f"{_ci(r['k'][max_rank-1], r['lo'][max_rank-1], r['hi'][max_rank-1]):>14s} "
                     f"{_ci(r['ceiling_k'], r['ceiling_ci'][0], r['ceiling_ci'][1]):>14s}")
    L.append("")
    L.append("MARGINAL PER-RANK HIT-RATE  (P[pose at rank k is valuable]); base = mean over 1–10")
    L.append("-" * 74)
    L.append(f"  {'category':26s} {'tool':16s} {'base':>6s} {'rank1':>8s} {'rank10':>8s}")
    for tc in tools:
        for c in cats:
            m = mrg[tc].get(c["ds"])
            if not m:
                continue
            def _fp(v):
                return "n/a" if v != v else f"{v*100:.0f}%"
            L.append(f"  {c['label']:26s} {tc:16s} {_fp(m['base_rate']):>6s} "
                     f"{_fp(m['frac'][0]):>8s} {_fp(m['frac'][max_rank-1]):>8s}")
    L.append("")
    L.append("NOTES / CAVEATS")
    L.append("-" * 74)
    L.append("· Exp. Ligands has only ~3 ligands × 4 frames (n=12); its Wilson bands are wide")
    L.append("  and overlap Benchmark at essentially every N — the datasets are statistically")
    L.append("  INDISTINGUISHABLE here; the overlap is the message, not a difference.")
    L.append("· Pseudoreplication: 4 receptor frames × ligand → complexes are correlated, so the")
    L.append("  independent-Bernoulli Wilson band understates the true width. A complex/frame-")
    L.append("  clustered bootstrap CI would be the honest interval; formal dataset contrasts")
    L.append("  (two-proportion / McNemar) are underpowered at n=12 and are deliberately NOT")
    L.append("  drawn on the panels.")
    L.append("· 'produced' = complexes the tool actually emitted a pose for; where it is below n")
    L.append("  the shortfall counts as recovery misses (production coverage, not ranking).")
    try:
        (out_dir / "fig_rank_comparison_stats.txt").write_text("\n".join(L).rstrip() + "\n")
        print(f"  rank-cmp stats → {out_dir/'fig_rank_comparison_stats.txt'}")
    except Exception as e:                                # pragma: no cover
        print(f"  [warn] writing rank-comparison stats txt failed: {e}")
    print(f"  rank-comparison figures → {out_dir/'fig_rank_recovery_comparison.png'}, "
          f"{out_dir/'fig_rank_enrichment_comparison.png'}")


# ---------------------------------------------------------------------------
# Per-dataset driver
# ---------------------------------------------------------------------------
def run_dataset(name: str, csv_path: Path, out_dir: Path, pad: float, workers: int,
                ligands: Optional[List[str]] = None,
                dd_variant: str = "gnina", eb_variant: str = "gnina"):
    print(f"\n=== {name} ===")
    print(f"reading {csv_path}")
    df = pd.read_csv(csv_path, low_memory=False)
    # essential columns
    for req in ("protein", "pose_file", "docking_method", "ligand"):
        if req not in df.columns:
            raise SystemExit(f"{csv_path}: missing required column '{req}'")

    # ── optional ligand restriction (case-insensitive substring match) ──
    if ligands:
        pats = [str(p).lower() for p in ligands]
        lig_lower = df["ligand"].astype(str).str.lower()
        mask = lig_lower.apply(lambda s: any(p in s for p in pats))
        matched = sorted(df.loc[mask, "ligand"].unique())
        if not matched:
            avail = sorted(df["ligand"].unique())
            raise SystemExit(f"--ligands {ligands} matched no ligand in {name}. "
                             f"Available ({len(avail)}): {avail}")
        print(f"  --ligands: keeping {int(mask.sum())} of {len(df)} poses across "
              f"{len(matched)} ligand(s): {matched}")
        df = df[mask].copy()

    n0 = len(df)
    df["pb_valid"] = _derive_pb_valid(df)

    # build a zone per frame from that frame's own receptor file
    out_dir.mkdir(parents=True, exist_ok=True)
    zones: Dict[str, Dict] = {}
    for frame, sub in df.groupby("protein"):
        recs = sub["protein_file_used"].dropna().unique() if "protein_file_used" in df.columns else []
        rec = None
        for cand in recs:
            if Path(cand).exists():
                rec = Path(cand); break
        if rec is None:  # fall back to canonical receptor location
            for cand in (Path(f"Data/Receptors/{frame}.pdb"),
                         Path(f"posebusters_results/_{name}_staging/receptors/{frame}.pdb")):
                if cand.exists():
                    rec = cand; break
        if rec is None:
            print(f"  !! no receptor file for frame {frame}; skipping its poses")
            continue
        z = build_zone(rec, pad)
        zones[str(frame)] = z
        _rout = f"{z['r_outer']:.1f}" if z.get("r_outer") is not None else "n/a"
        print(f"  {frame}: R91 z-end {np.round(z['p_R91'],2).tolist()}  "
              f"E106 z-end {np.round(z['p_E106'],2).tolist()}  L={z['length']:.1f}  "
              f"r_pore={z['r_pore_mean']:.1f}  r_outer={_rout}  "
              f"axis|PCA={z['axis_vs_pca_deg']}deg")
        write_zone_pdb(z, out_dir / f"{frame}_tm_zone.pdb")
        write_zone_pml(z, str(frame), out_dir / f"{frame}_tm_zone.pml")
    (out_dir / "tm_exclusion_zones.json").write_text(json.dumps(
        {"pad": pad, "definition": "axial slab between R91 and E106 CA rings; in_transmembrane = ligand centroid in slab",
         "zones": zones}, indent=2))

    # classify every pose (parallel)
    tasks = [(i, str(r.protein), str(r.pose_file)) for i, r in
             enumerate(df.itertuples(index=False))]
    results: List[Optional[Dict]] = [None] * len(tasks)
    print(f"  classifying {len(tasks)} poses on {workers} workers ...")
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(zones,)) as ex:
        futs = [ex.submit(_process_pose, t) for t in tasks]
        for fut in as_completed(futs):
            r = fut.result()
            results[r["_idx"]] = r
    geom = pd.DataFrame(results).drop(columns=["_idx"])

    carry = [c for c in _CARRY if c in df.columns]
    if "pb_valid" in df.columns and "pb_valid" not in carry:
        carry.append("pb_valid")
    out = pd.concat([df[carry].reset_index(drop=True), geom.reset_index(drop=True)], axis=1)
    out.insert(0, "dataset", name)

    # optimizer/refinement variant + display toolchain (asterisk marks gnina)
    _op = out["optimizer"] if "optimizer" in out.columns else pd.Series([None] * len(out))
    _rv = out["refine_variant"] if "refine_variant" in out.columns else pd.Series([None] * len(out))
    _pn = out["pose_name"] if "pose_name" in out.columns else pd.Series([None] * len(out))
    _ar = out["autodock_rank"] if "autodock_rank" in out.columns else pd.Series([None] * len(out))
    # from `df`, not `out`: row-aligned (geom is built positionally over df) and
    # immune to _CARRY dropping the column again.
    _or = (df["optimized_rank"].reset_index(drop=True) if "optimized_rank" in df.columns
           else pd.Series([None] * len(out)))
    out["variant"] = [_pose_variant(m, o, r, pn) for m, o, r, pn in
                      zip(out["docking_method"], _op, _rv, _pn)]
    out["toolchain"] = [_toolchain_label(_tool_key(m), v) for m, v in
                        zip(out["docking_method"], out["variant"])]

    # native pose rank (AutoDock Vina affinity rank, DiffDock confidence rank);
    # EquiBind has none, so it is (re)ranked by gnina energy just below.
    out["pose_rank"] = [_pose_rank(m, pn, pf, ar, orank)
                        for m, pn, pf, ar, orank in
                        zip(out["docking_method"], _pn, out["pose_file"], _ar, _or)]
    n_gr = _apply_equibind_gnina_rank(out)
    if n_gr:
        print(f"  ranked {n_gr} EquiBind poses by gnina energy (most-negative = rank 1)")

    # ── physical-placement gate (drop poses not actually on the protein) ──────
    if PLACEMENT_FLAG in df.columns:
        raw = df[PLACEMENT_FLAG].reset_index(drop=True)
        flag_map = {"True": True, "true": True, "False": False, "false": False,
                    True: True, False: False, 1: True, 0: False, 1.0: True, 0.0: False}
        pld = raw if raw.dtype == bool else raw.map(flag_map)
        pld_ok = pld.fillna(False).astype(bool)
    else:
        pld_ok = pd.Series(True, index=out.index)          # no PB flag → geometry only
    geo_sane = (out["radial_dist"] <= RADIAL_MAX) & out["axial_frac"].between(AXIAL_FRAC_LO, AXIAL_FRAC_HI)
    out["placed_on_protein"] = out["geom_ok"].astype(bool) & pld_ok & geo_sane.fillna(False)
    out["analysis_status"] = np.where(
        ~out["geom_ok"].astype(bool), "unevaluable",
        np.where(~out["placed_on_protein"], "invalid_placement",
                 np.where(out["in_transmembrane"] == True, "transmembrane", "kept")))  # noqa: E712

    n_bad = int((out["analysis_status"] == "unevaluable").sum())
    n_invalid = int((out["analysis_status"] == "invalid_placement").sum())
    if n_bad:
        print(f"  {n_bad} poses unevaluable (unreadable/corrupt/no-zone)")
    if n_invalid:
        print(f"  {n_invalid} poses dropped — not placed on protein (translated/exploded off the pocket)")

    class_csv = out_dir / "tm_pose_classification.csv"
    out.to_csv(class_csv, index=False)

    placed = out[out["analysis_status"].isin(["transmembrane", "kept"])].copy()
    excluded = out[out["analysis_status"] == "transmembrane"].copy()
    kept = out[out["analysis_status"] == "kept"].copy()
    invalid = out[out["analysis_status"] == "invalid_placement"].copy()
    _sort = ["docking_method", "protein", "ligand"]
    excluded.sort_values(_sort).to_csv(out_dir / "tm_excluded_poses.csv", index=False)
    kept.sort_values(_sort).to_csv(out_dir / "tm_kept_poses.csv", index=False)
    invalid.sort_values(_sort).to_csv(out_dir / "tm_invalid_poses.csv", index=False)

    # summary — computed over PLACED poses only (kept + transmembrane)
    def _by(col, frame=placed):
        res = {}
        for k, v in frame.groupby(col)["analysis_status"]:
            n = int(v.size); exc = int((v == "transmembrane").sum())
            res[str(k)] = {"n": n, "excluded": exc, "pct": round(100 * exc / max(n, 1), 1)}
        return res
    # reporting subset: the selected optimizer variant per tool (default gnina)
    _base = placed["docking_method"].map(_tool_key)
    sel = ((_base == "autodock") & placed["variant"].map(lambda v: _variant_matches(v, _AD_VARIANT))
           | ((_base == "diffdock") & placed["variant"].map(lambda v: _variant_matches(v, dd_variant)))
           | ((_base == "equibind") & placed["variant"].map(lambda v: _variant_matches(v, eb_variant))))
    reported = placed[sel]
    summary = {
        "dataset": name, "n_poses_total": n0,
        "n_evaluated": int(len(placed)),               # placed on protein (kept + TM), ALL variants
        "n_unevaluable": n_bad,
        "n_invalid_placement": n_invalid,              # parsed but not on the protein
        "n_excluded_transmembrane": int(len(excluded)),
        "pct_excluded": round(100 * len(excluded) / max(len(placed), 1), 2),
        "n_in_pore_lumen": int(placed["in_pore_lumen"].sum()),
        "n_excluded_but_pb_valid": int((excluded["pb_valid"] == True).sum())  # noqa: E712
            if "pb_valid" in excluded else None,
        "placement_gate": PLACEMENT_FLAG,
        "reported_diffdock_variant": dd_variant,
        "reported_equibind_variant": eb_variant,
        "by_toolchain": _by("toolchain", reported),    # selected variant, asterisk labels
        "by_method": _by("docking_method"),            # all variants, base tool
        "by_frame": _by("protein"),
        "pad_angstrom": pad,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  -> {len(excluded)}/{len(placed)} placed poses inside the TM layer "
          f"({summary['pct_excluded']}%); {summary['n_in_pore_lumen']} in the pore lumen")
    print(f"  wrote {class_csv}")
    print(f"        {out_dir/'tm_excluded_poses.csv'}  ({len(excluded)})")
    print(f"        {out_dir/'tm_kept_poses.csv'}  ({len(kept)})")
    print(f"        {out_dir/'tm_invalid_poses.csv'}  ({len(invalid)})")
    print(f"  reported toolchains (diffdock={dd_variant}, equibind={eb_variant}):")
    for m, s in summary["by_toolchain"].items():
        print(f"      {m:18s} {s['excluded']:5d}/{s['n']:<5d} ({s['pct']}%) in TM")

    try:
        make_figures(out, out_dir, name, dd_variant, eb_variant)
    except Exception as e:                          # never let plotting kill the run
        print(f"  [warn] figure generation failed: {e}")
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", choices=list(DATASETS) + ["all"], default="all")
    ap.add_argument("--pad", type=float, default=0.0,
                    help="Angstrom padding added beyond each residue ring plane (default 0).")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--ligands", nargs="+", default=None, metavar="NAME",
                    help="restrict the analysis to ligands whose name contains any of "
                         "these (case-insensitive substrings); default = all ligands.")
    ap.add_argument("--autodock-variant", default="all",
                    choices=["all", "none", "original", "raw", "gnina", "smina", "gnina_refinement"],
                    help="AutoDock optimizer variant used in the figures/table/summary "
                         "(default all). Set it whenever the run holds more than one "
                         "AutoDock variant, otherwise raw and optimised poses are pooled "
                         "and every AutoDock count is doubled. CSVs keep all variants.")
    ap.add_argument("--diffdock-variant", default="gnina",
                    choices=["gnina", "smina", "raw", "original", "all"],
                    help="DiffDock optimizer variant used in the figures/table/summary "
                         "(default gnina → labelled 'DiffDock*'). CSVs keep all variants.")
    ap.add_argument("--equibind-variant", default="gnina",
                    choices=["gnina", "smina", "raw", "all"],
                    help="EquiBind refine variant used in the figures/table/summary "
                         "(default gnina → 'EquiBind*', ranked by gnina energy).")
    ap.add_argument("--out-root", default="posebusters_results",
                    help="root under which <dataset>/transmembrane_filter is written")
    ap.add_argument("--combined-out", default=None,
                    help="dir for the cross-dataset (Exp. vs Benchmark) combined figures "
                         "(default <out-root>/transmembrane_filter_combined)")
    ap.add_argument("--no-combined", action="store_true",
                    help="skip the combined Exp.-vs-Benchmark figures")
    args = ap.parse_args(argv)

    global _AD_VARIANT
    _AD_VARIANT = args.autodock_variant

    names = list(DATASETS) if args.dataset == "all" else [args.dataset]
    summaries = []
    for name in names:
        csv_path = Path(DATASETS[name])
        if not csv_path.exists():
            print(f"skip {name}: {csv_path} not found")
            continue
        out_dir = Path(args.out_root) / name / "transmembrane_filter"
        summaries.append(run_dataset(name, csv_path, out_dir, args.pad, args.workers,
                                     ligands=args.ligands,
                                     dd_variant=args.diffdock_variant,
                                     eb_variant=args.equibind_variant))
    # cross-dataset combined figures (Exp. Ligands vs Benchmark × Orai). Reads each
    # dataset's persisted tm_pose_classification.csv, so it builds the full two-category
    # figure whenever both prior runs exist — not only when both were (re)run just now.
    if not args.no_combined:
        class_paths = {ds: Path(args.out_root) / ds / "transmembrane_filter" /
                       "tm_pose_classification.csv" for ds in DATASETS}
        combined_dir = (Path(args.combined_out) if args.combined_out
                        else Path(args.out_root) / "transmembrane_filter_combined")
        print(f"\n=== combined figures (Exp. Ligands vs Benchmark × Orai) → {combined_dir} ===")
        try:
            make_category_figures(class_paths, combined_dir,
                                  args.diffdock_variant, args.equibind_variant)
        except Exception as e:                        # never let combined plotting kill the run
            print(f"  [warn] combined figure generation failed: {e}")
        try:
            make_rank_comparison_figures(class_paths, combined_dir,
                                         args.diffdock_variant, args.equibind_variant)
        except Exception as e:                        # never let combined plotting kill the run
            print(f"  [warn] rank-comparison figure generation failed: {e}")

        # Cross-dataset PB-valid / TM-survival SHARE comparison (normalised so the two
        # ligand sets are directly comparable). Additive companion — needs BOTH datasets'
        # classification tables, so it is guarded on their presence and passed the same
        # refine-variant selection as this run so its toolchains match the summary.
        if all(p.exists() for p in class_paths.values()):
            share_out = Path(args.out_root) / "orai_pbvalid_tm_share_compare"
            print(f"\n=== PB-valid / TM-survival share comparison → {share_out} ===")
            try:
                from orai_pbvalid_tm_share_compare import main as _share_main
                # Pass this run's AutoDock variant through. Both arms now carry a gnina
                # AutoDock variant (the Orai x Benchmark gnina rescoring run landed
                # 2026-08-15), so pinning no longer empties the Benchmark side; the share
                # script still falls back per-dataset if a variant is absent.
                # --top-n-poses 10 is REQUIRED, not cosmetic: Experimental AutoDock and
                # EquiBind emit 30 poses per frame-ligand unit against the Benchmark's 10,
                # so an uncapped run silently triples the Experimental denominators and the
                # two panels stop being comparable.
                _share_argv = ["--results-root", args.out_root,
                               "--out-dir", str(share_out),
                               "--diffdock-variant", args.diffdock_variant,
                               "--equibind-variant", args.equibind_variant,
                               "--top-n-poses", "10"]
                # Only forward a CONCRETE AutoDock variant. This script's default is
                # 'all', which the share script cannot resolve against a dataset holding
                # both 'none' and 'gnina' — it refuses to pool them and raises. Leaving
                # the flag off lets the share script apply its own default instead.
                if args.autodock_variant not in ("all", None):
                    _share_argv += ["--autodock-variant", args.autodock_variant]
                _share_main(_share_argv)
            # SystemExit too: argparse raises it on an unknown flag, and it is not an
            # Exception subclass, so it would otherwise escape this guard and kill the
            # whole run after every figure had already been written.
            except (Exception, SystemExit) as e:      # never let it kill the run
                print(f"  [warn] PB-valid/TM share comparison failed: {e}")
        else:
            missing = [ds for ds, p in class_paths.items() if not p.exists()]
            print(f"  [skip] PB-valid/TM share comparison — need both datasets classified; "
                  f"missing: {missing}")

    print("\n==== DONE ====")
    for s in summaries:
        print(f"{s['dataset']}: {s['n_excluded_transmembrane']}/{s['n_evaluated']} "
              f"({s['pct_excluded']}%) excluded as transmembrane")


if __name__ == "__main__":
    main()
