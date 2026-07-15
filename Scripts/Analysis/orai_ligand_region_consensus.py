"""Consensus binding-region analysis of PB-valid poses on the Orai channel.

Question this answers (Orai × Benchmark cross-docking — benchmark ligands docked
into the Orai1 channel across 4 MD frames × 3 docking methods):

  (a) Are there clusters of PB-valid poses that recur ACROSS ALL DOCKED LIGANDS —
      i.e. binding regions of the channel that many different ligands hit the same
      way, regardless of method or frame? (ligand-agnostic consensus hotspots)

  (b) Do CLUSTERED SETS OF LIGANDS (chemo-types) dock in a similar region across
      all methods and all frames? i.e. "does a certain TYPE of ligand always end
      up in the same region of the channel, no matter the tool or the MD snapshot?"

Why a residue-contact "region", not a 3D centroid
--------------------------------------------------
The 4 Orai receptor frames (Orai1WT-START-Fr0 + MD snapshots Fr300/Fr400/Fr499)
are DIFFERENT conformations in DIFFERENT coordinate systems — a pose centroid is
NOT comparable across frames without superposing the receptors. So a region here is
defined by the SET OF RECEPTOR RESIDUES THE LIGAND CONTACTS (any heavy atom within
--contact-cutoff Å). Residue numbering (66-288, chains A-F) is identical across the
4 frames, so a contacted-residue fingerprint is directly comparable across frames,
methods AND ligands with no alignment. Orai1 is a symmetric homohexamer, so by
default the 6 subunits are folded onto the residue number (--no-fold-symmetry to
keep chains distinct): "the ligand sits against residue 106 of some subunit" is the
natural notion of "the same region" on a symmetric pore.

Pipeline
--------
  1. Load the Orai×Benchmark pose index, keep only PB-valid poses (pass every
     canonical PoseBusters test — same definition as the rest of the pipeline).
  2. Per pose, compute the contacted-residue fingerprint (cached; --force to redo).
  3. Cluster the PB-valid poses into binding REGIONS by fingerprint similarity
     (MiniBatchKMeans on L2-normalised residue-contact vectors; k by silhouette).
  4. Cluster the LIGANDS into chemo-types from the physicochemical-descriptor / PCA
     chemical space built by PoseBusters_DataSet_Analysis.ipynb (reused via
     ligand_descriptors.py): KMeans on the standardised descriptors.
  5. Cross-tabulate region × {ligand, ligand-cluster, method, frame} and score
     cross-ligand consensus (a) and per-ligand-type method/frame consistency (b).

Outputs (in --out-dir, default posebusters_results/orai_benchmark/region_consensus)
  pose_region_fingerprints.csv  — per PB-valid pose: method, frame, ligand and its
                                   contacted residues (the contact CACHE; reused on
                                   re-runs unless inputs change or --force)
  pose_region_assignments.csv   — the same poses annotated with their region label
                                   and ligand chemo-cluster (the analysis-ready table)
  region_summary.csv            — per region: n_poses / n_ligands / ligand coverage /
                                   #methods / #frames / signature residues / per-
                                   method & per-frame pose share / consensus flags
  ligand_clusters.csv           — per ligand: chemo-cluster + descriptors + dominant
                                   region + region purity + #methods/#frames that
                                   agree on the dominant region
  ligandcluster_region_matrix.csv — ligand-cluster × region occupancy fractions
  ligandcluster_consensus.csv   — per ligand-cluster: dominant region, purity, and
                                   how consistently its ligands dock there across
                                   methods/frames (the headline answer to (b))
  01_region_overview.png, 02_ligand_chemo_clusters.png,
  03_ligandcluster_region_consensus.png, 04_ligand_consistency.png
  05_snapshot3d_<frame>.png     — ONE per Orai snapshot: a 3D scatter of every
                                   PB-valid pose at its ligand centroid (colour =
                                   binding region / cluster, shape = docking tool)
                                   over a faint receptor envelope, plus a per-cluster
                                   tool-composition bar (how many poses in each
                                   cluster come from which tool). --no-3d to skip.
  pose_centroids.csv            — per-pose heavy-atom centroid cache for the 3D
                                   views (incremental, keyed by pose_file)
  poses_<frame>.pdb             — ONE per Orai snapshot: every PB-valid pose's heavy
                                   atoms combined into a single PDB in that frame's
                                   coordinate system. chain = docking tool, resName
                                   Rk = binding region/cluster, B-factor = region id,
                                   resSeq = pose index. Load WITH the receptor frame.
  view_<frame>.pml / .cxc       — PyMOL / ChimeraX scripts that load the receptor
                                   frame + poses_<frame>.pdb and apply the per-tool
                                   colour signatures (with a commented colour-by-
                                   cluster alternative). --no-pdb to skip all three.
  06a_jku_where.png             — (a) where the JKU ligands (GSK-7975A, Synta-66,
                                   2-APB-NH2) dock across the Orai frames: 3D centroids
                                   (colour = ligand, shape = tool) + ligand×region
                                   occupancy heatmap. JKU poses are projected onto the
                                   benchmark regions (same frames ⇒ same coords).
  06b_similar_where.png         — (b) where the chemically-similar benchmark ligands
                                   dock (ECFP4 Tanimoto ≥ --similar-threshold to a JKU
                                   ligand; top --top-n per tool by native score).
  06c_jku_vs_similar.png        — (c) JKU vs similar benchmark: per-tool region-occupancy
                                   bars + a region×tool commonality map (shared / JKU-only
                                   / similar-only). --no-jku to skip (a)/(b)/(c).
  chem_similarity.csv           — per benchmark ligand: max ECFP4 Tanimoto + nearest JKU
  jku_pose_region_assignments.csv, topn_{jku,similar}_per_tool.csv
  summary.json

Usage
-----
  python Scripts/Analysis/orai_ligand_region_consensus.py
  python Scripts/Analysis/orai_ligand_region_consensus.py --n-regions 8 --n-ligand-clusters 5
  python Scripts/Analysis/orai_ligand_region_consensus.py --limit-ligands 12   # smoke test
  python Scripts/Analysis/orai_ligand_region_consensus.py --force              # recompute contacts

Run in the `vina` conda env (RDKit + the posebusters_pose_comparison helpers).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import warnings
from pathlib import Path

import matplotlib
if os.environ.get("MPLBACKEND") is None:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from rdkit import Chem  # noqa: E402

# Reuse the single source of truth for pose loading, PDB parsing and the PB-valid
# definition (importing this in the parent also fixes the prolif/posebusters import
# order; forked workers inherit it and never re-import).
from posebusters_pose_comparison import (  # noqa: E402
    _build_pose_index, load_first_mol, load_protein_heavy_atoms, PB_CRITICAL_CHECKS,
)
from pocket_comparison_report import _label_panels  # noqa: E402
from ligand_descriptors import load_descriptors, add_pca, PCA_FEATURES, nice  # noqa: E402

DEFAULT_PB_CSV = Path("posebusters_results/orai_benchmark/dock/"
                      "posebusters_filtered_results.csv")
DEFAULT_OUT = Path("posebusters_results/orai_benchmark/region_consensus")
DEFAULT_FEATURES_CSV = "PoseBusters_Benchmark_Analysis/ligand_protein_features.csv"
DEFAULT_JKU_CSV = Path("posebusters_results/orai_jku/dock/"
                       "posebusters_filtered_results.csv")
DEFAULT_JKU_TEMPLATE_DIR = Path("posebusters_results/_orai_jku_staging/ligand_templates")
DEFAULT_LIGAND_SDF_DIR = Path("Data/Ligands/PoseBuster_Benchmark_Set")

CONTACT_CUTOFF = 4.0          # Å: a residue is "contacted" if any heavy atom is within
N_TOTAL_LIGANDS = 308          # benchmark ligands docked into Orai

METHOD_LABEL = {"autodock": "AutoDock Vina", "diffdock": "DiffDock",
                "equibind_guided": "EquiBind", "equibind": "EquiBind"}
METHOD_COLOR = {"autodock": "#1f77b4", "diffdock": "#ff7f0e",
                "equibind_guided": "#2ca02c", "equibind": "#2ca02c"}
# Marker shape per tool for the 3D snapshot views (colour there encodes region).
METHOD_MARKER = {"autodock": "o", "diffdock": "^",
                 "equibind_guided": "s", "equibind": "s"}
_MARKER_CYCLE = ["o", "^", "s", "D", "v", "P", "X", "*"]
# PDB chain ID per tool — lets PyMOL/ChimeraX colour each tool by `chain` selection.
TOOL_CHAIN = {"autodock": "A", "diffdock": "D",
              "equibind_guided": "E", "equibind": "E"}
_CHAIN_CYCLE = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")


def _mlabel(m: str) -> str:
    return METHOD_LABEL.get(str(m), str(m))


def _marker(m: str, fallback_idx: int = 0) -> str:
    return METHOD_MARKER.get(str(m), _MARKER_CYCLE[fallback_idx % len(_MARKER_CYCLE)])


def _safe(name: str) -> str:
    """Filesystem-safe token for a frame name used in output filenames."""
    return "".join(c if (c.isalnum() or c in "-._") else "_" for c in str(name))


# ───────────────────────────────────────────────────────────────────
# Contacted-residue fingerprints (the frame-invariant "region" descriptor)
# ───────────────────────────────────────────────────────────────────

def _contact_residues(mol: Chem.Mol, prot_xyz: np.ndarray,
                      prot_resid: list, cutoff: float) -> set:
    """Set of (chain, resseq) receptor residues with any heavy atom within *cutoff*
    Å of any ligand heavy atom. Vectorised, bounding-box pre-filtered."""
    mol = Chem.RemoveHs(mol)
    conf = mol.GetConformer()
    lig = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())],
                   dtype=np.float32)
    if len(lig) == 0 or len(prot_xyz) == 0:
        return set()
    lo = lig.min(0) - (cutoff + 1.0)
    hi = lig.max(0) + (cutoff + 1.0)
    mask = (prot_xyz >= lo).all(1) & (prot_xyz <= hi).all(1)
    p = prot_xyz[mask]
    pres = [r for r, m in zip(prot_resid, mask) if m]
    if len(p) == 0:
        return set()
    diff = lig[:, None, :] - p[None, :, :]
    d2 = (diff * diff).sum(-1)
    close = (d2 < cutoff * cutoff).any(0)
    return {pres[j] for j in np.where(close)[0]}


# Per-frame receptor arrays, set in each worker via the Pool initializer.
_W: dict = {}


def _init_worker(prot_xyz, prot_resid, cutoff):
    global _W
    _W = {"xyz": prot_xyz, "resid": prot_resid, "cutoff": cutoff}


def _fp_worker(pose_file: str):
    try:
        mol = load_first_mol(Path(pose_file))
    except Exception:
        mol = None
    if mol is None:
        return pose_file, None
    res = _contact_residues(mol, _W["xyz"], _W["resid"], _W["cutoff"])
    # store as "CHAIN:RESSEQ" tokens
    return pose_file, sorted(f"{c}:{int(r)}" for (c, r) in res)


def _heavy_atoms(mol):
    """(coords, element-symbols) over heavy atoms only. Filters H/D by atomic
    number — robust to mols loaded with sanitize=False, where Chem.RemoveHs can
    leave stray hydrogens behind. Returns None if no heavy atoms / no conformer."""
    try:
        mol = Chem.RemoveHs(mol)
        conf = mol.GetConformer()
    except Exception:
        return None
    coords, elems = [], []
    for i in range(mol.GetNumAtoms()):
        a = mol.GetAtomWithIdx(i)
        if a.GetAtomicNum() <= 1:
            continue
        p = conf.GetAtomPosition(i)
        coords.append((float(p.x), float(p.y), float(p.z)))
        elems.append(a.GetSymbol())
    return (coords, elems) if coords else None


def _centroid_worker(pose_file: str):
    """Heavy-atom centroid of a pose (no shared receptor state needed)."""
    try:
        mol = load_first_mol(Path(pose_file))
    except Exception:
        mol = None
    if mol is None:
        return pose_file, None
    h = _heavy_atoms(mol)
    if h is None:
        return pose_file, None
    c = np.asarray(h[0], dtype=np.float64).mean(0)
    return pose_file, (float(c[0]), float(c[1]), float(c[2]))


def _coords_worker(pose_file: str):
    """Heavy-atom (coords, element-symbols) of a pose for the combined PDB export."""
    try:
        mol = load_first_mol(Path(pose_file))
    except Exception:
        mol = None
    if mol is None:
        return pose_file, None
    return pose_file, _heavy_atoms(mol)


def _receptor_path_for_frame(idx_rows: pd.DataFrame, frame: str) -> Path | None:
    """Receptor PDB used to dock into a given frame (from protein_file_used)."""
    if "protein_file_used" in idx_rows.columns:
        v = idx_rows["protein_file_used"].dropna()
        if len(v):
            p = Path(str(v.iloc[0]))
            if p.exists():
                return p
    # fallback: common staging location
    for cand in (Path("posebusters_results/_orai_benchmark_staging/receptors") / f"{frame}.pdb",
                 Path("Data/Receptors") / f"{frame}.pdb"):
        if cand.exists():
            return cand
    return None


def _signature(args, n_poses: int, src_csv=None, limit: int | None = None) -> dict:
    pb = Path(src_csv) if src_csv is not None else Path(args.pb_csv)
    try:
        st = pb.stat()
        pb_sig = {"path": str(pb), "size": st.st_size, "mtime_ns": st.st_mtime_ns}
    except OSError:
        pb_sig = {"path": str(pb), "size": None, "mtime_ns": None}
    return {"schema": 1, "pb_csv": pb_sig, "contact_cutoff": float(args.contact_cutoff),
            "limit_ligands": int(args.limit_ligands if limit is None else limit),
            "n_pb_valid_poses": int(n_poses)}


def compute_fingerprints(args, raw_idx: pd.DataFrame, tag: str = "",
                         src_csv=None, apply_limit: bool = True) -> pd.DataFrame:
    """Per PB-valid pose: contacted-residue fingerprint. Cached + reused unless
    --force / inputs changed. Returns a frame with columns
    method, frame, ligand, pose_file, pose_name, residues (str), n_contacts.

    *tag* gives the cache its own filename so a second dataset (e.g. the JKU
    ligands, tag='jku') can be fingerprinted in the same out-dir without clobbering
    the benchmark cache. *apply_limit* honours --limit-ligands (benchmark smoke
    tests) but is turned off for the small JKU set."""
    valid = raw_idx[raw_idx["pb_valid"].astype(bool)].copy()
    valid = valid.rename(columns={"docking_method": "method", "protein": "frame"})
    if apply_limit and args.limit_ligands:
        keep = sorted(valid["ligand"].unique())[:args.limit_ligands]
        valid = valid[valid["ligand"].isin(keep)].copy()

    suffix = f".{tag}" if tag else ""
    cache = args.out_dir / f"pose_region_fingerprints{suffix}.csv"
    manifest = args.out_dir / f"pose_region_fingerprints{suffix}.manifest.json"
    sig = _signature(args, len(valid), src_csv=src_csv,
                     limit=(args.limit_ligands if apply_limit else 0))
    if not args.force and cache.exists() and manifest.exists():
        try:
            man = json.loads(manifest.read_text())
            if man.get("signature") == sig:
                df = pd.read_csv(cache, low_memory=False)
                if {"pose_file", "residues"}.issubset(df.columns):
                    print(f"Reusing cached fingerprints: {cache} ({len(df):,} poses). "
                          f"--force to recompute.")
                    fp_cols = ["method", "frame", "ligand", "pose_file",
                               "pose_name", "residues", "n_contacts"]
                    df = df[[c for c in fp_cols if c in df.columns]].copy()
                    df["residues"] = df["residues"].fillna("")
                    return df
        except Exception:
            pass

    print(f"Computing contact fingerprints for {len(valid):,} PB-valid poses "
          f"(cutoff {args.contact_cutoff} Å, {args.workers} workers) …")
    import multiprocessing as mp
    rows: list[dict] = []
    meta = valid.set_index("pose_file")
    for frame, g in valid.groupby("frame"):
        rec = _receptor_path_for_frame(g, frame)
        if rec is None:
            print(f"  ! no receptor PDB for frame {frame} — skipping its "
                  f"{len(g):,} poses"); continue
        prot_xyz, _elem, prot_resid = load_protein_heavy_atoms(rec)
        files = [f for f in g["pose_file"].tolist() if Path(f).exists()]
        print(f"  frame {frame}: {len(files):,} poses vs {len(prot_xyz):,} "
              f"receptor atoms ({rec.name})")
        if args.workers > 1:
            with mp.Pool(args.workers, initializer=_init_worker,
                         initargs=(prot_xyz, prot_resid, args.contact_cutoff)) as pool:
                results = pool.map(_fp_worker, files, chunksize=16)
        else:
            _init_worker(prot_xyz, prot_resid, args.contact_cutoff)
            results = [_fp_worker(f) for f in files]
        for pose_file, res in results:
            if res is None:
                continue
            r = meta.loc[pose_file]
            rows.append({
                "method": r["method"], "frame": frame, "ligand": r["ligand"],
                "pose_file": pose_file, "pose_name": r.get("pose_name", ""),
                "residues": " ".join(res), "n_contacts": len(res),
            })
    df = pd.DataFrame(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(cache, index=False)
    manifest.write_text(json.dumps({"signature": sig, "n_poses": len(df)}, indent=2))
    print(f"  wrote fingerprints → {cache} ({len(df):,} poses)")
    return df


def compute_pose_centroids(args, fp: pd.DataFrame) -> pd.DataFrame:
    """Add cx/cy/cz (ligand heavy-atom centroid) to *fp* for the 3D snapshot views.

    Cached incrementally by pose_file in pose_centroids.csv (centroids are
    deterministic per pose, so only poses not yet in the cache are recomputed —
    the expensive fingerprint cache is never touched). Centroids are only
    comparable WITHIN one Orai snapshot, since every frame is docked in its own
    coordinate system; the 3D figures are therefore drawn one frame at a time.
    """
    cache = args.out_dir / "pose_centroids.csv"
    have: dict[str, tuple] = {}
    if cache.exists():
        try:
            prev = pd.read_csv(cache)
            for r in prev.itertuples(index=False):
                have[str(r.pose_file)] = (float(r.cx), float(r.cy), float(r.cz))
        except Exception:
            have = {}
    todo = [f for f in fp["pose_file"].unique()
            if str(f) not in have and Path(str(f)).exists()]
    if todo:
        print(f"Computing heavy-atom centroids for {len(todo):,} poses "
              f"({args.workers} workers) …")
        import multiprocessing as mp
        if args.workers > 1:
            with mp.Pool(args.workers) as pool:
                results = pool.map(_centroid_worker, todo, chunksize=32)
        else:
            results = [_centroid_worker(f) for f in todo]
        for pose_file, c in results:
            if c is not None:
                have[str(pose_file)] = c
        args.out_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([{"pose_file": k, "cx": v[0], "cy": v[1], "cz": v[2]}
                      for k, v in have.items()]).to_csv(cache, index=False)
        print(f"  wrote centroids → {cache} ({len(have):,} poses)")
    out = fp.copy()
    coords = out["pose_file"].astype(str).map(have)
    out["cx"] = coords.map(lambda t: t[0] if isinstance(t, tuple) else np.nan)
    out["cy"] = coords.map(lambda t: t[1] if isinstance(t, tuple) else np.nan)
    out["cz"] = coords.map(lambda t: t[2] if isinstance(t, tuple) else np.nan)
    return out


# ───────────────────────────────────────────────────────────────────
# Region clustering (cluster poses by contacted-residue fingerprint)
# ───────────────────────────────────────────────────────────────────

def _fold(token: str, fold: bool) -> str:
    """'A:106' -> '106' (folded across symmetric subunits) or 'A:106' (kept)."""
    return token.split(":", 1)[1] if fold else token


def _fingerprint_matrix(df: pd.DataFrame, fold: bool):
    """(binary pose×feature matrix, feature list, per-pose token-sets)."""
    token_sets = [set(_fold(t, fold) for t in str(s).split() if t)
                  for s in df["residues"]]
    feats = sorted({t for s in token_sets for t in s},
                   key=lambda x: (len(x), x))
    fidx = {f: i for i, f in enumerate(feats)}
    M = np.zeros((len(token_sets), len(feats)), dtype=np.float32)
    for i, s in enumerate(token_sets):
        for t in s:
            M[i, fidx[t]] = 1.0
    return M, feats, token_sets


def cluster_regions(df: pd.DataFrame, fold: bool, n_regions: int | None,
                    seed: int = 0):
    """Cluster PB-valid poses into binding regions by contacted-residue fingerprint.

    Returns (labels aligned to df rows with -1 for empty-fingerprint poses,
    feature list, chosen k, silhouette, fitted KMeans, size-order remap). Uses
    MiniBatchKMeans on L2-normalised binary contact vectors (≈ spherical k-means /
    cosine), k by silhouette when not fixed via --n-regions. The model + remap are
    returned so a second pose set (JKU) can be PROJECTED onto these same regions
    via project_regions().
    """
    from sklearn.cluster import MiniBatchKMeans
    from sklearn.metrics import silhouette_score
    from sklearn.preprocessing import normalize

    M, feats, token_sets = _fingerprint_matrix(df, fold)
    nonempty = np.array([len(s) > 0 for s in token_sets])
    labels = np.full(len(df), -1, dtype=int)
    if nonempty.sum() < 4 or len(feats) < 2:
        labels[nonempty] = 0
        return labels, feats, 1, float("nan"), None, {}

    Xn = normalize(M[nonempty], norm="l2")
    rng = np.random.RandomState(seed)
    sub = (rng.choice(Xn.shape[0], 3000, replace=False)
           if Xn.shape[0] > 3000 else np.arange(Xn.shape[0]))

    if n_regions:
        ks = [n_regions]
    else:
        ks = list(range(4, min(15, Xn.shape[0] - 1)))
    best = None
    for k in ks:
        km = MiniBatchKMeans(n_clusters=k, random_state=seed, n_init=3,
                             batch_size=1024).fit(Xn)
        lab = km.labels_
        try:
            sil = silhouette_score(Xn[sub], lab[sub], metric="cosine")
        except Exception:
            sil = -1.0
        if best is None or sil > best[2]:
            best = (k, lab, sil, km)
    k, lab, sil, km = best

    # Relabel regions by descending size so region 0 is the most populated.
    order = pd.Series(lab).value_counts().index.tolist()
    remap = {old: new for new, old in enumerate(order)}
    labels[nonempty] = [remap[x] for x in lab]
    return labels, feats, k, float(sil), km, remap


def project_regions(df: pd.DataFrame, feats: list, fold: bool, km, remap: dict):
    """Assign poses in *df* to the benchmark regions defined by a fitted KMeans,
    using the SAME feature list (benchmark residue tokens). Residues a pose
    contacts that are not in *feats* are ignored; poses with no overlapping
    contact get region -1. Returns an int label array aligned to df rows.
    """
    from sklearn.preprocessing import normalize
    labels = np.full(len(df), -1, dtype=int)
    if km is None or not feats:
        token_sets = [set(_fold(t, fold) for t in str(s).split() if t)
                      for s in df["residues"]]
        labels[[i for i, s in enumerate(token_sets) if s]] = 0
        return labels
    fidx = {f: i for i, f in enumerate(feats)}
    token_sets = [set(_fold(t, fold) for t in str(s).split() if t)
                  for s in df["residues"]]
    M = np.zeros((len(token_sets), len(feats)), dtype=np.float32)
    for i, s in enumerate(token_sets):
        for t in s:
            j = fidx.get(t)
            if j is not None:
                M[i, j] = 1.0
    nonempty = M.sum(1) > 0
    if nonempty.any():
        pred = km.predict(normalize(M[nonempty], norm="l2"))
        labels[nonempty] = [remap.get(int(x), -1) for x in pred]
    return labels


def region_signature(df_region: pd.DataFrame, fold: bool, top: int = 8) -> str:
    """Top contacted residues of a region, by fraction of its poses contacting them."""
    toks = [set(_fold(t, fold) for t in str(s).split() if t)
            for s in df_region["residues"]]
    if not toks:
        return ""
    from collections import Counter
    c = Counter(t for s in toks for t in s)
    n = len(toks)
    ranked = sorted(c.items(), key=lambda kv: -kv[1])[:top]
    return ", ".join(f"{t}({100 * v / n:.0f}%)" for t, v in ranked)


# ───────────────────────────────────────────────────────────────────
# Ligand chemo-clustering (from the descriptor / PCA chemical space)
# ───────────────────────────────────────────────────────────────────

def cluster_ligands(ligands: list[str], features_csv: str,
                    n_clusters: int | None, seed: int = 0):
    """Cluster ligands into chemo-types from the physicochemical descriptors built
    by PoseBusters_DataSet_Analysis.ipynb (reused via ligand_descriptors).

    Returns (desc_df indexed by ligand with 'lig_cluster' + PCA scores,
    cluster_profile DataFrame, feats used, k). desc_df is None if the table is
    missing.
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    from sklearn.preprocessing import StandardScaler

    desc = load_descriptors(features_csv)
    if desc is None:
        return None, None, [], 0
    desc = desc[desc.index.isin(ligands)].copy()
    feats = [f for f in PCA_FEATURES if f in desc.columns]
    X = desc[feats].apply(pd.to_numeric, errors="coerce")
    ok = X.dropna()
    if len(ok) < 6:
        return desc, None, feats, 0
    Xz = StandardScaler().fit_transform(ok)

    if n_clusters:
        ks = [n_clusters]
    else:
        ks = list(range(3, min(9, len(ok) - 1)))
    best = None
    for k in ks:
        lab = KMeans(n_clusters=k, random_state=seed, n_init=10).fit_predict(Xz)
        try:
            sil = silhouette_score(Xz, lab)
        except Exception:
            sil = -1.0
        if best is None or sil > best[2]:
            best = (k, lab, sil)
    k, lab, _sil = best

    # Order clusters by median molecular weight so labels are stable/interpretable.
    tmp = pd.Series(lab, index=ok.index)
    mw = ok["mw"] if "mw" in ok else ok[feats[0]]
    order = tmp.groupby(tmp).apply(lambda s: mw.loc[s.index].median()).sort_values().index
    remap = {old: f"LC{new + 1}" for new, old in enumerate(order)}
    desc["lig_cluster"] = tmp.map(remap)
    desc, _load, _evr, _npc = add_pca(desc, features=feats)

    # Profile: median descriptors per cluster + an auto descriptor tag.
    prof = (desc.dropna(subset=["lig_cluster"])
            .groupby("lig_cluster")[feats].median())
    gmed = ok[feats].median()
    tags = {}
    for lc, row in prof.iterrows():
        hi = [nice(f).split(" (")[0] for f in feats if row[f] > 1.25 * gmed[f]]
        tags[lc] = (", ".join(hi[:3]) if hi else "near-median") + f" (n={int((desc['lig_cluster'] == lc).sum())})"
    prof["profile"] = pd.Series(tags)
    return desc, prof, feats, k


# ───────────────────────────────────────────────────────────────────
# Cross-tabulation & consensus scoring
# ───────────────────────────────────────────────────────────────────

def build_region_summary(fp: pd.DataFrame, fold: bool, n_lig_with_valid: int,
                         all_methods: list, all_frames: list) -> pd.DataFrame:
    rows = []
    for region, g in fp.groupby("region"):
        if region < 0:
            continue
        meths = sorted(g["method"].unique())
        frames = sorted(g["frame"].unique())
        rows.append({
            "region": f"R{region + 1}", "region_id": int(region),
            "n_poses": len(g), "n_ligands": g["ligand"].nunique(),
            "ligand_coverage_%": round(100 * g["ligand"].nunique() / n_lig_with_valid, 1),
            "n_methods": len(meths), "n_frames": len(frames),
            "methods": ", ".join(_mlabel(m) for m in meths),
            "frames": ", ".join(frames),
            "all_methods": len(meths) == len(all_methods),
            "all_frames": len(frames) == len(all_frames),
            "signature_residues": region_signature(g, fold),
        })
    out = pd.DataFrame(rows).sort_values("n_ligands", ascending=False)
    # Broad cross-ligand consensus = hit by every method & frame and by a large
    # share of distinct ligands.
    if not out.empty:
        out["broad_consensus"] = (out["all_methods"] & out["all_frames"]
                                  & (out["ligand_coverage_%"] >= 50))
    return out.reset_index(drop=True)


def per_ligand_consistency(fp: pd.DataFrame, all_methods: list,
                           all_frames: list, purity_thr: float) -> pd.DataFrame:
    rows = []
    for ligand, g in fp.groupby("ligand"):
        g = g[g["region"] >= 0]
        if g.empty:
            continue
        counts = g["region"].value_counts()
        dom = int(counts.index[0])
        purity = counts.iloc[0] / len(g)
        gd = g[g["region"] == dom]
        m_dom = sorted(gd["method"].unique())
        f_dom = sorted(gd["frame"].unique())
        m_tot = sorted(g["method"].unique())
        f_tot = sorted(g["frame"].unique())
        rows.append({
            "ligand": ligand, "n_valid_poses": len(g),
            "dominant_region": f"R{dom + 1}", "dominant_region_id": dom,
            "region_purity_%": round(100 * purity, 1),
            "n_methods_total": len(m_tot), "n_methods_in_dom": len(m_dom),
            "n_frames_total": len(f_tot), "n_frames_in_dom": len(f_dom),
            "all_methods_agree": len(m_dom) == len(all_methods),
            "all_frames_agree": len(f_dom) == len(all_frames),
            "consistent": (purity >= purity_thr
                           and len(m_dom) == len(all_methods)
                           and len(f_dom) == len(all_frames)),
        })
    return pd.DataFrame(rows).sort_values(
        ["consistent", "region_purity_%"], ascending=False).reset_index(drop=True)


def ligandcluster_region_matrix(fp: pd.DataFrame, regions: list) -> pd.DataFrame:
    sub = fp[(fp["region"] >= 0) & fp["lig_cluster"].notna()]
    if sub.empty:
        return pd.DataFrame()
    ct = pd.crosstab(sub["lig_cluster"], sub["region"])
    ct = ct.reindex(columns=[r for r in regions if r in ct.columns], fill_value=0)
    frac = ct.div(ct.sum(axis=1), axis=0)
    frac.columns = [f"R{c + 1}" for c in frac.columns]
    return frac


def ligandcluster_consensus(fp: pd.DataFrame, per_lig: pd.DataFrame,
                            prof: pd.DataFrame | None) -> pd.DataFrame:
    rows = []
    for lc, g in fp[fp["lig_cluster"].notna() & (fp["region"] >= 0)].groupby("lig_cluster"):
        counts = g["region"].value_counts()
        dom = int(counts.index[0])
        purity = counts.iloc[0] / len(g)
        ligs = per_lig[per_lig["ligand"].isin(g["ligand"].unique())]
        rows.append({
            "lig_cluster": lc,
            "profile": (prof.loc[lc, "profile"] if prof is not None
                        and lc in prof.index else ""),
            "n_ligands": g["ligand"].nunique(), "n_valid_poses": len(g),
            "dominant_region": f"R{dom + 1}", "dominant_region_id": dom,
            "pose_share_in_dom_%": round(100 * purity, 1),
            "ligands_consistent_%": round(100 * ligs["consistent"].mean(), 1) if len(ligs) else 0.0,
            "ligands_all_methods_%": round(100 * ligs["all_methods_agree"].mean(), 1) if len(ligs) else 0.0,
            "ligands_all_frames_%": round(100 * ligs["all_frames_agree"].mean(), 1) if len(ligs) else 0.0,
        })
    return pd.DataFrame(rows).sort_values("pose_share_in_dom_%",
                                          ascending=False).reset_index(drop=True)


# ───────────────────────────────────────────────────────────────────
# Plots
# ───────────────────────────────────────────────────────────────────

def plot_region_overview(fp, region_sum, all_methods, all_frames, out: Path):
    regions = list(region_sum["region_id"])
    rlabels = list(region_sum["region"])
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    ax = axes[0, 0]
    ax.bar(rlabels, region_sum["n_poses"], color="#4c72b0", edgecolor="black")
    ax.set_ylabel("Number of PB-valid poses")
    ax.set_title("Region size", fontweight="bold")
    ax.tick_params(axis="x", rotation=45)

    ax = axes[0, 1]
    ax.bar(rlabels, region_sum["ligand_coverage_%"], color="#55a868", edgecolor="black")
    ax.set_ylabel("% of ligands with a PB-valid pose here")
    ax.set_title("Cross-ligand coverage per region\n(how many distinct ligands hit it)",
                 fontweight="bold")
    ax.tick_params(axis="x", rotation=45)
    ax.axhline(50, color="grey", ls="--", alpha=0.6)

    # region × method and region × frame pose-share heatmaps
    for ax, key, order, title in (
            (axes[1, 0], "method", all_methods, "Region × method"),
            (axes[1, 1], "frame", all_frames, "Region × frame")):
        ct = pd.crosstab(fp[fp.region >= 0]["region"], fp[fp.region >= 0][key])
        ct = ct.reindex(index=regions, columns=order, fill_value=0)
        frac = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0).fillna(0).values
        im = ax.imshow(frac, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([_mlabel(o) if key == "method" else o for o in order],
                           rotation=30, ha="right")
        ax.set_yticks(range(len(regions)))
        ax.set_yticklabels(rlabels)
        ax.set_title(title + " — pose share within region", fontweight="bold")
        for i in range(frac.shape[0]):
            for j in range(frac.shape[1]):
                ax.text(j, i, f"{frac[i, j]:.2f}", ha="center", va="center",
                        color="white" if frac[i, j] < 0.6 else "black", fontsize=7)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    _label_panels(axes)
    fig.suptitle("Orai PB-valid binding regions — defined by contacted-residue "
                 "fingerprints (frame- & method-invariant)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def plot_ligand_clusters(desc, prof, feats, out: Path):
    if desc is None or "lig_cluster" not in desc.columns:
        return
    cl = sorted(desc["lig_cluster"].dropna().unique())
    cmap = plt.get_cmap("tab10")
    colors = {c: cmap(i % 10) for i, c in enumerate(cl)}
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    ax = axes[0]
    if {"PC1", "PC2"}.issubset(desc.columns):
        for c in cl:
            s = desc[desc["lig_cluster"] == c]
            ax.scatter(s["PC1"], s["PC2"], s=28, alpha=0.75, color=colors[c],
                       edgecolors="none", label=f"{c} (n={len(s)})")
        ax.set_xlabel("PC1 — ligand chemical space")
        ax.set_ylabel("PC2 — ligand chemical space")
        ax.legend(fontsize=8, title="Chemo-cluster")
    ax.set_title("Ligand chemo-clusters in descriptor PCA space", fontweight="bold")
    ax.grid(alpha=0.3)

    ax = axes[1]
    if prof is not None and feats:
        from sklearn.preprocessing import StandardScaler
        prof_feats = prof[feats]
        Z = StandardScaler().fit_transform(prof_feats.values)
        im = ax.imshow(Z, aspect="auto", cmap="coolwarm", vmin=-2, vmax=2)
        ax.set_yticks(range(len(prof_feats.index)))
        ax.set_yticklabels(prof_feats.index)
        ax.set_xticks(range(len(feats)))
        ax.set_xticklabels([nice(f).split(" (")[0] for f in feats],
                           rotation=45, ha="right", fontsize=8)
        ax.set_title("Cluster descriptor profile (z-scored medians)", fontweight="bold")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="z-score")

    _label_panels(axes)
    fig.suptitle("Ligand chemo-types — physicochemical descriptor clustering "
                 "(reused from PoseBusters_DataSet_Analysis)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def plot_ligandcluster_region(matrix, consensus, out: Path):
    if matrix is None or matrix.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.8),
                             gridspec_kw={"width_ratios": [1.25, 1]})

    ax = axes[0]
    im = ax.imshow(matrix.values, aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(matrix.shape[1])); ax.set_xticklabels(matrix.columns)
    ax.set_yticks(range(matrix.shape[0])); ax.set_yticklabels(matrix.index)
    ax.set_xlabel("Binding region"); ax.set_ylabel("Ligand chemo-cluster")
    ax.set_title("Where each ligand type docks\n(fraction of the cluster's PB-valid "
                 "poses in each region)", fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = matrix.values[i, j]
            if v > 0.02:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        color="white" if v < 0.6 else "black", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="pose fraction")

    ax = axes[1]
    if consensus is not None and not consensus.empty:
        y = np.arange(len(consensus)); h = 0.27
        ax.barh(y + h, consensus["pose_share_in_dom_%"], h, color="#b07aa1",
                edgecolor="black", label="poses in dominant region")
        ax.barh(y, consensus["ligands_all_methods_%"], h, color="#4c72b0",
                edgecolor="black", label="ligands: all methods agree")
        ax.barh(y - h, consensus["ligands_all_frames_%"], h, color="#dd8452",
                edgecolor="black", label="ligands: all frames agree")
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.lig_cluster}→{r.dominant_region}"
                            for r in consensus.itertuples()])
        ax.set_xlabel("%"); ax.set_xlim(0, 100)
        ax.legend(fontsize=8, loc="lower right")
        ax.set_title("Per chemo-cluster consistency\n(dominant region + method/frame "
                     "agreement)", fontweight="bold")
        ax.grid(axis="x", alpha=0.3)

    _label_panels(axes)
    fig.suptitle("Do ligand types dock in a consistent region across methods & frames?",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def plot_consistency(per_lig, region_sum, out: Path):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.4))

    ax = axes[0]
    ax.hist(per_lig["region_purity_%"], bins=20, color="#4c72b0", edgecolor="black")
    ax.axvline(per_lig["region_purity_%"].median(), color="red", ls="--",
               label=f"median {per_lig['region_purity_%'].median():.0f}%")
    ax.set_xlabel("Dominant-region purity\n(% of a ligand's PB-valid poses in its top region)")
    ax.set_ylabel("Number of ligands")
    ax.set_title("Dominant-region purity per ligand", fontsize=10, fontweight="bold")
    ax.legend(); ax.grid(alpha=0.3)

    ax = axes[1]
    cats = {
        "all methods\n& frames": int((per_lig["all_methods_agree"] & per_lig["all_frames_agree"]).sum()),
        "all methods\nonly": int((per_lig["all_methods_agree"] & ~per_lig["all_frames_agree"]).sum()),
        "all frames\nonly": int((~per_lig["all_methods_agree"] & per_lig["all_frames_agree"]).sum()),
        "neither": int((~per_lig["all_methods_agree"] & ~per_lig["all_frames_agree"]).sum()),
    }
    ax.bar(list(cats.keys()), list(cats.values()),
           color=["#55a868", "#4c72b0", "#dd8452", "#c44e52"], edgecolor="black")
    for i, v in enumerate(cats.values()):
        ax.text(i, v + 0.5, str(v), ha="center", fontsize=9, fontweight="bold")
    ax.set_ylabel("Number of ligands")
    ax.set_title("Dominant-region agreement", fontsize=10, fontweight="bold")

    ax = axes[2]
    rs = region_sum.head(10)
    ax.barh(range(len(rs)), rs["n_ligands"], color="#55a868", edgecolor="black")
    ax.set_yticks(range(len(rs)))
    ax.set_yticklabels([f"{r.region}: {r.signature_residues[:34]}"
                        for r in rs.itertuples()], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Number of distinct ligands")
    ax.set_title("Consensus regions (top residues)", fontsize=10, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)

    _label_panels(axes)
    fig.suptitle("Ligand-level docking consistency across methods & MD frames",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def plot_snapshot_clusters_3d(fp_xyz, frame, region_colors, out: Path,
                              receptor_path: Path | None = None,
                              max_ghost: int = 3000):
    """3D view of ONE Orai snapshot: every PB-valid pose plotted at its ligand
    heavy-atom centroid, COLOURED by binding region (cluster) and SHAPED by the
    docking tool that produced it, alongside a per-cluster tool-composition bar
    (how many poses in each cluster come from each tool).

    One figure per frame on purpose — pose coordinates are only comparable within
    a single snapshot (each MD frame is its own coordinate system), so coordinates
    are never mixed across frames. Returns the written path (or None if empty).
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    d = fp_xyz[(fp_xyz["frame"] == frame) & (fp_xyz["region"] >= 0)
               & fp_xyz["cx"].notna()].copy()
    if d.empty:
        return None
    methods = sorted(d["method"].unique())
    regions_here = sorted(int(r) for r in d["region"].unique())

    fig = plt.figure(figsize=(16.5, 7.5))
    ax = fig.add_subplot(1, 2, 1, projection="3d")

    # faint receptor heavy-atom envelope for spatial context (sub-sampled)
    if receptor_path is not None and Path(receptor_path).exists():
        try:
            pxyz, _e, _r = load_protein_heavy_atoms(Path(receptor_path))
            if len(pxyz) > max_ghost:
                sel = np.linspace(0, len(pxyz) - 1, max_ghost).astype(int)
                pxyz = pxyz[sel]
            ax.scatter(pxyz[:, 0], pxyz[:, 1], pxyz[:, 2], s=2, c="0.75",
                       alpha=0.10, marker=".", linewidths=0, depthshade=False)
        except Exception:
            pass

    # one scatter per (region, tool): colour = region, marker = tool
    for reg in regions_here:
        col = region_colors.get(reg, "#333333")
        for m in methods:
            s = d[(d["region"] == reg) & (d["method"] == m)]
            if s.empty:
                continue
            ax.scatter(s["cx"], s["cy"], s["cz"], s=46, color=col,
                       marker=_marker(m), edgecolors="black", linewidths=0.4,
                       alpha=0.92, depthshade=False)
    ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    ax.set_title("PB-valid pose centroids\n(colour = binding region · "
                 "shape = docking tool)", fontsize=10)

    # two legends: tool shapes (left) + region colours (right)
    tool_handles = [Line2D([0], [0], marker=_marker(m), color="0.3", ls="",
                           markerfacecolor="0.7", markeredgecolor="black",
                           markersize=9, label=_mlabel(m)) for m in methods]
    reg_handles = [Patch(facecolor=region_colors.get(r, "#333333"),
                         edgecolor="black", label=f"R{r + 1}")
                   for r in regions_here]
    leg1 = ax.legend(handles=tool_handles, title="Docking tool", loc="upper left",
                     fontsize=8, framealpha=0.9)
    ax.add_artist(leg1)
    ax.legend(handles=reg_handles, title="Region", loc="upper right",
              fontsize=8, framealpha=0.9, ncol=2)

    # right panel: per-cluster tool composition (stacked counts)
    ax2 = fig.add_subplot(1, 2, 2)
    ct = (pd.crosstab(d["region"], d["method"])
          .reindex(index=regions_here, columns=methods, fill_value=0))
    xpos = np.arange(len(regions_here))
    bottom = np.zeros(len(regions_here))
    for m in methods:
        vals = ct[m].to_numpy()
        ax2.bar(xpos, vals, bottom=bottom, color=METHOD_COLOR.get(m, "#888888"),
                edgecolor="black", linewidth=0.5, label=_mlabel(m))
        for xi, (v, b) in enumerate(zip(vals, bottom)):
            if v > 0:
                ax2.text(xi, b + v / 2, str(int(v)), ha="center", va="center",
                         fontsize=7, color="white")
        bottom += vals
    ax2.set_xticks(xpos)
    ax2.set_xticklabels([f"R{r + 1}" for r in regions_here])
    ax2.set_xlabel("Binding region (cluster)")
    ax2.set_ylabel("Number of PB-valid poses")
    ax2.set_title("Per-cluster tool composition\n(poses per tool in each cluster)",
                  fontsize=10)
    ax2.legend(fontsize=8, title="Docking tool")
    ax2.grid(axis="y", alpha=0.3)

    _label_panels([ax, ax2])
    fig.suptitle(f"Orai snapshot {frame}: clustered PB-valid poses by tool "
                 f"({len(d):,} poses · {len(regions_here)} regions · "
                 f"{len(methods)} tools)", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    return out


# ───────────────────────────────────────────────────────────────────
# Combined per-frame pose PDBs for PyMOL / ChimeraX
# ───────────────────────────────────────────────────────────────────

def _pdb_atom_name(elem: str, idx: int) -> str:
    """PDB atom-name field (cols 13-16, 4 wide). 1-char elements get a leading
    space per the PDB convention; longer names are clipped to fit."""
    nm = f"{elem}{idx}"
    if len(nm) > 4:
        nm = nm[:4]
    if len(elem) == 1 and len(nm) < 4:
        return f" {nm:<3s}"
    return f"{nm:<4s}"


def _hetatm_line(serial: int, atname: str, resname: str, chain: str, resseq: int,
                 x: float, y: float, z: float, occ: float, bfac: float,
                 elem: str) -> str:
    """One PDB HETATM record with strict column placement. Serial wraps at 99999
    and resSeq at 9999 (both are cosmetic here — selection is by chain / resName /
    B-factor, never by serial)."""
    s = (serial - 1) % 99999 + 1
    rs = (resseq - 1) % 9999 + 1
    return (f"HETATM{s:>5d} {atname:<4s} {resname:>3s} {chain:1s}{rs:>4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}          {elem:>2s}")


def _tool_chain(method: str, used: dict) -> str:
    """Stable single-char chain ID for a tool (known tools fixed; others cycle)."""
    c = TOOL_CHAIN.get(str(method))
    if c is None:
        taken = set(used.values())
        c = next((ch for ch in _CHAIN_CYCLE if ch not in taken), "Z")
    used[str(method)] = c
    return c


def export_frame_pose_pdbs(args, fp, region_colors, receptor_for_frame,
                           out_dir: Path):
    """Per Orai frame, write ONE combined PDB of every PB-valid pose's heavy atoms
    (chain = docking tool, resName = region/cluster Rk, B-factor = region id,
    resSeq = unique pose index) plus PyMOL (.pml) and ChimeraX (.cxc) scripts that
    load it together with the receptor frame and apply the per-tool colour
    signatures. Poses with no receptor contact (region < 0) are omitted.

    Returns the list of written PDB filenames.
    """
    import matplotlib.colors as mcolors
    written: list[str] = []
    sub = fp[(fp["region"] >= 0)].copy()
    if sub.empty:
        return written
    for fr in sorted(sub["frame"].unique()):
        d = sub[sub["frame"] == fr].sort_values(
            ["method", "region", "ligand", "pose_file"]).reset_index(drop=True)
        files = list(dict.fromkeys(d["pose_file"].tolist()))  # unique, order-stable
        files = [f for f in files if Path(str(f)).exists()]
        if not files:
            continue
        print(f"  PDB export {fr}: loading {len(files):,} pose geometries …")
        import multiprocessing as mp
        if args.workers > 1:
            with mp.Pool(args.workers) as pool:
                results = pool.map(_coords_worker, files, chunksize=32)
        else:
            results = [_coords_worker(f) for f in files]
        geom = {pf: g for pf, g in results if g is not None}

        methods = sorted(d["method"].unique())
        chain_used: dict = {}
        chain_for = {m: _tool_chain(m, chain_used) for m in methods}

        lines = [
            "REMARK 200 Combined PB-valid docking poses for Orai snapshot "
            f"{fr}",
            "REMARK 200 chain = docking tool | resName Rk = binding region/cluster "
            "| B-factor = region id | resSeq = pose index",
            "REMARK 200 tools: " + ", ".join(
                f"{chain_for[m]}={_mlabel(m)}" for m in methods),
        ]
        serial = 0
        pose_seq = 0
        n_atoms = 0
        for row in d.itertuples():
            g = geom.get(row.pose_file)
            if g is None:
                continue
            coords, elems = g
            pose_seq += 1
            chain = chain_for[row.method]
            reg = int(row.region)
            resname = f"R{reg + 1}"
            for ai, ((x, y, z), el) in enumerate(zip(coords, elems), start=1):
                serial += 1
                n_atoms += 1
                lines.append(_hetatm_line(
                    serial, _pdb_atom_name(el, ai), resname, chain, pose_seq,
                    x, y, z, 1.0, float(reg), el))
        lines.append("END")
        pdb_path = out_dir / f"poses_{_safe(fr)}.pdb"
        pdb_path.write_text("\n".join(lines) + "\n")
        written.append(pdb_path.name)
        print(f"    wrote {pdb_path.name}: {pose_seq:,} poses · {n_atoms:,} atoms")

        # Companion viewer scripts (absolute paths so they run from anywhere).
        rec = receptor_for_frame.get(fr)
        rec_abs = str(Path(rec).resolve()) if rec and Path(rec).exists() else None
        pdb_abs = str(pdb_path.resolve())
        regs_here = sorted(int(r) for r in d["region"].unique())
        _write_pymol_script(out_dir / f"view_{_safe(fr)}.pml", fr, rec_abs,
                            pdb_abs, methods, chain_for, regs_here, region_colors,
                            mcolors)
        _write_chimerax_script(out_dir / f"view_{_safe(fr)}.cxc", fr, rec_abs,
                               pdb_abs, methods, chain_for, regs_here,
                               region_colors, mcolors)
    return written


def _write_pymol_script(path, frame, rec_abs, pdb_abs, methods, chain_for,
                        regs_here, region_colors, mcolors):
    L = [f"# PyMOL: PB-valid poses on Orai snapshot {frame}",
         "# Colour signature = docking tool (chain). Run:  pymol "
         f"{path.name}", "bg_color white", "set ray_shadows, 0"]
    if rec_abs:
        L += [f'load {rec_abs}, receptor',
              "hide everything, receptor", "show cartoon, receptor",
              "color grey80, receptor", "set cartoon_transparency, 0.4, receptor"]
    L += [f'load {pdb_abs}, poses',
          "hide everything, poses", "show spheres, poses",
          "set sphere_scale, 0.30, poses",
          "# (or 'show sticks, poses' for ligand skeletons)"]
    L.append("")
    L.append("# ---- colour by DOCKING TOOL (default) ----")
    for m in methods:
        hx = METHOD_COLOR.get(m, "#888888").lstrip("#")
        L.append(f'color 0x{hx}, poses and chain {chain_for[m]}    '
                 f'# {_mlabel(m)}')
    L.append("")
    L.append("# ---- colour by REGION / CLUSTER instead (uncomment) ----")
    for r in regs_here:
        hx = mcolors.to_hex(region_colors.get(r, "#333333")).lstrip("#")
        L.append(f'# color 0x{hx}, poses and resn R{r + 1}')
    L.append("# or a continuous ramp over the region id stored in B-factor:")
    L.append("# spectrum b, rainbow, poses")
    L.append("orient")
    path.write_text("\n".join(L) + "\n")


def _write_chimerax_script(path, frame, rec_abs, pdb_abs, methods, chain_for,
                           regs_here, region_colors, mcolors):
    L = [f"# ChimeraX: PB-valid poses on Orai snapshot {frame}",
         f"# Run:  chimerax {path.name}", "set bgColor white"]
    rec_model, pose_model = None, "#1"
    if rec_abs:
        L.append(f'open {rec_abs}')
        rec_model, pose_model = "#1", "#2"
        L += [f'cartoon {rec_model}', f'color {rec_model} gray',
              f'transparency {rec_model} 60 target c']
    L += [f'open {pdb_abs}', f'style {pose_model} sphere',
          f'size {pose_model} atomRadius 0.45']
    L.append("")
    L.append("# ---- colour by DOCKING TOOL (default) ----")
    for m in methods:
        hx = METHOD_COLOR.get(m, "#888888")
        L.append(f'color {pose_model}/{chain_for[m]} {hx}    # {_mlabel(m)}')
    L.append("")
    L.append("# ---- colour by REGION / CLUSTER instead (uncomment) ----")
    for r in regs_here:
        hx = mcolors.to_hex(region_colors.get(r, "#333333"))
        L.append(f'# color {pose_model}:R{r + 1} {hx}')
    L.append('# or by the region id in the B-factor:  '
             f'color byattribute bfactor {pose_model} palette rainbow')
    L.append("view")
    path.write_text("\n".join(L) + "\n")


# ───────────────────────────────────────────────────────────────────
# JKU ligands × chemically-similar benchmark ligands cross-analysis
# ───────────────────────────────────────────────────────────────────

JKU_LABEL = {"2abp-nh2-OPT": "2-APB-NH2", "Synta-66-OPT-Singlet": "Synta-66",
             "gsk7975a-deprot-OPT": "GSK-7975A"}


def _jlabel(lig: str) -> str:
    return JKU_LABEL.get(str(lig), str(lig))


def _morgan_fp(sdf_path: Path, radius: int = 2, nbits: int = 2048):
    """ECFP4 (Morgan r2, 2048-bit) fingerprint from the first mol of an SDF."""
    from rdkit.Chem import AllChem
    try:
        for m in Chem.SDMolSupplier(str(sdf_path)):
            if m is not None:
                return AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=nbits)
    except Exception:
        return None
    return None


def compute_chem_similarity(jku_ligs, bench_ligs, jku_dir: Path, bench_dir: Path,
                            threshold: float) -> pd.DataFrame:
    """For each benchmark ligand, its max ECFP4 Tanimoto to any JKU ligand and the
    nearest JKU compound. is_similar = max Tanimoto ≥ threshold. Fingerprints are
    built from the clean start-conf / template SDFs (not the PDBQT round-trips)."""
    from rdkit import DataStructs
    jfp = {}
    for l in jku_ligs:
        fp = _morgan_fp(Path(jku_dir) / f"{l}.sdf")
        if fp is not None:
            jfp[l] = fp
    if not jfp:
        print("  ! no JKU template fingerprints could be built — similarity skipped.")
    rows = []
    for l in bench_ligs:
        fp = _morgan_fp(Path(bench_dir) / f"{l}_ligand_start_conf.sdf")
        if fp is None or not jfp:
            rows.append({"ligand": l, "max_tanimoto": np.nan,
                         "nearest_jku": None, "is_similar": False})
            continue
        sims = {jl: DataStructs.TanimotoSimilarity(fp, jf) for jl, jf in jfp.items()}
        nj = max(sims, key=sims.get)
        rows.append({"ligand": l, "max_tanimoto": round(float(sims[nj]), 3),
                     "nearest_jku": nj, "is_similar": bool(sims[nj] >= threshold)})
    return pd.DataFrame(rows).sort_values("max_tanimoto", ascending=False,
                                          na_position="last").reset_index(drop=True)


def _vina_affinities_from_pdbqt(pdbqt_path: Path) -> dict:
    """{model_number: Vina affinity (kcal/mol)} from a multi-model Vina output
    PDBQT (more-negative = better). This is the authoritative AutoDock score."""
    out, cur = {}, None
    try:
        for ln in Path(pdbqt_path).read_text(errors="ignore").splitlines():
            if ln.startswith("MODEL"):
                p = ln.split()
                cur = int(p[1]) if len(p) > 1 and p[1].isdigit() else None
            elif ln.startswith("REMARK VINA RESULT:") and cur is not None:
                m = re.search(r"(-?\d+\.\d+)", ln)
                if m:
                    out[cur] = float(m.group(1))
    except (OSError, ValueError):
        pass
    return out


def _vina_affinity_from_sdf(sdf_path: str):
    """Vina affinity from a converted SDF's <REMARK> tag, if the (meeko) conversion
    preserved it (the meeko-export path drops it — hence this is only a fallback)."""
    try:
        txt = Path(sdf_path).read_text(errors="ignore")
    except OSError:
        return None
    m = re.search(r"VINA RESULT:\s*(-?\d+\.\d+)", txt)
    return float(m.group(1)) if m else None


def _autodock_affinity_lookup(raw_df: pd.DataFrame) -> dict:
    """{pose_file: Vina affinity} for AutoDock poses. Prefers the 'autodock_affinity'
    CSV column (populated by run_posebusters from now on); otherwise recovers it from
    the source multi-model PDBQT (located via protein_file_used) and finally the
    converted SDF — so it also works on CSVs written before that column existed."""
    import re as _re
    if "docking_method" not in raw_df.columns:
        return {}
    ad = raw_df[raw_df["docking_method"] == "autodock"]
    if ad.empty:
        return {}
    if "autodock_affinity" in ad.columns and ad["autodock_affinity"].notna().any():
        return {str(pf): float(a) for pf, a in zip(ad["pose_file"], ad["autodock_affinity"])
                if pd.notna(a)}
    out, pdbqt_cache = {}, {}
    for r in ad.itertuples(index=False):
        pf = str(r.pose_file)
        mm = _re.search(r"_model(\d+)\.sdf$", pf)
        model = int(mm.group(1)) if mm else None
        aff = None
        pfu = getattr(r, "protein_file_used", None)
        if model is not None and isinstance(pfu, str) and pfu not in ("", "nan", "none"):
            base = _re.sub(r"_model\d+\.sdf$", ".pdbqt", Path(pf).name)
            src = Path(pfu).parent.parent / "autodock" / base
            if src not in pdbqt_cache:
                pdbqt_cache[src] = _vina_affinities_from_pdbqt(src) if src.exists() else {}
            aff = pdbqt_cache[src].get(model)
        if aff is None:
            aff = _vina_affinity_from_sdf(pf)
        if aff is not None:
            out[pf] = aff
    return out


def attach_native_scores(df: pd.DataFrame, score_lookup: dict,
                         autodock_aff: dict) -> np.ndarray:
    """native_score per pose, normalised so LOWER = better and comparable WITHIN a
    tool: AutoDock = Vina affinity, EquiBind = smina_affinity (or gnina_affinity for
    gnina-refined poses, see _score_lookup), DiffDock = -confidence."""
    out = np.full(len(df), np.nan)
    for i, (m, f) in enumerate(zip(df["method"].to_numpy(), df["pose_file"].to_numpy())):
        if m == "autodock":
            out[i] = autodock_aff.get(str(f), np.nan)
        elif m == "diffdock":
            conf = score_lookup.get(f, (np.nan, np.nan))[0]
            out[i] = -conf if conf == conf else np.nan  # higher conf → lower (better)
        else:  # equibind*
            out[i] = score_lookup.get(f, (np.nan, np.nan))[1]
    return out


def _score_lookup(raw_df: pd.DataFrame) -> dict:
    """pose_file → (diffdock_confidence, equibind_score) from a raw results CSV.

    The EquiBind score is ``smina_affinity``, falling back to ``gnina_affinity`` for
    gnina-refined poses that carry no smina re-score. This lets gnina-only EquiBind
    datasets — e.g. Orai×JKU, whose EquiBind was refined with gnina, leaving
    ``smina_affinity`` empty — still rank EquiBind by its minimisation energy instead
    of dropping to a NaN native_score (ranked last / by pose count). Both columns use
    the same LOWER = better (more-negative kcal/mol) convention, so coalescing is
    scale-consistent within EquiBind."""
    sub = raw_df.drop_duplicates("pose_file").set_index("pose_file")

    def _num(col):
        return (pd.to_numeric(sub[col], errors="coerce") if col in sub.columns
                else pd.Series(np.nan, index=sub.index))

    conf = _num("diffdock_confidence")
    smina = _num("smina_affinity")
    eq = smina.where(smina.notna(), _num("gnina_affinity"))   # smina, else gnina
    return {pf: (float(conf.get(pf, np.nan)), float(eq.get(pf, np.nan)))
            for pf in sub.index}


def top_n_ligands_per_tool(df_scored: pd.DataFrame, n: int):
    """Per tool, the n best ligands ranked by best (lowest) native_score, breaking
    ties — and ranking tools with no score at all (e.g. AutoDock on a CSV written
    before the autodock_affinity column existed) — by PB-valid pose count, so the
    selection is always deterministic and meaningful. Returns (dict method→[ligands],
    tidy ranking DataFrame)."""
    sel, rows = {}, []
    for m, g in df_scored.groupby("method"):
        order = pd.DataFrame({"score": g.groupby("ligand")["native_score"].min(),
                              "n_poses": g.groupby("ligand").size()})
        order = order.sort_values(["score", "n_poses"], ascending=[True, False],
                                  na_position="last")
        sel[m] = list(order.index[:n]) if n else list(order.index)
        for rank, (lig, r) in enumerate(order.iterrows(), 1):
            rows.append({"method": m, "ligand": lig, "best_native_score": r["score"],
                         "n_poses": int(r["n_poses"]), "rank": rank,
                         "selected": lig in sel[m]})
    return sel, pd.DataFrame(rows)


def _filter_top_n(df_scored: pd.DataFrame, sel: dict) -> pd.DataFrame:
    parts = [g[g["ligand"].isin(sel.get(m, []))] for m, g in df_scored.groupby("method")]
    return pd.concat(parts) if parts else df_scored.iloc[:0]


def _region_fraction(df: pd.DataFrame, regions: list) -> pd.Series:
    d = df[df["region"] >= 0]
    if d.empty:
        return pd.Series(0.0, index=regions)
    vc = d["region"].value_counts(normalize=True)
    return pd.Series([float(vc.get(r, 0.0)) for r in regions], index=regions)


def plot_pose_set_where(d, region_ids, receptor_for_frame, title, out,
                        lig_label=str, lig_sublabel=None, max_ghost=2500):
    """Generic 'where does this pose set dock' figure: a 3D centroid scatter per
    frame (colour = ligand, shape = tool) over the receptor envelope, plus a
    ligand×region occupancy heatmap and a per-ligand/tool pose-count bar. Used for
    both the JKU ligands (a) and the chemically-similar benchmark ligands (b)."""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from matplotlib.lines import Line2D
    d = d[d["cx"].notna() & (d["region"] >= 0)].copy()
    if d.empty:
        return None
    frames = sorted(d["frame"].unique())
    ligs = sorted(d["ligand"].unique())
    methods = sorted(d["method"].unique())
    cmap = plt.get_cmap("tab20" if len(ligs) > 10 else "tab10")
    lig_color = {l: cmap(i % cmap.N) for i, l in enumerate(ligs)}
    ncols = max(len(frames), 2)

    fig = plt.figure(figsize=(5.2 * ncols, 10))
    gs = fig.add_gridspec(2, ncols, height_ratios=[1.15, 1])
    axes3d = []
    for ci, fr in enumerate(frames):
        ax = fig.add_subplot(gs[0, ci], projection="3d"); axes3d.append(ax)
        rec = receptor_for_frame.get(fr)
        if rec is not None and Path(rec).exists():
            try:
                pxyz, _e, _r = load_protein_heavy_atoms(Path(rec))
                if len(pxyz) > max_ghost:
                    pxyz = pxyz[np.linspace(0, len(pxyz) - 1, max_ghost).astype(int)]
                ax.scatter(pxyz[:, 0], pxyz[:, 1], pxyz[:, 2], s=2, c="0.8",
                           alpha=0.08, marker=".", linewidths=0, depthshade=False)
            except Exception:
                pass
        s_fr = d[d["frame"] == fr]
        for l in ligs:
            for m in methods:
                s = s_fr[(s_fr["ligand"] == l) & (s_fr["method"] == m)]
                if s.empty:
                    continue
                ax.scatter(s.cx, s.cy, s.cz, s=40, color=lig_color[l],
                           marker=_marker(m), edgecolors="black", linewidths=0.4,
                           alpha=0.9, depthshade=False)
        ax.set_title(fr, fontsize=9)
        ax.set_xlabel("x (Å)"); ax.set_ylabel("y (Å)"); ax.set_zlabel("z (Å)")
    if axes3d:
        lh = [Line2D([0], [0], marker="o", color="w", markerfacecolor=lig_color[l],
                     markeredgecolor="black", markersize=8,
                     label=lig_label(l) + (f"  [{lig_sublabel[l]}]" if lig_sublabel else ""))
              for l in ligs]
        mh = [Line2D([0], [0], marker=_marker(m), color="0.3", ls="",
                     markerfacecolor="0.7", markeredgecolor="black", markersize=8,
                     label=_mlabel(m)) for m in methods]
        l1 = axes3d[0].legend(handles=lh, title="Ligand", loc="upper left", fontsize=7)
        axes3d[0].add_artist(l1)
        axes3d[-1].legend(handles=mh, title="Tool", loc="upper right", fontsize=7)

    # bottom-left: ligand × region occupancy heatmap (pooled over frames & tools)
    axh = fig.add_subplot(gs[1, :max(ncols - 1, 1)])
    mat = np.vstack([_region_fraction(d[d["ligand"] == l], region_ids).values for l in ligs])
    im = axh.imshow(mat, aspect="auto", cmap="magma", vmin=0, vmax=1)
    axh.set_xticks(range(len(region_ids)))
    axh.set_xticklabels([f"R{r + 1}" for r in region_ids])
    axh.set_yticks(range(len(ligs)))
    axh.set_yticklabels([lig_label(l) for l in ligs], fontsize=8)
    axh.set_xlabel("Binding region (benchmark consensus)")
    axh.set_title("Where each ligand docks — fraction of its PB-valid poses per "
                  "region\n(pooled over frames & tools)", fontsize=9)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if mat[i, j] > 0.02:
                axh.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center",
                         color="white" if mat[i, j] < 0.6 else "black", fontsize=7)
    fig.colorbar(im, ax=axh, fraction=0.046, pad=0.04, label="pose fraction")

    # bottom-right: PB-valid pose counts per ligand, stacked by tool
    axb = fig.add_subplot(gs[1, ncols - 1]) if ncols > 1 else None
    if axb is not None:
        ct = (pd.crosstab(d["ligand"], d["method"])
              .reindex(index=ligs, columns=methods, fill_value=0))
        xpos = np.arange(len(ligs)); bottom = np.zeros(len(ligs))
        for m in methods:
            axb.bar(xpos, ct[m].to_numpy(), bottom=bottom,
                    color=METHOD_COLOR.get(m, "#888888"), edgecolor="black",
                    linewidth=0.5, label=_mlabel(m))
            bottom += ct[m].to_numpy()
        axb.set_xticks(xpos)
        axb.set_xticklabels([lig_label(l) for l in ligs], rotation=25, ha="right",
                            fontsize=7)
        axb.set_ylabel("# PB-valid poses"); axb.set_title("Poses per ligand & tool",
                                                          fontsize=9)
        axb.legend(fontsize=7); axb.grid(axis="y", alpha=0.3)

    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    return out


def plot_jku_vs_similar(d_jku, d_sim, region_ids, reg_sig, out, top_n,
                        occ_thr: float = 0.10):
    """(c) Commonalities & differences across clusters across tools: per tool, the
    region-occupancy of JKU poses vs chemically-similar benchmark poses, and a
    region×tool map flagging where BOTH occupy (shared), only JKU, or only the
    similar set."""
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    methods = sorted(set(d_jku["method"].unique()) | set(d_sim["method"].unique()))
    regs = list(region_ids)
    fig = plt.figure(figsize=(max(5 * len(methods), 12), 9))
    gs = fig.add_gridspec(2, len(methods), height_ratios=[1, 0.95])

    top_axes = []
    for ci, m in enumerate(methods):
        ax = fig.add_subplot(gs[0, ci]); top_axes.append(ax)
        fj = _region_fraction(d_jku[d_jku.method == m], regs).values
        fs = _region_fraction(d_sim[d_sim.method == m], regs).values
        x = np.arange(len(regs)); w = 0.4
        ax.bar(x - w / 2, fj, w, color="#d62728", edgecolor="black", label="JKU")
        ax.bar(x + w / 2, fs, w, color="#1f77b4", edgecolor="black",
               label="similar benchmark")
        ax.set_xticks(x); ax.set_xticklabels([f"R{r + 1}" for r in regs], fontsize=7)
        ax.set_ylim(0, 1); ax.set_title(_mlabel(m), fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        if ci == 0:
            ax.set_ylabel("fraction of PB-valid poses")
        ax.legend(fontsize=7)

    # bottom: region × tool commonality map
    axc = fig.add_subplot(gs[1, :])
    codes = np.zeros((len(regs), len(methods)))
    for ci, m in enumerate(methods):
        fj = _region_fraction(d_jku[d_jku.method == m], regs)
        fs = _region_fraction(d_sim[d_sim.method == m], regs)
        for ri, r in enumerate(regs):
            j, s = fj[r] >= occ_thr, fs[r] >= occ_thr
            codes[ri, ci] = 3 if (j and s) else 1 if j else 2 if s else 0
    cmap = ListedColormap(["#eeeeee", "#f5b7b1", "#aed6f1", "#a9dfbf"])
    axc.imshow(codes, aspect="auto", cmap=cmap, vmin=0, vmax=3)
    axc.set_xticks(range(len(methods)))
    axc.set_xticklabels([_mlabel(m) for m in methods])
    axc.set_yticks(range(len(regs)))
    axc.set_yticklabels([f"R{r + 1}: {str(reg_sig.get(r, ''))[:30]}" for r in regs],
                        fontsize=7)
    txt = {0: "", 1: "JKU", 2: "sim", 3: "both"}
    for ri in range(len(regs)):
        for ci in range(len(methods)):
            t = txt[int(codes[ri, ci])]
            if t:
                axc.text(ci, ri, t, ha="center", va="center", fontsize=7)
    axc.set_title(f"Region occupied by ≥{int(occ_thr * 100)}% of a group's poses — "
                  "commonalities (both) vs differences (JKU-only / similar-only), per tool",
                  fontsize=9)
    leg = [Patch(facecolor="#a9dfbf", edgecolor="black", label="shared (both)"),
           Patch(facecolor="#f5b7b1", edgecolor="black", label="JKU only"),
           Patch(facecolor="#aed6f1", edgecolor="black", label="similar only"),
           Patch(facecolor="#eeeeee", edgecolor="black", label="neither")]
    axc.legend(handles=leg, fontsize=7, ncol=4, loc="upper center",
               bbox_to_anchor=(0.5, -0.07))
    _label_panels(top_axes + [axc])
    fig.suptitle("JKU vs chemically-similar benchmark ligands — where they cluster "
                 f"across tools (top {top_n}/tool, shared MD frames)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    return out


def run_jku_similarity_analysis(args, fp_bench, region_sum, feats, fold, km_reg,
                                remap_reg, region_colors, receptor_for_frame, raw_bench):
    """Load the JKU poses, project them onto the benchmark regions, find the
    chemically-similar benchmark ligands, rank the top-n per tool by native score,
    and emit figures (a)/(b)/(c) + CSVs. Returns a dict for summary.json."""
    jku_csv = Path(args.jku_csv)
    if not jku_csv.exists():
        print(f"  ! JKU results CSV not found ({jku_csv}) — skipping JKU analysis.")
        return {}
    print(f"\n=== JKU × chemically-similar benchmark analysis ===\nLoading JKU poses: {jku_csv}")
    raw_jku = pd.read_csv(jku_csv, low_memory=False)
    idx_jku = _build_pose_index(jku_csv, split_equibind=False)
    if "protein_file_used" in raw_jku.columns:
        idx_jku = idx_jku.merge(
            raw_jku[["pose_file", "protein_file_used"]].drop_duplicates("pose_file"),
            on="pose_file", how="left")
    fp_jku = compute_fingerprints(args, idx_jku, tag="jku", src_csv=jku_csv,
                                  apply_limit=False)
    if fp_jku.empty:
        print("  ! no JKU fingerprints — skipping."); return {}
    fp_jku["residues"] = fp_jku["residues"].fillna("")
    fp_jku["region"] = project_regions(fp_jku, feats, fold, km_reg, remap_reg)
    fp_jku = compute_pose_centroids(args, fp_jku)
    jku_ligs = sorted(fp_jku["ligand"].unique())
    shared_frames = sorted(set(fp_jku["frame"].unique()) & set(fp_bench["frame"].unique()))
    print(f"  JKU ligands: {[_jlabel(l) for l in jku_ligs]} · frames {sorted(fp_jku['frame'].unique())} "
          f"· shared with benchmark: {shared_frames}")

    # chemical similarity (benchmark ligand → nearest JKU)
    bench_ligs = sorted(fp_bench["ligand"].unique())
    sim = compute_chem_similarity(jku_ligs, bench_ligs, args.jku_template_dir,
                                  args.ligand_sdf_dir, args.similar_threshold)
    similar_ligs = sim.loc[sim["is_similar"], "ligand"].tolist()
    print(f"  chemically similar benchmark ligands (ECFP4 Tanimoto ≥ "
          f"{args.similar_threshold}): {len(similar_ligs)} → {similar_ligs[:12]}"
          f"{' …' if len(similar_ligs) > 12 else ''}")
    sim.to_csv(args.out_dir / "chem_similarity.csv", index=False)

    # native scores for ranking (only on the sets we rank: JKU + similar benchmark).
    # AutoDock affinity comes from the new CSV column when present, else recovered
    # from the source Vina PDBQT — see _autodock_affinity_lookup.
    lookup = {**_score_lookup(raw_bench), **_score_lookup(raw_jku)}
    autodock_aff = {**_autodock_affinity_lookup(raw_bench),
                    **_autodock_affinity_lookup(raw_jku)}
    fp_jku = fp_jku.copy()
    fp_jku["native_score"] = attach_native_scores(fp_jku, lookup, autodock_aff)
    sel_jku, rank_jku = top_n_ligands_per_tool(fp_jku, args.top_n)

    fp_sim = fp_bench[fp_bench["ligand"].isin(similar_ligs)
                      & fp_bench["frame"].isin(shared_frames)].copy()
    fp_sim = compute_pose_centroids(args, fp_sim)
    rank_sim = pd.DataFrame()
    if not fp_sim.empty:
        fp_sim["native_score"] = attach_native_scores(fp_sim, lookup, autodock_aff)
        sel_sim, rank_sim = top_n_ligands_per_tool(fp_sim, args.top_n)
        fp_sim_top = _filter_top_n(fp_sim, sel_sim)
    else:
        fp_sim_top = fp_sim
        print("  ! no chemically-similar benchmark ligands with PB-valid poses "
              "in the shared frames.")

    # persist tables
    region_ids = list(region_sum["region_id"])
    reg_sig = dict(zip(region_sum["region_id"], region_sum["signature_residues"]))
    rank_jku.assign(set="JKU").to_csv(args.out_dir / "topn_jku_per_tool.csv", index=False)
    if not rank_sim.empty:
        rank_sim.assign(set="similar_benchmark").to_csv(
            args.out_dir / "topn_similar_per_tool.csv", index=False)
    keep = ["method", "frame", "ligand", "region", "cx", "cy", "cz", "native_score",
            "pose_file"]
    fp_jku[[c for c in keep if c in fp_jku]].to_csv(
        args.out_dir / "jku_pose_region_assignments.csv", index=False)

    # figures (a) (b) (c)
    print("  generating JKU figures …")
    p_a = plot_pose_set_where(
        fp_jku, region_ids, receptor_for_frame,
        "(a) Where the JKU ligands dock in the Orai frames",
        args.out_dir / "06a_jku_where.png", lig_label=_jlabel)
    sub = {l: f"~{r.nearest_jku and _jlabel(r.nearest_jku)} {r.max_tanimoto:.2f}"
           for l, r in sim.set_index("ligand").iterrows()}
    p_b = None
    if not fp_sim_top.empty:
        p_b = plot_pose_set_where(
            fp_sim_top, region_ids, receptor_for_frame,
            f"(b) Where chemically-similar benchmark ligands dock "
            f"(top {args.top_n}/tool, ECFP4 Tanimoto ≥ {args.similar_threshold})",
            args.out_dir / "06b_similar_where.png",
            lig_label=str, lig_sublabel=sub)
    p_c = None
    if not fp_sim_top.empty:
        fp_jku_shared = fp_jku[fp_jku["frame"].isin(shared_frames)]
        p_c = plot_jku_vs_similar(fp_jku_shared, fp_sim_top, region_ids, reg_sig,
                                  args.out_dir / "06c_jku_vs_similar.png", args.top_n)

    made = [p.name for p in (p_a, p_b, p_c) if p is not None]
    if made:
        print(f"  JKU figures: {', '.join(made)}")
    return {
        "jku_ligands": jku_ligs, "jku_frames": sorted(fp_jku["frame"].unique()),
        "shared_frames": shared_frames,
        "n_similar_benchmark_ligands": len(similar_ligs),
        "similar_threshold": args.similar_threshold, "top_n": args.top_n,
        "jku_figures": made,
    }


# ───────────────────────────────────────────────────────────────────
# Main
# ───────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pb-csv", type=Path, default=DEFAULT_PB_CSV)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--features-csv", type=str, default=DEFAULT_FEATURES_CSV,
                    help="Per-ligand descriptor table from PoseBusters_DataSet_Analysis.")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 1))
    ap.add_argument("--contact-cutoff", type=float, default=CONTACT_CUTOFF,
                    help="Å within which a receptor residue counts as contacted "
                         "(default %(default)s).")
    ap.add_argument("--n-regions", type=int, default=0,
                    help="Fix the number of binding-region clusters (0 = pick by "
                         "silhouette).")
    ap.add_argument("--n-ligand-clusters", type=int, default=0,
                    help="Fix the number of ligand chemo-clusters (0 = silhouette).")
    ap.add_argument("--no-fold-symmetry", action="store_true",
                    help="Keep the 6 subunits (chains A-F) distinct instead of folding "
                         "them onto residue number (Orai is a symmetric homohexamer).")
    ap.add_argument("--purity-threshold", type=float, default=0.5,
                    help="Min dominant-region pose fraction for a ligand to count as "
                         "docking 'consistently' (default %(default)s).")
    ap.add_argument("--limit-ligands", type=int, default=0,
                    help="Process only the first N ligands (smoke test).")
    ap.add_argument("--no-3d", action="store_true",
                    help="Skip the per-snapshot 3D pose-cluster views "
                         "(05_snapshot3d_<frame>.png).")
    ap.add_argument("--ghost-atoms", type=int, default=3000,
                    help="Max receptor heavy atoms drawn as the faint channel "
                         "envelope in the 3D views (default %(default)s).")
    ap.add_argument("--no-pdb", action="store_true",
                    help="Skip the per-frame combined pose PDBs + PyMOL/ChimeraX "
                         "scripts (poses_<frame>.pdb, view_<frame>.pml/.cxc).")
    ap.add_argument("--jku-csv", type=Path, default=DEFAULT_JKU_CSV,
                    help="PoseBusters results CSV for the JKU ligands docked into "
                         "the same Orai frames (default %(default)s).")
    ap.add_argument("--jku-template-dir", type=Path, default=DEFAULT_JKU_TEMPLATE_DIR,
                    help="Dir with the clean JKU ligand SDFs (<name>.sdf) used for "
                         "ECFP4 fingerprints.")
    ap.add_argument("--ligand-sdf-dir", type=Path, default=DEFAULT_LIGAND_SDF_DIR,
                    help="Dir with benchmark <ligand>_ligand_start_conf.sdf files "
                         "(for ECFP4 fingerprints).")
    ap.add_argument("--top-n", type=int, default=5,
                    help="Top N ligands per tool (by native docking score) to use in "
                         "the JKU/similar analysis (default %(default)s).")
    ap.add_argument("--similar-threshold", type=float, default=0.20,
                    help="ECFP4 Tanimoto ≥ this to count a benchmark ligand as "
                         "chemically similar to a JKU ligand (default %(default)s). "
                         "The JKU drugs are structurally distinct from the PDB "
                         "benchmark — max observed Tanimoto is ~0.30, so this is a "
                         "deliberately loose 'nearest-analog' cut.")
    ap.add_argument("--no-jku", action="store_true",
                    help="Skip the JKU × chemically-similar benchmark cross-analysis "
                         "(figures 06a/06b/06c).")
    ap.add_argument("--force", action="store_true",
                    help="Recompute contact fingerprints even if a matching cache exists.")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fold = not args.no_fold_symmetry

    print(f"Loading Orai×Benchmark pose index: {args.pb_csv}")
    raw = pd.read_csv(args.pb_csv, low_memory=False)
    idx = _build_pose_index(args.pb_csv, split_equibind=False)
    # carry protein_file_used for receptor lookup
    if "protein_file_used" in raw.columns:
        idx = idx.merge(raw[["pose_file", "protein_file_used"]].drop_duplicates("pose_file"),
                        on="pose_file", how="left")
    n_valid = int(idx["pb_valid"].astype(bool).sum())
    print(f"  poses: {len(idx):,}  PB-valid: {n_valid:,} "
          f"({100 * n_valid / len(idx):.1f}%)  ligands: {idx['ligand'].nunique()}  "
          f"frames: {sorted(idx['protein'].unique())}")

    fp = compute_fingerprints(args, idx)
    if fp.empty:
        print("No fingerprints produced — aborting."); return
    fp["residues"] = fp["residues"].fillna("")
    # Drop any stale annotation columns so re-runs from cache don't collide on merge.
    fp = fp.drop(columns=[c for c in ("region", "lig_cluster") if c in fp.columns])
    all_methods = sorted(fp["method"].unique())
    all_frames = sorted(fp["frame"].unique())
    n_lig_valid = fp["ligand"].nunique()
    print(f"  fingerprinted {len(fp):,} PB-valid poses · {n_lig_valid} ligands · "
          f"methods={[_mlabel(m) for m in all_methods]} · {len(all_frames)} frames")

    # ── (1) cluster poses into binding regions ──────────────────────
    labels, feats, k_reg, sil, km_reg, remap_reg = cluster_regions(
        fp, fold, args.n_regions or None)
    fp["region"] = labels
    n_empty = int((labels < 0).sum())
    print(f"\nBinding regions: k={k_reg} (silhouette={sil:.3f}); "
          f"{n_empty} poses had no receptor contact and were left unclustered.")

    # ── (2) cluster ligands into chemo-types ────────────────────────
    desc, prof, lfeats, k_lig = cluster_ligands(
        sorted(fp["ligand"].unique()), args.features_csv,
        args.n_ligand_clusters or None)
    if desc is not None and "lig_cluster" in desc.columns:
        fp = fp.merge(desc["lig_cluster"], left_on="ligand", right_index=True, how="left")
        print(f"Ligand chemo-clusters: k={k_lig} "
              f"({desc['lig_cluster'].notna().sum()} ligands classified)")
    else:
        fp["lig_cluster"] = np.nan
        print("Ligand descriptor table unavailable — skipping chemo-cluster analysis.")

    # ── (3) cross-tabulation & consensus scoring ────────────────────
    region_sum = build_region_summary(fp, fold, n_lig_valid, all_methods, all_frames)
    per_lig = per_ligand_consistency(fp, all_methods, all_frames, args.purity_threshold)
    regions_order = list(region_sum["region_id"])
    lcr = ligandcluster_region_matrix(fp, regions_order)
    lcc = ligandcluster_consensus(fp, per_lig, prof)

    # Annotated per-pose table (kept SEPARATE from the fingerprint cache so the
    # cache stays pure fingerprints and re-runs don't re-merge stale columns).
    fp.to_csv(args.out_dir / "pose_region_assignments.csv", index=False)
    region_sum.to_csv(args.out_dir / "region_summary.csv", index=False)
    if desc is not None:
        cols = ["lig_cluster"] + [c for c in (lfeats or []) if c in desc.columns]
        lig_out = desc[cols].join(
            per_lig.set_index("ligand")[["dominant_region", "region_purity_%",
                                         "n_methods_in_dom", "n_frames_in_dom",
                                         "all_methods_agree", "all_frames_agree",
                                         "consistent"]], how="left")
        lig_out.to_csv(args.out_dir / "ligand_clusters.csv")
    per_lig.to_csv(args.out_dir / "per_ligand_consistency.csv", index=False)
    if not lcr.empty:
        lcr.to_csv(args.out_dir / "ligandcluster_region_matrix.csv")
    if not lcc.empty:
        lcc.to_csv(args.out_dir / "ligandcluster_consensus.csv", index=False)

    # ── (4) figures ─────────────────────────────────────────────────
    print("\nGenerating figures …")
    plot_region_overview(fp, region_sum, all_methods, all_frames,
                         args.out_dir / "01_region_overview.png")
    plot_ligand_clusters(desc, prof, lfeats, args.out_dir / "02_ligand_chemo_clusters.png")
    if not lcr.empty:
        plot_ligandcluster_region(lcr, lcc, args.out_dir / "03_ligandcluster_region_consensus.png")
    if not per_lig.empty:
        plot_consistency(per_lig, region_sum, args.out_dir / "04_ligand_consistency.png")

    # Per-snapshot 3D views + combined pose PDBs share a global region→colour map
    # (so a region keeps its colour across every snapshot, figure and viewer) and a
    # receptor-PDB-per-frame lookup. Both are per-frame because pose coordinates are
    # only comparable within a single Orai snapshot.
    reg_ids = sorted(int(r) for r in fp["region"].unique() if r >= 0)
    cmap = plt.get_cmap("tab20" if len(reg_ids) > 10 else "tab10")
    region_colors = {r: cmap(i % cmap.N) for i, r in enumerate(reg_ids)}
    receptor_for_frame = {fr: _receptor_path_for_frame(g, fr)
                          for fr, g in idx.groupby("protein")}

    snapshot3d_files: list[str] = []
    if not args.no_3d:
        fp = compute_pose_centroids(args, fp)
        for fr in all_frames:
            p = plot_snapshot_clusters_3d(
                fp, fr, region_colors,
                args.out_dir / f"05_snapshot3d_{_safe(fr)}.png",
                receptor_for_frame.get(fr), max_ghost=args.ghost_atoms)
            if p is not None:
                snapshot3d_files.append(p.name)
        if snapshot3d_files:
            print(f"  3D snapshot views: {', '.join(snapshot3d_files)}")

    # Combined per-frame pose PDBs + PyMOL/ChimeraX scripts (chain = tool,
    # resName/B-factor = cluster) to inspect the clusters on the real receptor.
    pose_pdb_files: list[str] = []
    if not args.no_pdb:
        print("Writing combined pose PDBs + viewer scripts …")
        pose_pdb_files = export_frame_pose_pdbs(
            args, fp, region_colors, receptor_for_frame, args.out_dir)

    # JKU ligands × chemically-similar benchmark ligands cross-analysis (a/b/c).
    jku_summary: dict = {}
    if not args.no_jku:
        if "cx" not in fp.columns:                 # need benchmark centroids for (b)
            fp = compute_pose_centroids(args, fp)
        jku_summary = run_jku_similarity_analysis(
            args, fp, region_sum, feats, fold, km_reg, remap_reg,
            region_colors, receptor_for_frame, raw)

    # ── (5) printed answer + summary.json ───────────────────────────
    print("\n=== (a) Consensus regions shared across ligands ===")
    print(region_sum[["region", "n_poses", "n_ligands", "ligand_coverage_%",
                      "n_methods", "n_frames", "broad_consensus",
                      "signature_residues"]].to_string(index=False))
    if not lcc.empty:
        print("\n=== (b) Does each ligand type dock in a consistent region? ===")
        print(lcc.to_string(index=False))
    n_consistent = int(per_lig["consistent"].sum()) if not per_lig.empty else 0
    print(f"\n{n_consistent}/{len(per_lig)} ligands dock 'consistently' "
          f"(≥{100 * args.purity_threshold:.0f}% of PB-valid poses in one region, "
          f"agreed by all methods AND all frames).")

    summary = {
        "n_pb_valid_poses": int(len(fp)),
        "n_ligands": int(n_lig_valid), "n_total_ligands": N_TOTAL_LIGANDS,
        "methods": all_methods, "frames": all_frames,
        "n_regions": int(k_reg), "region_silhouette": sil,
        "n_ligand_clusters": int(k_lig),
        "contact_cutoff_A": args.contact_cutoff, "fold_symmetry": fold,
        "n_consistent_ligands": n_consistent,
        "broad_consensus_regions": region_sum.loc[
            region_sum.get("broad_consensus", pd.Series(dtype=bool)) == True,
            "region"].tolist() if "broad_consensus" in region_sum else [],
        "snapshot3d_figures": snapshot3d_files,
        "pose_pdb_files": pose_pdb_files,
        "jku_analysis": jku_summary,
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nAll figures + CSVs → {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
