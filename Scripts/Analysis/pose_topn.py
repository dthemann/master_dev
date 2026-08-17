#!/usr/bin/env python3
"""Top-N pose cap shared by the Orai comparison analyses.

"Top-N poses of each tool" = each tool's N best-ranked poses per (protein, ligand),
using the ranking native to each emitted variant: raw AutoDock = Vina mode / affinity
rank, optimized AutoDock = ``optimized_rank`` from ``optimization_log.csv``, tiled
Uni-Dock = merged MODEL / affinity order, DiffDock = confidence rank (from the pose
filename), and EquiBind = gnina-minimisation-energy rank (most-negative
``gnina_affinity`` = 1, since EquiBind emits no native confidence score).

The cap is applied on the FULL (transmembrane-inclusive) produced-pose set, so it is
"top-N first, then any downstream TM / PB-valid filtering" — never "filter first, then
top-N of the survivors". Only groups with more than N rankable poses are trimmed; smaller
groups and poses with no derivable rank are kept untouched.

Lightweight (stdlib + pandas/numpy only) so it can be imported from analysis scripts and
notebook cells alike.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Set

import numpy as np
import pandas as pd


def _tool_key(method) -> str:
    m = str(method).strip().lower()
    if m.startswith("autodock") or m == "vina":
        return "autodock"
    if m.startswith("diffdock"):
        return "diffdock"
    if m.startswith("unidock") or m == "uni-dock":
        return "unidock"
    if m.startswith("equibind"):
        return "equibind"
    return m


def _native_rank(method, pose_name, pose_file, autodock_rank=None,
                 optimized_rank=None, unidock_rank=None, optimizer=None):
    """Native rank of a pose (1 = top). None if it cannot be derived from the columns."""
    m = _tool_key(method)
    pf = str(pose_file)
    name = Path(pf).name if pf and pf != "nan" else str(pose_name)
    if m == "autodock":
        opt = str(optimizer or "original").strip().lower()
        if opt not in {"", "original", "raw", "native", "none", "nan"}:
            if optimized_rank is not None and pd.notna(optimized_rank):
                try:
                    return int(float(optimized_rank))
                except (TypeError, ValueError):
                    return None
            # An optimized pose without its optimizer rank is not safely
            # rankable; falling back to the inherited Vina rank would silently
            # defeat re-ranking.
            return None
        if autodock_rank is not None and pd.notna(autodock_rank):
            try:
                return int(float(autodock_rank))
            except (TypeError, ValueError):
                pass
        mt = re.search(r"_(?:model|pose)(\d+)", name)
        return int(mt.group(1)) if mt else None
    if m == "unidock":
        if unidock_rank is not None and pd.notna(unidock_rank):
            try:
                return int(float(unidock_rank))
            except (TypeError, ValueError):
                pass
        mt = re.search(r"_(?:model|pose|rank)(\d+)", name)
        return int(mt.group(1)) if mt else None
    if m == "diffdock":
        mt = re.search(r"rank(\d+)", name)
        return int(mt.group(1)) if mt else None
    mt = re.search(r"(?:unguided|guided|pose|rank|model)_?(\d+)", name)
    return int(mt.group(1)) if mt else None


def _variant_key(row) -> str:
    """Isolate a tool's refinement variant so ranking stays within one pose set
    (raw / smina / gnina are distinct pose files that must not be pooled when ranking)."""
    m = _tool_key(row.get("docking_method"))
    if m == "autodock":
        opt = str(row.get("optimizer", "original") or "original").lower()
        if opt in {"", "raw", "native", "none", "nan"}:
            opt = "original"
        return f"autodock_{opt}"
    if m == "diffdock":
        return f"diffdock_{str(row.get('optimizer', '')).lower()}"
    if m == "equibind":
        rv = str(row.get("refine_variant", "") or row.get("optimizer", "")).lower()
        return f"equibind_{rv}"
    return m


def _group_rank(sub: pd.DataFrame) -> pd.Series:
    m = _tool_key(sub["docking_method"].iloc[0])
    if m == "equibind":
        gaff = pd.to_numeric(sub.get("gnina_affinity"), errors="coerce")
        order = gaff.rank(method="first", ascending=True)     # most-negative energy = 1
        order[gaff.isna()] = np.nan
        return order
    def _column(name):
        return (sub[name] if name in sub.columns
                else pd.Series([None] * len(sub), index=sub.index))

    ar = _column("autodock_rank")
    opt_r = _column("optimized_rank")
    uni_r = _column("unidock_rank")
    optimizer = _column("optimizer")
    return pd.Series(
        [_native_rank(m, pn, pf, a, o, u, v) for m, pn, pf, a, o, u, v in
         zip(sub["docking_method"], _column("pose_name"), _column("pose_file"),
             ar, opt_r, uni_r, optimizer)],
        index=sub.index, dtype="float64")


def cap_frame(df: pd.DataFrame, top_n: int = 10) -> pd.DataFrame:
    """Return *df* with only each tool-variant's top-N poses per (protein, ligand)."""
    if not top_n or top_n <= 0:
        return df
    df = df.copy()
    df["_vk"] = df.apply(_variant_key, axis=1)
    keep = []
    for _, sub in df.groupby(["_vk", "protein", "ligand"], sort=False):
        if len(sub) <= top_n:
            keep.extend(sub.index)
            continue
        r = _group_rank(sub)
        sel = sub.index[(r <= top_n).fillna(False).values]
        keep.extend(sel if len(sel) else sub.index[:top_n])   # nothing rankable -> file order
    return df.loc[sorted(keep)].drop(columns=["_vk"])


def top_n_allowlist(csv, top_n: int = 10) -> Set[str]:
    """Set of ``pose_file`` strings surviving a top-N cap of a per-pose CSV."""
    df = pd.read_csv(csv, low_memory=False)
    return set(cap_frame(df, top_n)["pose_file"].astype(str))


def top_n_keys(csv, top_n: int = 10) -> Set[tuple]:
    """Set of ``(protein, ligand, pose-file-basename)`` identity keys surviving a top-N
    cap. Path-independent, so it matches poses across pipelines that store pose_file with
    different absolute paths (e.g. the PandaMap summary vs the PoseBusters CSV)."""
    df = pd.read_csv(csv, low_memory=False)
    cap = cap_frame(df, top_n)
    return set(zip(cap["protein"].astype(str), cap["ligand"].astype(str),
                   cap["pose_file"].astype(str).map(lambda p: Path(str(p)).name)))


def cap_csv(in_csv, out_csv, top_n: int = 10) -> pd.DataFrame:
    """Write a top-N-capped copy of a per-pose CSV; returns the capped frame."""
    capped = cap_frame(pd.read_csv(in_csv, low_memory=False), top_n)
    capped.to_csv(out_csv, index=False)
    return capped
