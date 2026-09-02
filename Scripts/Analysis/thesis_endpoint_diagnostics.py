#!/usr/bin/env python3
"""Regenerate the endpoint diagnostics quoted in the Results chapter.

Two blocks of numbers in the thesis had no committed generator until now:

1. VALIDITY GATE COST — what the PoseBusters conjunct costs the headline endpoint.
   Counts complexes recovered under "within 2 A" alone and under "within 2 A AND
   PB-valid", per variant and ranking depth. The gap is what the gate discards.

2. BOUNDED CLAIM — the AutoDock Vina + gnina vs DiffDock + smina contrast reported
   as a paired difference with a Newcombe interval, a power floor and a TOST
   equivalence verdict, rather than as a failed significance test.

Run:
    python Scripts/Analysis/thesis_endpoint_diagnostics.py [--csv OUTDIR]

Conventions that matter and are easy to get wrong:
  * Near-nativeness uses the `rmsd` column, NOT `pb_rmsd`. Only `rmsd` reproduces
    the pose counts in the pose-production table.
  * A complex counts when ANY pose in the top-k pool meets the criterion
    (the "existence gate"), matching the recovery tables.
  * Ranking source differs per arm: AutoDock + gnina by `optimized_rank`, DiffDock
    by its own confidence (`rank`), EquiBind by gnina affinity ascending, since
    EquiBind writes a 999 sentinel into `rank`.
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import method_filter as mf  # noqa: E402  (shared single-point method exclusion)
from stats_utils import (  # noqa: E402
    newcombe_paired_diff_ci, mcnemar_power, tost_paired_proportions, mcnemar_exact,
)

DEFAULT_METRICS = ('posebusters_results/benchmark_full_protein_vina_scoring/dock/'
                   'pose_comparison_report/per_pose_metrics.csv')
N_SET = 303
FULL_BENCHMARK = 428
DEPTHS = (1, 5, 15, 30)

ARMS = [
    ('AutoDock Vina (raw)', 'autodock_mgltools_exh128', 'rank'),
    ('AutoDock Vina + gnina', 'autodock_mgltools_exh128_gnina', 'optimized_rank'),
    ('DiffDock (raw)', 'diffdock', 'rank'),
    ('DiffDock + smina', 'diffdock_smina', 'rank'),
    ('DiffDock + gnina', 'diffdock_gnina', 'rank'),
    ('EquiBind + gnina', 'equibind_unguided_gnina', 'rank'),
]
HEADLINE = ('AutoDock Vina + gnina', 'DiffDock + smina', 'EquiBind + gnina')


def load(path):
    df = pd.read_csv(path, low_memory=False)
    df['cid'] = df.protein.astype(str) + '_' + df.ligand.astype(str)
    return df


def _effective_rank(x, rankcol):
    """EquiBind has no native ranking (999 sentinel), so order by gnina affinity,
    most negative first. Mirrors _effective_rank in posebusters_pose_comparison.py."""
    r = pd.to_numeric(x[rankcol], errors='coerce') if rankcol in x else None
    if r is None or r.notna().sum() == 0 or (r.dropna() == 999).all():
        aff = pd.to_numeric(x['gnina_affinity'], errors='coerce')
        return aff.groupby(x['cid']).rank(method='first', ascending=True)
    return r


def recovered(df, meth, depth, rankcol, require_valid):
    x = df[df.method == meth].copy()
    x['_rk'] = _effective_rank(x, rankcol)
    x = x[x._rk <= depth]
    ok = x.rmsd <= 2.0
    if require_valid:
        ok = ok & x.pb_valid.astype('boolean').fillna(False)
    return x.assign(ok=ok).groupby('cid').ok.any()


def gate_cost(df):
    rows = []
    for label, meth, rc in ARMS:
        for d in DEPTHS:
            near = int(recovered(df, meth, d, rc, False).sum())
            gated = int(recovered(df, meth, d, rc, True).sum())
            rows.append({'arm': label, 'k': d, 'near_native': near, 'gated': gated,
                         'cost_complexes': near - gated,
                         'cost_pp': round(100 * (near - gated) / N_SET, 2)})
    return pd.DataFrame(rows)


def bounded_claim(df, arm_a=('autodock_mgltools_exh128_gnina', 'optimized_rank'),
                  arm_b=('diffdock_smina', 'rank')):
    # arm_a must track the dominant AutoDock arm. It was left at the superseded
    # Meeko 'autodock_gnina' after the exh128 migration, which silently published
    # the old arm's contrast into the appendix. Reported figures are AutoDock
    # minus DiffDock, i.e. the NEGATION of this function's diff_pp/lo_pp/hi_pp,
    # which are p_b - p_a.
    rows = []
    for d in DEPTHS:
        A = recovered(df, arm_a[0], d, arm_a[1], True)
        B = recovered(df, arm_b[0], d, arm_b[1], True)
        idx = A.index.union(B.index)
        a = A.reindex(idx, fill_value=False).astype(int).values
        b = B.reindex(idx, fill_value=False).astype(int).values
        ci = newcombe_paired_diff_ci(a, b)
        pw = mcnemar_power(a, b)
        _, _, p = mcnemar_exact(a, b)
        row = {'k': d, 'n': ci['n'], 'A': int(a.sum()), 'B': int(b.sum()),
               'both': ci['n11'], 'A_only': ci['n10'], 'B_only': ci['n01'],
               'neither': ci['n00'], 'n_discordant': ci['n_discordant'],
               'diff_pp': round(100 * ci['diff'], 2),
               'lo_pp': round(100 * ci['lo'], 2), 'hi_pp': round(100 * ci['hi'], 2),
               'mcnemar_p': round(p, 4),
               'power': round(pw['observed_power'], 3),
               'mde_pp': round(100 * pw['mde'], 2),
               'n_needed': None if pd.isna(pw['n_needed']) else int(round(pw['n_needed']))}
        for m in (5, 10, 12, 15):
            row['equiv_%dpp' % m] = tost_paired_proportions(a, b, margin=m / 100)['equivalent']
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--metrics', default=DEFAULT_METRICS)
    ap.add_argument('--csv', default=None, help='directory to write both tables into')
    mf.add_method_filter_args(ap)
    args = ap.parse_args()

    df = load(args.metrics)
    df = mf.apply_method_filter(df, 'method', args, label='endpoint-diagnostics',
                                out_dir=args.csv)
    # ARMS drives the pivot's row order via reindex; an arm excluded from the
    # frame but left in the constant would print as a row of NaN rather than
    # being absent.
    global ARMS, HEADLINE
    _present = set(df['method'].astype(str))
    ARMS = [a for a in ARMS if a[1] in _present]
    if not ARMS:
        raise SystemExit('endpoint-diagnostics: every arm was excluded.')
    HEADLINE = tuple(h for h in HEADLINE if h in {a[0] for a in ARMS})

    gc = gate_cost(df)
    print('=' * 78)
    print('1. VALIDITY GATE COST  (complexes of %d)' % N_SET)
    print('=' * 78)
    piv = gc.pivot(index='arm', columns='k', values='cost_complexes')
    piv = piv.reindex([a[0] for a in ARMS])
    print(piv.to_string())
    head = gc[gc.arm.isin(HEADLINE)]
    raw = gc[gc.arm == 'DiffDock (raw)']
    print('\nHeadlined arms : %.2f to %.2f pp  (%d to %d complexes)'
          % (head.cost_pp.min(), head.cost_pp.max(),
             head.cost_complexes.min(), head.cost_complexes.max()))
    print('Raw DiffDock   : %.2f to %.2f pp  (%d to %d complexes)'
          % (raw.cost_pp.min(), raw.cost_pp.max(),
             raw.cost_complexes.min(), raw.cost_complexes.max()))

    bc = bounded_claim(df)
    print()
    print('=' * 78)
    print('2. BOUNDED CLAIM  AutoDock Vina + gnina (A) vs DiffDock + smina (B)')
    print('=' * 78)
    cols = ['k', 'A', 'B', 'both', 'A_only', 'B_only', 'neither', 'n_discordant',
            'diff_pp', 'lo_pp', 'hi_pp', 'mcnemar_p', 'power', 'mde_pp', 'n_needed']
    print(bc[cols].to_string(index=False))
    print()
    print('Equivalence (TOST, alpha = 0.05):')
    print(bc[['k'] + [c for c in bc.columns if c.startswith('equiv_')]].to_string(index=False))
    r1 = bc[bc.k == 1].iloc[0]
    scale = (N_SET / FULL_BENCHMARK) ** 0.5
    print('\nAt rank-1: %d of %d complexes discordant (%.1f%%). MDE %.1f pp at 80%% power; '
          'observed power %.2f.' % (r1.n_discordant, N_SET,
                                    100 * r1.n_discordant / N_SET, r1.mde_pp, r1.power))
    print('Docking the full %d-complex benchmark would still leave the MDE near %.1f pp.'
          % (FULL_BENCHMARK, r1.mde_pp * scale))

    if args.csv:
        os.makedirs(args.csv, exist_ok=True)
        gc.to_csv(os.path.join(args.csv, 'validity_gate_cost.csv'), index=False)
        bc.to_csv(os.path.join(args.csv, 'bounded_claim.csv'), index=False)
        print('\nwritten to %s' % args.csv)


if __name__ == '__main__':
    main()
