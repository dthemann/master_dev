"""Shared statistical helpers for the docking-analysis suite.

The `vina` conda env ships scipy / numpy / pandas / sklearn but **not**
statsmodels / scikit-posthocs / pingouin, so every test here is hand-rolled from
scipy primitives. The implementations of Holm, Cochran's Q, exact McNemar,
Wilcoxon rank-biserial, Dunn, Jonckheere-Terpstra, the co-occurrence permutation
test, Wilson CI and Cliff's delta were consolidated verbatim from the trusted
copies that already lived in ``pose_cluster_crystal_pocket_report.py`` and
``pocket_comparison_report.py`` (see STATISTICAL_VALIDATION_PLAN.md). New tests
(G-test / chi-square of independence with Cramer's V and adjusted residuals,
Cochran-Armitage trend, Mann-Whitney + Cliff's delta bundle, bootstrap CIs,
cluster bootstrap, and the paired-across-tools bundles) follow the same
scipy-only style.

Design notes that callers must respect (from the validation plan):
  * Unit of analysis: aggregate correlated poses to one value per complex
    (Benchmark) or per (frame, ligand) pair (Orai x JKU) BEFORE calling these,
    or pass ``clusters=`` to the bootstrap helpers. Feeding raw poses inflates n.
  * The design is paired (same complexes across tools) -> prefer the paired
    entry points (``paired_proportions``, ``paired_continuous``, ``mcnemar_exact``,
    ``wilcoxon_rankbiserial``) over independent-sample tests.
  * Correct for multiplicity with ``holm`` (small pairwise families) or
    ``bh_fdr`` (large families: the 20 PB checks, 14 descriptors, per-residue).
  * Report an effect size + CI alongside every p-value.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    # formatting / multiplicity
    "fmt_p", "p_stars", "holm", "bh_fdr",
    # proportions
    "wilson_ci", "cochran_q", "mcnemar_exact", "paired_proportions",
    "paired_2x2_association", "two_proportion_test", "cochran_armitage",
    "newcombe_paired_diff_ci", "mcnemar_power", "tost_paired_proportions",
    # continuous
    "wilcoxon_rankbiserial", "friedman_kendall_w", "paired_continuous",
    "kruskal_dunn", "dunn", "jonckheere", "mannwhitney_cliffs",
    "cliffs_delta", "cliffs_delta_ci", "hodges_lehmann", "median_diff_ci",
    # correlation
    "spearman", "kendall",
    # contingency
    "chi2_independence", "gtest_independence", "cramers_v",
    # resampling
    "permutation_test_paired", "bootstrap_ci", "cluster_bootstrap_ci",
    "consensus_perm_test",
]


# --------------------------------------------------------------------------- #
# formatting & multiplicity correction
# --------------------------------------------------------------------------- #
def fmt_p(p) -> str:
    """Compact p-value string for figure annotations."""
    if p is None or p != p:
        return "n/a"
    if p < 1e-4:
        return "<1e-4"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3f}"


def p_stars(p) -> str:
    """Significance stars: *** <.001, ** <.01, * <.05, else 'ns'."""
    if p is None or p != p:
        return ""
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def holm(pvals):
    """Holm-Bonferroni step-down adjusted p-values (input order preserved).

    NaN inputs pass through unchanged (mirrors bh_fdr): a non-finite raw p is not a
    hypothesis to correct. A degenerate pairwise test — e.g. two identical conditions
    whose paired differences are all zero, where wilcoxon_rankbiserial returns NaN —
    must never inherit the running step-down maximum and surface as a finite (let alone
    significant) adjusted value. The family size m counts only the finite p-values."""
    p = np.asarray(pvals, float)
    out = np.full(p.shape, np.nan)
    idx_ok = np.flatnonzero(~np.isnan(p))
    m = len(idx_ok)
    if m == 0:
        return out
    vals = p[idx_ok]
    order = np.argsort(vals)
    running = 0.0
    for rank, j in enumerate(order):
        running = max(running, (m - rank) * vals[j])
        out[idx_ok[j]] = min(running, 1.0)
    return out


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR-adjusted p-values (q-values), order preserved.

    Uses scipy.stats.false_discovery_control when available (scipy>=1.11);
    falls back to a manual step-up. NaNs are passed through unchanged.
    """
    p = np.asarray(pvals, float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if not ok.any():
        return out
    vals = p[ok]
    try:
        from scipy.stats import false_discovery_control
        out[ok] = false_discovery_control(vals, method="bh")
        return out
    except Exception:
        m = len(vals)
        order = np.argsort(vals)
        ranked = vals[order]
        adj = ranked * m / (np.arange(1, m + 1))
        adj = np.minimum.accumulate(adj[::-1])[::-1]
        res = np.empty(m)
        res[order] = np.clip(adj, 0, 1)
        out[ok] = res
        return out


# --------------------------------------------------------------------------- #
# proportions
# --------------------------------------------------------------------------- #
def wilson_ci(k, n, z=1.96):
    """Wilson score CI for a binomial rate k/n -> (lo, hi). Default 95%."""
    if not n:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def cochran_q(X):
    """Cochran's Q omnibus for k paired binary conditions.

    X : (n x k) array of 0/1 (n subjects, k conditions). Returns (Q, p, df).
    """
    from scipy.stats import chi2
    X = np.asarray(X, float)
    n, k = X.shape
    Cj = X.sum(axis=0)
    Ri = X.sum(axis=1)
    T = X.sum()
    denom = k * T - float(np.sum(Ri ** 2))
    if denom == 0:
        return np.nan, np.nan, k - 1
    Q = (k - 1) * (k * float(np.sum(Cj ** 2)) - T ** 2) / denom
    return float(Q), float(chi2.sf(Q, k - 1)), k - 1


def mcnemar_exact(a, b):
    """Exact (binomial) McNemar for two paired binary vectors.

    Returns (n10, n01, p) where n10 = #(a=1,b=0), n01 = #(a=0,b=1).
    """
    from scipy.stats import binomtest
    a = np.asarray(a); b = np.asarray(b)
    n10 = int(np.sum((a == 1) & (b == 0)))
    n01 = int(np.sum((a == 0) & (b == 1)))
    disc = n10 + n01
    if disc == 0:
        return n10, n01, 1.0
    p = binomtest(min(n10, n01), disc, 0.5).pvalue
    return n10, n01, float(min(p, 1.0))


def _counts_2x2(a, b):
    """Paired 2x2 cell counts for two binary vectors -> (n11, n10, n01, n00)."""
    a = np.asarray(a).astype(bool)
    b = np.asarray(b).astype(bool)
    if a.shape != b.shape:
        raise ValueError("paired vectors must have the same length")
    return (int(np.sum(a & b)), int(np.sum(a & ~b)),
            int(np.sum(~a & b)), int(np.sum(~a & ~b)))


def newcombe_paired_diff_ci(a, b, z=1.96):
    """Newcombe (1998) method 10 score interval for the difference of two
    CORRELATED proportions, p_b - p_a, on paired binary data.

    This is the interval to quote alongside McNemar. McNemar tests whether the
    difference is zero using the discordant pairs only; this bounds how large the
    difference could be, using the marginals and their correlation.

    a, b : paired 0/1 vectors of equal length (a = reference arm, b = comparator).
    Returns dict with the 2x2 cells, both marginal rates, the difference
    (p_b - p_a) and its (lo, hi). Default z gives a 95% interval.

    Reference: Newcombe RG, Stat Med 1998;17:2635-2650, method 10.
    """
    n11, n10, n01, n00 = _counts_2x2(a, b)
    n = n11 + n10 + n01 + n00
    if n == 0:
        nan = float("nan")
        return {"n": 0, "n11": 0, "n10": 0, "n01": 0, "n00": 0,
                "p_a": nan, "p_b": nan, "diff": nan, "lo": nan, "hi": nan,
                "n_discordant": 0, "phi": nan}
    p_a = (n11 + n10) / n
    p_b = (n11 + n01) / n
    lo_a, hi_a = wilson_ci(n11 + n10, n, z)
    lo_b, hi_b = wilson_ci(n11 + n01, n, z)

    # phi: correlation between the two binary variables; 0 when a margin is degenerate
    denom = (n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00)
    phi = ((n11 * n00 - n10 * n01) / np.sqrt(denom)) if denom > 0 else 0.0

    diff = p_b - p_a
    lo = diff - np.sqrt(max(0.0, (p_b - lo_b) ** 2
                            - 2 * phi * (p_b - lo_b) * (hi_a - p_a)
                            + (hi_a - p_a) ** 2))
    hi = diff + np.sqrt(max(0.0, (hi_b - p_b) ** 2
                            - 2 * phi * (hi_b - p_b) * (p_a - lo_a)
                            + (p_a - lo_a) ** 2))
    return {"n": n, "n11": n11, "n10": n10, "n01": n01, "n00": n00,
            "p_a": float(p_a), "p_b": float(p_b), "diff": float(diff),
            "lo": float(max(-1.0, lo)), "hi": float(min(1.0, hi)),
            "n_discordant": n10 + n01, "phi": float(phi)}


def mcnemar_power(a, b, alpha=0.05, power=0.80):
    """Observed power and minimum detectable effect for an exact McNemar contrast.

    A paired binary test conditions on the DISCORDANT pairs, so its power is
    governed by how often the two arms disagree, not by the sample size. This
    returns both the power actually achieved against the observed effect and the
    smallest difference the design could have detected at the requested power.

    a, b : paired 0/1 vectors. Returns dict with n_discordant, the discordant
    proportion, observed power, the minimum detectable difference in proportion
    units, and the number of pairs that would be needed for the observed effect.
    """
    from scipy.stats import norm
    n11, n10, n01, n00 = _counts_2x2(a, b)
    n = n11 + n10 + n01 + n00
    m = n10 + n01
    nan = float("nan")
    if n == 0 or m == 0:
        return {"n": n, "n_discordant": m, "p_discordant": nan,
                "observed_power": nan, "mde": nan, "n_needed": nan}
    psi = m / n                      # discordance rate
    diff = (n01 - n10) / n           # observed difference p_b - p_a
    z_a = norm.ppf(1 - alpha / 2)
    z_b = norm.ppf(power)

    # normal approximation to the conditional binomial on m discordant pairs
    if diff == 0:
        obs_power = alpha
    else:
        lam = abs(diff) * n / np.sqrt(m)      # noncentrality in SD units
        obs_power = float(norm.sf(z_a - lam) + norm.cdf(-z_a - lam))

    mde = (z_a + z_b) * np.sqrt(psi / n)      # detectable |p_b - p_a| at `power`
    n_needed = (((z_a + z_b) ** 2) * psi / diff ** 2) if diff != 0 else nan
    return {"n": n, "n_discordant": m, "p_discordant": float(psi),
            "diff": float(diff), "observed_power": float(min(1.0, obs_power)),
            "mde": float(mde), "n_needed": float(n_needed)}


def tost_paired_proportions(a, b, margin, alpha=0.05):
    """Two one-sided tests for EQUIVALENCE of two correlated proportions.

    Declares equivalence when the (1 - 2*alpha) Newcombe interval for p_b - p_a
    lies entirely inside (-margin, +margin). Use this instead of reading a
    non-significant McNemar as "no difference": a null result bounds nothing on
    its own, whereas this states the largest difference the data still allow.

    margin : equivalence bound in proportion units (e.g. 0.10 for 10 points).
    Returns dict with the interval used, the margin and an `equivalent` flag.
    """
    from scipy.stats import norm
    z = norm.ppf(1 - alpha)          # one-sided z -> (1-2*alpha) two-sided interval
    res = newcombe_paired_diff_ci(a, b, z=z)
    lo, hi = res["lo"], res["hi"]
    res.update({"margin": float(margin), "alpha": float(alpha),
                "conf_level": float(1 - 2 * alpha),
                "equivalent": bool(lo > -margin and hi < margin)})
    return res


def paired_proportions(data, labels=None, z=1.96):
    """Full paired-binary comparison across k conditions (the validity-report /
    reach-figure workhorse).

    data   : (n x k) 0/1 array, or dict{label: 1d 0/1}, or pandas DataFrame.
             Rows are the SAME units (complexes) across all k columns; drop rows
             with any missing condition before calling (listwise complete).
    Returns a dict with:
      rates     : {label: (rate, lo, hi)} Wilson CIs
      omnibus   : {Q, p, df}  Cochran's Q  (k>=2; identical to McNemar for k=2)
      pairwise  : [{a, b, a_wins, b_wins, n, p_raw, p_holm, star}, ...]
                  oriented so 'a' is the higher-rate condition of the pair.
      n, k
    """
    import pandas as pd
    if isinstance(data, dict):
        df = pd.DataFrame(data)
    elif isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        arr = np.asarray(data)
        cols = labels if labels is not None else [f"c{i}" for i in range(arr.shape[1])]
        df = pd.DataFrame(arr, columns=cols)
    if labels is None:
        labels = list(df.columns)
    df = df[labels].dropna()
    M = df.to_numpy(float)
    n, k = M.shape
    rates = {lab: (float(df[lab].mean()),) + wilson_ci(int(df[lab].sum()), n, z)
             for lab in labels}
    Q, pQ, dfQ = cochran_q(M) if n and k >= 2 else (np.nan, np.nan, max(k - 1, 0))
    pairwise = []
    for i in range(k):
        for j in range(i + 1, k):
            a, b = labels[i], labels[j]
            if rates[b][0] > rates[a][0]:
                a, b = b, a
            n10, n01, p = mcnemar_exact(df[a].to_numpy(), df[b].to_numpy())
            pairwise.append({"a": a, "b": b, "a_wins": n10, "b_wins": n01,
                             "n": n, "p_raw": float(p)})
    if pairwise:
        for pr, pa in zip(pairwise, holm([pp["p_raw"] for pp in pairwise])):
            pr["p_holm"] = float(pa)
            pr["star"] = p_stars(pa)
    return {"rates": rates, "omnibus": {"Q": Q, "p": pQ, "df": dfQ},
            "pairwise": pairwise, "n": n, "k": k}


def paired_2x2_association(a, b, min_margin=20, z=1.96):
    """Association between two paired binary vectors (same units): do the two
    conditions succeed on the SAME units more (or less) than independence predicts?

    This is the co-reach / co-occurrence test. Unlike ``mcnemar_exact`` -- whose
    null is marginal homogeneity ('do the two conditions succeed at equal RATES'),
    which uses only the discordant cells and is blind to the joint-success cell --
    this tests INDEPENDENCE of the two indicators and reports a SIGNED effect, so
    redundancy (both tend to succeed together, phi > 0) is distinguished from
    complementarity (one tends to succeed where the other fails, phi < 0).

    a, b : 1d 0/1 arrays over the SAME units (drop unpaired rows before calling).
    min_margin : an association is only estimable when every margin of the 2x2 has
        at least this many units; when a condition is near-constant (e.g. a tool at
        a success ceiling) the 'fails' margin collapses and phi/OR are undefined in
        practice -> returns estimable=False with the descriptive counts still filled.

    Returns dict:
      n, n11, n10, n01, n00   the 2x2 counts (1 = success; n11 = both succeed)
      rate_a, rate_b          marginal success rates of a and b
      observed, expected      joint-success count vs n*rate_a*rate_b (independence)
      phi                     SIGNED phi = sign(n11*n00 - n10*n01) * sqrt(chi2/n);
                              +1 perfect concordance, -1 perfect complementarity
      odds_ratio, or_ci       Haldane-Anscombe (+0.5) odds ratio + Wald log-CI
      p                       two-sided Fisher exact p (independence null)
      estimable, reason       whether every margin >= min_margin
    """
    from scipy.stats import fisher_exact
    a = np.asarray(a).astype(int); b = np.asarray(b).astype(int)
    n = int(a.size)
    n11 = int(np.sum((a == 1) & (b == 1)))
    n10 = int(np.sum((a == 1) & (b == 0)))
    n01 = int(np.sum((a == 0) & (b == 1)))
    n00 = int(np.sum((a == 0) & (b == 0)))
    rate_a = (n11 + n10) / n if n else np.nan
    rate_b = (n11 + n01) / n if n else np.nan
    expected = n * rate_a * rate_b if n else np.nan
    margins = [n11 + n10, n01 + n00, n11 + n01, n10 + n00]      # row0,row1,col0,col1
    estimable = n > 0 and min(margins) >= min_margin
    tab = np.array([[n11, n10], [n01, n00]], float)             # signed phi via chi2
    rs = tab.sum(1); cs = tab.sum(0); tot = tab.sum()
    if tot > 0 and np.all(rs > 0) and np.all(cs > 0):
        E = np.outer(rs, cs) / tot
        chi2 = float(np.sum((tab - E) ** 2 / E))
        phi = float(np.sign(n11 * n00 - n10 * n01) * np.sqrt(chi2 / tot))
    else:
        phi = np.nan
    a_, b_, c_, d_ = n11 + 0.5, n10 + 0.5, n01 + 0.5, n00 + 0.5  # Haldane-Anscombe OR
    OR = (a_ * d_) / (b_ * c_)
    se = np.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
    or_ci = (float(np.exp(np.log(OR) - z * se)), float(np.exp(np.log(OR) + z * se)))
    try:
        _, p = fisher_exact([[n11, n10], [n01, n00]], alternative="two-sided")
        p = float(p)
    except Exception:
        p = np.nan
    reason = "" if estimable else (
        f"min margin {min(margins)} < {min_margin}: a condition is near-constant, "
        "association not estimable")
    return {"n": n, "n11": n11, "n10": n10, "n01": n01, "n00": n00,
            "rate_a": float(rate_a), "rate_b": float(rate_b),
            "observed": n11, "expected": float(expected),
            "phi": phi, "odds_ratio": float(OR), "or_ci": or_ci,
            "p": p, "estimable": bool(estimable), "reason": reason}


def two_proportion_test(k1, n1, k2, n2):
    """Independent two-group proportion comparison (Fisher exact on the 2x2).

    Returns dict with Fisher p, risk difference, and Wilson CIs for each rate.
    Use for genuinely disjoint groups (not the paired same-complex design).
    """
    from scipy.stats import fisher_exact
    table = [[k1, n1 - k1], [k2, n2 - k2]]
    _, p = fisher_exact(table)
    p1 = k1 / n1 if n1 else np.nan
    p2 = k2 / n2 if n2 else np.nan
    return {"p": float(p), "risk_diff": float(p1 - p2),
            "rate1": (p1, *wilson_ci(k1, n1)), "rate2": (p2, *wilson_ci(k2, n2))}


def cochran_armitage(successes, totals, scores=None):
    """Cochran-Armitage test for a linear trend in proportions across ordered
    groups (e.g. success rate vs flexibility tertiles, TM-loss vs MD frame).

    successes, totals : per-group counts (ordered). scores : group scores
    (default 0,1,2,...). Returns (z, p, slope_sign). Normal approximation.
    """
    from scipy.stats import norm
    r = np.asarray(successes, float)
    n = np.asarray(totals, float)
    x = np.arange(len(r), dtype=float) if scores is None else np.asarray(scores, float)
    N = n.sum()
    R = r.sum()
    if N == 0 or R == 0 or R == N:
        return np.nan, np.nan, 0
    p = R / N
    xbar = np.sum(n * x) / N
    T = np.sum(r * (x - xbar))
    var = p * (1 - p) * np.sum(n * (x - xbar) ** 2)
    if var <= 0:
        return np.nan, np.nan, 0
    z = T / np.sqrt(var)
    return float(z), float(2 * norm.sf(abs(z))), int(np.sign(T))


# --------------------------------------------------------------------------- #
# continuous
# --------------------------------------------------------------------------- #
def wilcoxon_rankbiserial(x, y):
    """Wilcoxon signed-rank p + matched-pairs rank-biserial effect size.

    Returns (rank_biserial, p, n_nonzero_pairs). rank_biserial in [-1,1];
    positive => x tends to exceed y.
    """
    from scipy.stats import rankdata, wilcoxon
    x = np.asarray(x, float); y = np.asarray(y, float)
    d = x - y
    d = d[d != 0]
    n = len(d)
    if n < 1:
        return np.nan, np.nan, 0
    r = rankdata(np.abs(d))
    r_plus = float(np.sum(r[d > 0])); r_minus = float(np.sum(r[d < 0]))
    total = r_plus + r_minus
    rb = (r_plus - r_minus) / total if total else np.nan
    try:
        p = float(wilcoxon(x, y, zero_method="wilcox").pvalue)
    except Exception:
        p = np.nan
    return float(rb), p, n


def friedman_kendall_w(*cols):
    """Friedman omnibus + Kendall's W for k paired continuous conditions.

    Each col is a 1d array of the SAME n units. Returns (chi2, p, W, n, df).
    """
    from scipy.stats import friedmanchisquare
    cols = [np.asarray(c, float) for c in cols]
    n = len(cols[0]); k = len(cols)
    chi2, p = friedmanchisquare(*cols)
    W = float(chi2) / (n * (k - 1)) if n and k > 1 else np.nan
    return float(chi2), float(p), W, n, k - 1


def paired_continuous(data, labels=None):
    """Full paired-continuous comparison across k conditions (oracle-RMSD /
    interaction-count / wall-clock workhorse).

    data : (n x k) array, dict{label: 1d}, or DataFrame; rows are the same units
           across conditions (drop incomplete rows before calling).
    Returns dict:
      omnibus  : {chi2, p, df, kendall_w, n}   Friedman (+ Kendall's W)
      pairwise : [{a, b, rank_biserial, p_raw, p_holm, star, n}, ...]
      medians  : {label: median}
    """
    import pandas as pd
    if isinstance(data, dict):
        df = pd.DataFrame(data)
    elif isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        arr = np.asarray(data)
        cols = labels if labels is not None else [f"c{i}" for i in range(arr.shape[1])]
        df = pd.DataFrame(arr, columns=cols)
    if labels is None:
        labels = list(df.columns)
    df = df[labels].dropna()
    cols = [df[l].to_numpy(float) for l in labels]
    n, k = len(df), len(labels)
    chi2, p, W, _, dfree = friedman_kendall_w(*cols) if n and k > 1 else (np.nan, np.nan, np.nan, n, k - 1)
    pairwise = []
    for i in range(k):
        for j in range(i + 1, k):
            rb, praw, npair = wilcoxon_rankbiserial(cols[i], cols[j])
            pairwise.append({"a": labels[i], "b": labels[j], "rank_biserial": rb,
                             "p_raw": praw, "n": npair})
    if pairwise:
        for pr, pa in zip(pairwise, holm([pp["p_raw"] for pp in pairwise])):
            pr["p_holm"] = float(pa)
            pr["star"] = p_stars(pa)
    return {"omnibus": {"chi2": chi2, "p": p, "df": dfree, "kendall_w": W, "n": n},
            "pairwise": pairwise,
            "medians": {l: float(np.median(c)) for l, c in zip(labels, cols)}}


def dunn(groups):
    """Dunn's post-hoc (tie-corrected z, uncorrected two-sided p) for k groups.
    Returns [(i, j, z, p), ...]."""
    from scipy.stats import rankdata, norm
    data = np.concatenate([np.asarray(g, float) for g in groups])
    N = len(data)
    ranks = rankdata(data)
    _, counts = np.unique(data, return_counts=True)
    ties = float(np.sum(counts ** 3 - counts))
    sigma2 = (N * (N + 1) / 12.0) - ties / (12.0 * (N - 1))
    sizes, meanranks, idx = [], [], 0
    for g in groups:
        m = len(g)
        meanranks.append(float(np.mean(ranks[idx:idx + m])))
        sizes.append(m); idx += m
    out = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            se = np.sqrt(sigma2 * (1.0 / sizes[i] + 1.0 / sizes[j]))
            z = (meanranks[i] - meanranks[j]) / se if se else np.nan
            out.append((i, j, float(z), float(2 * norm.sf(abs(z)))))
    return out


def kruskal_dunn(groups, labels=None):
    """Kruskal-Wallis omnibus + epsilon^2 effect size + Dunn post-hoc (Holm).

    For k INDEPENDENT groups (e.g. Ro5 pass/fail, per-tool tightness with
    informative missingness). Returns dict{H, p, df, epsilon_sq, n, medians,
    posthoc:[{pair, z, p_raw, p_holm, star}]}.
    """
    from scipy.stats import kruskal
    groups = [np.asarray(g, float) for g in groups]
    groups = [g[~np.isnan(g)] for g in groups]
    k = len(groups)
    labels = labels or [f"g{i}" for i in range(k)]
    N = sum(len(g) for g in groups)
    if any(len(g) == 0 for g in groups) or k < 2:
        return None
    H, p = kruskal(*groups)
    eps2 = (float(H) - k + 1) / (N - k) if N > k else np.nan
    posthoc = []
    for (i, j, z, praw) in dunn(groups):
        posthoc.append({"pair": f"{labels[i]} vs {labels[j]}", "z": round(z, 3),
                        "p_raw": praw})
    if posthoc:
        for h, pa in zip(posthoc, holm([h["p_raw"] for h in posthoc])):
            h["p_holm"] = float(pa)
            h["star"] = p_stars(pa)
    return {"H": float(H), "p": float(p), "df": k - 1,
            "epsilon_sq": float(eps2),
            "n": {labels[i]: int(len(groups[i])) for i in range(k)},
            "medians": {labels[i]: float(np.median(groups[i])) for i in range(k)},
            "posthoc": posthoc}


def jonckheere(groups):
    """Jonckheere-Terpstra trend test across ORDERED groups (normal approx).
    Returns (JT, z, p)."""
    from scipy.stats import norm
    sizes = [len(g) for g in groups]
    N = sum(sizes)
    JT = 0.0
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            gj = np.asarray(groups[j], float)
            for a in np.asarray(groups[i], float):
                JT += float(np.sum(gj > a) + 0.5 * np.sum(gj == a))
    mean = (N ** 2 - sum(s ** 2 for s in sizes)) / 4.0
    var = (N ** 2 * (2 * N + 3) - sum(s ** 2 * (2 * s + 3) for s in sizes)) / 72.0
    z = (JT - mean) / np.sqrt(var) if var > 0 else np.nan
    return float(JT), float(z), float(2 * norm.sf(abs(z)))


def cliffs_delta(a, b) -> float:
    """Cliff's delta effect size: P(a>b) - P(a<b) in [-1, 1]."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if a.size == 0 or b.size == 0:
        return float("nan")
    gt = sum(np.sum(a > y) for y in b)
    lt = sum(np.sum(a < y) for y in b)
    return float((gt - lt) / (a.size * b.size))


def cliffs_delta_ci(a, b, n_boot=2000, seed=0, alpha=0.05):
    """Percentile bootstrap CI for Cliff's delta. Returns (delta, lo, hi)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    d = cliffs_delta(a, b)
    if a.size == 0 or b.size == 0:
        return d, np.nan, np.nan
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        aa = rng.choice(a, a.size, replace=True)
        bb = rng.choice(b, b.size, replace=True)
        boots[i] = cliffs_delta(aa, bb)
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return d, float(lo), float(hi)


def mannwhitney_cliffs(a, b):
    """Independent two-group bundle: Mann-Whitney U p + Cliff's delta + CI.

    Returns dict{U, p, cliffs_delta, delta_ci:(lo,hi), n_a, n_b, medians}.
    """
    from scipy.stats import mannwhitneyu
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if a.size == 0 or b.size == 0:
        return None
    U, p = mannwhitneyu(a, b, alternative="two-sided")
    d, lo, hi = cliffs_delta_ci(a, b)
    return {"U": float(U), "p": float(p), "cliffs_delta": d, "delta_ci": (lo, hi),
            "n_a": int(a.size), "n_b": int(b.size),
            "medians": (float(np.median(a)), float(np.median(b)))}


def hodges_lehmann(x, y=None):
    """Hodges-Lehmann estimator. One-sample (y=None) -> median of Walsh averages
    of x. Paired -> pass the differences as x. Two-sample -> median of pairwise
    differences x_i - y_j."""
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    if y is None:
        walsh = (x[:, None] + x[None, :]) / 2.0
        iu = np.triu_indices(len(x), k=0)
        return float(np.median(walsh[iu]))
    y = np.asarray(y, float); y = y[~np.isnan(y)]
    diffs = x[:, None] - y[None, :]
    return float(np.median(diffs))


def median_diff_ci(x, y, n_boot=2000, seed=0, alpha=0.05, paired=True):
    """Bootstrap CI for the (Hodges-Lehmann) median difference between paired or
    independent samples. Returns (estimate, lo, hi)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    rng = np.random.default_rng(seed)
    if paired:
        d = x - y
        d = d[~np.isnan(d)]
        est = hodges_lehmann(d)
        n = len(d)
        boots = np.array([hodges_lehmann(rng.choice(d, n, replace=True))
                          for _ in range(n_boot)])
    else:
        xm = x[~np.isnan(x)]; ym = y[~np.isnan(y)]
        est = hodges_lehmann(xm, ym)
        boots = np.array([hodges_lehmann(rng.choice(xm, len(xm), replace=True),
                                         rng.choice(ym, len(ym), replace=True))
                          for _ in range(n_boot)])
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return est, float(lo), float(hi)


# --------------------------------------------------------------------------- #
# correlation
# --------------------------------------------------------------------------- #
def spearman(x, y):
    """Spearman rho with p-value (nan_policy=omit). Returns (rho, p, n)."""
    from scipy.stats import spearmanr
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 3:
        return np.nan, np.nan, int(m.sum())
    r = spearmanr(x[m], y[m])
    return float(r.statistic), float(r.pvalue), int(m.sum())


def kendall(x, y):
    """Kendall tau-b with p-value. Returns (tau, p, n)."""
    from scipy.stats import kendalltau
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = ~(np.isnan(x) | np.isnan(y))
    if m.sum() < 3:
        return np.nan, np.nan, int(m.sum())
    r = kendalltau(x[m], y[m])
    return float(r.statistic), float(r.pvalue), int(m.sum())


# --------------------------------------------------------------------------- #
# contingency
# --------------------------------------------------------------------------- #
def cramers_v(table):
    """Bias-corrected Cramer's V for an r x c contingency table."""
    from scipy.stats import chi2_contingency
    table = np.asarray(table, float)
    chi2, _, _, _ = chi2_contingency(table, correction=False)
    n = table.sum()
    r, c = table.shape
    if n == 0:
        return np.nan
    phi2 = chi2 / n
    phi2corr = max(0.0, phi2 - (r - 1) * (c - 1) / (n - 1))
    rcorr = r - (r - 1) ** 2 / (n - 1)
    ccorr = c - (c - 1) ** 2 / (n - 1)
    denom = min(rcorr - 1, ccorr - 1)
    return float(np.sqrt(phi2corr / denom)) if denom > 0 else np.nan


def _adjusted_residuals(table):
    """Adjusted standardized residuals (O-E)/sqrt(E(1-p_row)(1-p_col)); |z|>~2
    flags cells that drive an association."""
    table = np.asarray(table, float)
    n = table.sum()
    rowsum = table.sum(axis=1, keepdims=True)
    colsum = table.sum(axis=0, keepdims=True)
    E = rowsum @ colsum / n
    with np.errstate(divide="ignore", invalid="ignore"):
        adj = (table - E) / np.sqrt(E * (1 - rowsum / n) * (1 - colsum / n))
    return adj


def chi2_independence(table):
    """Pearson chi-square test of independence + Cramer's V + adjusted residuals.

    Returns dict{chi2, p, df, cramers_v, residuals, expected, n_low_expected}.
    Note: chi-square approximation is unreliable when expected counts <5 are
    common (n_low_expected) -> prefer gtest_independence or a permutation there.
    """
    from scipy.stats import chi2_contingency
    table = np.asarray(table, float)
    chi2, p, dof, exp = chi2_contingency(table, correction=False)
    return {"chi2": float(chi2), "p": float(p), "df": int(dof),
            "cramers_v": cramers_v(table),
            "residuals": _adjusted_residuals(table), "expected": exp,
            "n_low_expected": int(np.sum(exp < 5))}


def gtest_independence(table):
    """G-test (likelihood-ratio) of independence + Cramer's V + adjusted
    residuals. More robust than Pearson chi-square for sparse tables.

    Returns dict{G, p, df, cramers_v, residuals, expected, n_low_expected}.
    """
    from scipy.stats import chi2 as _chi2, chi2_contingency
    table = np.asarray(table, float)
    _, _, dof, exp = chi2_contingency(table, correction=False)
    mask = table > 0
    G = 2.0 * np.sum(table[mask] * np.log(table[mask] / exp[mask]))
    p = float(_chi2.sf(G, dof))
    return {"G": float(G), "p": p, "df": int(dof), "cramers_v": cramers_v(table),
            "residuals": _adjusted_residuals(table), "expected": exp,
            "n_low_expected": int(np.sum(exp < 5))}


# --------------------------------------------------------------------------- #
# resampling
# --------------------------------------------------------------------------- #
def permutation_test_paired(a, b, statistic=None, n_perm=10000, seed=0,
                            alternative="two-sided"):
    """Paired permutation test: under H0 the two conditions are exchangeable
    WITHIN each unit, so each pair (a_i, b_i) is independently swapped.

    a, b     : paired 1d arrays (same length; NaN pairs dropped).
    statistic: callable(a, b) -> scalar; default mean(a) - mean(b).
    Returns (observed, p).
    """
    a = np.asarray(a, float); b = np.asarray(b, float)
    m = ~(np.isnan(a) | np.isnan(b))
    a, b = a[m], b[m]
    n = len(a)
    if statistic is None:
        statistic = lambda x, y: float(np.mean(x) - np.mean(y))
    obs = statistic(a, b)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        swap = rng.random(n) < 0.5
        aa = np.where(swap, b, a)
        bb = np.where(swap, a, b)
        s = statistic(aa, bb)
        if alternative == "two-sided":
            count += abs(s) >= abs(obs) - 1e-12
        elif alternative == "greater":
            count += s >= obs - 1e-12
        else:
            count += s <= obs + 1e-12
    return float(obs), (count + 1) / (n_perm + 1)


def bootstrap_ci(values, statistic=np.mean, n_boot=2000, seed=0, alpha=0.05):
    """Percentile bootstrap CI of a statistic over 1d values. Returns
    (estimate, lo, hi)."""
    v = np.asarray(values, float); v = v[~np.isnan(v)]
    if v.size == 0:
        return np.nan, np.nan, np.nan
    est = float(statistic(v))
    rng = np.random.default_rng(seed)
    boots = np.array([statistic(rng.choice(v, v.size, replace=True))
                      for _ in range(n_boot)])
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return est, float(lo), float(hi)


def cluster_bootstrap_ci(values, clusters, statistic=np.mean, n_boot=2000,
                         seed=0, alpha=0.05):
    """Cluster (block) bootstrap CI: resample whole CLUSTERS with replacement,
    not individual rows. Use when rows are correlated within a cluster
    (poses within a complex / (frame,ligand) pair) so a naive bootstrap would
    understate the CI. Returns (estimate, lo, hi).

    values   : 1d array of the quantity.
    clusters : same-length array of cluster ids (complex / pair).
    statistic: callable(values_subset) -> scalar (default mean).
    """
    values = np.asarray(values, float)
    clusters = np.asarray(clusters)
    m = ~np.isnan(values)
    values, clusters = values[m], clusters[m]
    uniq = np.unique(clusters)
    if uniq.size == 0:
        return np.nan, np.nan, np.nan
    by = {c: values[clusters == c] for c in uniq}
    est = float(statistic(values))
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.choice(uniq, uniq.size, replace=True)
        sample = np.concatenate([by[c] for c in pick])
        boots[i] = statistic(sample)
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return est, float(lo), float(hi)


def consensus_perm_test(R, n_perm=10000, seed=42):
    """Independence-null permutation test for co-occurrence of per-unit binary
    successes across k conditions: do the conditions succeed on the SAME units
    more than if their per-unit successes were independent? Each column is
    shuffled across units independently (marginals preserved). R is (n x k) 0/1.

    Returns dict{observed, expected, p_more_consensus} over #-conditions-hit bins.
    """
    R = np.asarray(R, int)
    n, k = R.shape
    rng = np.random.default_rng(seed)

    def dist(mat):
        s = mat.sum(axis=1)
        return np.array([np.sum(s == c) for c in range(k + 1)], float)

    obs = dist(R)
    cols = [R[:, j].copy() for j in range(k)]
    null = np.zeros((n_perm, k + 1))
    for t in range(n_perm):
        perm = np.column_stack([rng.permutation(c) for c in cols])
        null[t] = dist(perm)
    expected = null.mean(axis=0)
    # one-sided p that the all-k-agree bin is higher than chance
    p_all = (np.sum(null[:, k] >= obs[k]) + 1) / (n_perm + 1)
    # two-sided p (deviation from expected in EITHER direction): unlike the one-
    # sided version it can register complementarity (fewer joint successes than
    # chance), not only redundancy/inflation.
    dev = np.abs(null[:, k] - expected[k])
    p_all_two = (np.sum(dev >= abs(obs[k] - expected[k])) + 1) / (n_perm + 1)
    return {"observed": obs.tolist(), "expected": expected.tolist(),
            "p_all_agree": float(p_all), "p_all_agree_two_sided": float(p_all_two),
            "k": k, "n": n}


# --------------------------------------------------------------------------- #
# self-test
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    import numpy as _np
    rng = _np.random.default_rng(0)
    fails = []

    def check(name, cond):
        print(f"  {'ok ' if cond else 'FAIL'}  {name}")
        if not cond:
            fails.append(name)

    print("stats_utils self-test")

    # Holm vs known
    hp = holm([0.01, 0.04, 0.03])
    check("holm monotone & <=1", hp.max() <= 1 and hp[0] <= hp[2])

    # BH matches manual on a simple set
    q = bh_fdr([0.001, 0.5, 0.02, 0.7])
    check("bh_fdr finite in [0,1]", _np.all((q >= 0) & (q <= 1)))

    # Wilson CI brackets the point estimate
    lo, hi = wilson_ci(8, 10)
    check("wilson brackets 0.8", lo < 0.8 < hi)

    # McNemar: strong asymmetry is significant
    a = _np.array([1] * 20 + [0] * 20)
    b = _np.array([0] * 20 + [0] * 20)
    _, _, pmn = mcnemar_exact(a, b)
    check("mcnemar detects asymmetry", pmn < 0.001)

    # Cochran's Q reduces to something sane for k=2 (matches McNemar direction)
    X = _np.column_stack([a, b])
    Q, pQ, dfQ = cochran_q(X)
    check("cochran_q k=2 significant", pQ < 0.001 and dfQ == 1)

    # paired_proportions bundle
    pp = paired_proportions({"t1": a, "t2": b})
    check("paired_proportions omnibus sig", pp["omnibus"]["p"] < 0.001)
    check("paired_proportions pairwise holm present", "p_holm" in pp["pairwise"][0])

    # Friedman + Kendall W on separated conditions
    x1 = rng.normal(0, 1, 40); x2 = x1 + 1.5; x3 = x1 + 3.0
    chi2, pf, W, nF, dfF = friedman_kendall_w(x1, x2, x3)
    check("friedman separates", pf < 1e-4 and 0 <= W <= 1)

    pc = paired_continuous({"a": x1, "b": x2, "c": x3})
    check("paired_continuous omnibus sig", pc["omnibus"]["p"] < 1e-4)
    check("paired_continuous kendall_w in [0,1]", 0 <= pc["omnibus"]["kendall_w"] <= 1)

    # Wilcoxon rank-biserial sign
    rb, pw, npair = wilcoxon_rankbiserial(x2, x1)
    check("wilcoxon rank-biserial positive", rb > 0 and pw < 1e-4)

    # Mann-Whitney + Cliff's delta on separated groups
    g1 = rng.normal(0, 1, 60); g2 = rng.normal(2, 1, 60)
    mw = mannwhitney_cliffs(g2, g1)
    check("mannwhitney sig & delta>0", mw["p"] < 1e-6 and mw["cliffs_delta"] > 0.5)
    check("cliffs CI excludes 0", mw["delta_ci"][0] > 0)

    # Kruskal + Dunn
    kd = kruskal_dunn([rng.normal(0, 1, 30), rng.normal(1, 1, 30), rng.normal(3, 1, 30)],
                      ["lo", "mid", "hi"])
    check("kruskal sig & eps2 in [0,1]", kd["p"] < 1e-4 and 0 <= kd["epsilon_sq"] <= 1)

    # Jonckheere trend on ordered groups
    JT, zJ, pJ = jonckheere([rng.normal(0, 1, 20), rng.normal(1, 1, 20), rng.normal(2, 1, 20)])
    check("jonckheere detects trend", pJ < 1e-3)

    # Cochran-Armitage trend
    zc, pca, sgn = cochran_armitage([5, 10, 20], [40, 40, 40])
    check("cochran-armitage trend sig & positive", pca < 0.01 and sgn > 0)

    # chi2 / G-test of independence on an associated table
    table = _np.array([[40, 10], [10, 40]])
    ci = chi2_independence(table); gi = gtest_independence(table)
    check("chi2 independence sig", ci["p"] < 1e-4 and ci["cramers_v"] > 0.4)
    check("gtest independence sig", gi["p"] < 1e-4)
    check("adjusted residuals shape", ci["residuals"].shape == (2, 2))

    # spearman surfaces a p-value
    xx = rng.normal(0, 1, 50); yy = xx + rng.normal(0, 0.3, 50)
    rho, prho, nn = spearman(xx, yy)
    check("spearman sig positive", rho > 0.7 and prho < 1e-6)

    # paired permutation test
    obs, pperm = permutation_test_paired(x2, x1, n_perm=2000, seed=1)
    check("paired permutation sig", pperm < 0.01 and obs > 0)

    # cluster bootstrap CI wider than naive on clustered data
    clusters = _np.repeat(_np.arange(20), 10)
    vals = _np.repeat(rng.normal(0, 1, 20), 10) + rng.normal(0, 0.01, 200)
    _, clo, chi = cluster_bootstrap_ci(vals, clusters)
    _, nlo, nhi = bootstrap_ci(vals)
    check("cluster bootstrap wider than naive", (chi - clo) > (nhi - nlo))

    # consensus permutation test returns bins
    R = (rng.random((50, 3)) < 0.5).astype(int)
    cpt = consensus_perm_test(R, n_perm=1000)
    check("consensus perm returns k+1 bins", len(cpt["observed"]) == 4)

    print(f"\n{'ALL PASS' if not fails else 'FAILURES: ' + ', '.join(fails)}")
    raise SystemExit(1 if fails else 0)
