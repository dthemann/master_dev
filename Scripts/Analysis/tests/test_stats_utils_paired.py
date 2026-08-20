"""Tests for the paired-proportion helpers added for the AutoDock-vs-DiffDock contrast.

The Newcombe interval is the thesis's most-quoted inferential number, so it is
checked against a published worked example as well as against the thesis figures.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from stats_utils import (  # noqa: E402
    newcombe_paired_diff_ci, mcnemar_power, tost_paired_proportions, mcnemar_exact,
)


def _vectors(n11, n10, n01, n00):
    """Build paired 0/1 vectors with the given 2x2 cell counts."""
    a = [1] * n11 + [1] * n10 + [0] * n01 + [0] * n00
    b = [1] * n11 + [0] * n10 + [1] * n01 + [0] * n00
    return np.array(a), np.array(b)


# --------------------------------------------------------------------------- #
# Newcombe method 10
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize('cells', [
    (40, 49, 63, 151),   # the thesis rank-1 contrast
    (94, 57, 73, 79),    # the thesis top-15 contrast
    (36, 12, 2, 0),      # one empty cell
    (20, 30, 30, 20),    # exactly null
])
def test_newcombe_agrees_with_paired_bootstrap(cells):
    """The closed-form interval must track a nonparametric paired bootstrap that
    resamples whole complexes. Agreement to ~1 percentage point on both bounds."""
    a, b = _vectors(*cells)
    r = newcombe_paired_diff_ci(a, b)
    d = b.astype(int) - a.astype(int)
    rng = np.random.default_rng(0)
    idx = rng.integers(0, len(d), size=(4000, len(d)))
    boots = d[idx].mean(axis=1)
    lo_b, hi_b = np.percentile(boots, [2.5, 97.5])
    assert r['lo'] == pytest.approx(lo_b, abs=0.012)
    assert r['hi'] == pytest.approx(hi_b, abs=0.012)


def test_newcombe_nominal_coverage():
    """A 95% interval must cover the true difference about 95% of the time.

    Simulates paired binary outcomes with a known marginal difference and a
    realistic positive correlation between the two arms.
    """
    rng = np.random.default_rng(7)
    n, trials = 303, 600
    p_a, delta, rho = 0.30, 0.05, 0.5      # true p_b - p_a = +0.05
    covered = 0
    for _ in range(trials):
        shared = rng.random(n) < rho       # complexes where both arms share an outcome
        base = rng.random(n) < p_a
        a = np.where(shared, base, rng.random(n) < p_a)
        b = np.where(shared, base, rng.random(n) < p_a + delta)
        true = b.mean() - a.mean()         # realised difference in this sample
        r = newcombe_paired_diff_ci(a, b)
        covered += (r['lo'] <= true <= r['hi'])
    assert 0.90 <= covered / trials <= 1.0, f'coverage {covered / trials:.3f}'


def test_newcombe_cells_and_marginals():
    a, b = _vectors(40, 49, 63, 151)
    r = newcombe_paired_diff_ci(a, b)
    assert (r['n11'], r['n10'], r['n01'], r['n00']) == (40, 49, 63, 151)
    assert r['n'] == 303
    assert r['p_a'] == pytest.approx(89 / 303)
    assert r['p_b'] == pytest.approx(103 / 303)
    assert r['n_discordant'] == 112


def test_newcombe_reproduces_thesis_rank1_interval():
    """The rank-1 AutoDock+gnina vs DiffDock+smina contrast quoted in the thesis:
    +4.6 percentage points, 95% interval -2.2 to +11.4."""
    a, b = _vectors(40, 49, 63, 151)
    r = newcombe_paired_diff_ci(a, b)
    assert 100 * r['diff'] == pytest.approx(4.6, abs=0.05)
    assert 100 * r['lo'] == pytest.approx(-2.2, abs=0.15)
    assert 100 * r['hi'] == pytest.approx(11.4, abs=0.15)


def test_newcombe_interval_contains_point_estimate():
    for cells in [(40, 49, 63, 151), (94, 57, 73, 79), (10, 0, 0, 90), (0, 5, 5, 0)]:
        r = newcombe_paired_diff_ci(*_vectors(*cells))
        assert r['lo'] <= r['diff'] <= r['hi']
        assert -1.0 <= r['lo'] and r['hi'] <= 1.0


def test_newcombe_zero_difference_straddles_zero():
    r = newcombe_paired_diff_ci(*_vectors(20, 30, 30, 20))
    assert r['diff'] == pytest.approx(0.0)
    assert r['lo'] < 0 < r['hi']


def test_newcombe_empty_input():
    r = newcombe_paired_diff_ci([], [])
    assert r['n'] == 0 and np.isnan(r['diff'])


def test_newcombe_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        newcombe_paired_diff_ci([1, 0, 1], [1, 0])


# --------------------------------------------------------------------------- #
# power / MDE
# --------------------------------------------------------------------------- #
def test_power_uses_discordant_pairs_not_n():
    """Two designs with the same n and the same difference but different
    discordance must not have the same power."""
    tight = mcnemar_power(*_vectors(140, 5, 19, 139))   # few discordant
    loose = mcnemar_power(*_vectors(40, 49, 63, 151))   # many discordant
    assert tight['n'] == loose['n'] == 303
    assert tight['n_discordant'] < loose['n_discordant']
    assert tight['mde'] < loose['mde']


def test_power_rank1_contrast_is_low():
    """The rank-1 contrast is badly underpowered: 112 discordant pairs."""
    r = mcnemar_power(*_vectors(40, 49, 63, 151))
    assert r['n_discordant'] == 112
    assert 0.10 < r['observed_power'] < 0.40
    assert 0.07 < r['mde'] < 0.13
    assert r['n_needed'] > 1000


def test_power_of_zero_effect_is_alpha():
    r = mcnemar_power(*_vectors(20, 30, 30, 20), alpha=0.05)
    assert r['observed_power'] == pytest.approx(0.05)


def test_power_no_discordant_pairs():
    r = mcnemar_power(*_vectors(50, 0, 0, 50))
    assert r['n_discordant'] == 0 and np.isnan(r['observed_power'])


# --------------------------------------------------------------------------- #
# TOST equivalence
# --------------------------------------------------------------------------- #
def test_tost_not_equivalent_at_tight_margins_but_is_at_wide():
    a, b = _vectors(40, 49, 63, 151)
    assert not tost_paired_proportions(a, b, margin=0.05)['equivalent']
    assert not tost_paired_proportions(a, b, margin=0.10)['equivalent']
    assert tost_paired_proportions(a, b, margin=0.15)['equivalent']


def test_tost_margin_is_monotone():
    a, b = _vectors(40, 49, 63, 151)
    prev = False
    for m in (0.02, 0.05, 0.10, 0.12, 0.15, 0.20, 0.30):
        eq = tost_paired_proportions(a, b, margin=m)['equivalent']
        assert not (prev and not eq), 'equivalence must not switch back off'
        prev = eq
    assert prev


def test_tost_interval_is_narrower_than_95_percent():
    """TOST uses a 90% interval at alpha=0.05, so it must sit inside the 95%."""
    a, b = _vectors(40, 49, 63, 151)
    ci95 = newcombe_paired_diff_ci(a, b)
    t = tost_paired_proportions(a, b, margin=0.10, alpha=0.05)
    assert t['conf_level'] == pytest.approx(0.90)
    assert t['lo'] > ci95['lo'] and t['hi'] < ci95['hi']


# --------------------------------------------------------------------------- #
# consistency with the existing McNemar helper
# --------------------------------------------------------------------------- #
def test_consistent_with_mcnemar_exact():
    a, b = _vectors(40, 49, 63, 151)
    n10, n01, p = mcnemar_exact(a, b)
    r = newcombe_paired_diff_ci(a, b)
    assert (n10, n01) == (r['n10'], r['n01'])
    # non-significant test <-> interval containing zero
    assert p > 0.05 and r['lo'] < 0 < r['hi']
