"""Tests for the shared single-point method-exclusion filter.

The behaviour that matters here is not "does it drop rows" but the two ways the
suite has previously got this wrong: a bare prefix swallowing sibling arms, and a
selector that matches nothing yet reports success. Both are asserted explicitly.
"""
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import method_filter as mf  # noqa: E402


BENCHMARK_METHODS = [
    "autodock", "autodock_gnina", "autodock_gnina_refinement",
    "autodock_mgltools", "autodock_mgltools_gnina",
    "autodock_mgltools_exh18", "autodock_mgltools_exh64",
    "autodock_mgltools_exh64_gnina", "autodock_mgltools_exh92",
    "autodock_mgltools_exh128", "autodock_mgltools_exh128_gnina",
    "autodock_vinardo", "autodock_vinardo_gnina",
    "autodock_vinardo_gnina_refinement",
    "diffdock", "diffdock_smina", "diffdock_gnina",
    "equibind_unguided_raw", "equibind_unguided_gnina",
    "unidock2",
]


def _frame(methods=None, column="method", per_method=3):
    methods = methods or BENCHMARK_METHODS
    return pd.DataFrame({column: [m for m in methods for _ in range(per_method)],
                         "rmsd": 1.0})


class _Args:
    def __init__(self, preset=None, methods=None, strict=False):
        self.exclude_preset = preset
        self.exclude_methods = methods
        self.exclude_strict = strict


# --------------------------------------------------------------------------- #
# matching semantics
# --------------------------------------------------------------------------- #
def test_bare_key_is_exact_never_a_prefix():
    """'autodock' must not swallow autodock_gnina / autodock_mgltools*."""
    matched, unmatched = mf.match_methods(BENCHMARK_METHODS, ["autodock"])
    assert matched == ["autodock"]
    assert unmatched == []


def test_glob_is_opt_in_and_bounded_to_its_family():
    matched, _ = mf.match_methods(BENCHMARK_METHODS, ["autodock_vinardo*"])
    assert set(matched) == {"autodock_vinardo", "autodock_vinardo_gnina",
                            "autodock_vinardo_gnina_refinement"}


def test_unmatched_patterns_are_reported_not_swallowed():
    matched, unmatched = mf.match_methods(BENCHMARK_METHODS, ["autodock", "nope_*"])
    assert matched == ["autodock"]
    assert unmatched == ["nope_*"]


def test_matches_are_deduplicated_across_overlapping_patterns():
    matched, _ = mf.match_methods(BENCHMARK_METHODS,
                                  ["autodock_vinardo", "autodock_vinardo*"])
    assert len(matched) == len(set(matched))


# --------------------------------------------------------------------------- #
# the meeko preset
# --------------------------------------------------------------------------- #
def test_meeko_preset_drops_exactly_the_meeko_ligand_arms():
    df = _frame()
    out = mf.apply_method_filter(df, "method", _Args(preset="meeko"))
    kept = set(out["method"])
    assert kept == {
        "autodock_mgltools", "autodock_mgltools_gnina",
        "autodock_mgltools_exh18", "autodock_mgltools_exh64",
        "autodock_mgltools_exh64_gnina", "autodock_mgltools_exh92",
        "autodock_mgltools_exh128", "autodock_mgltools_exh128_gnina",
        "diffdock", "diffdock_smina", "diffdock_gnina",
        "equibind_unguided_raw", "equibind_unguided_gnina", "unidock2",
    }


def test_meeko_preset_spares_every_adfrsuite_arm():
    out = mf.apply_method_filter(_frame(), "method", _Args(preset="meeko"))
    kept = set(out["method"])
    assert all(m in kept for m in BENCHMARK_METHODS if "mgltools" in m)


def test_meeko_preset_leaves_the_learned_tools_untouched():
    out = mf.apply_method_filter(_frame(), "method", _Args(preset="meeko"))
    kept = set(out["method"])
    assert {"diffdock", "diffdock_smina", "diffdock_gnina",
            "equibind_unguided_raw", "equibind_unguided_gnina",
            "unidock2"} <= kept


def test_preset_works_on_a_presplit_frame_of_base_keys_only():
    """Frames read before a script's optimizer split carry base tree keys.

    The preset must still remove the whole Meeko trees there, which is what makes
    one preset usable in both posebusters_pose_comparison (post-split, variant
    keys) and posebusters_validity_report (pre-split, base keys).
    """
    base = ["autodock", "autodock_vinardo", "autodock_mgltools",
            "autodock_mgltools_exh128", "diffdock", "equibind_guided"]
    out = mf.apply_patterns(_frame(base, column="docking_method"),
                            "docking_method", mf.PRESETS["meeko"][0])
    assert set(out["docking_method"]) == {
        "autodock_mgltools", "autodock_mgltools_exh128", "diffdock",
        "equibind_guided"}


# --------------------------------------------------------------------------- #
# guard rails
# --------------------------------------------------------------------------- #
def test_no_flags_is_a_byte_identical_no_op():
    df = _frame()
    out = mf.apply_method_filter(df, "method", _Args())
    assert out is df


def test_strict_makes_an_unmatched_pattern_fatal():
    with pytest.raises(SystemExit):
        mf.apply_method_filter(_frame(), "method",
                               _Args(methods="not_a_method", strict=True))


def test_emptying_the_frame_is_fatal_even_without_strict():
    with pytest.raises(SystemExit):
        mf.apply_method_filter(_frame(["autodock"]), "method", _Args(preset="meeko"))


def test_missing_column_is_fatal():
    with pytest.raises(SystemExit):
        mf.apply_method_filter(_frame(), "docking_method", _Args(preset="meeko"))


def test_preset_and_explicit_methods_union():
    out = mf.apply_method_filter(
        _frame(), "method", _Args(preset="meeko", methods="unidock2,diffdock_gnina"))
    kept = set(out["method"])
    assert "unidock2" not in kept and "diffdock_gnina" not in kept
    assert "diffdock" in kept


def test_provenance_sidecar_records_dropped_and_kept(tmp_path):
    import json
    mf.apply_method_filter(_frame(), "method", _Args(preset="meeko"),
                           out_dir=tmp_path)
    payload = json.loads((tmp_path / "method_filter.json").read_text())
    assert payload["preset"] == "meeko"
    assert "autodock" in payload["dropped_methods"]
    assert "autodock_mgltools" in payload["kept_methods"]
    assert payload["preset_provenance"]


# --------------------------------------------------------------------------- #
# ordering: the filter must run before a key-folding step
# --------------------------------------------------------------------------- #
def test_filtering_after_the_autodock_fold_would_destroy_the_adfrsuite_arms():
    """Why load_and_score filters BEFORE _apply_autodock_split, not after.

    _autodock_scoring_base rewrites every autodock_mgltools* key onto plain
    'autodock', so once the split has run the ADFRsuite arms are indistinguishable
    from the Meeko ones and the meeko preset removes both. Asserting the wrong
    order actively fails keeps that from being "simplified" back later.
    """
    V = pytest.importorskip("posebusters_validity_report")
    raw = pd.DataFrame({
        "docking_method": ["autodock", "autodock_mgltools",
                           "autodock_mgltools_exh128", "diffdock"],
        "optimizer": ["original", "original", "gnina", "original"],
    })

    # WRONG order — fold first, then filter.
    folded = V._apply_autodock_split(raw.copy())
    assert set(folded["docking_method"]) == {"autodock", "autodock_gnina", "diffdock"}, \
        "the fold is what makes ordering matter; if this changes, revisit the wiring"
    late = mf.apply_patterns(folded, "docking_method", mf.PRESETS["meeko"][0])
    assert set(late["docking_method"]) == {"diffdock"}, \
        "filtering after the fold silently takes the ADFRsuite arms with it"

    # RIGHT order — filter first, as load_and_score now does.
    early = mf.apply_patterns(raw.copy(), "docking_method", mf.PRESETS["meeko"][0])
    early = V._apply_autodock_split(early)
    assert set(early["docking_method"]) == {"autodock", "autodock_gnina", "diffdock"}
    assert len(early) == 3, "both ADFRsuite arms plus DiffDock survive"


def test_row_counts_drop_by_exactly_the_excluded_methods():
    df = _frame(per_method=7)
    out = mf.apply_method_filter(df, "method", _Args(preset="meeko"))
    dropped_methods = 6           # autodock, _gnina, _gnina_refinement + 3 vinardo
    assert len(df) - len(out) == dropped_methods * 7
