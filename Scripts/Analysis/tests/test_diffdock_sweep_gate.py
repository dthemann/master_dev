"""Pin the three-stage gate used by DiffDock_Parameter_Sweep.ipynb.

The gate is a conjunction over a SINGLE pose: PoseBusters-valid AND rmsd <= 2 A
(in-place, no superposition) AND bestfit_rmsd <= 1 A with rmsd < 1000 (the
exploded-pose guard). The two ways it has gone wrong before are asking each
stage of the pool separately and then AND-ing the answers, which counts
complexes no single pose satisfies, and dropping the exploded-pose guard, which
lets DiffDock coordinate explosions pass Form on a near-perfect internal
conformation while sitting thousands of angstrom off the receptor.

These tests lift gate_flags() and endpoints() out of the notebook itself rather
than re-implementing them, so they fail if the notebook's definitions drift.

Run with an interpreter that has pytest, e.g.
    /home/manndo/anaconda3/envs/sgai/bin/python -m pytest \
        Scripts/Analysis/tests/test_diffdock_sweep_gate.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[3]
NOTEBOOK = REPO / "DiffDock_Parameter_Sweep.ipynb"

# The cell that defines the gate. Located by content, not by index, so inserting
# a cell above it does not silently test the wrong code.
GATE_MARKER = "def gate_flags("


def _load_gate():
    nb = json.loads(NOTEBOOK.read_text())
    sources = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]
    cells = [s for s in sources if GATE_MARKER in s]
    assert len(cells) == 1, f"expected exactly one gate cell, found {len(cells)}"
    # Everything up to the first use of the definitions — the definitions alone,
    # without the notebook state they would otherwise need.
    head = cells[0].split("EP = endpoints(")[0]
    ns: dict = {"pd": pd, "np": np}
    exec(compile(head, "gate_cell", "exec"), ns)
    return ns


@pytest.fixture(scope="module")
def gate():
    return _load_gate()


def _pose(rank, rmsd, kabsch, valid):
    return {"arm": "a", "protein": "P", "variant": "smina", "rank": rank,
            "rmsd": rmsd, "bestfit_rmsd": kabsch, "pb_valid": valid}


def test_thresholds_are_the_documented_ones(gate):
    assert gate["RMSD_OK"] == 2.0
    assert gate["KABSCH_OK"] == 1.0
    assert gate["EXPLODED"] == 1000.0
    assert gate["DEPTH_CAP"] == 30


def test_three_poses_each_passing_one_stage_do_not_count(gate):
    """The failure this recompute exists to prevent."""
    df = pd.DataFrame([
        _pose(1, 50.0, 5.0, True),    # valid only
        _pose(2, 1.0, 5.0, False),    # near-native only
        _pose(3, 50.0, 0.5, False),   # form only
    ])
    flags = gate["gate_flags"](df)
    assert flags.s_valid.any() and flags.s_near.any() and flags.s_form.any()
    assert not flags.s_triple.any(), "no single pose satisfies all three"

    ep = gate["endpoints"](df, "smina").iloc[0]
    assert not ep.gate_top30
    assert ep.leak == "validity_form"


def test_one_pose_passing_all_three_counts(gate):
    df = pd.DataFrame([
        _pose(1, 50.0, 5.0, True),
        _pose(2, 1.0, 0.5, True),
    ])
    ep = gate["endpoints"](df, "smina").iloc[0]
    assert ep.gate_top30
    assert not ep.gate_top1, "the qualifying pose is at rank 2, not rank 1"
    assert ep.leak == "ranking"


def test_exploded_pose_fails_form_despite_perfect_conformation(gate):
    """Kabsch 0.07 A but 2.6e6 A from the crystal — Form is the only stage it
    could otherwise pass, and the rmsd < 1000 guard is what stops it."""
    df = pd.DataFrame([_pose(1, 2.6e6, 0.07, False)])
    flags = gate["gate_flags"](df)
    assert not flags.s_form.any()
    assert not flags.s_triple.any()
    assert gate["endpoints"](df, "smina").iloc[0].leak == "sampling"


def test_a_thirty_pose_arm_is_left_in_rank_order(gate):
    """No subsampling when the pool is already at the cap: depth is the
    confidence rank, and a rank-1 hit counts at top-1."""
    df = pd.DataFrame([_pose(1, 0.9, 0.4, True)]
                      + [_pose(r, 50.0, 5.0, False) for r in range(2, 31)])
    ep = gate["endpoints"](df, "smina").iloc[0]
    assert ep.n_poses == 30
    assert ep.gate_top1 and ep.leak == "passes"


def test_over_sampled_pool_is_reduced_to_the_cap(gate):
    """A 60-pose arm is scored on 30, so it never competes on a bigger pool."""
    df = pd.DataFrame([_pose(r, 50.0, 5.0, False) for r in range(1, 61)])
    assert gate["endpoints"](df, "smina").iloc[0].n_poses == 30


def test_over_sampled_pool_is_drawn_at_random_not_confidence_best(gate):
    """The correction that matters. Taking rank <= 30 of a 60-pose pool would be
    the confidence-best 30, a selection a real 30-sample arm never gets. Put the
    only qualifying pose deep in the pool: under a rank cap it could never count,
    under a random draw it counts about half the time."""
    # endpoints() reseeds from a fixed seed on every call, which is what makes a
    # run reproducible. The draw therefore varies between COMPLEXES inside one
    # call, not between calls, so the distribution has to be measured that way.
    rows = []
    for c in range(40):
        rows += [_pose(r, 50.0, 5.0, False) | {"protein": f"P{c}"} for r in range(1, 60)]
        rows.append(_pose(60, 0.9, 0.4, True) | {"protein": f"P{c}"})
    ep = gate["endpoints"](pd.DataFrame(rows), "smina")

    assert (ep.n_poses == 30).all(), "every 60-pose complex must be reduced to 30"
    counted = int(ep.gate_top30.sum())
    assert 5 < counted < 35, (
        f"the deep pose was drawn for {counted}/40 complexes; 0 would mean a rank "
        f"cap is still applied and 40 would mean the pool was not reduced. A random "
        f"30 of 60 should include it about half the time.")


def test_leak_strata_are_exhaustive_and_ordered(gate):
    """Every complex lands in exactly one stratum, and a complex that fails
    sampling is never attributed to ranking."""
    frames = {
        "sampling": pd.DataFrame([_pose(1, 40.0, 5.0, False)]),
        "validity_form": pd.DataFrame([_pose(1, 1.0, 5.0, False)]),
        "ranking": pd.DataFrame([_pose(1, 40.0, 5.0, False), _pose(2, 1.0, 0.4, True)]),
        "passes": pd.DataFrame([_pose(1, 1.0, 0.4, True)]),
    }
    for expected, df in frames.items():
        assert gate["endpoints"](df, "smina").iloc[0].leak == expected


def test_variant_column_selects_the_geometry(gate):
    """The notebook keys on `variant`, not `optimizer`; a mismatched geometry
    must contribute nothing rather than silently leaking in."""
    df = pd.DataFrame([
        _pose(1, 1.0, 0.4, True) | {"variant": "raw"},
        _pose(1, 40.0, 5.0, False) | {"variant": "smina"},
    ])
    ep = gate["endpoints"](df, "smina").iloc[0]
    assert not ep.gate_top30, "the raw-geometry pose must not count for the smina arm"
