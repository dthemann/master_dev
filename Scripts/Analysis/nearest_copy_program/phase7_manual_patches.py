#!/usr/bin/env python
"""Phase 7 manual transcriptions: the 19 harness checks whose expected value is a yaml LIST (no scalar token to
trace). Values come from the harness copy's own recomputation on the _nearest tree (phase7_manual_ledger.csv,
2026-09-08). Exact-match replacements on thesis_expected_values_nearest.yaml only."""
from pathlib import Path
HERE = Path(__file__).resolve().parent
p = HERE / "thesis_expected_values_nearest.yaml"; s = p.read_text()
REPS = [
 # table_3 reach_pct [rate, wilson_lo, wilson_hi] per depth 1/5/10/15, any-copy rule
 ("    AutoDock:  [[59, 53, 64], [79, 74, 83], [81, 76, 85], [83, 79, 87]]",
  "    AutoDock:  [[78, 73, 82], [86, 81, 89], [88, 84, 91], [90, 86, 93]]"),
 ("    DiffDock:  [[54, 49, 60], [76, 70, 80], [79, 74, 83], [81, 76, 85]]",
  "    DiffDock:  [[73, 68, 78], [83, 78, 86], [85, 80, 88], [86, 82, 90]]"),
 ("    EquiBind:  [[43, 37, 49], [46, 40, 51], [47, 41, 52], [47, 42, 53]]",
  "    EquiBind:  [[46, 40, 51], [49, 43, 54], [50, 45, 56], [51, 46, 57]]"),
 ("  all_three_observed: [89, 114, 121, 125]", "  all_three_observed: [103, 125, 132, 137]"),
 ("  all_three_expected: [42, 82, 90, 96]", "  all_three_expected: [79, 105, 113, 120]"),
 # table_22 DiffDock (raw): accurate-but-invalid share at k = 1 and 5 (k = 10 and 15 unchanged)
 ('    "DiffDock (raw)":        [[15.5, 11.9, 20.0], [19.1, 15.1, 23.9], [19.5, 15.4, 24.3], [18.8, 14.8, 23.6]]',
  '    "DiffDock (raw)":        [[21.5, 17.2, 26.4], [19.8, 15.7, 24.7], [19.5, 15.4, 24.3], [18.8, 14.8, 23.6]]'),
 # appendix H.9: the 10 pp TOST equivalence at rank-1 no longer holds under nearest-copy scoring (plan v2 predicted this)
 ("      equivalent_within_10pp: [true, false, false, false]", "      equivalent_within_10pp: [false, false, false, false]"),
]
for old, new in REPS:
    assert s.count(old) == 1, old
    s = s.replace(old, new)
p.write_text(s); print(f"{len(REPS)} manual replacements applied to {p.name}")
