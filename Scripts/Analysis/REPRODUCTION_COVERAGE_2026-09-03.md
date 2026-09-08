# Reproduction coverage of the short build

What `Thesis_Reproduction.ipynb` asserts about `thesis_latex/Thesis_short.tex`, what it
did not assert before 2026-09-03, and the defects the gap was hiding.

Written 2026-09-03. The companion document `FINDINGS_2026-09-03.md` covers the
Table 21 baseline defect, which was found the same night from the other direction.

---

## 1. Why this was worth doing

The short build prints 27 tables and 39 figures. Before this pass the harness read
8 of those tables and md5-checked 24 of the figure assets. Everything else was
unchecked, and an unchecked float is where the pipeline drifts: Table 21 had been
tested against the wrong AutoDock arm since the first commit of its test list, and
nothing noticed because no assertion read that table.

This pass closes the gap on every float that a script can reach, and records the
ones it cannot together with the reason.

## 2. Coverage

| | before | after |
|---|---|---|
| assertions | 459 | 1,282 (1,362 after the Table 12 fix later the same day) |
| tables read | 8 of 27 | 23 of 27 |
| figure assets tied to a pipeline file | 24 of 39 | 37 of 39, 2 declared |
| completeness assertion over printed figures | none | present |

Tables now asserted: 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 14, 16, 17, 18, 19, 20,
21, 22, 24, 25, 26, 27.

Not asserted, with the reason:

- **Table 7**, physics-based and AI-based docking studies. Authored from the cited
  literature. No artifact in this repository produces it.
- **Table 13**, the PoseBusters check list. A description of the battery, not a
  measurement. Its one countable claim, that twenty-two checks were applied, is
  asserted inside Table 24 where the three groups must partition exactly that many.
- **Table 15**, reported potencies and proposed sites of the Orai1 modulators.
  Authored from the cited pharmacology.
- **Table 23**, the mapping from design to test and effect size. Prose in a grid;
  it prints no numbers.

## 3. What the coverage gap was hiding

### 3.1 Table 12 does not reproduce, and neither does the Figure 26 caption

The AutoDock column of Table 12 cannot be produced from anything in the repository.

| cell | thesis | canonical run, 2026-09-01 | older run, 2026-08-30 |
|---|---|---|---|
| AutoDock hydrogen bond, benchmark | 3.85 | 4.15 | 3.97 |
| AutoDock hydrogen bond, experimental | 2.42 | 2.53 | 2.53 |
| AutoDock Cliff's delta | +0.42 | +0.45 | +0.40 |
| DiffDock hydrogen bond, benchmark | 2.47 | 2.47 | — |
| DiffDock hydrogen bond, experimental | 1.31 | 1.31 | — |

Eighteen of the table's thirty-five cells differ and every one of them is in the
AutoDock column. All fourteen DiffDock cells reproduce exactly. That pattern is the
signature of a superseded AutoDock arm rather than of a different aggregation.

Four things establish that the artifact is current and the table is not:

1. The shipped `image20` and `image22` are byte-identical to the 2026-09-01 run.
2. Re-running `orai_pandamap_interaction_compare.py` into a scratch directory
   reproduces that run's stats file exactly, modulo its own timestamp line.
3. The main body already carries the current values. It quotes Cliff's delta of
   +0.45 and +0.54 for hydrogen bonds and −0.59 and −0.58 for hydrophobic contacts,
   and says fourteen of twenty-six cells clear Benjamini-Hochberg correction, which
   is what the current run gives as 8 AutoDock plus 6 DiffDock. The table's own
   footnote says eleven. Body and table contradict each other and the body is right.
4. The Orai PandaMap trees carry only `autodock_gnina`, so no arm on disk could
   produce 3.85.

The aggregation is not the explanation. The caption says per-complex mean counts per
pose, and `_per_complex_means` in the compare script computes exactly that, so the
table and the stats file are on the same basis.

The caption of Figure 26 has the same problem. It says the residue axis is ranked
over 7,231 poses and runs from Asp110 at 30.2 per cent down to Asp114 at 16.7. The
current run pools 7,458 poses, puts Asp110 at 28.8 per cent, and ends the axis at
Pro146 rather than Asp114. The figure itself is current; only its caption is not.

**Status: FIXED 2026-09-03.** Table 12's AutoDock column and the Figure 26 caption were
regenerated onto the canonical run, and the harness now expects them to pass. The
`expected_to_fail` marker is gone and the recorded-defect count is zero.

Two corrections to what this section originally claimed, both found while fixing it:

1. **"cannot be produced from anything in the repository" was wrong.** All twenty-one
   printed AutoDock cells reproduce character-for-character from
   `obsolete/pandamap_results/orai_interaction_compare_PRE_FR0/orai_pandamap_interaction_compare_stats.txt`,
   and a second file, `obsolete/posebusters_results/orai_pbvalid_tm_share_compare_PRE_FR0/orai_ligand_level_contrasts.txt`,
   carries the same values from a different script. The search that produced this section
   looked in `pandamap_results/` but not in `obsolete/`. The table was two generations
   stale, not unattributable, which is a housekeeping defect rather than an integrity one.
   The Figure 26 caption traces to the same PRE_FR0 generation.

2. **Fixing the count alone would have left a false sentence.** The footnote also claims
   every clearing cell appears in the table. Canonically AutoDock clears carbon-pi and
   pi-pi stacking, neither of which was a printed row, so two rows were added. The
   fourteen clearing cells span exactly the nine types now printed.

Three further locations carried the same PRE_FR0 substitution and were not listed here:
the total-contacts effect sizes in the appendix, the AutoDock residue sentence after the
hot-spot figure, and the fifteen-versus-thirteen worked example. The last two asserted a
direction and a significance verdict that the canonical data contradict, so they would
have survived any purely numeric assertion.

### 3.2 The ligand-contrasts stage read superseded inputs

`orai_ligand_level_contrasts.py` takes nine input paths. The notebook stage passed
one of them. The other eight defaulted to the plain `orai_jku`, `orai_benchmark` and
`orai_interaction_compare` trees, which are one generation behind the `_matched`
ones the thesis reports.

Re-running the stage as it was registered produced materially different numbers:
EquiBind's cluster-quality contrasts moved from six frame-ligand units to four,
Compactness for EquiBind flipped from −0.09 ns to +0.39, and adjusted p-values moved
throughout the file. Nothing errored, and the output landed in the canonical
directory, which is what made it invisible.

This is the same defect already documented for the `orai_pandamap_compare` stage,
fixed there and left standing here.

**Status: fixed.** The stage now names all nine paths. Verified: with them the script
reproduces the shipped sidecar byte for byte apart from its generated-on line.

### 3.3 A body sentence traces to the superseded run

`body_main_short.tex:683` says that at the ligand unit no interaction-type cell
passes correction, "where the smallest adjusted values are 0.086 and 0.106".

Isolating the ligand-unit rows of the interaction-type family, per tool:

| run | AutoDock | DiffDock |
|---|---|---|
| canonical, `_matched` trees | 0.089 | 0.106 |
| the stage as previously registered, superseded trees | 0.086 | 0.106 |
| Table 12's footnote | 0.142 | 0.106 |

The printed 0.086 is the superseded run. The canonical value is 0.089. The footnote's
0.142 is a third generation again. DiffDock agrees across all three, which is why the
discrepancy is invisible unless the AutoDock arm is isolated.

**Status: FIXED 2026-09-03.** The body now prints 0.089. The conclusion the sentence draws
was unchanged either way, since neither value passes correction.

One refinement to the diagnosis above. The sentence was not wholly stale but internally
mixed. Its "fourteen of twenty-six" was already the canonical count while the 0.086 beside
it was one generation back, so a previous pass had updated the count and left the p-value.

### 3.4 Table 14 had no generator

The receptor geometry across the four frames existed only in the document. No script
wrote it and no stage produced it.

It is now measured in the harness directly from the four receptor PDBs, which are
copy-only inputs. The pore axis is the first principal component of the CA cloud and
the origin is that cloud's centroid; each diagnostic atom's perpendicular distance is
averaged over the six subunits. All twenty printed cells reproduce, and so does the
surrounding prose: 0.14 and 0.08 angstrom of spread at the gate and the filter, 0.96
at Asp114, and Fr0 sitting 3.5 to 5.2 angstrom more compact at the three aspartates
against 0.3 to 0.8 at the gate and filter.

### 3.5 A latent trap in the endpoint diagnostics

`thesis_endpoint_diagnostics.py` defaults its metrics path to
`benchmark_full_protein_vina_scoring`, which is not the canonical tree, and the two
per-pose tables differ. They agree on the AutoDock and DiffDock arms that Table 22
and Appendix H.9 use, so nothing printed is wrong today. The default would carry a
future divergence into the appendix in silence, so the spec pins the canonical path
rather than relying on it.

## 4. Figures

Thirteen shipped assets had no provenance pair. Each was located by hashing every PNG
in the repository and matching on md5 rather than by name, and every one is
byte-identical to a live pipeline output:

- `image20` and `image22`, the Orai total-contacts and residue-hotspot panels. These
  were the consequential pair: `validate_regeneration.py` re-renders and byte-compares
  them, but they were absent from the provenance layer, so the two layers disagreed
  about what was covered.
- `image1`, the docking flow chart. Author-drawn, and the regeneration check rightly
  declines to re-render it, but it has a single source file and its provenance is
  checkable even though its regeneration is not.
- `image33` to `image42`, the chemical-space chapter. A byte comparison of a
  *re-render* is excluded because matplotlib moves 1 to 2.5 per cent of pixels between
  runs. That is an argument against re-rendering them, not against checking that the
  shipped asset is the file the notebook last wrote.

Two have no source anywhere in the repository and are now recorded rather than left as
a silent hole in the count. `image12` and `image14` are one-off interactive PyMOL
renders. Both captions disclaim measurement in the document itself, so nothing numeric
depends on them, and the harness asserts that they still have no generator: if one
acquires a source, the check fails and the entry should be promoted.

A completeness assertion now reads `\includegraphics` from the two body files and
fails if the thesis prints an image that is neither paired nor declared. Without it,
adding a figure silently adds an unchecked asset while the count still reads as full
coverage.

## 5. How the values were validated

Each value went through three independent passes.

1. **Recomputation from the primary data**, by hand, before any of it was encoded.
   Tables 3, 4, 10, 14, 20, 24, 25, 26 and 27 and the Appendix H.9 block were each
   derived from the canonical CSVs or the receptor files and compared cell by cell to
   the printed table.
2. **An independent inventory.** Eighteen agents transcribed every printed cell of all
   27 tables, every number quoted in a figure caption or its prose, and 1,184 prose
   numbers, and located the producing artifact for each without seeing pass 1. The two
   passes agreed on every table's source and on the classification of all 39 figures.
3. **Adversarial mutation of the encoded assertions.** Each new block was perturbed
   value by value in a scratch copy of the spec to confirm the checker rejects the
   perturbed value. A check that still passes after its expected value is changed is
   vacuous and worse than no check, and this pass exists to find those.

## 6. Traps encoded in the spec

Recorded at the point of use so they are enforced rather than remembered.

- **Table 24's join.** The PoseBusters table's `docking_method` is the engine, not the
  variant: `autodock_mgltools_exh128` holds 17,952 rows because the raw and rescored
  arms share it. Joining on it pools them and the pose count reads double. The join is
  on `pose_file`, unique in both tables.
- **Table 27's input.** Its sidecar carries depths 1, 5 and 15. The sibling
  `filmstrip_stats__per_tool_depth.csv` has the same columns at depths 1, 3 and 5, so
  reading it supplies no top-15 row and lines depth 3 up against top-5 without erroring.
- **Table 22's derivation.** Computed from counts, not by subtracting the printed cells
  of Table 18. Differencing two values already rounded to one decimal moves three of the
  twenty entries by a tenth.
- **Appendix H.9's orientation.** `bounded_claim()` computes p_b − p_a with AutoDock as
  a, so its difference and interval are the negation of the printed ones.
- **Appendix H.9's equivalence bounds.** TOST at alpha 0.05 reads a 90 per cent
  interval, not the 95 per cent one quoted beside it, and the smallest passing margin is
  the larger absolute bound rounded *up*. Rounding down names a margin that fails.
- **Table 17's GSK file.** The row is the neutral species and must come from
  `gsk7975a-deprot-as-neutral-OPT.sdf`. The sibling `gsk7975a-deprot-OPT.sdf` is the −1
  phenolate, a different species, and reading it changes four cells plausibly. The
  harness asserts the two stay distinct.
- **Table 17's hydrogen basis.** Rotatable bonds must be counted hydrogen-suppressed.
  On the explicit basis the three ligands read 6, 7 and 4 instead of 5, 5 and 4.
- **Table 4's rounding.** Compared unrounded. The CSV carries 0.5035 and 0.4635, which
  the thesis prints as 0.504 and 0.463, and Python's `round` takes both the other way.
- **Table 11's arm keys.** The Orai control PoseBusters table names the engine and
  carries one arm per engine. The composed variant keys used on the benchmark tree are
  not in that file.
