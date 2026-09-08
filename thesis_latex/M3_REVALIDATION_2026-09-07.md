# Re-validation of examiner finding M3 (2026-09-07)

**Finding under test** — "Appendix D's Fr0 receptor description does not match the receptor the
reported AutoDock arms docked, and the caveat built on it is stale." (Finding I-1, upheld at major.)
**Source** — `thesis_latex/EXAMINER_REVIEW_2026-09-06_multiagent.md:269-282`
**Method** — 15-agent workflow: six independent lenses (residue alignment, run provenance, origin of
the numbers, downstream dependencies, M1 replacement numbers, steelman), three adversarial refuters
(evidence, logic, proportionality), three drafted fix options each checked by an accuracy-and-style
judge. Key counts were then re-measured by hand.
**Adjudicated against** — the live short build, both reported docking trees
(`Dockings/Orai_Benchmark_MGLTools_exh128`, `Dockings/Orai_JKU_MGLTools_exh128`), the retired
`Dockings/vina_results/JKU_meeko` arm and its PoseBusters backup.

---

## Verdict: CONFIRMED in substance — but the finding's rationale needs one correction

The paragraph at `body_appendix_short.tex:560` is false of every reported arm and should go. It was
not, however, fabricated or unexecuted. It is a bit-accurate description of a retired arm that went
stale when the panels were repointed on 2026-08-30.

| # | Sub-claim | Verdict |
|---|---|---|
| 1 | Reported arms docked a 1,338-residue / 10,368-heavy-atom Fr0 file (md5 9fd4928d) | **CONFIRMED** — byte-identical in both trees; 3,079 + 120 provenance JSONs carry its sha256 6b3ab996 |
| 2 | The 1,333 / 10,321 / five-residue / 41-atom figures reproduce from the Meeko file | **CONFIRMED** — set difference is exactly Phe253 A (11), Ser93 B (6), Leu248 C (8), His256 C (10), Ser93 D (6) |
| 3 | No residue was dropped from the reported search receptor | **CONFIRMED** — strict superset of the delivered 10,362 atoms; only change is +1 OXT on His288 of each chain |
| 4 | The caveat "attaches to nothing" | **CONFIRMED** — search and validation receptors for every reported Fr0 AutoDock row are coordinate-identical (10,368/10,368); reported Fr0 AutoDock validity is 100.0 % (JKU) and 98.90 % (control, failures at Ala112 D, Glu210 C, Val37 C, none a named residue) |
| 5 | No Results or Discussion sentence cites the caveat | **CONFIRMED** — zero hits in the short build, the reproduction apparatus (`thesis_expected_values.yaml`, `thesis_assertions.py`, `REGENERATE.md`, notebooks, glossary) or `Future Research/`; line 560 carries no `\label`, `\footnote` or `\ref` |
| 6 | "No docking run located in the repository ever staged that Meeko file, so the 1,333-residue search may never have been executed" | **REFUTED** — `Dockings/vina_results/JKU_meeko/` is intact: Vina batch log "Rigid receptor: .../Orai1WT-START-Fr0_meeko_retry1.pdbqt", docking log `prot_num_residues 1333 / prot_num_atoms 12509`, exhaustiveness 32, 120 gnina SDFs with `receptor_sha256 b418b635` (the Meeko file). No staging copy was needed, which is why a staging search missed it |
| 7 | (implicit) the numbers were never supported | **REFUTED** — `posebusters_results/orai_jku_UFFON_backup_20260830` holds exactly 240 Fr0 AutoDock rows (4 ligands × 30 modes × {raw, gnina}); 42/240 reproduce within 8 Å on an all-atom basis (35 heavy-only); exactly two Synta-66 rows fail only min-distance + overlap, both at protein atom 1421 = Phe253 A CB. They are the same Vina pose (model 4) counted raw and gnina |
| 8 | "Fr0 gained six repaired atoms rather than losing 41, so no conservative correction is owed" | **CONFIRMED but not load-bearing** — the validation receptor carries the same six OXT, so there is no composition mismatch of either sign; receptor identity is the decisive fact |

### Two defects the finding under-states

- **The paragraph contradicts the thesis's own Appendix C.10.** `:560` says PoseBusters screened Fr0
  "against the intact heavy-atom structure"; `:489` correctly says the START-Fr0 AutoDock rows were
  rescreened against "the energy-minimised one the search actually read" (10,368 atoms). Deleting
  `:560` removes a live contradiction.
- **The paragraph was over-scoped from the day it was written.** It says "Its AutoDock search
  receptor alone", but the Orai × Benchmark control panel never docked the Meeko file. Even on
  2026-08-20 it was true of one of the two Orai AutoDock panels.

### Timeline

| Date | Event |
|---|---|
| 2026-07-24 | Retired Orai × JKU Meeko arm docks `Orai1WT-START-Fr0_meeko_retry1.pdbqt` (exh 32) |
| 2026-08-18 20:04 | 42-of-240 and Synta-66 numbers computed (audit session 52095a0c) on that arm |
| 2026-08-20 08:21 | Paragraph enters the long build (commit 40e5c318) — true of the then-reported JKU arm |
| 2026-08-30 ~20:12 | Panels repointed to the MGLTools exh128 trees — paragraph goes stale |
| 2026-08-30 22:26 | Paragraph carried verbatim into the short build (commit 5c150c08) |

---

## What the reported Fr0 AutoDock receptor actually is (verified numbers for any replacement)

| Quantity | Value | Basis |
|---|---|---|
| Residues / heavy atoms / OXT | 1,338 / 10,368 / 6 (His288 A–F) | staged PDBQT, md5 9fd4928d |
| Heavy-atom RMSD vs delivered PDB | 0.819 Å (same after Kabsch; centroid shift 0.029 Å) | 10,362 matched atoms with CHARMM Ile CD→CD1 alias |
| Cα RMSD | 0.636 Å over 1,338 | same |
| Largest displacement | 3.631 Å, Lys203 NZ chain A | same |
| Fr300 / Fr400 / Fr499 staged vs delivered | 0.000 Å over 10,362 | MGLTools moves no heavy atoms |
| DiffDock and EquiBind Fr0 inputs vs delivered | 0.000 Å, 10,362/10,362 | `prepared_proteins`, `_prep_Orai1WT-START-Fr0` |
| Parent of the staged file | `Dockings/Orai_Benchmark_old/_staging/receptors/pdbqt/Orai1WT-START-Fr0_minimised.pdb` (OpenMM 8.3.1, 2026-06-22) | coordinate-multiset identical |
| Recipe | PDBFixer findMissing*/addMissingAtoms/addMissingHydrogens(7.0), amber14-all + GBn2, ≤500 steps | `Scripts/Utilities/prep_docking.py:239-286`, fallback at `:993-1006` |
| Atom-type trace of pH-7 placement | N +11 / NA −11 vs Fr300 (Asn/Gln amide N), OA +6 (OXT) | histidines identical in both frames |

**Caveat:** the fallback firing is inferred from the code path and the Fr0-only `_minimised.pdb`
files. No log records "hydrogen addition failed". Do not write that a log shows it.

---

## Options

All three edit the short build only. `thesis_latex/obsolete/body_appendix.tex:378` and
`obsolete/body_main.tex:906` carry the same text and are reported, not edited.

### Option A — delete the paragraph (lowest risk, −128 words)

Remove lines 560 and 561 (the paragraph and one of its two flanking blank lines):

```bash
sed -n '560p' thesis_latex/body_appendix_short.tex | cut -c1-60   # must print "One qualification attaches to the Fr0 snapshot."
sed -i '560,561d' thesis_latex/body_appendix_short.tex
```

No bridging clause is needed. "All 1,338 Cα atoms" at `:583` and "10,362 atoms" at `:558` are
measured from the delivered PDBs and survive independently. Closes M3 only; leaves the Fr0
minimisation disclosed solely at `:489`.

### Option B — replace with the correct disclosure (closes M1 and M3 together; net about +30 words)

Replacement for line 560:

```latex
One qualification attaches to the Fr0 snapshot. Its AutoDock search receptor alone was passed through a PDBFixer repair and an OpenMM energy minimisation before conversion, the fallback the preparation pipeline applies when the converter cannot add hydrogens to a frame. That second minimisation retained all 1,338 residues and added six OXT atoms, one at the C-terminal His288 of each chain, giving 10,368 heavy atoms. The minimised coordinates lie 0.82\angstrom{} from the delivered ones over all heavy atoms and 0.64\angstrom{} over the C\ensuremath{\alpha} atoms, with no rigid-body shift. The largest single displacement is 3.63\angstrom{} at the Lys203 side chain of chain A. DiffDock and EquiBind read the delivered coordinates without that repair. The Fr0 frame is therefore the one place where the three tools did not search the same receptor geometry. The Fr0 AutoDock poses were rescreened against their own search receptor, as recorded in Appendix~\ref{posebusters-validity-checks}.
```

Line 558, last two sentences (the pH clause is otherwise contradicted by `addMissingHydrogens(7.0)`):

```latex
The preparation step removes every hydrogen and renames HSD, HSE and HSP to plain HIS. The delivered receptor is therefore a purely heavy-atom structure of 10,362 atoms. Polar hydrogens are then re-added by the PDBQT converter under its own geometric and atom-type rules rather than by a titration model. No pH-based protonation was applied at any stage, apart from the Fr0 repair described below, which placed hydrogens at pH 7.0 before the converter rebuilt them.
```

`body_main_short.tex:97`, first sentence:

```latex
The Orai1 molecular-dynamics snapshots, Fr0, Fr300, Fr400 and Fr499, represented limited receptor flexibility \cite{ref038}, \cite{ref072} and received no pH-based protonation, apart from one repair step applied to Fr0 alone.
```

`body_main_short.tex:156`, first sentence:

```latex
All three tools dock into the same receptor coordinates for a given frame, except that the Fr0 AutoDock file is an energy-minimised copy displaced by 0.82\angstrom{} over heavy atoms without a rigid-body shift. Their poses therefore share one coordinate system, and the within-frame site metrics apply no symmetry folding.
```

Checks after applying B: `:163` ("The pH-based hydrogen placement PDBFixer offers is left disabled on
this path") is scoped to the benchmark path and stays true. "Energy-minimised geometry" at `:95` and
`:609` keeps its MD-origin sense; the new paragraph says "second minimisation" to keep the two apart.

**B-lite:** use the new `:560` paragraph but stay silent on protonation (drop the `:558` and `:97`
changes). This leaves the "no pH-based protonation at any stage" clause technically contradicted by
the fallback's pH-7 placement, which M1 flags.

### Option C — one-to-three-sentence corrected caveat (−68 words)

```latex
The AutoDock search receptor for the Fr0 snapshot carried all 1,338 residues and 10,368 heavy atoms, the delivered set plus six C-terminal OXT atoms added during its preparation for docking. An earlier converter pass on that frame dropped five residues under a permissive residue filter. That file was searched only by a retired arm whose results appear nowhere in this work.
```

The judge ranked sentence one alone above the full three, because sentences two and three describe a
retired arm no reader can locate and Appendix D nowhere else discusses superseded preparations.
C does not disclose M1.

### Recommendation

Option B. It removes the stale paragraph, resolves the `:560`-vs-`:489` contradiction, and puts the
real Fr0 asymmetry (M1) in the chapter an examiner reads for it, with every number verified on the
exact file both panels docked. If M1 is deferred, take Option A now rather than keeping `:560`; the
deletion must not wait on M1.

---

## Housekeeping found on the way

- `posebusters_results/orai_jku/dock/converted_pdbqt/` still holds 124 stale `*_meeko_retry1__*` SDFs
  beside the 120 live Fr0 ones. The current CSV never references them, but a glob-based re-bust would.
- The 08-31 validation report (`obsolete/EXAMINER_ASSESSMENT_VALIDATION_2026-08-31.md:114`) called
  these numbers unsupported by "any committed receptor or result file". That charge is wrong and should
  not be repeated.
