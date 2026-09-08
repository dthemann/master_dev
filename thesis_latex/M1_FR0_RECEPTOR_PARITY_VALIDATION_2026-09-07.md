# M1 / E1-1 validation: Fr0 AutoDock receptor minimised, parity with DiffDock and EquiBind broken

Date 2026-09-07. Validation of the examiner finding M1 (E1-1) from `EXAMINER_REVIEW_2026-09-06_multiagent.md`.
Method: inline measurement, then a nine-agent workflow (six adversarial verifiers with distinct lenses, two independent
LaTeX drafters, one judge). Workflow run `wf_2c74e282-814`; scratch outputs under
`/tmp/claude-1000/-home-manndo-master-dev/1e9da683-e8b2-4212-a04f-69a12a5fddec/scratchpad/{parity,protonation,fallback-why,fr0-contact,m1_downstream}`.
The drop-in edits in section 8 were APPLIED on 2026-09-07 (eleven line edits plus the new Appendix D paragraph at `body_appendix_short.tex:562`); `Thesis_short.pdf` rebuilt at 142 pages with no errors or undefined references.

## 1. Verdict

**The finding is CONFIRMED on every measurable sub-claim, and it is understated in one respect.** The AutoDock Fr0 search
receptor in both reported arms is an OpenMM-minimised, PDBFixer-repaired copy that DiffDock and EquiBind never read. The
"no pH-based protonation" sentences are false as written for that file. Two refinements: the examiner's "bounded below
critical by magnitude" argument does not survive a contact-level test, and the examiner's edit list is incomplete
(six further lines, one cited line is blank).

| sub-claim | verdict | evidence |
|---|---|---|
| Staged Fr0 PDBQT (sha256 `6b3ab996`, byte-identical in `Orai_Benchmark_MGLTools_exh128` and `Orai_JKU_MGLTools_exh128`) is minimised | CONFIRMED | 0.819 A heavy-atom, 0.636 A C-alpha, 3.631 A max (Lys203 chain A NZ) from `Data/Receptors/Orai1WT-START-Fr0.pdb`; 0.000 A from the OpenMM-written `Orai_Benchmark_old/_staging/receptors/pdbqt/Orai1WT-START-Fr0_minimised.pdb` |
| Six OXT added on Fr0 alone; 10,368 vs 10,362 heavy atoms | CONFIRMED | PDBFixer `findMissingAtoms` reports exactly six missing C-terminal OXT on the delivered Fr0 |
| Fr300 / Fr400 / Fr499 staged receptors are 0.000 A from delivered | CONFIRMED | 10,362 / 10,362 coordinate multiset matches on each; `prepare_receptor` on the delivered files reproduces the staged PDBQTs bit-for-bit |
| Both AutoDock variants (raw Vina and gnina) on Fr0 read the minimised file | CONFIRMED | 3,079 + 120 gnina provenance JSONs record `receptor_sha256 6b3ab996`; both raw-Vina batch logs print `Rigid receptor:` pointing at that file; `docking_log` `prot_num_atoms` 12,636 = its ATOM count |
| DiffDock read the unminimised Fr0 | CONFIRMED | `prepared_proteins/Orai1WT-START-Fr0_prepared.pdb` (both trees, sha `60a95388`) 0.000 A from delivered, 10,362 atoms; no minimisation code path in `run_diffdock.py` |
| EquiBind read the unminimised Fr0 | CONFIRMED | all reported EquiBind trees stage a file byte-identical (`d9890e05`) to the delivered PDB; pipeline only runs `reduce` or copies |
| Fallback fired because `prepare_receptor` hydrogen addition crashes on delivered Fr0 | CONFIRMED by re-execution | exit 1, `PyBabel/addh.py:303 IndexError` on Fr0; exit 0 on Fr300/400/499. Trigger: Ser93 C-alpha to C-beta bond stretched to 1.85 to 2.06 A in all six chains, above MolKit's 1.791 A C-C bond-perception cutoff |
| No run log records the firing | CONFIRMED, but closed | `prep_docking.py` logs only to stdout, no notebook output captured it. Converting the Jun-22 `_minimised.pdb` with the pipeline flags reproduces sha `6b3ab996` byte-for-byte, which is stronger evidence than a log |
| pH-7 hydrogen model ran on Fr0 | CONFIRMED | 2,208 of the 2,268 polar hydrogens in the staged Fr0 PDBQT are coordinate-identical to PDBFixer `addMissingHydrogens(7.0)` output and carry OpenMM names (`H`, `HH11`, `H2/H3`); Fr300 carries converter names (`HN`, `1HH1`, `HN1`) |
| All 60 histidines doubly protonated exactly as in Fr300 | CONFIRMED | PDBFixer wrote all 60 as HID (HD1 only); the converter then added HE2 to each, giving the same HD1+HE2 typing it builds from bare heavy atoms on the other frames. Asp, Glu, Lys, Arg, Tyr and N-termini typed identically on all four frames |
| 11 NA lost, 6 OA gained | CONFIRMED, mechanism reinterpreted | the 6 OA are the OXT; the 11 NA on Fr300 (6 on Fr400, 10 on Fr499) are Asn/Gln amide nitrogens whose hydrogens the converter failed to place. Hydrogen-completeness effects, not titration effects. Fr0 alone also carries 18 Cys thiol hydrogens |
| Three `_minimised.pdb` artefacts, all Fr0 | UNDERCOUNT | five `.pdb` (four distinct conformers: Apr 19, Jun 22, Aug 6, Aug 28) plus two `.pdbqt`, all Fr0. The reported arms use the Jun-22 conformer |
| `:489` already discloses the minimisation | CONFIRMED | but it is framed as a screening-file matter, does not say the AutoDock preparation performed it, does not say DiffDock/EquiBind read the delivered file, and contradicts `:560` |

## 2. Where the finding is wrong or incomplete

1. **`body_main_short.tex:150` is a blank line.** The preparation sentences the examiner meant are `body_appendix_short.tex:150`
   (ligand-scoped, "no pK_a model ran at any stage") and `body_appendix_short.tex:163` ("The pH-based hydrogen placement PDBFixer
   offers is left disabled on this path" and "Chain termini are taken as deposited and are not repaired"). Both are violated by the Fr0 path.
2. **The edit list misses six lines.** `main:95` ("contrasts one energy-minimised geometry"), `app:163`, `app:572` ("Every statement ...
   rests on the coordinates as delivered"), `app:583` (the 5.88 to 6.02 A spread was computed from the delivered Fr0), `app:609`
   ("controlled contrast"), and the `app:560` vs `app:489` internal contradiction on what PoseBusters screened Fr0 against
   (the current `_orai_matched_root` CSVs prove `:489`: AutoDock Fr0 rows are screened against `_orai_benchmark_staging_fr0corrected`,
   DiffDock and EquiBind rows against `Data/Receptors`).
3. **"Table 20 / App C.8 / C.10" hold no Fr0 statements in the current build.** The screening statements are `app:489` (C.9) and `app:560` (App D).
4. **Minor number drift.** 149 atoms above 2 A and 0.819 A heavy-atom under a keying that maps ILE CD to CD1 and excludes the OXT
   (the finding's 145 / 0.815 are the same measurement under a slightly different atom set).
5. **The "0.394 A between two preps" figure in `REGENERATE.md:664` and `run_orai_mgltools_arm.py:8` is not reproduced** by any pair of the
   four surviving minimised conformers (0.21 to 0.30 A heavy-atom). Do not quote it.
6. **The `:156` consequence clause is true.** The minimisation was in place (centroid shift 0.024 A, superposed vs unsuperposed C-alpha RMSD
   0.635 vs 0.636 A), so "their poses already share one coordinate system" and the within-frame site metrics are geometrically sound.
   Only the premise "same receptor file" is false.

## 3. Magnitude: the C-alpha number understates what docking sees

The examiner bounds severity by 0.64 A C-alpha against the 5.9 to 6.0 A frame separation. At the contact level the two Fr0 receptors are
not interchangeable in either direction:

| test | result |
|---|---|
| Fr0 AutoDock poses (control panel, 3,079 raw Vina) with any receptor atom within 4 A that moved > 1 A | 99.3 % (rank-1: 99.7 %) |
| ... moved > 2 A | 65.8 % (rank-1: 70.5 %) |
| median fraction of each pose's 4 A shell that moved > 1 A | 36 % |
| residues with any atom moved > 2 A | 97 of 1,338; 25 inside the Arg91 to Glu106 slab, including Glu106 (2 chains), Lys87 (4), Tyr80 (6) |
| Vina `--score_only` of rank-1 poses against a PDBQT carrying the delivered coordinates | mean +3.5 kcal/mol worse (median +2.7); 95.5 % worsen by > 1 kcal/mol |
| within-complex rank correlation, staged vs delivered scores (92 complexes) | median Spearman -0.04 (the Fr0 rank order is receptor-specific) |
| approximate PoseBusters clash (0.75 x vdW) against delivered atoms | 23.5 % of poses (rank-1: 31.8 %) vs 0.1 % against staged |
| reverse test (quarantined CSVs under `obsolete/`): DiffDock Fr0 pass-all against the minimised receptor | 79.1 % to 39.6 % |
| reverse test: EquiBind Fr0 | 58.2 % to 15.8 % (minimised arm 68.3 % to 14.7 %) |

gnina-optimised poses inherit exactly the same shell statistics, so gnina minimisation does not decouple them from the receptor.
This is consistent with the 733 spurious minimum-distance failures at `:489`: the receptor difference is contact-relevant, and the
correction paragraph should not call it negligible. It remains one frame of four and the minimisation is in place, which is why
"major" rather than "critical" is still the right grade.

A physical cause exists and is worth one sentence in Appendix D: the delivered Fr0 is geometrically defective (nine heavy-atom bonds
above 1.9 A, Ser93 C-alpha to C-beta 1.85 to 2.06 A in every chain, His256 chain C C-alpha to C 2.21 A; Fr300/400/499 have only normal
sulphur bonds above 1.9 A). This matches the memory that Fr0 is a restrained pre-dynamics frame. DiffDock and EquiBind therefore docked
into a receptor with stretched covalent bonds, AutoDock into the repaired one.

## 4. Downstream: Fr0 is where the between-tool picture changes

The Results chapter contains no Fr0-resolved between-tool number, and Fr0 sits inside every pooled between-tool Orai statement.
Recomputed from the current `_orai_matched_root` CSVs (reproduces Table `tab:results-orai-yield` exactly under the top-10 cap):

| control panel | Fr0 | Fr300 | Fr400 | Fr499 |
|---|---|---|---|---|
| AutoDock Vina usable yield | **67.3 %** | 86.3 % | 89.4 % | 88.3 % |
| AutoDock valid poses lost to the slab | **32.0 %** | 12.9 % | 9.8 % | 10.9 % |
| DiffDock usable yield | **76.2 %** | 68.3 % | 67.9 % | 71.8 % |
| DiffDock valid poses lost to the slab | 3.7 % | 13.1 % | 13.3 % | 8.5 % |
| AutoDock to DiffDock agreement at 5 A | **27.0 %** | 7.6 % | 12.6 % | 10.9 % |

Pooled with vs without Fr0: AutoDock to DiffDock agreement 14.3 % to 10.4 % (text says 14 %); all-pairs 9.1 % to 6.8 % (text 9 %);
all-pairs unit flag 8.1 % to 4.8 % (Fr0 supplies 53 of the 97 units); control slab loss AutoDock 16.4 % to 11.2 %, DiffDock 9.6 % to 11.6 %.
The `:635` sentence "DiffDock loses less on both panels" holds on the control panel only because of Fr0. The tool ordering inverts on Fr0
alone (AutoDock 67.3 % below DiffDock 76.2 %). The two headline conclusions survive dropping Fr0 (usable-yield reversal across panels;
AutoDock cross-panel Mann-Whitney p 1.9e-4 to 6.2e-4).

Whether the Fr0 outlier is caused by the preparation or by Fr0's conformation cannot be separated from existing data, because no
unminimised AutoDock Fr0 arm exists. The wide cytosolic vestibule where AutoDock's excess Fr0 poses sit is present in the delivered Fr0
and changes by at most 1.6 A under minimisation, which weakly favours conformation. EquiBind, with no preparation difference, is also a
Fr0 outlier on validity (68.3 % vs 40.6 to 51.1 %). The transmembrane slab for Fr0 was built from the delivered file for all three tools
(`tm_exclusion_zones.json`), so the placement rule is common even though the AutoDock screening receptor is not.

## 5. Options

| option | what it does | cost | what it settles |
|---|---|---|---|
| **A. Documentation only** (the examiner's proportionate correction) | apply the 11 line edits in section 8 and the new Appendix D paragraph | about one hour plus a LaTeX rebuild; no compute | states the parity break, the numbers, and who read what. Does not resolve whether Fr0's between-tool outlier is preparation or conformation |
| **B. A plus a Fr0-resolved sentence in Results 4.2** | add the per-frame usable-yield and agreement numbers above, and the with/without-Fr0 agreement rates, with a small reproducible script or notebook cell | about half a day | tells the examiner the pooled between-tool numbers that Fr0 carries; the numbers exist but currently only in scratch scripts |
| **C. Restore parity by re-docking AutoDock Fr0 on the delivered coordinates** | a hydrogen-only minimisation (heavy atoms fixed, 2 s on CUDA) yields a PDBQT that `prepare_receptor -A hydrogens` accepts with heavy atoms identical to the delivered file, the pipeline's own hydrogen naming and double-His typing, no OXT. Verified feasible in scratch. Then re-run the Fr0 slice of both arms | search 4.2 h serial for the 308 control ligands plus 2 min for the modulators; gnina rescoring about 2.1 h; PoseBusters re-bust; the 26 Orai `REGENERATE.md` stages; every Orai number in 4.2 and the appendices moves; the sha-256 pin in `run_orai_mgltools_arm.py` and the notebook must be changed deliberately. One to two days | the only option that answers "preparation or conformation" and makes Fr0 receptor-matched. Caveats: six Arg/Val atom-type and 49 charge differences relative to the current typing; the receptor keeps its stretched bonds, which Vina ignores; the naive PDBFixer-hydrogens route silently produces 221 chargeless hydrogens because the vina env writes X-H at 1.18 A, so the heavy-atoms-fixed minimisation step is required |
| **D. Re-dock DiffDock and EquiBind into the minimised Fr0** | parity in the other direction, all three tools on the repaired receptor | DiffDock 34 GPU-hours for the control Fr0 slice plus 40 min for the modulators and 3.7 h of smina/gnina optimisation; EquiBind about 9 h; regeneration as in C; known CUDA OOM history on Orai DiffDock | same as C but at roughly ten times the compute |
| **E. Drop Fr0 from between-tool statements** | keep Fr0 as an AutoDock-internal sensitivity arm, report cross-tool Orai on Fr300/400/499 | 26 regeneration stages with a frame filter; loses the `:609` minimised-vs-thermal contrast; discards data | avoids the confound but the examiner asked for disclosure, not removal, and the pooled numbers move anyway (14 % to 10.4 %) |

## 6. Recommendation

Do A now and B with it. The drop-in text is ready and every number in it is measured from the reported trees. Add B because the
downstream check shows the pooled agreement rate and the `:635` sentence are carried by Fr0, and an examiner who reads the new appendix
paragraph will ask what Fr0 does to the Results numbers. Offer C in the defence as the experiment that would settle attribution, and run it
only if there is time for a full Orai regeneration, since it moves every Orai number. Do not run D or E.

Two housekeeping items while you are there: correct or delete the "0.394 A" claim in `REGENERATE.md:664` and `run_orai_mgltools_arm.py:8`,
and note that the 733 figure at `:489` comes from the superseded exh32 control tree (the current exh128 control validates Fr0 at 98.9 %
against its own receptor), so the sentence should say which run it describes.

## 7. Placement notes from the judge

The new paragraph goes directly after the replaced `:560` paragraph and before the residue-numbering paragraph at `:562`, so that
"That fallback" refers to the PDBFixer/OpenMM fallback just described and "described above" at `:609` resolves. Apply the two line-97
edits in the listed order. Net length about +200 words in the appendix and +30 in Methods. Left untouched on purpose: `main:95` still
reads at design level and is now qualified by `:609`; `app:167` and `:402` (the Fr0 search box was fitted to the delivered file, not the
prepared receptor, exact match 98.3 x 91.8 x 92.6 A) is a minor inconsistency; `app:361` and `:363` stay true between the two AutoDock
panels since both stage sha `6b3ab996`. Every `current_excerpt` below was grep-verified verbatim and unique at its line, and the three
`\ref` labels exist.

## 8. Drop-in edits (verbatim, applied 2026-09-07; four lines later corrected, see section 10)

### 8.1 New Appendix D paragraph, insert after the replaced `body_appendix_short.tex:560`

```latex
That fallback breaks receptor parity on Fr0. DiffDock and EquiBind read the delivered coordinates on every frame. The AutoDock Fr0 search receptor lies 0.82\angstrom{} from them over the 10,362 shared heavy atoms and 0.64\angstrom{} over the 1,338 C\ensuremath{\alpha} atoms, with no rigid-body component. The shift exceeds 2\angstrom{} for 149 atoms in 97 residues, 25 of them inside the Arg91 to Glu106 pore slab, and peaks at 3.63\angstrom{} at the side-chain nitrogen of Lys203 in chain A. The displaced residues include Tyr80 and Lys87, the two that occur most often within 4\angstrom{} of a rank-1 Fr0 AutoDock pose. Of the 2,268 polar hydrogens in that file, 2,208 were placed by PDBFixer and the converter added only the second imidazole hydrogen on each of the 60 histidines. The histidine, aspartate, glutamate, lysine and arginine typing therefore matches the other frames, while Fr0 alone carries 18 cysteine thiol hydrogens and six carboxylate termini. The Fr0 cross-tool comparison is consequently not receptor-matched, and the difference is felt at contact level. Of the 3,079 control-panel Fr0 AutoDock poses, 99.3\% touch a receptor atom within 4\angstrom{} that moved by more than 1\angstrom{}, which is the difference behind the spurious minimum-distance failures recorded in Appendix~\ref{posebusters-validity-checks}. The shift is nonetheless in place and about a tenth of the 5.88 to 6.02\angstrom{} that separate Fr0 from the later frames after superposition. The within-frame coordinate system the site metrics rely on is therefore unchanged.
```

### 8.2 Line edits

#### 1. `thesis_latex/body_appendix_short.tex:558`

**Current (verbatim):**

```latex
Polar hydrogens are then re-added by the PDBQT converter under its own geometric and atom-type rules rather than by a titration model, and no pH-based protonation was applied at any stage.
```

**Replacement:**

```latex
Polar hydrogens are then re-added by the PDBQT converter under its own geometric and atom-type rules rather than by a titration model. The AutoDock receptor for Fr0 alone took a different path, which the next paragraph sets out.
```

*Why:* Draft A base. The clause 'no pH-based protonation was applied at any stage' is false for the AutoDock Fr0 receptor, whose 2,208 of 2,268 polar hydrogens are coordinate-identical to the PDBFixer addMissingHydrogens(7.0) output. Splitting into two sentences avoids the exception clause of Draft B and keeps the converter statement true for Fr300, Fr400, Fr499 and for every DiffDock and EquiBind receptor.

#### 2. `thesis_latex/body_appendix_short.tex:560`

**Current (verbatim):**

```latex
One qualification attaches to the Fr0 snapshot. Its AutoDock search receptor alone required the permissive residue filter of the PDBQT converter, which dropped five residues carrying 41 heavy atoms, namely Ser93 in chains B and D, Leu248 in chain C, Phe253 in chain A and His256 in chain C. That search ran against 1,333 residues and 10,321 heavy atoms rather than the full assembly. Every other frame carried the complete 1,338 residues, and PoseBusters screened all frames including Fr0 against the intact heavy-atom structure. Forty-two of the 240 Fr0 AutoDock poses lie within 8\angstrom{} of a dropped residue, and two Synta-66 poses fail only the protein-distance and protein-overlap checks against an atom of the deleted Phe253. The reported Fr0 AutoDock validity is therefore marginally conservative for those poses.
```

**Replacement:**

```latex
One qualification attaches to the Fr0 snapshot. ADFRsuite \texttt{prepare\_receptor} aborts on the delivered Fr0 coordinates, because its hydrogen builder fails at the bond between C\ensuremath{\alpha} and C\ensuremath{\beta} of Ser93, which is stretched to 1.85 to 2.06\angstrom{} in every chain of that frame. The AutoDock preparation therefore fell back to PDBFixer and OpenMM for that frame alone. PDBFixer added the six missing C-terminal OXT atoms and placed hydrogens at a nominal pH of 7.0, and OpenMM then energy-minimised the complete structure in the amber14-all force field with GBn2 implicit solvent. The converter was run on that output. The AutoDock Fr0 search receptor therefore carries the complete 1,338 residues and 10,368 heavy atoms at minimised coordinates, and PoseBusters screened the Fr0 AutoDock poses against those same coordinates. The Fr300, Fr400 and Fr499 search receptors are coordinate-identical to the delivered files, and DiffDock and EquiBind read the delivered heavy-atom coordinates on every frame including Fr0.
```

*Why:* Draft B base for the trigger (chain A Ser93 CB, CA-CB 1.873 A above the 1.791 A MolKit cutoff, 1.850 to 2.062 A across chains A-F, re-executed IndexError at PyBabel/addh.py:303) grafted with Draft A's two closing facts that Draft B dropped: what PoseBusters screened Fr0 AutoDock against (current _orai_matched_root CSVs use _orai_benchmark_staging_fr0corrected, 0.000 A from the staged PDBQT) and that DiffDock and EquiBind read the delivered coordinates (DiffDock prepared Fr0 0.000 A; EquiBind staging digests byte-identical to Data/Receptors on all four frames, re-verified). The old paragraph described the never-docked Meeko retry1 file (1,333 residues, 10,321 atoms) and a 240-pose slice that no longer exists (current JKU Fr0 AutoDock slice is 90 rows). Both reported arms stage sha256 6b3ab996 with 1,338 residues, 10,368 heavy atoms and 6 OXT, byte-reproduced by prepare_receptor on the OpenMM-written _minimised.pdb. 'Those same coordinates' replaces Draft A's 'the same file' because the screening receptor is a PDB coordinate-identical to the PDBQT, not the PDBQT itself.

#### 3. `thesis_latex/body_main_short.tex:97`

**Current (verbatim):**

```latex
The Orai1 molecular-dynamics snapshots, Fr0, Fr300, Fr400 and Fr499, represented limited receptor flexibility \cite{ref038}, \cite{ref072} and received no pH-based protonation. The PDBQT converter added polar hydrogens by its geometric rules rather than a titration model (Appendix~\ref{orai1-receptor-model}).
```

**Replacement:**

```latex
The Orai1 molecular-dynamics snapshots, Fr0, Fr300, Fr400 and Fr499, represented limited receptor flexibility \cite{ref038}, \cite{ref072}. The PDBQT converter added polar hydrogens by its geometric rules rather than a titration model, except on the AutoDock copy of Fr0, which the preparation fallback repaired, protonated at pH 7.0 and energy-minimised before conversion (Appendix~\ref{orai1-receptor-model}).
```

*Why:* Draft A's fuller excerpt (unambiguous match) with Draft B's verb 'protonated'. 'Received no pH-based protonation' is false for the AutoDock Fr0 receptor, which carries OpenMM hydrogen names (H, H2/H3, HH11) and 6 added OXT from the PDBFixer fallback at Scripts/Utilities/prep_docking.py:236-286. The exception is named in the Methods with the appendix cross-reference kept.

#### 4. `thesis_latex/body_main_short.tex:97`

**Current (verbatim):**

```latex
Because no titration model was applied,
```

**Replacement:**

```latex
Because no titration model was applied to the ligands,
```

*Why:* Draft B only. Once the receptor exception is admitted two sentences earlier, the unscoped clause would contradict it. The sentence explains the neutral modulator forms, so the ligand scope is what it means. Excerpt verified unique at line 97.

#### 5. `thesis_latex/body_main_short.tex:156`

**Current (verbatim):**

```latex
All three tools dock into the same receptor file for a given frame, so their poses already share one coordinate system and the within-frame site metrics apply no symmetry folding.
```

**Replacement:**

```latex
All three tools dock into the same receptor coordinates for a given frame, the AutoDock copy of Fr0 excepted, which was energy-minimised in place by 0.64\angstrom{} over the C\ensuremath{\alpha} atoms (Appendix~\ref{orai1-receptor-model}). Their poses therefore already share one coordinate system and the within-frame site metrics apply no symmetry folding.
```

*Why:* Merged. Draft B's 'coordinates' replaces 'file' (DiffDock and EquiBind never read the same file, only the same coordinates) and Draft A's 0.64 A CA figure (0.636 A) states the size of the exception. The consequence survives because the minimisation was in place (centroid shift 0.024 A, superposed vs unsuperposed CA RMSD 0.635 vs 0.636 A), and 'in place' carries that justification. The ', so ...' chain both drafts kept is broken into two sentences to match the house style.

#### 6. `thesis_latex/body_appendix_short.tex:150`

**Current (verbatim):**

```latex
no protomer or tautomer enumeration and no pK\textsubscript{a} model ran at any stage.
```

**Replacement:**

```latex
no protomer or tautomer enumeration and no pK\textsubscript{a} model ran on any ligand at any stage.
```

*Why:* Identical in both drafts. The paragraph is ligand-scoped but the clause reads as global and is the line the examiner's ':150' most plausibly meant (body_main_short.tex:150 is blank). Scoping to ligands keeps it true alongside the Fr0 receptor exception.

#### 7. `thesis_latex/body_appendix_short.tex:163`

**Current (verbatim):**

```latex
The pH-based hydrogen placement PDBFixer offers is left disabled on this path. No benchmark receptor was therefore titrated at a nominal pH.
```

**Replacement:**

```latex
The pH-based hydrogen placement PDBFixer offers is left disabled on this path. No benchmark receptor was therefore titrated at a nominal pH. A fallback to that placement exists for receptors the converter cannot hydrogenate, and it fired on the Orai1 Fr0 receptor alone (Appendix~\ref{orai1-receptor-model}).
```

*Why:* Draft A's three-sentence structure, with the fallback introduced before it is referred to (neither draft had introduced the word 'fallback' at this point in the file). prep_docking.py:992-1012 fires the PDBFixer/OpenMM branch on an addh.py IndexError, reproducible on Fr0 and on no other reported receptor, and every surviving *_minimised* artefact is Fr0, so the benchmark-only conclusion stands. The 'chain termini are taken as deposited' sentence later in the same paragraph is benchmark-scoped and stays true.

#### 8. `thesis_latex/body_appendix_short.tex:489`

**Current (verbatim):**

```latex
They ran under the shared-receptor policy, so every pose was screened against a single file per frame and the exactness check never applied. On the START-Fr0 frame that file was the deposited structure rather than the energy-minimised one the search actually read, which cost the control panel 733 spurious minimum-distance failures until the frame was rescreened against its own search receptor.
```

**Replacement:**

```latex
They ran under the shared-receptor policy, so every pose was screened against a single file per frame and tool and the exactness check never applied. On the START-Fr0 frame that file was the deposited structure rather than the energy-minimised one that the AutoDock search alone read (Appendix~\ref{orai1-receptor-model}), which cost the control panel 733 spurious minimum-distance failures until the AutoDock poses of that frame were rescreened against their own search receptor.
```

*Why:* Merged. Draft A's 'per frame and tool' removes the contradiction with the next sentence, which already says DiffDock and EquiBind were screened against their own files. Draft B's 'the AutoDock poses of that frame' is the accurate subject, since the current _orai_matched_root CSVs screen AutoDock Fr0 rows against _orai_benchmark_staging_fr0corrected and DiffDock/EquiBind rows against Data/Receptors. Both drafts add the cross-reference that resolves the former contradiction with :560.

#### 9. `thesis_latex/body_appendix_short.tex:572`

**Current (verbatim):**

```latex
Every statement in this work that depends on the receptor therefore rests on the coordinates as delivered, and the four frames are treated as four unlabelled conformational samples of one membrane simulation rather than as timed points along a characterised trajectory.
```

**Replacement:**

```latex
Every statement in this work that depends on the receptor therefore rests on the coordinates as delivered, the AutoDock copy of Fr0 excepted, and the four frames are treated as four unlabelled conformational samples of one membrane simulation rather than as timed points along a characterised trajectory.
```

*Why:* Draft A wording ('the AutoDock copy of Fr0 excepted' reads better than Draft B's 'the AutoDock Fr0 minimisation excepted' and matches the phrase now used at main:156). Without the exception the sentence is false for every AutoDock Fr0 number.

#### 10. `thesis_latex/body_appendix_short.tex:583`

**Current (verbatim):**

```latex
whereas Fr0 lies 5.88 to 6.02\angstrom{} from each of them.
```

**Replacement:**

```latex
whereas the delivered Fr0 lies 5.88 to 6.02\angstrom{} from each of them and its minimised AutoDock copy 5.84 to 5.98\angstrom{}.
```

*Why:* Draft A. The frame-spread figures and Table tab:appendix-frame-geometry were computed from Data/Receptors/Orai1WT-START-Fr0.pdb (thesis_expected_values.yaml maps Fr0 to that file), the DiffDock and EquiBind receptor. Recomputed Kabsch CA against the staged AutoDock Fr0 gives 5.86/5.84/5.98 A, so naming the file and giving the second range shows the conclusion is file-independent.

#### 11. `thesis_latex/body_appendix_short.tex:609`

**Current (verbatim):**

```latex
Fr0 against the three later frames is accordingly a controlled contrast between a minimised and a thermally sampled receptor geometry.
```

**Replacement:**

```latex
Fr0 against the three later frames is accordingly a contrast between a minimised and a thermally sampled receptor geometry, controlled for DiffDock and EquiBind and confounded for AutoDock by the preparation minimisation described above.
```

*Why:* Draft B only, with 'its second Fr0 minimisation' replaced by 'the preparation minimisation described above'. The factsheet does not verify that the delivered Fr0 is itself a minimised structure (it carries stretched covalent bonds up to 2.21 A), so 'second' would rest on the thesis's own unverified label. For the AutoDock arm the Fr0 receptor is 0.82 A heavy-atom from the delivered file, so the contrast is controlled only for the two learned tools. Table 14 and the 3.5 to 5.2 A compaction figures remain those of the delivered Fr0 and need no change.

## 9. Verifier corrections to the orchestrator's own inline claims

- The all-atom `_minimised.pdb` histidines are HID (HD1 only), not HIE. The converter adds HE2. Outcome unchanged.
- Atoms above 2 A: 149 (ILE CD mapped to CD1, OXT excluded), not 145.
- `Data/Receptors/pdbqt/Orai1WT-START-Fr0_mgl_tools.pdbqt` (Aug 28) is a further minimised prep, not an unminimised one. The only
  unminimised Fr0 PDBQT in the repository is the Meeko `retry1` file (1,333 residues, 10,321 heavy atoms, coordinate-identical to the
  delivered file), which no reported arm docked. This is the file `:560` currently describes (finding M3), so replacing `:560` fixes M3 too.

## 10. Post-application verification (2026-09-07, workflow `wf_27e0a488-7ac`)

Three verifiers re-checked the applied text (numbers against the files, house style and LaTeX, residual contradictions and M3 closure).
The residual-contradiction lens passed: no "at any stage" clause without ligand scope, no "same receptor file" claim, no 1,333 / 10,321 /
five-residue / 41-atom / "Forty-two of the 240" / "marginally conservative" text survives anywhere in the three `.tex` files, so **M3 is closed**.
The other two lenses required the following corrections, all applied and rebuilt (142 pages, no errors, no undefined references):

| line | problem | fix applied |
|---|---|---|
| `app:560` | Ser93 alone is not the trigger. Re-execution showed a copy with all six Ser93 C-beta atoms repaired still crashes; only a copy that also repairs stretched bonds in Lys203 (chains C, D, F), Leu248 (C), Ile251 (F), Phe253 (A, E) and His256 (C) converts. Ser93 chain A is merely the first offending atom in file order | "aborts ... because several covalent bonds in that frame lie beyond the bond-perception cutoff of its hydrogen builder. The first is the bond between C-alpha and C-beta of Ser93, stretched to 1.85 to 2.06 A in every chain, and some chains carry further stretched bonds in Lys203, Leu248, Ile251, Phe253 and His256." Opener "One qualification attaches to the Fr0 snapshot." dropped because `:558` now announces the paragraph |
| `app:489` | 733 is the superseded exh32 gnina control tree (3,030 poses, `obsolete/.../orai_benchmark_gnina/.../PRE_FR0.csv`). The reported exh128 control arm (3,079 Fr0 AutoDock poses, gnina variant) has **893** minimum-distance failures against the delivered Fr0 (`obsolete/posebusters_results/orai_benchmark_mgltools_exh128/dock/`) and **3** against its own receptor (current `_orai_matched_root`). Verified per variant, not pooled | "That mismatch produced 893 spurious minimum-distance failures among the 3,079 Fr0 AutoDock poses of the control panel. Rescreening those poses against their own search receptor left three." The carried-over ", so" chain and the 55-word sentence were split |
| `app:562` | "the two that occur most often within 4 A" is true only among the 97 displaced residues (overall the most frequent rank-1 contacts are Arg83, Leu95, Tyr80, Phe99, Arg91) | "Among the displaced residues, Tyr80 and Lys87 are the two that occur most often within 4 A of a rank-1 Fr0 AutoDock pose." |
| `app:562` | second sentence repeated the last sentence of `:560`; "felt at contact level" | sentence removed, "visible at the contact level" |
| `main:156` | "the AutoDock copy of Fr0 excepted, which was energy-minimised" attaches the relative clause to Fr0 and conflates the operation with the displacement | "The one exception is the AutoDock copy of Fr0, which was energy-minimised in place and whose C-alpha atoms moved by 0.64 A ... The poses of all three tools therefore already share one coordinate system" |
| `main:97` | "the preparation fallback" has no antecedent in the main text | "a preparation fallback" |

Section 8 above records the text as first applied; the lines in this table supersede it. Verified facts that were NOT in sections 1 to 4 but are
now in the thesis: the amber14-all / GBn2 force-field pair (from `prep_docking.py:264`), and the overall rank-1 contact ranking above.

### 10.1 Second verification round (single verifier over the five corrected lines)

Factual lens PASSED (893 and 3 recounted per variant; Ser93 1.850 to 2.062 A across chains; Lys203, Leu248, Ile251, Phe253 and His256 are
the only residues in the delivered Fr0 with an intra-residue C-C or C-N bond above 1.791 A). Wording items applied and rebuilt (142 pages, clean):

| line | change |
|---|---|
| `app:489` | "deposited structure" became "delivered structure" (the Orai frames were never deposited; Appendix D says "delivered" throughout). The failure sentence now reads "Against the delivered structure, 893 of the 3,079 Fr0 AutoDock poses of the control panel failed the minimum-distance check. Rescreening them against their own search receptor left three, and the reported verdicts are those of the rescreen." |
| `app:560` | "The first is the bond" became "The first in file order is the bond" |
| `app:562` | opener "That fallback" became "The PDBFixer and OpenMM fallback"; second sentence subject became "The minimised search receptor" to avoid repeating :560; "the difference behind" became "the displacement behind" |
| `main:95` | **beyond the eleven:** "contrasts one energy-minimised geometry with three thermally sampled ones" became "contrasts one pre-dynamics starting geometry with three thermally sampled ones" |
| `app:611` | **beyond the eleven:** "a contrast between a minimised and a thermally sampled receptor geometry" became "a contrast between a starting and a thermally sampled receptor geometry" |

The last two were changed because the new `:560` text documents C-alpha to C-beta bonds of 1.85 to 2.06 A in every chain and a 2.21 A
C-alpha to C bond in His256, which no energy minimisation would leave in place, so the thesis can no longer call the delivered Fr0
"energy-minimised". "Starting geometry" is what the file's START label and the recovered MD provenance support. The general sentence at
`:611` about energy-minimised models being closer to the training distribution is untouched. Left as pre-existing text outside the M1
scope: the ", so" chain in the last sentence of `:489` and the "Float too large" warning at input line 480.
