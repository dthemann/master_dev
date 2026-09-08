# Answers to the four MD-frame questions, from the delivered files

**Date:** 2026-09-03
**Source files:** `/mnt/c/.../Orai1-Docking-MScProjekt-DMann/Orai/MDFiles/`
(`hORAI1_WT_ARN.psf`, `Orai1_WT_bar3.dcd`, `Orai1_WT_fix1.dcd`, the four snapshot PDBs, `FH Wien.7z`)
**Scripts and data:** `Scripts/Analysis/md_provenance_probe/` (four scripts, three JSON result files)
**Figure:** `Scripts/Analysis/MD_FRAME_PROVENANCE_2026-09-03.png`

---

## Headline

Three of the four questions can now be answered from the delivered files alone. The fourth
still needs Heinrich. The most important result is a correction, not a confirmation:

> **The delivered production-stage trajectory spans roughly a third of a nanosecond, not tens of
> nanoseconds.** Frame 300 is at approximately 0.2 nanoseconds. Writing "Frame 300 = 30 ns"
> would have been wrong by a factor of about 150.

The DCD headers cannot tell you this because they were rewritten by VMD and carry placeholder
timing. The interval was measured instead from how far water molecules move between saved frames.

---

## Question 1 — Trajectory length and frame interval

### What the headers say, and why you cannot use them

Both DCD files were re-written by the VMD molfile plugin on 23 May 2025. The header carries:

| field | value | meaning |
|---|---|---|
| `remarks` | `Created by DCD plugin ... Created 23 May, 2025` | not the original NAMD/CHARMM header |
| `istart` | 0 | placeholder |
| `nsavc` | 1 | placeholder |
| `delta` | 1.0 | placeholder, exactly the float `1.0` |

`istart=0, nsavc=1, delta=1.0` is what VMD writes when it has no timestep information to pass on.
The earlier note that the header reports `dt = 0.049 ps` is simply `delta` multiplied by the
AKMA-to-picosecond factor; it is an artefact of the placeholder, not a measurement. **The original
timestep and save frequency are not present in the delivered files.**

### What the trajectory itself says

The interval between saved frames was measured from the mean squared displacement of water oxygen
atoms, using the well-characterised self-diffusion coefficient of the CHARMM TIP3 water model
(D approximately 5.5e-5 cm^2/s). In the diffusive regime, mean squared displacement = 6 D t.

**Water near the membrane must be excluded.** Interfacial water diffuses more slowly than bulk
water, and including it biases the interval upward. The slope was therefore measured against
distance from the bilayer, which in `bar3` spans z from -35.5 to +24.2 Angstrom:

| water selection | molecules | slope (square Angstrom per frame) | implied interval |
|---|---|---|---|
| all sampled | 12,000 | 1.880 | 0.570 ps |
| more than 30 Angstrom out | 8,945 | 2.103 | 0.637 ps |
| more than 45 Angstrom out | 5,024 | 2.185 | **0.662 ps** |
| more than 52 Angstrom out | 3,077 | 2.180 | 0.661 ps |
| more than 58 Angstrom out | 1,369 | 2.166 | 0.657 ps |

The slope plateaus beyond 45 Angstrom, which is genuine bulk. **Take 0.66 ps per frame.**

The same measurement repeated on three separate windows of `bar3` and two disjoint sets of water
molecules agrees to within 4%, and on `fix1` it gives 0.074 to 0.077 ps per frame once the initial
heating ramp is past.

**Consequently:**

| quantity | best estimate | plausible range |
|---|---|---|
| `bar3` interval per frame | **0.66 ps** | 0.55 to 0.80 ps |
| `bar3` total span (500 frames) | **0.33 ns** | 0.28 to 0.40 ns |
| `fix1` interval per frame | 0.075 ps | |
| `fix1` total span (501 frames) | 0.038 ns | |
| Fr300 | approximately 0.20 ns into `bar3` | |
| Fr400 | approximately 0.26 ns into `bar3` | |
| Fr499 | approximately 0.33 ns into `bar3`, the last frame | |
| Fr300 to Fr400 separation | approximately 66 ps | |
| Fr300 to Fr499 separation | approximately 131 ps | |

The range comes from the spread of published self-diffusion coefficients for CHARMM TIP3 water,
about 5.3 to 6.0e-5 cm^2/s, plus the fact that the run temperature is not recorded.

### Why this is trustworthy

- **It cannot be an artefact of striding.** The measurement is of the displacement between the
  frames you actually have. If a long run had been subsampled to 500 frames, consecutive frames
  would be hundreds of picoseconds apart and water would have moved tens of Angstrom, aliasing
  against the 118 Angstrom box. It moves 1.5 Angstrom.
- **A second tracer agrees.** Potassium and chloride ions give a mean squared displacement 0.55
  times that of water at the same lag, which is the expected ratio for these ions in TIP3 water.
- **The error bar is small compared with the claim.** Even a fivefold error in the assumed water
  diffusion coefficient would put the run at 1.7 ns. To reach "Frame 300 = 30 ns" the water would
  have to be diffusing about 150 times slower than the TIP3 model does.
- **The protein motion is consistent.** Consecutive frames differ by 0.45 Angstrom alpha-carbon
  root-mean-square deviation and frames 100 apart by 1.0 Angstrom, which is the growth expected
  over tens of picoseconds, not tens of nanoseconds.
- **It was independently reproduced.** Three separate re-derivations from the raw trajectory files,
  written without sight of the code above, returned 0.54 to 0.56 picoseconds per frame using all
  sampled water. That matches the 0.57 obtained here on the same basis. The bulk-water correction
  then raises the figure to 0.66.

**Assumption to state if you use this:** the run is at ordinary simulation temperature, near 300 K.
Nothing in the coordinates fixes the temperature directly, though the amplitude of the atomic
fluctuations is consistent with it.

### The one alternative reading, and why it does not rescue the nanosecond story

A reasonable objection is that `bar3` might be a **500-frame window cut out of a much longer
production run** at its native save frequency, rather than a short run in its own right. The
measurement cannot distinguish those two directly, because it measures the spacing between the
frames you have, not their offset within any parent run.

It does not matter for the question you asked, and here is why.

- **Under either reading the three snapshots span about 130 picoseconds.** Fr300 to Fr499 is 199
  frames at 0.66 ps. So they are not independent conformational samples either way.
- **The window reading is the less likely of the two.** If `bar3` were cut from deep inside an
  equilibrated production run, the protein root-mean-square deviation against the window's own
  first frame would jump quickly to a thermal floor near 0.5 Angstrom and then fluctuate flat.
  What it actually does is climb monotonically along a relaxation curve to 1.5 Angstrom over the
  whole file. That is a system still settling, which is what an equilibration stage looks like.
- **The save frequency argues the same way.** Saving roughly every 0.7 ps is an equilibration-stage
  setting. This group's own deposited production input for a comparable system uses
  `dcdfreq 10000` with a 2 fs timestep, that is one frame per 20 picoseconds, roughly 40 times
  coarser. A 300 ns production run saved at 0.7 ps would be well over 400,000 frames.

For context, published all-atom Orai hexamer production runs cluster between 100 and 550
nanoseconds per replica. Whatever `bar3` is, it is not one of those.

### Is Frame 300 after the equilibration?

Yes, in the sense that matters, and this can be shown without knowing the timestep at all.

- The alpha-carbon root-mean-square deviation of the protein against frame 0 of `bar3` rises
  steeply, then flattens from about frame 250 onwards at 1.44 to 1.59 Angstrom. Fr300, Fr400 and
  Fr499 all lie in that plateau. See panel A of the figure.
- The box volume in `bar3` is already stationary at its own frame 0: mean 1750.2 cubic nanometres
  over the first 50 frames against 1749.8 over the last 200, with a standard deviation of 1.7. The
  density had finished equilibrating before this file starts.

**But be careful how you word it.** "After the equilibration" is true relative to the start of
this file. It is not a claim that the system had been equilibrated for nanoseconds beforehand,
and the file itself does not show how much simulation preceded it.

---

## Question 2 — How the frames were selected

You do not need Heinrich to answer this, and the systematic reading is now essentially certain.

- `bar3` contains exactly **500 frames**, indices 0 to 499. **Fr499 is the final frame of the file.**
- Fr300 and Fr400 are round numbers exactly 100 frames apart, and Fr400 to Fr499 is 99 frames.
- Fr0 is not from this file at all. It is frame 0 of `fix1`, a different and earlier stage.

A random draw producing 300, 400 and the last index is not a credible account. The honest and
defensible description is a **systematic selection: the final frame of the trajectory and two
earlier frames at equal spacing**, plus a separate pre-dynamics starting structure.

Suggested wording for the methods section:

> Drei Konformationen wurden systematisch aus dem gelieferten Produktionsabschnitt entnommen,
> naemlich der letzte Frame sowie zwei um jeweils rund 100 Frames frueher liegende Frames. Die
> vierte Struktur stammt nicht aus diesem Abschnitt, sondern ist die minimierte Ausgangsstruktur
> vor Beginn der Dynamik.

This is methodologically just as defensible as a random draw, and it is what the files show.

### The selection is good; the trajectory is short

This distinction is worth getting right, because it lets you defend the choice rather than
apologise for it. Every pair of the 100 frames sampled at stride five was superposed and compared:

| set | mean pairwise deviation | maximum |
|---|---|---|
| all of `bar3` | 1.07 Angstrom | 1.62 Angstrom |
| the plateau only, frames 250 onwards | 0.88 Angstrom | 1.19 Angstrom |
| **the three chosen frames** | **1.03 Angstrom** | **1.17 Angstrom** |

**The three frames span 98.7% of the conformational spread present in the equilibrated portion of
the trajectory**, and 72.5% of the spread in the file as a whole. Within what was available, the
selection is close to optimal. You can say that.

**What you must nonetheless disclose** is the other half. The frames are about 66 picoseconds
apart, their pairwise deviations sit barely above the 0.45 Angstrom adjacent-frame thermal floor,
and the plateau they sample is itself only about 165 picoseconds long. So the receptor ensemble is
narrow, and it is narrow because the delivered trajectory is short, not because the frames were
badly chosen. That is a limitation of the input data, and stating it that way is both accurate and
better for you than the alternative.

---

## Question 3 — The homology model and what Fr0 is

### What the delivered model is, measured

| property | value |
|---|---|
| subunits | 6, segments MONA to MONF, identical sequence |
| residues per subunit | 66 to 288, **223 residues, no gaps** |
| termini | acetyl patch at the N-terminus, methyl patch at the C-terminus |
| numbering | human Orai1; Arg91, Val102, Glu106, His113, Asp110/112/114, Tyr115 all sit at their human positions |
| histidine states | HSE and HSP present, so the states were assigned, not defaulted |
| system | CHARMM-GUI v3.7, built 6 April 2022, "assemble previously generated components" |

### Is Fr0 the model unchanged, or minimised?

**Neither exactly. Fr0 is the minimised, membrane-embedded structure at the instant dynamics
began, before the protein was ever allowed to move.** Three measurements say so.

1. **The protein is held fixed throughout `fix1`.** Its alpha-carbon root-mean-square deviation
   against frame 0 reaches 0.029 Angstrom and stays there for all 501 frames. That is not a
   restrained protein, that is a fixed one.
2. **The box is exactly constant** at 129.200 x 128.700 x 136.700 Angstrom for every frame, with
   a standard deviation of zero. Constant volume, no barostat.
3. **`fix1` starts cold and warms up.** The per-frame water displacement ramps monotonically from
   0.012 to 0.19 square Angstrom over the first fifteen frames, then settles at 0.23. That ramp is
   a heating curve from near-zero kinetic energy. A structure at the start of a heating run is a
   minimised structure.

Supporting: Fr0's backbone heavy-atom bonds are systematically 0.019 Angstrom shorter than in the
three production frames (N-CA 1.4391 against 1.457), which is what a minimised structure looks
like next to a thermally expanded one.

**So Fr0 is not "the equilibrated starting structure"** as the thesis currently calls it at
`body_main_short.tex:591`. It is the pre-equilibration structure. This finally explains, with a
mechanism, why Fr0 sits 5.88 to 6.02 Angstrom away from the other three.

### Which publication the model comes from

**Your guess was right. The model is from Frischauf et al. 2015, which is already your ref038.**

> I. Frischauf, V. Zayats, M. Deix, A. Hochreiter, I. Jardin, M. Muik, B. Lackner, B. Svobodova,
> T. Pammer, M. Litvinukova, A. A. Sridhar, I. Derler, I. Bogeski, C. Romanin, R. H. Ettrich,
> R. Schindl. "A calcium-accumulating region, CAR, in the channel Orai1 enhances Ca2+ permeation
> and SOCE-induced gene transcription." *Sci. Signal.* 8(408):ra131, 2015.
> doi 10.1126/scisignal.aab1901

That paper builds the model rather than citing an earlier one. From its methods:

> "The published *D. melanogaster* Orai1 crystal structure (PDB ID: 4HKR) was used to construct a
> full atomistic model of human Orai1. Two extracellular loops and one intracellular loop are
> unresolved in the *Drosophila* crystal structure and were modeled by retrieving loop
> conformations bridging two sets of anchor residues from the Protein Data Bank."

> "Homology modelling was performed in Yasara, using the standard macro hm_built.mcr with the
> FixModelRes option to constrain conserved residues."

**The model is deposited, and the identification was confirmed by superposition rather than
inferred.** ModelArchive accession **`ma-akdjp`** (doi 10.5452/ma-akdjp), "Structure of human Orai1
based on homology with dOrai", deposited by R. H. Ettrich, written by YASARA on 11 November 2015.

The deposited coordinates were downloaded and superposed on Fr0:

| test | result |
|---|---|
| chains and residue range | six chains, residues 66 to 288, in both |
| residue-name mismatches | **zero** |
| heavy atoms per chain | 1,727 in both |
| alpha-carbon superposition, 1,338 pairs | **0.122 Angstrom** |
| all common heavy atoms, 10,296 | 1.62 Angstrom, side chains rebuilt by CHARMM-GUI |

The deposit also carries an idiosyncrasy that no generic 4HKR-based model would share. It breaks
sixfold symmetry at the M4 extension, with subunits A, C and E differing from B, D and F by a mean
of 8.5 Angstrom in that region. **Fr0 reproduces that same asymmetry.** The sequence is separately
a byte-perfect match to UniProt Q96D31 over residues 66 to 288.

The match is to Fr0, which is the pre-dynamics structure. The three later frames descend from it
through the simulation and sit about 5.9 Angstrom away, as expected.

So the citation chain is:

| role | paper | your ref |
|---|---|---|
| template crystal structure, *Drosophila* Orai, PDB 4HKR | Hou et al., Science 2012 | ref037 |
| **the human Orai1 homology model itself** | **Frischauf et al., Sci. Signal. 2015** | **ref038** |
| later re-use of the same model | Frischauf 2017, Tiffner 2021, Hopl 2024, Najjar 2025 | ref072 for Najjar |

**Do not cite Hou et al. 2012 for the model.** It is the template, and it is *Drosophila*. There are
still no experimental coordinates for the human Orai1 channel in the Protein Data Bank.

One precision worth keeping. That the template was 4HKR specifically rests on the Frischauf methods
text, which names it. The coordinates cannot establish it on their own, because the model sits
0.20 Angstrom from 4HKR and 0.23 Angstrom from its companion entry 4HKS while those two are only
0.15 Angstrom apart. Cite the paper's statement, not a structural inference.

Two caveats worth carrying into the appendix.

- **A substantial part of the model has no template.** 4HKR leaves two loops unresolved, and
  Frischauf et al. rebuilt them from Protein Data Bank fragments. Independent estimates of the
  templated fraction ranged from 59% to 71% of each subunit depending on how the alignment is
  registered, so quote it qualitatively rather than as a number. What matters is which parts are
  de novo: the calcium-accumulating region loop and most of the M3-M4 loop, which are two of the
  regions the Orai1 chapter reasons about.
- **The 73% identity figure is the authors' own and is not reproducible as stated.** It comes
  verbatim from Hou et al. 2012, which never defines the region it applies to. Recomputed over the
  four annotated transmembrane helices it is 64 to 67%, and it is strongly graded: about 94% over
  M1, but only 33% over M4. Quote it as their figure, and do not lean on it to justify treating the
  fly structure as a proxy in an M4-relevant context.

---

## Question 4 — A citable reference for the simulation setup

**Najjar et al. 2025 is definitively the wrong paper, and no substitute was found.** The Najjar
deposit was downloaded from Zenodo (doi 10.5281/zenodo.15021906) and its topology compared with
yours:

| | your `hORAI1_WT_ARN.psf` | Najjar et al. 2025 |
|---|---|---|
| CHARMM-GUI version | **v3.7, 6 April 2022** | **V2.0, 24 April 2020** |
| atoms | 182,692 | 304,917 |
| segment names | MONA to MONF | PROA to PROF |
| residue numbering | 66 to 288, native | 2 to 237, offset by 60 |
| lipids | 426: cholesterol 213, DSM 71, DLPC 71, DLPE 71 | 694: DDPC 260, DLPS 128, DLPE 106, PSM 97, DMPI 62, DLPG 41 |
| cholesterol | **50 mol%** | **none** |
| ions | 109 K+, 150 Cl-, no calcium | 115 Ca2+, 5 Cl- |
| box | 129 x 129 x 137 Angstrom | 156 x 156 x 121 Angstrom |

Eight independent mismatches. Citing Najjar et al. 2025 for your setup would be a misattribution
that an examiner with the paper open could catch.

**No other published paper was found to match.** Every Orai membrane simulation in this lineage
was checked against its primary text: Frischauf 2015 and 2017 use POPC, Hopl 2024 uses 394 POPC,
Najjar 2025 uses the Dawaliby composition. A wider literature search found no simulation paper
combining these lipids at this ratio.

**State this as a bounded result, not a universal negative.** The search covered the nine
Orai and CRAC papers co-authored by Heinrich Krobath between 2020 and 2026 and a general literature
sweep. It cannot exclude a paper by someone else in the wider Linz group, or an unpublished
manuscript. The defensible sentence is that no published description of this system was found, so
it is treated as unpublished pending confirmation.

**What your membrane actually is**, confirmed from the atom composition of one residue of each type:

| residue | formula | identity | count |
|---|---|---|---|
| CHL1 | C27H46O | cholesterol | 213 (50.0 mol%) |
| DSM | C47H97N2O6P | N-lignoceroyl dihydrosphingomyelin, d18:0/24:0 | 71 |
| DLPC | C32H64NO8P | 1,2-dilauroyl phosphatidylcholine, 12:0/12:0 | 71 |
| DLPE | C29H58NO8P | 1,2-dilauroyl phosphatidylethanolamine, 12:0/12:0 | 71 |

The system is exactly charge-neutral. The DSM assignment was confirmed by parsing the CHARMM-GUI
conformer libraries directly: DSM has 153 atoms with a fully saturated sphingoid base, against 151
for the singly unsaturated LSM and 127 for the common PSM.

A membrane of 50 mol% cholesterol with a very-long-chain saturated dihydro-sphingomyelin alongside
two 12-carbon phospholipids is unusual, and the chain-length mismatch between the 24-carbon
sphingomyelin and the 12-carbon phospholipids is not a standard raft mimic. Do not assume it was
designed that way. An equally good reading is that it reflects CHARMM-GUI list order and defaults
at build time. Either way, Heinrich should be able to say which.

**How to cite it until he answers.** Describe the system explicitly from the topology file and
attribute it as unpublished data provided by H. Krobath, Johannes Kepler University Linz, personal
communication. Cite Najjar et al. 2025 only for the group's Orai1 gating findings, which is what
Appendix D already uses it for.

**One thing that did check out:** Heinrich Krobath is author 23 of 25 on Najjar et al. 2025, and
its author contributions read "F.H., H.K., and T. Renger performed and supervised MD simulations."
He is the right person to ask.

---

## What this changes in the thesis

Five items, in descending order of how much they matter.

**1. `body_main_short.tex:591` can be read as sourcing the frames to the wrong publications.** It
reads "The contributing Orai group supplied all four from its structural and molecular-dynamics
work \cite{ref038}, \cite{ref072}." That is defensible if it means only that the group whose work
these papers represent supplied the files. It is false under the more natural reading that the
frames come from the systems those papers describe. The *model* is indeed from ref038, and that
half is now confirmed. The *trajectory* is from neither: ref038 simulated 512 POPC lipids in
GROMACS with OPLS-AA, and ref072 simulated a six-component cholesterol-free bilayer in NAMD.
Split the sentence so that ref038 carries the model and the trajectory is attributed explicitly as
unpublished data from the contributing group.

**2. The same line calls Fr0 "the equilibrated starting structure".** It is the *pre*-equilibration
structure: minimised, membrane-embedded, taken at the instant dynamics began. Note that
`body_main_short.tex:95` already gets this right, calling it "one energy-minimised geometry"
contrasted with "three thermally sampled ones". That sentence is now measurement-backed. Only
line 591 needs the word changed.

**3. A new limitation belongs in the Limitations chapter.** The three production snapshots span
about 130 picoseconds. Their pairwise alpha-carbon deviations of 0.93 to 1.19 Angstrom sit barely
above the 0.45 Angstrom adjacent-frame thermal floor. The thesis argues at
`body_main_short.tex:95` that four receptor states are used because "conformational variation can
therefore enter only through separate receptor structures". That argument is sound, but the
variation actually supplied is far narrower than "four molecular-dynamics snapshots" implies, and
the honest version of the claim says so.

**4. Appendix D's provenance paragraph at `body_appendix_short.tex:572` is now largely answerable.**
It currently declares unrecorded: the modelling program, the template identifiers, the force field,
the lipid composition, the solvent and ion treatment, the trajectory length, the replicate identity
and the frame times. Of those, the modelling program (YASARA), the template (4HKR), the builder
(CHARMM-GUI v3.7, 6 April 2022), the force field family (CHARMM with CMAP), the full lipid
composition, the solvent and ion content, and the frame spacing are now all established. Replicate
identity and the absolute position of `bar3` within any longer run remain genuinely unknown.

**5. The 73% identity figure is load-bearing in a place it should not be.**
`body_appendix_short.tex:581` locates the modulator binding region at the M3 interface with the
non-pore-lining helices. Sequence conservation between the fly template and human Orai1 is graded:
about 94% over M1, but only about 33% over M4. The model is also roughly 29% de novo per subunit,
and that de novo fraction includes the calcium-accumulating region loop and most of the M3-M4 loop.
Both of the regions the Orai1 chapter reasons about are therefore the least well determined parts
of the model. The appendix already says loop conformations are the least well determined part of a
homology model; it can now say so quantitatively.

### Other sentences a full sweep flagged

A separate pass read the appendix and methods against the measured facts. These are the further
hits worth acting on, in addition to the five above.

| location | the sentence | why it now fails |
|---|---|---|
| `body_appendix_short.tex:572` | "Nothing further was supplied." | The topology file was supplied, and it answers most of the list that follows. |
| `body_appendix_short.tex:583` | "Fr0 is thus not a member of the ensemble the later frames sample, which is consistent with its role as the starting structure" | Correct, but the reason is now known and mechanical rather than inferred from a label. |
| `body_appendix_short.tex:611` | "The frames may also be statistically dependent, and the simulation times that would settle the question are among the metadata the section above records as unavailable." | The times are now estimated, and they settle it. The frames are dependent. |
| `body_main_short.tex:541` | "the twelve experimental units represent only three chemotypes and one receptor trajectory" | Two stages, not one trajectory. Fr0 comes from a different run than the other three. |
| `body_main_short.tex:812` | "A stronger study would add independent trajectories, systematic frame selection..." | Frame selection already was systematic. What is missing is independent trajectories and a longer one. |

Two flags from that sweep are **not** supported and should be ignored. Nothing measured here bears
on the ligand optimisation provenance at `body_appendix_short.tex:793`. And "All four states derive
from one simulation of the wild-type channel without STIM1" at `:611` remains true at the level of
the simulated system, even though the frames come from two stages of it.

---

## What the files could not tell us

| item | status |
|---|---|
| integration timestep and `dcdfreq` as actually configured | **not recoverable**, only inferable from the measured interval |
| what ran between `fix1` and `bar3` | **not delivered**. The protein moves 5.76 Angstrom and the box shrinks from 2273 to 1750 cubic nanometres between them, so at least one unrestrained equilibration stage is missing. The name `bar3` implies `bar1` and `bar2`. |
| whether a longer production run exists beyond `bar3` | **unknown, and this is now the key question for Heinrich** |
| the simulation temperature | not in the coordinates |
| replicate identity | not in the coordinates |
| homology-model template identifiers | **resolved**: 4HKR, via Frischauf 2015 and ModelArchive `ma-akdjp` |
| which publication describes this membrane and these parameters | **none exists**; treat as unpublished |

The archive `FH Wien.7z` was opened and contains only the same three files. A search of the whole
`Master Project` tree found exactly six trajectory-related files: the three delivered ones and an
identical backup copy of them. **There are no NAMD configuration files, logs, restart files or
extended-system files anywhere in the delivered material.** Nothing further can be recovered
without asking.

---

## Appendix A — Drop-in replacement text

Nothing in the thesis has been edited. These are drafts for when you decide to act, written in the
thesis register with no mid-sentence colons, semicolons or clause dashes.

**For `body_main_short.tex:591`**, replacing the two sentences beginning "These are the equilibrated
starting structure Fr0" and "The contributing Orai group supplied all four":

```latex
These are the energy-minimised starting structure Fr0 and the later frames Fr300, Fr400 and Fr499.
The homology model is the one published by the contributing Orai group \cite{ref038}. The
molecular-dynamics frames come from an unpublished membrane simulation of that model supplied by
the same group.
```

**For `body_appendix_short.tex:572`**, replacing the provenance paragraph. This version states what
is known and keeps the two genuine gaps explicit:

```latex
The provenance of the model is recorded here as far as the delivered files establish it. The
homology model is the human Orai1 model of Frischauf et al., built in YASARA on the
\emph{Drosophila} Orai crystal structure 4HKR and deposited in ModelArchive under accession
ma-akdjp \cite{ref038}, \cite{ref037}. The delivered coordinates match that deposit, namely six
chains of residues 66 to 288. The simulation system was assembled with CHARMM-GUI version 3.7 on
6 April 2022 and carries a CHARMM topology with CMAP. The bilayer holds 426 lipids, of which 213
are cholesterol, 71 a long-chain saturated sphingomyelin, 71 dilauroyl phosphatidylcholine and 71
dilauroyl phosphatidylethanolamine, so cholesterol accounts for half the membrane. The system is
solvated with 40,123 water molecules and neutralised with 109 potassium and 150 chloride ions, and
it contains no calcium. Two trajectory stages were supplied. Fr0 is the first frame of a
fixed-protein constant-volume stage, and Fr300, Fr400 and Fr499 are frames of a later
constant-pressure stage whose density is already stationary at its first frame. No publication
describes this particular system, and it is therefore cited as unpublished data from the
contributing group. Two items remain unrecorded. The integration timestep and the trajectory save
frequency were removed when the trajectory files were rewritten, and the position of the supplied
stage within any longer simulation is not known. The interval between saved frames was instead
estimated from the displacement of bulk water at approximately 0.7 picoseconds, which places the
three later frames within a window of roughly 130 picoseconds.
```

**For the Limitations chapter**, a new sentence:

```latex
The receptor ensemble is also narrower than the term ensemble suggests. The three
molecular-dynamics frames lie within about 130 picoseconds of one another, and their pairwise
backbone deviations of 0.93 to 1.19\angstrom{} sit close to the 0.45\angstrom{} deviation between
adjacent frames. They do span 98.7\% of the pairwise spread present in the equilibrated part of the
supplied trajectory, so the limitation is the length of that trajectory rather than the choice of
frames.
```

---

## Appendix B — How to reproduce

The DCD files are plain Fortran-record binaries. Each frame is
`56 + 3 * (8 + 4 * natoms)` bytes with `natoms = 182692`, so the protein block (the first 21,026
atoms) and the water block can be read by seeking, without loading 1.1 GB per pass.

```
frame_size = 56 + 3 * (8 + 4 * 182692) = 2192384
header     = filesize - n_frames * frame_size        # 276 bytes for both files
unit cell  = 6 doubles at (header + frame*frame_size + 4), order a, gamma, b, beta, alpha, c
x block    = header + frame*frame_size + 56 + 4      # then y, z at +730776, +1461552
```

Segment ranges in `hORAI1_WT_ARN.psf`:

| segment | atoms |
|---|---|
| MONA-MONF (protein) | 1 to 21,026 |
| MEMB (426 lipids) | 21,027 to 62,064 |
| TIP3 (40,123 waters) | 62,065 to 182,433 |
| IONS (109 K+, 150 Cl-) | 182,434 to 182,692 |
