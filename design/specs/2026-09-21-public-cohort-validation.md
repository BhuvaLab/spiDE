# External validation of spiDE on four public single-cell spatial cohorts

**Status:** draft spec, 2026-09-21
**Scope:** calibration and power/recall for spiDE on public data, for the paper
**Code:** `research/public/`  ·  **Data:** `/scratch/project_mnt/S0249/R_projects/spiDE_public_data/`

---

## 0. State as of 2026-09-21

| step | state |
|---|---|
| download (`download.sh`) | **done** — ~8 GB, only the files the pipeline reads |
| sample sheets from GEO SOFT | **done** — `research/public/metadata/*.csv`, tracked |
| schema + CosMx reader (`00_common.R`) | **done** |
| builders `01`–`04` | **done**; running as array `28723330` |
| GSE282639 built | **done** — 974 genes x 104,680 cells, 28 subjects |
| annotation (`05_annotate.R`) | written, unrun |
| QC / bandwidth grid (`06_qc_report.R`) | written, unrun |
| sizing pilot (§10 step 0) | not started — **the next thing to do** |
| any spiDE run | not started |

Two implementation facts worth carrying forward. The interactive allocation is
**32 GB / 4 cores**, and running the builds in it OOM-killed all three large
ones; every step from here goes through `sbatch` with measured memory. And
`.readCosMx()` now builds the count matrix in cell chunks — densifying a CosMx
6K slide costs 12.8 GB for the matrix plus as much again for the `data.table`
it is copied from, which is what the OOM actually was.

**Pilot observation, 2026-09-21 (GSE282639, 28 patients).** The Satterthwaite
reference df on some `ResponseNiche` columns come out at 25,388 — cell-level,
not patient-level. This is the *documented* behaviour (`CLAUDE.md`:
satterthwaite "separates between-sample contrasts from within-sample ones —
`ResponseNiche`: larger df, more power", traded against "a mild liberal drift
at larger S"), not a defect, and whether the drift bites on these cohorts is
exactly what the null arms measure. Noted here because a reader of the results
will ask. Separately: `research/fdr-ordering/pkgfinal13_all_fit.sbatch` claims
its arm "carries the finished package: every condition-bearing column bounded
at its own compartment's patient count, and `df.method` defaulting to
`between`". **No such bounding exists** in 0.99.21, in the live working tree,
or in `FINDINGS.md`; the comment describes something unmerged. This study runs
the shipped defaults, which is the point of it.

One open item against §5: the live YTMA arm (`pkgfinal13`) is running against a
working tree whose `df.method` bounding is in flight. Every arm here must pin
`SPIDE_PKG` to a frozen snapshot taken at submission, and the report must name
it — otherwise this analysis inherits a mid-flight change it did not intend to
test.

---

## 1. What this analysis has to establish

Every calibration and power claim spiDE currently makes rests on one cohort
(YTMA CosMx WTA, 55 patients) plus a synthetic simulator. That is a single
point of failure for the paper: the shuffle-null result, the nested-intercept
fix, the QL scale and the per-gene calibration were all measured on the same
tissue, the same platform, the same panel and the same annotation. A reviewer
is entitled to ask whether any of it generalises.

This analysis must support four claims, in decreasing order of importance.

**C1 — Calibration generalises.** Under a null where the tested effect is zero
by construction, spiDE's test statistic has unit spread and its FDR control is
honest, on tissues, platforms, panel sizes and cohort sizes the estimator was
not tuned on. This is the claim the paper cannot do without.

**C2 — The nested (sample × cell type) intercept is necessary off YTMA.**
CLAUDE.md records that the block was justified on YTMA shuffle grids and that
the synthetic benchmark *cannot* show its benefit, because the simulator plants
no between-sample composition effect. Four independent real cohorts either
reproduce that finding or they do not.

**C3 — The estimator has power on real data, quantified.** Real cohorts have no
ground truth, so power is measured by injection into real covariance structure
and corroborated by cross-dataset and cross-platform replication.

**C4 — ~~The panel-size regime~~ WITHDRAWN 2026-09-21.** This claim needed
GSE282193's 6,175-gene panel and that cohort is parked (§2.1). The three
remaining cohorts read 343 / 950 / 974 genes, too narrow a range to say
anything about panel depth. Not claimed.

Anything this analysis **does not** establish is listed in §12 and must not be
claimed.

---

## 2. The four cohorts, and what each one is for

| | GSE250346 | GSE282639 | GSE289194 | ~~GSE282193~~ |
|---|---|---|---|---|
| tissue | lung | lung | lymph node / extranodal | kidney allograft |
| platform | Xenium | CosMx | CosMx | CosMx |
| panel | 343 custom | 1,000-plex (974 after QC) | 1,000-plex | **6,175** |
| patients (S) | **35** (26 PF / 9 control) | **28** (15 IPF / 13 control) | **74** (37 nodal / 37 extranodal) | 23 pieces, subjects unresolved |
| sections | 45 TMA cores | 43 FOVs | 74 TMA cores | 23 |
| cells (post-QC) | ~1.6 M | 104,680 | ~700 k (est.) | ~700 k (est.) |
| second condition | — | **MUC5B GG 16 / GT 12** | — | — |
| cell types | 47, **deposited** | transferred | 19, authors' own | InSituType |
| published niches | CNiche ×12, TNiche ×12 | — | 7 niches | — |

Each is in the set for a reason that no other member covers.

- **GSE250346** is the only cohort whose cell types, niches and pathologist
  annotations are all deposited, so it is the one where a power claim can be
  checked against something external to spiDE. It is also the largest by cells
  and the only Xenium member.
- **GSE282639** is the *same disease* as GSE250346 on the *other platform*,
  which makes cross-platform replication of PF calls possible; and it carries a
  second, architecture-independent binary (MUC5B rs35705950 genotype) that
  crosses the disease factor almost evenly (8/5/8/7). A three-way call driven
  by tissue architecture cannot be driven by a germline genotype, so the two
  contrasts on identical cells separate "the estimator found structure" from
  "the estimator found the condition".
- **GSE289194** is the scale member: 74 deposited cases with a 37/37 split, of
  which **44 are analysed** after the slide restriction below, which still puts
  spiDE in the regime CLAUDE.md records as its best (raw power 0.660 at S = 30,
  FDP 0.041 by S = 16). Nothing else public reaches it. The table above gives
  the deposited counts; the measured, analysed ones are 358,912 cells over 44
  patients (not the ~700 k estimated).
- **GSE282193** is the panel-depth member: 6,175 genes, the only public ≥5K-plex
  cohort with a balanced binary condition, and therefore the only test of C4.

### 2.1 GSE282193 is PARKED (2026-09-21)

Both available sources were exhausted, on instruction, before parking it.

**GEO metadata fields.** The finest mapping GEO carries is
`!Sample_description = raw data file: slideN.tar`, which places each GSM on a
*slide*, not on a FOV or a tissue piece:

| deposited tar | subjects declared |
|---|---|
| slide1 | IGF_001, IGF_003, DGF_002, DGF_003 (4) |
| **slide2 — NOT DEPOSITED** | IGF_002, DGF_001 (2) |
| slide3 | IGF_004–007, DGF_004–006 (7) |
| slide4 | IGF_008–010, DGF_007–009 (6) |

Two consequences. **`slide2.tar` is absent from the deposit**, so 2 of the 19
subjects have no data at all — the cohort is at most 17. And a slide carries
4–7 subjects, so slide-level assignment does not identify a patient.

**The geometric reconstruction does not close the gap.** Single-linkage on the
FOV origins at a 0.60 mm gap gives 7 / 7 / 9 pieces on slides 1 / 3 / 4 against
4 / 7 / 6 declared subjects. Only slide 3 matches. Slide 1 matches at a 1.00 mm
gap instead. A threshold that has to be re-tuned per slide to reach a known
answer is fitting, not recovery, and even a correct piece count would not say
which piece is `IGF_001`.

**Paper supplementary material.** There is none: no PubMed id is linked to the
series, and a literature search for the depositors (Wang, Lin, Zhang, Jen,
Than, Wang; UC Davis) on CosMx kidney delayed graft function returns no
matching publication.

**Decision: parked.** GSE282193 is removed from the run matrix. The built
object (`GSE282193_nocond_unannotated.rds`, 6,175 genes x 772,699 cells x 23
pieces) and `04_build_GSE282193.R` stay on disk and in the tree, and the
condition arms remain gated behind `SPIDE_GSE282193_MAP`, so the cohort is one
CSV away from running if the depositors ever supply the map.

**What parking costs.** Claim C4 — whether the per-gene heavy tail is a
function of panel depth — loses its only instrument. The three remaining
cohorts read 343 / 950 / 974 genes, all within a factor of three of each other
and all far below YTMA's 13,348. C4 is therefore **out of scope** for this
analysis and must not be claimed; §12 records it.

### 2.2 GSE289194 is restricted to three of its six slides (2026-09-21)

Its six slides are not one assay. All six carry 1,000 genes, but three carry
**208 negative probes against 10**, and measuring signal against background
rather than reading the headers separates them cleanly: per-gene signal divided
by per-probe background is 15.5 / 17.5 / 15.8 on F1 / F5 / F6 and
**0.9 / 3.8 / 5.6** on F2 / F3 / F4. On F2 the mean per-gene count *equals* the
mean per-negative-probe count, and F2 alone contributed 666,188 of 1.15 M cells
-- so the majority of the cohort came from its worst slide, and slide is nearly
nested in case.

The restriction costs 27 cases and leaves **44 patients and 358,912 cells** on
one configuration. `SPIDE_SLIDES` overrides it. Every number measured on this
cohort is at S = 44; quoting 74 anywhere downstream is a bug, and was one in
the first draft of `research/public/FINDINGS.md`.

---

## 3. Data processing — decisions already taken

Code: `research/public/R/00_common.R` (schema) and `01`–`04` (builders). Sample
sheets derived from the GEO SOFT records are tracked in
`research/public/metadata/*.csv`. Only the files the pipeline reads were
downloaded (`download.sh`); transcript-level and image files are ignored, which
is the difference between ~8 GB and ~450 GB.

Three harmonisations are substantive, not cosmetic.

**GSE289194 is restricted to slides F1/F5/F6, and this is the largest single
decision in the data work.** The six slides look like two panel revisions — 1,010
vs 1,208 features — and the first plan was to build on the 950-gene intersection.
Measuring the slides rather than reading their headers showed they are two
different assay configurations. Both carry 1,000 genes; the second group has 208
negative probes against 10, and a different signal-to-background regime entirely
(first 40k cells of each slide, 2026-09-21):

| slide | features | neg probes | median counts/cell | per-gene signal ÷ per-probe background |
|---|---|---|---|---|
| F1 | 1,010 | 10 | 557 | 15.5 |
| F5 | 1,010 | 10 | 667 | 17.5 |
| F6 | 1,010 | 10 | 492 | 15.8 |
| F2 | 1,208 | 208 | **22** | **0.9** |
| F3 | 1,208 | 208 | 184 | 3.8 |
| F4 | 1,208 | 208 | 126 | 5.6 |

On F2 the mean per-gene count *equals* the mean per-negative-probe count — the
cells carry no usable signal — and F2 alone contributed 666,188 of the 1.15 M
cells in the unrestricted build, i.e. the majority of the cohort came from its
worst slide. Since slide is nearly nested in case, keeping them mixes a
slide-sized depth artefact into the between-case term the nested intercept
exists to absorb, on a cohort whose entire purpose is calibration.

The restriction costs 27 cases and leaves **46 (25 nodal / 21 extranodal)** on one
configuration with identical feature sets, ~500–670 counts per cell and signal at
15× background. S = 46 is still the largest cohort in the study. `SPIDE_SLIDES`
overrides it, and the builder asserts that the retained slides share a feature set
so a future override that re-mixes configurations fails at build time.

A related bug is fixed with it: `nCount` and `neg_frac` are measured on the
slide's **full** gene set before the matrix is subset to the analysed panel.
Using the subset as the denominator inflated `neg_frac` wherever the panel was
trimmed, which is a slide-differential QC filter.

**Coordinates are microns.** `buildNiches()`'s `sigma` is a bandwidth in
coordinate units. CosMx reports global pixels, Xenium microns; mixing them would
compare a 30 µm niche against a 3.6 µm one. `PX_TO_UM = 0.12028` is *checked*
against the FOV pitch (expected ~509 µm) rather than trusted, and a mismatch is
an error, not a warning.

**`sample_id` and `section_id` are different fields with different jobs.**
`section_id` is a physical piece of tissue and is what `buildNiches()` splits on:
two cores of one patient sitting centimetres apart on a TMA are not one spatial
field, and a KDE evaluated across the gap would put a patient's own cells in each
other's neighbourhoods. `sample_id` is the patient and is what `fitSpiDE()` puts
the random intercept on, because the condition is a patient-level contrast.
GSE250346 is the case that makes this concrete: 45 sections nested in 35 donors.
Per-cohort consequence: GSE282639 uses **local** FOV pixels (each FOV is its own
piece of tissue) while GSE289194 uses **global** pixels (a core's FOVs tile one
contiguous piece).

**`assay("counts")` is raw integer counts.** `polishSpiDE()` refuses non-integer
input. Negative-control features are removed from the gene set — they would enter
`edgeR`'s cross-gene dispersion moderation — and their per-cell total is kept as
`neg_frac`, available as a covariate.

QC is deliberately loose (`nCount ≥ 20`, `nFeature ≥ 5`, `neg_frac ≤ 0.1`): the
purpose is to remove cells the model cannot use, not to curate. On GSE282639 it
keeps 90.5% of 115,623 cells.

`.checkSchema()` refuses any object whose condition is not a 2-level factor, whose
condition varies within a patient, whose sections are not nested in patients,
whose counts are non-integer, or that still carries negative-control features.
Each check corresponds to a failure that is otherwise found hours into a fit.

---

## 4. Annotation and the analysis partition

Code: `research/public/R/05_annotate.R`. Two separate jobs, and conflating them
is the usual mistake.

**Labelling.** Labels are assigned to **clusters**, not cells: a per-cell call
fragments a type across a neighbourhood and the KDE then reads a density that
flickers between types. Clustering runs on log-normalised counts with the sample
effect left in — a batch correction would move cells across type boundaries in a
way that can correlate with the condition, and a cell-type assignment that
depends on the condition is a direct route to a false niche-dependent call.
Cluster purity is recorded and a cluster below 50% purity is reported.

Reference per cohort, each closer than a public atlas would be:

**Are author annotations available?** Checked per cohort, since a deposited
label always beats a transferred one:

| cohort | author annotation | what we use | bridge |
|---|---|---|---|
| GSE250346 | **yes, complete** — `final_CT`, 47 types, all 1.63 M cells, in the deposited Seurat | the authors' labels, unchanged | none |
| GSE289194 | **partial** — 19 states on 79,364 cells (5.2% of the retained slides) via the authors' GitHub | those labels as a within-dataset SingleR reference | none — same cells, same panel |
| GSE282639 | **no** — GEO has only `RAW.tar` and a logCPM matrix; the paper says *"Code is available upon request"* | SingleR from Habermann et al. 2020 lung scRNA (GSE135893) | scRNA ↔ CosMx, over ~500+ shared genes |

**The GSE282639 reference was changed, and the first choice was a real defect.**
The original plan transferred from GSE250346 — same tissue, same disease, and
convenient. It is also **circular**: §8.2 scores cross-platform replication
*between* these two cohorts, and deriving one's cell types from the other builds
part of that concordance in before any model is fitted. It would have inflated
the headline replication result.

The fix is also the better transfer on its own terms. GSE250346 is a 343-gene
panel sharing only **177** genes with this 974-gene CosMx panel, and over 177
genes the CD4+ T / CD8+ T / Treg / NK profiles are nearly indistinguishable: the
transfer collapsed T into NK (of 4,384 cells with ≥ 2 of six canonical T markers,
1,514 landed in "NK" and 526 in "T cell", while only 866 cells in the whole
cohort are NK-marker positive), forcing a T/NK merge. Habermann is
whole-transcriptome, so the reference is built on the full panel, and it is what
the depositors themselves used ("a previously published 10X single cell dataset
from IPF and controls"). Its 31 labels are mapped onto the **same 13-compartment
vocabulary** GSE250346 uses, so the two lung cohorts stay comparable while their
references stay independent. The T/NK merge is now conditional on the data
rather than assumed: it applies only if either compartment comes back under 500
cells, and the outcome is recorded.

The GSE289194 join deserves a note because it was not obvious: the authors'
released metadata carries no cell id, so labels are keyed onto our cells by
`(slide, FOV, integer local centroid)`. That key is unique on both sides and
matches **100%** of the 6,375 labelled F1 cells.

**Partitioning.** The label set is *not* the analysis partition. Each index type
carries `2 × n_niche` parameters and a gene can only support them where it is
expressed; CLAUDE.md records that at 22 × 22 on YTMA, 58 of 769 genes died on
singular per-gene Grams. Every cohort is therefore collapsed to an **index side
≤ 15 and a niche side ≤ 10**, targeting ~150 tested three-way columns — the same
testing-family size as the YTMA v11 arm (156), so the FDR cascade operates in a
regime this project has already calibrated. The maps are written out as data in
`MAPS`, reviewable and never inferred. Compartments under 500 cells are **merged
into `Other`, not dropped**: removing cells changes the niche density every
remaining cell sees.

**The annotation is guarded by markers, because a lost lineage is invisible in
every summary of what was assigned.** Clustering resolution had to be raised
twice: `k = 15`, resolution 1 gave 12 clusters for a 47-type reference, and
resolution 2 gave 22 — at which point GSE282639 came back with **no T-cell and no
B-cell compartment at all**, T having merged into NK/NKT and B into Plasma. The
markers were never missing (CD3D/E/G, CD2, CD8A, CD4, IL7R, MS4A1, CD79A, CD19
are all in the 177 shared genes); ~4,800 cells per cluster was the problem, and
a table of assigned compartment sizes looks entirely reasonable when a lineage
has been absorbed rather than dropped.

`.markerCheck()` therefore counts cells carrying ≥ 2 canonical markers of each
lineage and compares that to the size of the compartment they were assigned to.
A lineage with ≥ 1,000 marker-positive cells and under 25% of them in its own
compartment **stops the build**. Default resolution is now 4.

**Annotation robustness is an arm, not an assumption.** For each cohort the
calibration null is re-run under a *second, independent* partition (a coarser
map; for GSE282639 also an HLCA-based SingleR transfer instead of the GSE250346
one). Calibration must be invariant to this. If it is not, the calibration claim
is about a partition, not about the estimator, and the paper must say so.

---

## 5. Model configuration — fixed across every arm

The point of fixing these is that the external validation tests the **shipped
defaults**, not a per-cohort tuning. Any deviation is an arm with its own token.

```r
fitSpiDE(spe, condition = <cond>, sigma = <grid>,
         sample_id  = "sample_id",      # patient — random-effect group
         covariates = "loglib",
         random = "intercept", re.celltype = TRUE,
         df.method = "satterthwaite",
         re.prop = 1, re.maxit = 2L, re.maxit.psi = 1L,
         winsor = 4, lambda.a = 0)
polishSpiDE(fit, spe, psi = "profile", tau2 = TRUE)
testSpiDE(fit, spe, fdr = 1, dispersion = "ql", combine = "cauchy")
```

Niches are built **per section**: `buildNiches(spe, sigma = grid, sample_id = "section_id")`.

`fdr = 1` throughout. CLAUDE.md records that gating the stored q-values at a
level reproduces `testSpiDE(fdr = level)` exactly, so every nominal level is
scored from one stored table without refitting.

`loglib` is mandatory, not optional: "There is no library-size term, and that is
a real gap" (CLAUDE.md) measured the best-calibrated configuration anywhere in
the investigation to be raw counts **plus** a depth term, and the slope the data
want is 0.979 — a conventional offset. It enters through `covariates` as a
centred `log(nCount)`.

Deviating arms, each one token:

| token | change | why |
|---|---|---|
| `_nonest` | `re.celltype = FALSE` | **C2** — the whole point |
| `_psimod` | `psi = "moderated"` | cost/benefit; CLAUDE.md 2026-09-14 says not a free default |
| `_pearson` | `dispersion = "pearson"` | the 0.99.17 scale, for continuity |
| `_brown` | `combine = "brown"` | combiner sensitivity |
| `_cov` | `+ area_um2, aspect, neg_frac` | imaging artefacts track local density by construction |
| `_coarse` | second annotation | §4 robustness |

---

## 6. The bandwidth grid

σ is a physical length in microns and must be set from the data, not copied from
YTMA. Two constraints bound it.

*Lower:* a bandwidth below the mean nearest-neighbour cell spacing estimates the
cell itself, not its neighbourhood. *Upper:* a bandwidth approaching the short
side of the smallest section makes the KDE nearly constant within a section, and
the niche covariate then aliases the section intercept. **The binding case is
GSE282639, where a section is a single 509 µm FOV.**

Rule: σ ∈ {15, 30, 60, 120} µm, subject to `σ_max ≤ (1/4) × min section short
side`, computed per cohort by `06_qc_report.R` and recorded. A cohort whose
sections cannot support 120 µm runs the truncated grid and the report says which.
The grid is shared across cohorts so that a bandwidth means the same
neighbourhood everywhere — which is the entire reason coordinates were
harmonised to microns.

---

## 7. Calibration — the primary result

Three nulls, because they kill different axes and only the pair brackets the
three-way term. Implementation is ported unchanged from
`research/fdr-ordering/R/package_fixed_design.R` (`SPIDE_SHUF`).

| null | what is permuted | what becomes null | what it tests |
|---|---|---|---|
| `free` | niche matrix rows within (section × cell type) | `CellType:niche`, `CellType:condition:niche` | the niche axis |
| `block` | toroidal shift of the niche field within a section | the same | the niche axis **with spatial smoothness preserved** |
| `perm` | the condition label, one draw per patient | `CellType:condition`, `CellType:condition:niche` | the condition axis |

`block` is the one that matters most. CLAUDE.md records that on YTMA the free
shuffle calibrates flat at 0.96–0.98 while block keeps 0.988 → 1.206 in the
brightest genes, i.e. the nested intercept does **not** remove the
spatially-smooth-covariate component. Whether that residual is a property of the
estimator or of that cohort is exactly what four cohorts answer.

**Do not shuffle expression across cells.** It breaks each cell's pairing with
its own library size and tests the offset construction at the same time.

**Replication:** 5 seeds per null per cohort per bandwidth. Five is what the YTMA
block grids used and what the exceedance ranges quoted in CLAUDE.md are built on.

**Metrics**, computed per gene and per index cell type:

1. `sd(t)` overall and by expression band (the per-gene granularity CLAUDE.md
   establishes as the right one — per-index pools over genes and under-calibrates
   the bright ones).
2. `P(|t| > 1.96)`, target 0.05.
3. Exceedances of `|t| > 4.89` against the null range.
4. **Empirical FDP of the FDR cascade** at nominal 0.01 / 0.05 / 0.10 / 0.20:
   discoveries on the complete null are all false by construction.
5. λ (the calibration slope) via `.claude/skills/calibration-check`.

**Pass criteria, pre-registered.** Per cohort, on `free` at the primary
bandwidth: median per-gene `sd(t)` in **[0.95, 1.05]**, no expression band
outside **[0.90, 1.10]**, and cascade FDP at nominal 0.05 **≤ 0.10**. On `block`:
the same, with a documented exemption for the top expression band, whose YTMA
value is 1.10–1.21. A cohort failing `free` is a finding and is reported as one —
it is not re-tuned until it passes.

---

## 8. Power and recall

Real data has no ground truth, so power is measured three ways and only the
first is quantitative.

**8.1 Injection (the quantitative instrument).** Reuse
`research/fdr-ordering/R/inject.R`: plant a known three-way `celltype ×
condition × niche` effect into the real counts of a chosen (index, niche) pair,
at a grid of effect sizes, then run the full pipeline and score TPR against the
planted set and FDP against everything else. Because the host is real data, the
covariance structure, the depth distribution and the per-gene tail are all real —
which is precisely what the synthetic benchmark cannot supply.

Grid: effect size β ∈ {0.25, 0.5, 1, 2} on the log scale, 200 planted genes per
replicate, 5 replicates per cohort per β. Report TPR at realised FDP, and the
ROC over the stored `fdr = 1` table.

Injection targets are chosen to span the cell-count range, because CLAUDE.md
records that the residual inflation tracks the number of cells in the index type:
one populous index type, one mid, one near the 500-cell floor.

**8.2 Cross-platform replication (PF).** GSE250346 and GSE282639 are the same
disease contrast on two platforms, two labs, two panels. Score the concordance of
called (gene, index, niche) triplets over the **shared gene and compartment
space**, against a null concordance from label-permuted calls. This is the
strongest evidence available that calls are biology and not platform artefact,
and no previous spiDE analysis could run it.

**The shared space is 177 genes** (measured 2026-09-21: CosMx 974 ∩ Xenium 343).
Each cohort is fit on its own full panel and concordance is scored on the
intersection only, so the estimate is unbiased but its precision is bounded by
177 genes × ~90 index/niche slots. That is enough to detect substantial
concordance and not enough to bound it tightly; the report quotes the interval,
not a point estimate, and a null result here is uninformative rather than
negative.

**8.3 External biological comparator.** GSE250346 ships the authors' CNiche and
TNiche assignments and pathologist geojson regions; GSE289194 ships seven
published niches and a nodal-vs-extranodal T-cell phenotype result. Ask whether
spiDE's calls in an index type land where the published niche structure says the
relevant neighbours are. Qualitative, reported as such, never as a validation
statistic.

---

## 9. Scoring and multiplicity

**Quote per-gene calibrated statistics on real data, never the cascade's raw
counts.** CLAUDE.md is explicit: on the real transcriptome the cascade returns
1,079 triplets against a null mean of 311, an empirical FDP of 0.29. Every
discovery count in this analysis is a per-gene calibrated count with its
empirical FDP from the matched null attached, and the cascade's count is reported
beside it as a *diagnostic of the cascade*, which is a separate and legitimate
result (C1 extends to "and here is how the multiplicity step behaves off YTMA").

The gene filter (`sd(shuffle t) > 1.3`) is applied **only** with the hot list
built from grids other than the one scored. CLAUDE.md records that defining it on
the scored grid reproduces a withdrawn 6× figure; the honest effect is 1.5×.

`spiGSEA()` is out of scope for this analysis except to record the stored `rho`
per cohort, because the set-level count scales with it and no set-level null has
ever been run. A gene-set claim is not made.

---

## 10. Arms and compute

**Step 0 — sizing pilot, before anything is submitted.** One bandwidth, 200
genes, one null seed, per cohort, on CPU, timed. The four cohorts differ by 200×
in cells × genes (GSE282639 at 974 × 105 k against GSE282193 at 6,175 × ~700 k)
and nothing about YTMA's timings transfers. The pilot fixes cores, memory,
walltime and whether a cohort needs the GPU fit path. Its output is a table in
the report; no full arm is submitted without it.

**Staging.** Every cohort runs as the chained `fit → polish → test` arrays the
cohort driver already implements (`SPIDE_STAGE`), with checkpoints keyed to the
arm, because the polish is CPU work of hours and must not hold a GPU idle.

**Two operational rules that are not negotiable**, both from losses recorded in
CLAUDE.md:

1. **Pin `SPIDE_PKG` to a frozen snapshot** (`freeze_snapshot.sh`). A sweep whose
   tasks load the live working tree is not one experiment.
2. **Never edit a driver while an array is live.** R parses a script through an
   8 KB buffer; an in-place edit killed 255 benchmark tasks and eight 6–7 h cohort
   tasks on 2026-09-07. Each sbatch copies its driver to a task-private temp file
   at task start.

Additionally: **one BLAS thread per forked polish worker**
(`.bplapplySingleBLAS`), and check `OMP_NUM_THREADS × SPIDE_CPUS` in every
sbatch — set-threads-then-fork cost a 9× penalty on the 0.99.19 runs.

**Queue discipline.** A chained YTMA cohort run (`pkgfinal13`, 11 tasks × 64
cores) is live as of writing. These arms queue behind it; `mem` requests are
sized from `sstat` on the pilot, not guessed, because the 16 TB per-user memory
cap turns over-requesting into halved concurrency.

**Arm matrix** (per cohort; `×5` = seeds):

| arm | GSE250346 | GSE282639 | GSE289194 |
|---|---|---|---|
| real, condition mode | ✓ | ✓ (disease) + ✓ (MUC5B) | ✓ |
| real, niche mode | ✓ | ✓ | ✓ |
| `free` null ×5 | ✓ | ✓ | ✓ |
| `block` null ×5 | ✓ | ✓ | ✓ |
| `perm` null ×5 | ✓ | ✓ | ✓ |
| `_nonest` (C2) `block` ×5 | ✓ | ✓ | ✓ |
| injection, 4 β × 5 reps | ✓ | ✓ | ✓ |
| `_coarse`, `block` ×5 | ✓ | ✓ | ✓ |

GSE282193 is parked (§2.1); its rows are removed rather than left as "gated",
so the matrix is the thing that is actually being run.

Sensitivity arms (`_psimod`, `_pearson`, `_brown`, `_cov`) run on **GSE282639
only** — smallest cohort, fastest turnaround — unless it shows a difference, in
which case the differing arm is repeated on GSE289194.

---

## 11. Deliverables

1. `research/reports/benchmarks/spiDE-public-cohorts.Rmd` — the report, in the
   reading order the other six use, rendered into the research site.
2. `research/reports/benchmarks/tables/public_cohorts.rds` — **one canonical
   table**, method arms as extra *rows*. Never a parallel file.
3. A cohort-summary table (§2) with the measured numbers, generated not typed.
4. Figures: per-gene `sd(t)` by expression band per cohort per null; cascade FDP
   vs nominal; injection TPR/FDP curves; the PF cross-platform concordance.
5. `research/fdr-ordering/FINDINGS.md` entry per measured result, via the
   `record-finding` skill, with date, job ids, snapshot and seed.
6. CLAUDE.md gains a short section pointing at the report; numbers stay in the
   report and `vignettes/spiDE-calibration.Rmd`, per the existing rule that the
   calibration vignette is the only place in the package that quotes benchmark
   numbers.

---

## 12. What this analysis does not establish

State these in the paper rather than leaving them to be found.

- **Not a transcriptome-scale result.** The deepest public panel is 6,175 genes
  against YTMA's 13,348, and no public CosMx WTA multi-patient cohort exists. The
  per-gene tail claim remains YTMA's.
- **Not a test of `spiGSEA()`.** No set-level null has been run anywhere, and the
  stored `rho` moves the set count by 3× (618 vs 183 sets on identical fits).
- **Not a biological discovery.** Any calls reported are calibration and
  replication evidence. A biological claim needs independent validation,
  which is out of scope here.
- **Not a test of `random = "slope"`.** Every arm is `random = "intercept"`; the
  slope model needs `re.maxit = 10` and is a separate question.
- **Nothing about panel depth.** GSE282193 is parked (§2.1), so the deepest
  panel here is 974 genes against YTMA's 13,348. C4 is withdrawn.
- **Annotation is inferred for three of four cohorts.** §4's robustness arm
  bounds the exposure; it does not remove it.

---

## 13. Risks

| risk | mitigation |
|---|---|
| GSE250346's 1.6 M cells make the design matrix intractable | sizing pilot first; section-stratified subsample as a fallback arm, with the subsample fraction reported |
| cross-platform PF gene overlap too small for §8.2 | `05_annotate.R` stops below 150 shared genes; concordance is scored on the shared space only and the size is reported |
| a cohort fails calibration | reported as a finding, not tuned away; §7's criteria are pre-registered here for exactly this reason |
| annotation drives the result | §4 robustness arm; calibration must be invariant |
| queue contention with the live YTMA run | pilot-derived sizing; arms queue behind |

---

## 14. Provenance

Every object records its source accession, the build date, the spiDE version,
the pixel scale, the compartment map and the annotation route in
`metadata(spe)`. Every run records the snapshot path, the seed, the shuffle mode
and the arm token in its output, as the YTMA driver already does. A number that
reaches the paper must be traceable to a job id and a frozen package.
