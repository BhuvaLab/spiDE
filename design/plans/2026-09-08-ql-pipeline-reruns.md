# What the 0.99.18 pipeline invalidates, and the rerun plan

Date: 2026-09-08. Status: **plan, awaiting go-ahead**; nothing below has been
submitted.

The package defaults moved twice in two releases: 0.99.17 added the nested
(sample × cell type) intercept and the per-gene convergence step (profile-ML
dispersion, Pearson scale), and 0.99.18 keeps the moderated dispersion at the
converged mean and scales the standard errors by the quasi-likelihood
dispersion. Every number quoted as "the shipped pipeline" was produced under
one of three pipelines, and only the most recent switch arms were produced
under the one that ships now. This audit says which is which, and what a
rerun costs.

## 1. Audit: what was produced under which pipeline

Canonical tables (`research/reports/benchmarks/tables/`), by the `design`
column of their rows, and the cohort grids (`research/fdr-ordering/out_convnull/`),
by arm token. Pipelines: **A** = design term only (0.99.16 flags: no nested
block, no convergence), **B** = 0.99.17 defaults (nested + converged, profile
ψ, Pearson), **C** = 0.99.18 defaults (nested + converged, moderated ψ, QL).

| evidence | pipeline | rows / tasks | verdict |
|---|---|---|---|
| null type-I, all random modes (`fixed`/`intercept`/`slope`) and both df methods, both layouts, full S sweep — `null_type1`, `null_qq` (`celltype-response`) | A | 400 jobs | **rerun**: the random-mode and df comparisons the calibration vignette opens with are pipeline-A numbers |
| power: samples / cells / effect / layout sweeps, `power_fdr`, `power_points`, `power_by_category` (`celltype-response`) | A | 720 jobs | **rerun**: only the samples and effect sweeps exist under B/C |
| combiner (Cauchy vs Brown), `combiner` | A | 80 jobs | **rerun** with the study (the t distribution feeding the combiners changed) |
| `ctresp` (the `CellType:condition` design comparison) | A vs legacy design | 120 + 120 jobs | not needed: it compares designs, not pipelines; the result is a property of the term |
| gene-set calibration (`gsea_*`, self-contained vs competitive) | A | 228 jobs | **rerun**: spiGSEA averages per-gene t, whose scale changed |
| `pql_*` (the `re.prop` subsampling study) | A | 200 jobs | not needed: a rejected option, kept as history |
| switch arms `nested-converged`, `nested-only`, `converged-only` | B (and A/B mixes) | 1,160 jobs | keep: they are the record of *why* the defaults moved |
| switch arms `psi-moderated`, `psi-moderated-ql`, `ql-only` | C for `psi-moderated-ql` | 640 jobs, null at S = 4/10/16/30 + power samples/effect | keep: already pipeline C, intercept mode only |
| two-stage rows everywhere, `adversarial_*`, `lsstruct_nullrate` | (no GLM inference) | — | unaffected |
| `shuffle_null`, `shuffle_decile`, `spatial`, `timing` | A | — | shuffle: keep (the pre-fix record of the defect); `spatial`/`timing`: rerun with the study |
| cohort base arm bw 30 (real + 5 block + 3 free) | B | 9 tasks | keep as the B record; pipeline C already measured on the same panel (`_psimod_ql`: same calls) |
| cohort artefact arms `_cov5`, `_plasma` (bw 30), base bw 10/50/70 | B | 18 tasks | optional: the identity-instability conclusion is about the perturbations, and C was shown not to change calls |
| cohort **final configuration** `_cov5_plasma` bw 30 (real + 5 block) and bw 10/50/70 (real + 2 block each) | B | 15 tasks | **rerun**: this is the reported configuration |
| cohort **transcriptome** `_cov5_plasma_all` (in flight, job 28132880, 14 h) | B | 6 tasks at 512G | **rerun** after it lands; the in-flight run becomes the pipeline-B comparator at scale |
| `compositionTest()` results | (limma pseudobulk) | — | unaffected |
| long tests (`longtests/testthat/`) | run under current defaults | 4 files | **rerun** as part of the release gate |
| toy numbers quoted in CLAUDE.md (planted t 1.68 → 10.19; niche-mode 7.82 / 5.63 → 14.27) | B | minutes | **re-measure** and requote |

## 2. What the reruns need from the harness

1. **A shipped-pipeline arm.** `design_switch_arms` gains
   `"nested-converged-ql" = list(re.celltype = TRUE, converge = TRUE,
   fit = list(polish.psi = "moderated"), test = list(dispersion = "ql"),
   families = <all but pql>, modes = "all")`, on a frozen 0.99.18 library
   (`research/libs/0.99.18`, `SPIDE_LIB_0918`). `.arm_modes()` currently
   restricts switch arms to intercept mode; the new `modes` field lifts that
   so the random-mode and df comparisons are produced. Check that the
   combiner, ctresp and gsea families route through `.arm_fit_args()` /
   `.arm_test_args()` (the gsea family has its own `.gsea_fit()`).
2. **A `design_shipped` config value** the reports read instead of the
   literal `"celltype-response"`, so "the shipped arm" in every figure moves
   with the release; the palette gains the level.
3. **Cohort driver defaults.** `package_fixed_design.R` still defaults
   `SPIDE_PSI=profile` and `SPIDE_DISP=pearson` and suffixes the *new*
   pipeline (`_psimod`, `_ql`). Invert: the package defaults need no
   suffix; `_profile` / `_pearson` mark the old ones. Existing pipeline-B
   grids of the rerun arms are archived with a `v17_` prefix (as the
   `prereview_` grids were), which the scorer and the site summaries ignore
   by pattern, so the reruns take the unsuffixed names.
4. **GPU gate.** The H100 validation (job 28159434) ran the QL pre-pass and
   the inference on the toy; the earlier inference failure was on the
   cohort's design. One real grid at bw 30 with `backend = "gpu"` against
   the existing CPU `_psimod_ql` grid (same configuration) decides whether
   the cohort and transcriptome reruns use the accelerator (fit 4× faster).

## 3. The plan, in order

| phase | what | size | wall |
|---|---|---|---|
| 0 gate | finish the release (suites, R CMD check, merges, branch deletions); freeze a 0.99.18 snapshot and build its library; long tests; re-measure the toy numbers; GPU gate on one cohort grid | 1 GPU task + local | half a day |
| 1 cohort | archive B grids of the final configuration; rerun bw 30 (real + 5 block) and bw 10/50/70 (real + 2 block each) under C; rescore; site summaries; cohort report | 15 tasks, 64G, ~3 h each (moderated ψ) | a day, QOS-bound |
| 2 benchmark | the `nested-converged-ql` arm, all families but pql, all modes, 40 reps; install rows; reports and the calibration vignette switch to the new shipped arm | ~1,350 jobs, 20–70 CPU-min each | a day at ~200 concurrent |
| 3 transcriptome | after 28132880 lands (kept as the B comparator), rerun 6 tasks under C, GPU fit if the gate passed | 6 tasks, 512G or GPU | one to two days |
| 4 optional | the artefact arms under C, only if the cohort report must be single-pipeline | 18 tasks | a day |

Phases 1 and 2 are independent of each other and of phase 3; 3 waits on the
memory the in-flight run holds. Nothing in phase 2 changes a published
conclusion by expectation — pipeline C was measured against B on the
families that carry the decisions — but the vignette's opening numbers
(random modes, df methods, the small-sample regime) are pipeline A and should
be replaced, not annotated.
