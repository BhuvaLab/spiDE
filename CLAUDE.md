# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this package does

spiDE identifies **context-specific, neighbourhood-dependent differential expression** in spatial
transcriptomics data. Within an *index* cell type, it tests how gene expression changes with an
experimental *condition* as a function of the local density (the *niche*) of surrounding cell types.
It is a from-scratch R/Bioconductor reimplementation of the analysis previously done in flat scripts
(`batch_nichede_v9.R` / `YTMA_nicheDE_v9.md`), built on top of the author's own **SpaNorm** package,
which supplies the generic per-gene negative binomial GLM fitting engine (`SpaNorm::fitNB`,
`SpaNorm::calculateMu`, `SpaNorm::invert_mat`).

## Commands

This is a standard R/Bioconductor package (roxygen2-documented, testthat edition 3). Run everything
from an R session in the package root, e.g. `Rscript -e '<command>'`.

```r
devtools::load_all()                 # load package for interactive dev
devtools::document()                 # regenerate NAMESPACE + man/*.Rd from roxygen comments
devtools::test()                     # run the full testthat suite
devtools::test_active_file()         # run only the currently open test file
testthat::test_file("tests/testthat/test-fitSpiDE.R")  # run one test file
devtools::check()                    # full R CMD check
BiocCheck::BiocCheck()               # Bioconductor-specific compliance checks
devtools::build_vignettes()          # render vignettes/spiDE.Rmd
```

Regenerate the shipped example dataset (`data/toySpiDE.rda`) with `source("data-raw/make_toySpiDE.R")`
(loads the package via `devtools::load_all()` first, since it calls the internal `.toySPE()`).

`longtests/testthat/` holds the slow numerical checks (mixed-effects numerics, GSEA numerics,
niche-only calibration). `devtools::test()` does **not** run them and neither does CI — run one
explicitly, e.g. `testthat::test_file("longtests/testthat/test-mixed-numerics.R")`.

Project automations live in `.claude/` (allowlisted in `.gitignore`, so they are shared):
two hooks (`r-parse-check.sh` parses every edited `.R` file; `protect-canonical-tables.sh` blocks
direct writes to `research/reports/benchmarks/tables/*.rds`), two agents
(`numerical-robustness-reviewer`, `evidence-auditor`) and two skills (`run-benchmark-arm`,
`calibration-check`). `.mcp.json` adds a GitHub server that reads `${GITHUB_PAT}` from the
environment — no token is committed.

CI is `.github/workflows/check-bioc.yml` (R CMD check + `BiocCheck::BiocCheck()` across four
R/Bioconductor configurations on push/PR to `main`) plus `pkgdown.yaml`. There is still no lint
config (`.lintr`) — lintr diagnostics surfaced by the editor reflect default rules, not a
project-specific config. Dot-separated argument names like
`lambda.a`, `winsor`, `maxit.psi` are intentional: they mirror `SpaNorm::fitNB`'s own argument names
so they can be forwarded via `...` without renaming.

## Architecture

**Entry points.** `spiDE()` chains the three stages below; `buildNiches()` → `fitSpiDE()` →
`testSpiDE()` is the same pipeline unrolled. `twoStageSpiDE()` is a *separate estimator* over the
same niche `reducedDims` (see "The two-stage estimator"), and `spiGSEA()` runs on an already-fitted
object.

### Pipeline (three stages, chained by `spiDE()`)

1. **Niche construction** — `buildNiches()` (`R/buildNiches.R`). Per sample and per bandwidth `sigma`,
   a Gaussian KDE of every cell type is evaluated at each cell's location (`.effectiveNiche`, via
   `spatstat.explore::densityfun`), producing a cells × cellTypes density matrix stored as
   `reducedDim(spe, "Niche<sigma>")`. Samples are processed with `.speApply()`
   (`BiocParallel::bplapply`, replacing the original scripts' `foreach`/`doParallel`). Cell types
   absent from a given sample are zero-filled via `.fillMissingDims()`. `mergeNiches()` can coarsen
   niche columns into groups afterwards; `computeSizeFactors()` derives a per-sample library-size
   `colData` offset from counts/area.

2. **Design + fit** — `fitSpiDE()` (`R/fitSpiDE.R`, design assembly in `R/design.R`). For each
   bandwidth, `.buildNicheDesign()` builds
   `~ 0 + <covariates> + CellType + CellType:condition + CellType:(niche cols) +
   CellType:condition:(niche cols) + niche cols`,
   drops symmetric self-interactions (an index cell type against its own niche density), and tags
   every column via `.tagCovtype()` as one of `CellType` / `Niche` / `Response` /
   `ResponseCellType` / `ResponseNiche` / `Other`. Note the **cell-means coding**: there is no bare
   `condition` main effect — the condition contrast is carried by the `CellType:condition`
   (`ResponseCellType`) columns, one per cell type. Code that looks up a single `"Response"` column
   will find nothing; match on `ResponseCellType` instead.

   **Two modes.** With `condition = NULL` the design drops the condition terms *and* the
   bare niche main effects (`~ 0 + <covariates> + CellType + CellType:(niche cols)`), and the
   two-way `CellType:niche` columns — tagged `Niche` in *both* designs — become the tested
   effects, cell-means coded. The main effects must go: `niche_n = sum_c CellType_c:niche_n`,
   so keeping them makes `model.matrix()` alias away one interaction per niche (always the
   alphabetically first cell type's), which is harmless in condition mode but leaves that cell
   type untestable in niche mode. Which tag is the tested tag is decided in exactly one place —
   `.testedCols()` / `.nicheTestCols()` in `R/design.R`, keyed off the `mode` slot
   (`"condition"` / `"niche"`) carried on `SpiDEFit`/`SpiDEResults` and read via `.fitMode()`.
   **Never match a covtype literal downstream** — go through those predicates, or a new mode
   silently skips your code path. In niche mode `results(type = "celltype")` /
   `results(type = "patient")` are empty, `spiGSEA(type = "celltype")` errors, and
   `df.method = "between"` means the cell-level residual df under `random = "intercept"` (the
   tested slope is a within-sample contrast) or `S - 1` under `random = "slope"`. Niche mode is
   mildly anti-conservative because the niche covariate is spatially autocorrelated and spiDE
   does not model that; see `longtests/testthat/test-nicheOnly-mixed.R`.

   **`ResponseNiche` columns are the scientifically important ones** — the three-way
   `celltype:condition:niche` interactions. The whole gene set is fit in a single
   `SpaNorm::fitNB(Y, W, ...)` call per bandwidth — **never slice genes before fitting**, because
   `fitNB` moderates dispersion across the full gene set via `edgeR::estimateDisp(robust=TRUE,
   tagwise=TRUE)`; blocking genes at this stage would change the fit. Produces one `SpiDEFit` per
   bandwidth (S4 class in `R/AllClasses.R`), collected into a `SpiDEResults`.

3. **Inference + combination + FDR** — `testSpiDE()` (`R/testSpiDE.R`), which chains three internal
   stages:
   - `.blockedInference()` (`R/inference.R`) — per gene *block* (genes are independent post-fit, so
     this stage IS blockable, dispatched via `BiocParallel`), computes working weights from
     `SpaNorm::calculateMu`, Wald t-statistics/SEs for the `Response`/`ResponseNiche` columns
     (`SpaNorm::invert_mat` for the covariance inverse), then combines correlated `ResponseNiche`
     p-values within a gene — separately for up/down directions and at both gene-level and
     per-index-cell-type granularity — using the combiner chosen by the `combine` argument of
     `testSpiDE()`/`spiDE()`: **`"cauchy"`** (the default; the correlation-agnostic tan-transform
     Cauchy/ACAT test via `.cauchyCombine()`, with a `1e-15` clamp on the one-sided p-values) or
     **`"brown"`** (Brown's method, `poolr::mvnconv` + `poolr::fisher`, which consumes the
     coefficient correlation matrix). Brown's method mirrors the per-gene loop in the original
     `batch_nichede_v9.R` almost line-for-line; Cauchy combines **two-sided** p-values while Brown keeps
     **one-sided** ones — `tan((0.5 - p)pi)` diverges to `-Inf` as `p -> 1`, so under one-sided input
     a gene up in one niche and down in another cancels exactly, whereas Brown's `-2log(p)` is
     bounded at 0 and is safe one-sided. Under Cauchy both `p.combined.pos` and `p.combined.neg`
     therefore carry the same combined value and `SpiDEFit@two.sided` tells the FDR cascade not to
     apply the direction split. Cauchy was made the default after a calibration/power
     study (`research/reports/benchmarks/spiDE-cauchy-vs-brown.Rmd`) showed it matches or beats Brown
while controlling
     type-I error under correlation without estimating `R`. The two combiners populate the same
     `p.combined.pos`/`p.combined.neg` slots (see below), so downstream code is combiner-agnostic.
   - `.cauchyCombine()` / `.geneWeights()` / `.combineBandwidths()` (`R/combine.R`) — combines
     p-values **across bandwidths** with a tan-transform Cauchy combination test, weighted by each
     gene's relative log-likelihood across bandwidths (`exp(loglik - rowMax)`, thresholded).
   - `.hierarchicalFDR()` (`R/fdr.R`) — a three-level nested Benjamini-Hochberg cascade: gene level →
     per-index-cell-type level (both gated at `fdr/2` per direction, then merged into `"Both"` when
     both directions pass) → per-(gene, index, niche-cell-type) level (gated at `fdr`). Produces the
     tidy `results()` table keyed by `(gene, ct_index, ct_niche, bandwidth.max)`.

### Key invariant: fit whole, infer blocked

The counts matrix `Y` may be dense, sparse, or a `DelayedArray` and is never eagerly densified. The
**fit** stage (`fitSpiDE`) always sees the entire gene set in one `fitNB` call. The **inference**
stage (`.blockedInference`) is the only place genes are chunked (`.chunkGenes()`,
user-controlled via `block.size`) and parallelised (`BPPARAM`) — because per-gene Wald +
combination results depend only on that gene's own `alpha`/`psi`, this is exact, not an
approximation. When modifying either stage, preserve this split.

Both stages take `backend = c("auto", "cpu", "gpu")`, forwarded to `fitNB` for the fit and used by
`.blockedInference()` for the batched per-gene Wald covariance. The batching helpers live in
`R/inference-batch.R` and must behave identically on a base R matrix and a torch tensor (`.rowsOf()`,
`.gramBatch()`); two independent memory budgets bound them (`.inferenceBlockSize()` for the gene
block, `.covBatchSize()` for the covariance sub-batch — the latter applies on **both** backends), and
`gpu.mem.budget` overrides the GPU one. GPU is opt-in via `SpaNorm::checkGPU()`; `torch` is only in
`Suggests`, so nothing here may hard-depend on it.

### S4 classes

- `SpiDEFit` (`R/AllClasses.R`) — one bandwidth's fit + inference (design `W`, per-column `covtype`
  tags, `coefmap`, per-gene `alpha`/`psi`/`loglik`, and once inferred: `t_stat`, `se`,
  `p.combined.pos`, `p.combined.neg` — the last two hold the within-gene combined p-values from
  whichever combiner (`"cauchy"` default / `"brown"`) was used; they are a
  `genes × (1 + n_index)` matrix, column `"Gene"` then one per index cell type).
- `SpiDEResults` (`R/AllClasses.R`) — container for a list of `SpiDEFit` (one per bandwidth) plus
  cross-bandwidth combined p-values and the final tidy `results` data.frame (with
  `results.celltype`/`results.patient` behind `results(type = )`, and `diagnostics` for the two-stage
  path). Both classes expose a `$` accessor (`slot(x, name)`) and a `show` method; validity is
  enforced via `setValidity()` — note `validSpiDEResults()` requires `length(fits) == length(sigma)`
  **only when `fits` is non-empty**, so `twoStageSpiDE()` can return a bandwidth with no GLM fit.
- Generics live in `R/AllGenerics.R`; methods are implemented per-file (`buildNiches` in
  `buildNiches.R`, `fitSpiDE` in `fitSpiDE.R`, etc.) — when adding a new exported function, add the
  generic there, not inline in the implementation file.

### Toy fixture

`.toySPE()` (`R/toydata.R`, internal, `@noRd`) generates a seeded synthetic `SpatialExperiment` with a
deliberately planted, recoverable effect: gene `G1` is up-regulated in index cell type `A`, in
`Responder`s, in proportion to the local density of niche cell type `B` (B cells cluster at high `x`).
Exported examples and the vignette use the pre-baked `data(toySpiDE)` instead (built by
`data-raw/make_toySpiDE.R`), since exported-function `@examples` cannot call internal helpers under
`R CMD check`. `field`/`n_per` were tuned (500 units, 80 cells/sample) specifically so all four default
bandwidths (10/30/50/70) fit without IRLS collinearity failures — don't shrink the field without
re-checking every bandwidth still converges.

### The sample-level correction (and the default)

`random = "none"` uses reduced-design Wald standard errors formed from cell-level information. That
treats every cell as an independent replicate of a patient-level contrast, and is badly
anti-conservative on data with few samples but many cells: on a null with per-sample intercepts it
rejects at **~0.71** against a nominal 0.05 (worse as cell counts become imbalanced), where a random
intercept holds **~0.04** against a calibrated pseudobulk reference of ~0.04.

**The default is therefore `random = "intercept"`.** `"none"` remains available — it reproduces the
original `batch_nichede_v9.R` behaviour and is the back-compatible path — but it should not be used
for inference. Note `nicheDesign()` deliberately keeps `"none"` as *its* default: the random-effect
columns are collinear with the cell-type block and identified only by the penalty applied at fit
time, so a design returned with them included is rank-deficient (rank 23 of 25 on the toy), which is
correct for fitting and surprising from a constructor.
`tests/testthat/test-spiDE-e2e.R` checks niche-*specificity* of the planted G1/A/B signal (is B the
strongest niche association for G1 in index A) rather than asserting G1 has the single largest test
statistic genome-wide, since null genes can outrank it on raw |t|.

`fitSpiDE()`/`spiDE()` take a `random = c("intercept", "none", "slope")` argument that selects the
**mixed-effects correction** for the pseudo-replication (see `vignettes/spiDE-model.Rmd` and the plan
in the PR). Random effects are implemented via the ridge = random-effects equivalence, reusing
`SpaNorm::fitNB`'s per-column `lambda.a` penalty (no SpaNorm change): patient-level random effects are
added as ridge-penalised design columns (tagged `"Random"` in `covtype`), targeting *only* the
response-related fixed effects — a random intercept per sample (counterpart of the `Response` main
effect) and, under `"slope"`, per-sample random slopes on the `CellType:niche` bases (counterparts of
the `ResponseNiche` β terms). `.fitNBmixed()` (`R/mixed.R`) estimates the variance components
`tau2` with a shared-across-genes Schall/PQL loop. `.fitNBmixed()` and the df machinery live in
`R/mixed.R`; the batched/GPU helpers live in `R/inference-batch.R`. That loop is the mixed fit's
dominant cost (it re-fits every gene per iteration), so it is sped up the same
way `fitNB` subsamples cells for dispersion: the inner iterations fit on a
stratified cell subsample (`re.prop`, sampled per cell type × sample with a
`re.min.cells` floor) with a single dispersion iteration (`re.maxit.psi`), then a
**final fit on all cells with full dispersion** supplies the coefficients/`psi`
inference uses — so subsampling only perturbs the shared `tau2`, not the per-gene
effects. Defaults (`re.prop=1`, i.e. subsampling OFF, and `re.maxit.psi=1L`) speed up the mixed fit on
CPU while keeping the response-niche t-stats highly correlated with the full fit
(see `research/reports/benchmarks/spiDE-mixed-benchmark.Rmd`); `re.prop=1` restores the
reproducible, all-cell path. No seed is set internally (set one externally).
`.blockedInference()` (`R/inference.R`) then uses
the **full** penalised covariance `(X'WX + Λ)⁻¹`, the working **Pearson** dispersion (not the NB
`psi`), and a reference df from `SpiDEFit@df` (see `df.method` below) — the three together
are what restore calibration (see `tests/testthat/test-mixedEffects.R`). New `SpiDEFit` slots:
`re_group`, `tau2`, `penalty`, `df` (all `NULL` for a fixed-effects fit). Because a per-sample random
intercept absorbs all between-sample effects, `checkSample()` rejects sample-constant covariates when
`random != "none"`.

### The reference df (`df.method`, default `"satterthwaite"`)

`fitSpiDE()`/`spiDE()` take `df.method = c("satterthwaite", "between")`, used only when
`random != "none"`. **`"satterthwaite"` is the default**; `SpiDEFit@df` is then a *named
per-tested-column vector* (aligned to the columns of `t_stat`/`se`), whereas under `"between"` it is
a *scalar*. Anything reading `@df` must handle both shapes.

- `"between"` tests every `ResponseCellType`/`ResponseNiche` coefficient against the same scalar
  between-sample df `S − 2` — the original back-compatible behaviour, and a misnomer in the
  condition-free case (see the niche-mode note above, and `.fitNBmixed()`'s comments in `R/mixed.R`
  for the per-mode values).
- `"satterthwaite"` derives a df per tested column from the shared variance-component fit
  (`.varParamCov()` / `.satterthwaiteDF()` in `R/mixed.R`), separating between-sample contrasts
  (`Response`: small df, close to `S − 2`) from within-sample ones (`ResponseNiche`: larger df, more
  power).

The default changed on measurement (`research/`, and
`research/reports/benchmarks/spiDE-simulation.Rmd`): `"between"`
is severely over-conservative when samples are few (null type-I ≈ 0.001 at `S = 4` against a nominal
0.05, with near-zero power), while `"satterthwaite"` holds type-I in 0.042–0.065 over the whole
sampled range and gains ≈ 0.10 mean TPR. The trade is a mild liberal drift at larger `S` (worst
measured ≈ 0.065); a Kenward–Roger correction is the indicated next step. An lmerTest oracle check
lives in `tests/testthat/test-satterthwaite.R`.

### `re.maxit`, and a documented failed experiment

`re.maxit` defaults to **2**, lowered from 10 on measurement: for `random = "intercept"` one
iteration is indistinguishable from ten on null type-I (to three decimal places) and on `tau2` (to
two), because the loop converges in a couple of steps. It also largely dissolves a hazard of the
larger cap — `tau2` can enter a 2-cycle, so the answer depends on the **parity** of `re.maxit`.
**That evidence covers the intercept model only**: under `random = "slope"` the slope variance
component decays monotonically across all ten iterations without meeting `re.tol`, so slope fits
should pass `re.maxit = 10`.

A fixed-effects alternative was built and removed: samples coded as `contr.sum` contrasts nested
within condition, each coefficient tested against its own split-plot error stratum. It is calibrated
at **one of eighteen** measured design points and collapses to zero rejections as cells per sample
grow, because the between-sample mean square is the wrong scale for a cell-means condition contrast
— the inflation the contrast needs is constant while that mean square grows as `sqrt(cells per
stratum)`. Do not rebuild it without reading `research/reports/between-sample-stratum.html`, which
records the measurements and the two intermediate findings that *were* correct.

### The two-stage estimator (`twoStageSpiDE()`)

`twoStageSpiDE()` (`R/twostage.R`, stages in `R/twostage-stage1.R` / `R/twostage-stage2.R`) is a
**different estimator**, not a fourth `random` mode — `random = "none"/"intercept"/"slope"` all fit
the *same* design and differ only in which columns are ridge-penalised, whereas this never fits the
niche design at all. `condition` is assigned per **patient**, so patients are the experimental units:
stage 1 estimates a niche slope per (sample, index cell type), stage 2 pools those slopes by
precision within each patient and contrasts the pooled patient slopes. None of `W`, `alpha`, `psi`,
`penalty`, `tau2` or the Satterthwaite df machinery applies, which is why it is not slotted under
`random`. Design spec: `design/specs/2026-08-11-twostage-fixes-design.md`.

**Stage 1** (`.sampleSlopes()`) runs one *joint* weighted fit of the working response on **all** niche
columns per (sample, index) subset — so restricting `niche` changes every remaining slope, not just
which rows are reported — and drops the index type's own niche column, matching the GLM design's
symmetric self-interaction rule. `stage1` selects the working response:
- `"spanorm"` (the default, but **not the best** — see "Which method to use") reads the stored
  `SpaNorm::SpaNorm()` fit from `metadata(spe)$SpaNorm` (`.spanormComponents()` errors clearly when
  absent) and passes its **library-size/batch** components to `fitNB()` as a fixed **offset**
  (`.stage1Offset()`), so spiDE re-models all biology itself and assumes only that the LS effects are
  right. The former `epsilon` argument (`"addback"`/`"residual"`) is **deprecated and ignored**; it
  warns if supplied. Any numbers quoted from an `epsilon`-era run are superseded.
- `"ols"` regresses log-CPM with unit weights — no stored fit and no dispersion needed, so it is the
  path the toy examples/tests use, **and it is the best-performing stage-1 path on measurement**.
- `"nb"` fits a fresh `fitNB` per (sample, index) subset; cost is per-*gene*, so it is a
  small-restriction reference path, not a default.

**Stage 2** pools per-core slopes by `1/v` within patient (`.poolPatientSlopes()`), estimates one
between-patient variance `tau2` per (index, niche) by DerSimonian–Laird pooled as the median over
genes (`.tau2DL()`), and contrasts patients with `limma::lmFit(weights = 1/(v + tau2))` +
`eBayes(robust = TRUE)` (`.limmaStage2()`). `patient.covariates` adjusts at the patient level —
exactly the sample-constant covariates `fitSpiDE()` rejects under `random != "none"`. FDR here is a
**plain BH over triplets**, not the hierarchical cascade.

The returned `SpiDEResults` keeps the tidy `results` schema but has an **empty `@fits`** with `@sigma`
still set; `validSpiDEResults()` was relaxed to key that on the empty list rather than on a third
`mode` value, deliberately — a new `mode` would change what `.testedCols()`/`.nicheTestCols()` treat
as tested. Diagnostics land in the new `@diagnostics` slot: `r2` (niche columns against the SpaNorm
biology basis, `"spanorm"` only), `inclusion` (per-index patient inclusion under `min.cells`, which
*warns* when dropout is associated with `condition`), and `tau2`.

**Three caveats, all load-bearing.** (1) **The old permutation numbers are withdrawn.** The
often-quoted pair (two-stage raw type-I 0.036 with zero false calls, against ~576 for
`fitSpiDE(random = "intercept")`) was measured on this function's *predecessor* and has now been
superseded by a direct measurement that points the **other way** — see "What the niche-shuffle null
showed" below. Do not cite them. (2) **`stage1 = "spanorm"` is the weakest of the three stage-1
paths**, and it is the default only for historical reasons; `"ols"` beats it on calibration, power
and cost (see below). (3) It does **not** solve multiplicity: a full-panel space of ~1.8M triplets
buries real signal, and ACAT over a gene's ~137 mostly-null triplets is no better than Bonferroni.
Restrict `index`, `niche` and the gene set to a pre-specified hypothesis (~4,000 tests is the order
at which a `p ≈ 1e-5` effect survives).
Benchmarked against the published simulation study in
`research/reports/benchmarks/spiDE-twostage-benchmark.Rmd`; its arm lives as extra **rows**
(`method`/`df.method == "twostage"`, labelled by `stage1` and `ls.model`) in the one canonical table
per scenario under `research/reports/benchmarks/tables/`, never a parallel file.

### What the niche-shuffle null showed (2026-08-21, real YTMA cohort)

The sharpest calibration test run on this project so far, and the one that should be repeated before
any future claim of niche-dependent DE. **Permuting the patient label — the older null — cannot
detect a spurious slope**, because it leaves every patient's real niche slopes intact and randomises
only who is a Responder; it tests stage 2 only. The niche shuffle destroys the niche↔expression
association itself, so the true slope is **zero by construction** and every call is a false positive.

Implementation (`research/plasmode/niche_shuffle.R`, `niche_shuffle_glm.R`): permute the **rows of
the niche matrix within (cell type × sample)**. Do *not* shuffle expression across cells — that
breaks each cell's pairing with its SpaNorm offset, whose LS component is position-dependent, and so
tests the offset construction at the same time. Two modes: `free` (random rows) and `block` (toroidal
shift, which preserves local spatial smoothness so architecture/FOV confounding survives). Both were
run; they agree, which rules out unmodelled spatial structure as the driver.

| run (bw 30, valid = Tumor+Fibroblast) | all sd(t) | all frac \|t\|>1.96 | valid sd(t) |
|---|---|---|---|
| REAL two-stage | 1.36 | 0.145 | 1.08 |
| two-stage shuffles (n = 6) | 1.35–1.38 | 0.141–0.149 | 1.06–1.08 |
| REAL intercept GLM | 1.05 | 0.062 | 1.10 |
| intercept GLM shuffles (n = 4) | 0.93–1.04 | 0.036–0.059 | 1.01–1.10 |

Two conclusions, both important:

1. **Neither estimator's real data is distinguishable from its own shuffled null.** Two-stage: 912
   discoveries on real against a null range of 595–986 (package `fdr.niche`); sd(t) 1.36 vs
   1.35–1.38. Intercept: sd(t) 1.100 real against a shuffle maximum of 1.101. On this cohort at this
   resolution there is **no detectable niche-dependent differential expression** for either method.
   A better estimator buys a better-calibrated null, not a finding.
2. **The intercept GLM is far better calibrated than the two-stage estimator across all index
   types** (sd(t) 1.04 vs 1.36; 5.9% vs 14.5% exceeding \|t\|>1.96). It fits all cells jointly, so
   thin cell types borrow strength, where two-stage fits each (patient, index) subset independently.
   In the *valid* subset the two are equally mildly liberal (~1.06–1.10), so the GLM's advantage is
   entirely in the thin index types.

### Which method to use

Measured, not assumed. Simulation numbers are from the structured-LS sweep at 85 of 99 parts
(`research/plasmode/summary/twostage_*.csv`), so they may still move slightly.

**Default to `fitSpiDE(random = "intercept")` when samples are plentiful (S ≥ 16).** Raw power 0.660
at S = 30 against 0.451 (`ols`), 0.293 (`nb`), 0.274 (`spanorm`); FDP 0.041 at a nominal 0.05 by
S = 16; and the best real-data shuffle calibration above.

**It fails badly at small S**: FDP **0.35** against a nominal 0.05 at S = 4 (0.33 against a nominal
0.01), recovering to 0.04 by S = 16. Note this does *not* show up in null type-I (0.049 at S = 4,
respectable) — the two measure different things, and reading only the null table would have declared
it fine. Any sub-analysis that thins the patient count re-enters this regime.

**If using `twoStageSpiDE()`, prefer `stage1 = "ols"` over the `"spanorm"` default**: better
calibrated (0.043–0.054 vs 0.058–0.071), **1.6× the raw power** (0.451 vs 0.274 at S = 30), much
better precision (FDP 0.032 vs 0.173 at S = 16), and it needs no stored SpaNorm fit so it is the
cheapest. `"nb"` is the most conservative and the least powerful of the three.

**Keep `pool.psi = TRUE`.** A paired ablation on identical datasets
(`research/plasmode/poolpsi_ablation.R`, 30 pairs) gives type-I 0.0673 pooled vs 0.0729 unpooled,
consistently lower at every S, paired p < 0.0001. Small but unambiguous.

### Two-stage calibration is governed by cells per subset

Measured on the real cohort at bw 30: `sd(t)` per index type tracks **median cells per (patient,
index) subset** at r = **−0.838**, against r = −0.659 for patient count. Tumor (388 cells/subset)
is calibrated at sd(t) 1.05; Mast (38 cells) is inflated at 1.65. Ruled out as drivers: **bandwidth**
(sd(t) 1.363 at σ = 30 vs 1.364 at σ = 70 — essentially identical) and **patient count** (restricting
to ≥ 30 patients moved 1.363 only to 1.277). The mechanism is ~12 niche columns fit to as few as
30–40 cells: a near-singular design whose Fisher variances understate uncertainty. The lever is
therefore `min.cells` or niche resolution (p/n), not σ.

A flat `min.cells` that reaches calibration is expensive on this cohort: 30 → 11 index types, 75 → 7,
100 → 5, 150 → 4, 200 → 3. Reaching 100 deletes the entire T cell / DC / Monocyte / Mast compartment.

**Separately, dropout can be confounded with condition.** `@diagnostics$inclusion` warns on this and
it fires here: B cell (Fisher p = 0.027), DC (0.006), Monocyte (0.016) — which patients contribute
depends on their outcome. No threshold or variance correction fixes informative missingness. Only
Tumor and Fibroblast have complete inclusion (55/55) **and** calibrated variance.

### Supplying `psi` changes the objective (winsorisation), it does not just save time

`fitNB()`'s `winsor` caps counts against the **current fitted mu**, so supplying `psi` changes the
IRLS weights, hence mu, hence *which cells are capped* — the estimating and supplied-psi paths
maximise **different** objectives. On spiDE's near-singular stage-1 subsets the gap reaches 64
log-likelihood units at `winsor = 4` and closes to a coin flip at `winsor = Inf` (median +0.04). On a
well-conditioned synthetic design it does **not** close, so conditioning is also involved and the
mechanism is not fully characterised. `research/notes/fitnb-offset-psi-disagreement.R` prints both
regimes including its own counter-example.

Consequences: (a) neither path "under-converges" — a two-pass refit was built and **reverted** on
measurement (1.70× cost, median −0.41 loglik, better in only 2 of 6 subsets); (b) validate any
shared-dispersion speedup on **type-I/power**, never on likelihood, which cannot rank two different
objectives.

Related: `fitNB()` subsamples cells for dispersion above a size threshold and sets **no seed
internally**. At 60 cells repeated fits agree exactly; at 600 cells two identical calls differ by
max|d psi| = 0.245. Any script whose numbers will be cited must `set.seed()` and record it.

### Every per-gene matrix inversion is guarded

A singular per-gene information matrix once aborted an 80-minute real-cohort run: the matrix is built
from each gene's own weights, so it is singular for one gene of 13,348 (rcond 2.7e-20) and fine for
the rest — the same design on a 200-gene subset completed. Six sites were unguarded, including the
default Cauchy path in `.blockedInference()` that every `fitSpiDE()` call reaches. Treatment is
matched to the site, not uniform:

- per-gene sites (`R/inference.R`) → the gene drops out as `NA`; `p.adjust()` already ignores `NA`.
- batched site (`invert_mat_batched`, which a grep for `invert_mat(` **misses**) → diagnostic error
  naming `cov.batch`/`backend` as the levers.
- shared sites (`R/mixed.R`, weights averaged across genes) → `tau2` stops with a diagnosis, since a
  bad variance component would corrupt every gene's inference; the Satterthwaite df degrades to the
  conservative `between` df with a warning, failing toward validity.

### GPU: usable for the fit, broken for the inference

On this cluster's torch build, `fitSpiDE(backend = "gpu")` completes and is ~4× faster than CPU
(16.3 min vs 69–100 min on the real cohort), but `testSpiDE()` then dies in a torch CUDA kernel JIT
(nvrtc), consistent with the known nvrtc-builtins soname gap. `research/plasmode/gpu_smoke.R` passes
because it exercises specific fp64 kernels, not the blocked-inference path — so **the gate passing
does not imply inference will run**. h100 only regardless: a100 fails fp64 `digamma` *inside* NB
dispersion, and l40s runs fp64 at ~1:64.

### Gene-set inference (`spiGSEA()`)

`spiGSEA()` (`R/spiGSEA.R`) averages per-gene statistics over a set, converting **t to z before
averaging** and inflating the variance by `sqrt((1 + rho(m - 1))/m)` for inter-gene correlation
(`rho` estimated from the counts, or reused from the fit). `test = "competitive"` (vs genes outside
the set, camera-style) is the **default**: the `"self-contained"` form is not calibrated — on a null
benchmark it called 20.6 of 208 sets per replicate, all false (FDP 1.00), against 0.05 for
competitive, because it assumes the averaged z have unit spread and they do not under signal. It is
kept only to reproduce the flat script's `fry_res`. `type = "celltype"` errors in niche mode.

### The toy fixture's effect size

`.toySPE()` plants `log_effect = beta * (x / field)`, so `beta` is the **maximum log-fold-change
across the field**, not a linear signal knob. G1's dynamic range inflates its own estimated NB
dispersion, hence its standard error, so the response is **non-monotonic** — measured t of 2.22,
2.87, **5.28**, 2.49, 1.42, 0.25 at beta 1, 1.5, **2**, 2.5, 3, 4. The default sits at the peak.
Raising it makes the planted effect *harder* to recover, not easier.

### Where the evidence lives

Most non-obvious defaults in this package were chosen on measurement, and the measurement is written
down. Before changing one, read the corresponding record:

- `vignettes/spiDE-model.Rmd` — the model, the pseudo-replication problem, the random-effects and df
  machinery.
- `research/reports/benchmarks/` — the five validation reports (simulation study, combiner,
  mixed-fit speedups, two-stage arm, spiGSEA calibration), rendered to a static site at
  `research/docs/` by `build_site.R` (https://bhuvalab.github.io/spiDE-research/). They moved out of
  `vignettes/` because the eight built vignettes alone exceeded Bioconductor's 10 MB tarball cap;
  the package keeps only the quickstart, model, and calibration vignettes.
- `research/reports/benchmarks/tables/*.rds` — the canonical benchmark tables the reports read
  (they render without the HPC runs; refreshed by `research/R/install_results.R` and
  `research/plasmode/install_twostage.R`). **One canonical table per scenario**: a new method arm is
  extra *rows*, not a parallel file that would carry a stale copy of the others.
- `research/` — a git submodule (`BhuvaLab/spiDE-research`) holding the benchmark harness and the
  written-up negative results (e.g. `research/reports/between-sample-stratum.html`).
- `design/specs/` — design specs and implementation plans for larger changes.
- `research/notes/fitnb-offset-psi-disagreement.R` — runnable, self-contained; prints both regimes
  of the winsorisation/psi finding *including its own counter-example*.
- `research/plasmode/niche_shuffle.R` + `niche_shuffle_glm.R` — the shuffle null for both estimators
  (`SPIDE_SHUF` = `free` | `block` | `none`, the last being the real-data comparator).
- `research/plasmode/poolpsi_ablation.R` — the paired `pool.psi` ablation.
- `.claude/skills/calibration-check/` — scores any results table for calibration, per-index
  breakdown, dropout-vs-condition confounding, and re-FDRs the valid subset. **Run it before quoting
  any discovery count.**
- `.claude/agents/` — `numerical-robustness-reviewer` (unguarded inversions, p/n conditioning,
  unseeded stochastic paths) and `evidence-auditor` (numeric claims vs the canonical tables).

**A warning about the benchmark harness.** A sweep whose tasks load `SPIDE_HOME` from the **live
working tree** is not one experiment: tasks start at different times, so edits mid-run give different
tasks different code. A previously reported result ("two-stage null inflation grows with S, 0.078 →
0.127") came from such a run and **did not survive re-measurement on a frozen snapshot** — the paired
ablation shows no trend with S in either configuration (p = 0.82 / 0.62). Always pin `SPIDE_PKG` to a
snapshot (`.claude/skills/run-benchmark-arm/scripts/freeze_snapshot.sh`).

### Checkers

`R/checkers.R` holds shared input validation (`checkSPE`, `checkCondition`, `checkCovariates`,
`checkNiche`, `checkCounts`) called at the top of every exported entry point — extend these rather than
duplicating validation logic in a new function.
