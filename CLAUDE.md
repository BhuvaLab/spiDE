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
- `"spanorm"` (default) linearises the one-step response at the stored `SpaNorm::SpaNorm()` fit read
  from `metadata(spe)$SpaNorm` (`.spanormComponents()` errors clearly when absent), and under
  `epsilon = "addback"` (the default) adds the fitted biology component back so only the
  library-size/batch part is removed; `epsilon = "residual"` leaves the bare working residual, which
  under-states the slope when a niche column overlaps the biology basis.
- `"ols"` regresses log-CPM with unit weights — no stored fit and no dispersion needed, so it is the
  path the toy examples/tests use.
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

**Two caveats, both load-bearing.** (1) The often-quoted permutation numbers (raw type-I 0.036, zero
false calls where `fitSpiDE(random = "intercept")` returned ~576) were measured on this function's
**predecessor** (the nbresid/Welch estimator, since replaced by the SpaNorm-anchored joint one);
re-measurement on the current implementation is pending — don't quote them as current. (2) It does
**not** solve multiplicity: a full-panel space of ~1.8M triplets buries real signal, and ACAT over a
gene's ~137 mostly-null triplets is no better than Bonferroni. Restrict `index`, `niche` and the gene
set to a pre-specified hypothesis (~4,000 tests is the order at which a `p ≈ 1e-5` effect survives).
Benchmarked against the published simulation study in
`research/reports/benchmarks/spiDE-twostage-benchmark.Rmd`; its arm lives as extra **rows**
(`method`/`df.method == "twostage"`, labelled by `stage1` and `ls.model`) in the one canonical table
per scenario under `research/reports/benchmarks/tables/`, never a parallel file.

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

### Checkers

`R/checkers.R` holds shared input validation (`checkSPE`, `checkCondition`, `checkCovariates`,
`checkNiche`, `checkCounts`) called at the top of every exported entry point — extend these rather than
duplicating validation logic in a new function.
