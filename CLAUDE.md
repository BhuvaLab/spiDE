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

There is no lint config (`.lintr`) or CI workflow in the repo yet — lintr diagnostics surfaced by the
editor reflect default rules, not a project-specific config. Dot-separated argument names like
`lambda.a`, `winsor`, `maxit.psi` are intentional: they mirror `SpaNorm::fitNB`'s own argument names
so they can be forwarded via `...` without renaming.

## Architecture

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
     study (`vignettes/spiDE-cauchy-vs-brown.Rmd`) showed it matches or beats Brown while controlling
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

### S4 classes

- `SpiDEFit` (`R/AllClasses.R`) — one bandwidth's fit + inference (design `W`, per-column `covtype`
  tags, `coefmap`, per-gene `alpha`/`psi`/`loglik`, and once inferred: `t_stat`, `se`,
  `p.combined.pos`, `p.combined.neg` — the last two hold the within-gene combined p-values from
  whichever combiner (`"cauchy"` default / `"brown"`) was used; they are a
  `genes × (1 + n_index)` matrix, column `"Gene"` then one per index cell type).
- `SpiDEResults` (`R/AllClasses.R`) — container for a list of `SpiDEFit` (one per bandwidth) plus
  cross-bandwidth combined p-values and the final tidy `results` data.frame. Both classes expose a
  `$` accessor (`slot(x, name)`) and a `show` method; validity is enforced via `setValidity()`.
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

### Known statistical property (not a bug) and the mixed-effects correction

The **default** (`random = "none"`) reduced-design Wald standard errors make the method
anti-conservative on data with few samples but many cells (cell-level pseudo-replication) — this
faithfully reproduces the original `batch_nichede_v9.R` behaviour and is the back-compatible path.
`tests/testthat/test-spiDE-e2e.R` checks niche-*specificity* of the planted G1/A/B signal (is B the
strongest niche association for G1 in index A) rather than asserting G1 has the single largest test
statistic genome-wide, since null genes can outrank it on raw |t|.

`fitSpiDE()`/`spiDE()` take a `random = c("none", "intercept", "slope")` argument that turns on a
**mixed-effects correction** for the pseudo-replication (see `vignettes/spiDE-model.Rmd` and the plan
in the PR). Random effects are implemented via the ridge = random-effects equivalence, reusing
`SpaNorm::fitNB`'s per-column `lambda.a` penalty (no SpaNorm change): patient-level random effects are
added as ridge-penalised design columns (tagged `"Random"` in `covtype`), targeting *only* the
response-related fixed effects — a random intercept per sample (counterpart of the `Response` main
effect) and, under `"slope"`, per-sample random slopes on the `CellType:niche` bases (counterparts of
the `ResponseNiche` β terms). `.fitNBmixed()` (`R/fitSpiDE.R`) estimates the variance components
`tau2` with a shared-across-genes Schall/PQL loop. That loop is the mixed fit's
dominant cost (it re-fits every gene per iteration), so it is sped up the same
way `fitNB` subsamples cells for dispersion: the inner iterations fit on a
stratified cell subsample (`re.prop`, sampled per cell type × sample with a
`re.min.cells` floor) with a single dispersion iteration (`re.maxit.psi`), then a
**final fit on all cells with full dispersion** supplies the coefficients/`psi`
inference uses — so subsampling only perturbs the shared `tau2`, not the per-gene
effects. Defaults (`re.prop=0.2`, `re.maxit.psi=1L`) speed up the mixed fit on
CPU while keeping the response-niche t-stats highly correlated with the full fit
(see `vignettes/spiDE-mixed-benchmark.Rmd`); `re.prop=1` restores the
reproducible, all-cell path. No seed is set internally (set one externally).
`.blockedInference()` (`R/inference.R`) then uses
the **full** penalised covariance `(X'WX + Λ)⁻¹`, the working **Pearson** dispersion (not the NB
`psi`), and a **between-patient** reference df (`S − 2`, stored in `SpiDEFit@df`) — the three together
are what restore calibration (see `tests/testthat/test-mixedEffects.R`). New `SpiDEFit` slots:
`re_group`, `tau2`, `penalty`, `df` (all `NULL` for a fixed-effects fit). Because a per-sample random
intercept absorbs all between-sample effects, `checkSample()` rejects sample-constant covariates when
`random != "none"`.

### Checkers

`R/checkers.R` holds shared input validation (`checkSPE`, `checkCondition`, `checkCovariates`,
`checkNiche`, `checkCounts`) called at the top of every exported entry point — extend these rather than
duplicating validation logic in a new function.
