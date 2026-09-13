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

Build the pkgdown site with the article sync first, never `pkgdown::build_site()` alone:

```sh
Rscript vignettes/articles/pkgdown-sync-articles.R && Rscript -e 'pkgdown::build_site()'
```

The validation reports are copied from the research submodule into `vignettes/articles/`
(gitignored) by that script, which takes the list from the articles index in `_pkgdown.yml`
and deletes any copy not in it. pkgdown refuses an article on disk that the index does not
list (`vignette missing from index`) and one the index lists that is not on disk, and it has
no ignore mechanism; a copy synced before a report was retired (the two-stage benchmark,
0.99.19) survives every `git pull` on that clone until the sync runs. `vignettes/articles/README.md`
has the error text and the remedy.

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

**Entry points.** The pipeline is **fit → polish → test → gsea**: `buildNiches()` →
`fitSpiDE()` → `polishSpiDE()` → `testSpiDE()` → `spiGSEA()`; `spiDE()` chains the first four
(`polish = TRUE` by default). `polishSpiDE()` (`R/polish.R`) is a *stage*, not an option of the
fit: it converges each gene, sets its dispersion (`psi = "profile"` by default) and, for a mixed
fit, re-estimates the variance components from the converged coefficients (`tau2 = TRUE`),
refreshing `@tau2`, `@penalty` and the Satterthwaite `@df`; it clears the inference slots. The
user runs it or skips it according to their data (it needs integer counts). The former two-stage
estimator (`twoStageSpiDE()`) left the package in 0.99.19 and is archived as the standalone
research package `spiDEtwostage` (`research/twostage/`); the mixed-effects model is the
recommended approach.
`compositionTest()` (`R/composition.R`) is a *different question*: the between-sample
composition association that `fitSpiDE()`'s nested intercept deliberately absorbs, tested at the
patient level — pseudobulk per (sample, index type), the sample's mean niche density around those
cells, `limma` across samples, with `"niche"` (pooled) and `"condition:niche"` (the patient-level
counterpart of the three-way term). It is real signal and it is not niche-dependent DE; never
report one as the other.

### Pipeline (four stages; `spiDE()` chains fit, polish and test)

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
   bandwidth (S4 class in `R/AllClasses.R`), collected into a `SpiDEResults`. The fit converges
   nothing per gene: that is the next stage.

2b. **Polish** — `polishSpiDE()` (`R/polish.R`). Per gene, damped Newton on the gene's own
   penalised NB log-likelihood from a sane start (blocked and dispatched like inference), the
   dispersion by profile ML at the converged mean (or `psi = "moderated"` to keep `fitNB`'s), then
   for a mixed fit the variance components re-estimated from the converged fit: one Schall step
   (`.schallStep()`, shared with the fit's loop) with the polished coefficients and the
   gene-averaged weights at the polished mean and dispersion, a re-polish at the new penalty,
   iterated to `tau2.tol`, and `.satterthwaiteDF()` recomputed. **Why the re-estimate exists**
   (`research/fdr-ordering/FINDINGS.md`, 2026-09-08): the fit's loop reads the shared fit's own
   unconverged coefficients, and on `.toyClustered()` that reported a between-sample variance of
   10 against a planted 0.49 (sample intercepts three times too wide, dispersion fifteen times too
   large) under every release since 2026-08-08, while calibration survived because the QL scale is
   robust to the dispersion and an inflated component only weakens the ridge. The long test
   `longtests/testthat/test-mixed-numerics.R` and `tests/testthat/test-polish-stage.R` hold the
   stage to the planted window.

   **Its cost is one cold pass, and the rest was measured away (2026-09-10).** The 0.99.19
   cohort runs spent 620-700 min per task on this stage (769 genes, four workers): a cold pass
   of 245-395 min and three re-polish passes of 110-168 min each. Three causes. (1) **BLAS
   oversubscription**: the cohort driver sets OpenBLAS to `OMP_NUM_THREADS` for the fit and
   then forks `NCPU` polish workers, each inheriting that count -- 4 x 8 threads on 8 cores.
   At the cohort's design shape, 4 workers x 4 threads on 4 cores take 2.3-2.6 s per Newton
   step against 0.24-0.28 s at one thread each (9x), and one worker gains only 1.5x from four
   threads (the gram is memory-bound). Both blocked stages (`.polishFit()` and `.blockedInference()`) now dispatch through
   `.bplapplySingleBLAS()`, which sets one BLAS thread inside each worker (RhpcBLASctl in
   Suggests); any driver that forks workers must do the same for code that predates it. (2) **The re-polish is warm** (`.polishGene(warm =
   TRUE)`): a few Newton steps at the held dispersion, no psi search, no restart check (that
   check threw ~100 of 769 genes back to the sane start on every re-polish at bandwidth 10).
   (3) **The loop converges**: Schall's map is linear and slow for the nested block (6-12% above
   its limit after the old cap of three, in 13 of 15 runs), so `.tau2Iterate()` is
   Steffensen-accelerated with `tau2.maxit = 10`, `tau2.tol = 1e-2`, and a final warm pass at
   the reported penalty. The tolerance is set where the estimate stops mattering: a 5% change
   in a component moves individual `t` by at most 0.05, a tenfold error in a near-zero nested
   one by 0.13. A fixed-point stop does not exist for a component that is truly zero (the
   default `.toyClustered()`), where the map creeps toward the floor sublinearly and the cap
   decides; `.toyClustered(sd_nested = 0.3)` gives an interior one. The gene-subset estimate
   of `tau2` (one scalar per component, a mean over genes) is the next lever at transcriptome
   scale and is not built. Record: `research/fdr-ordering/FINDINGS.md`, 2026-09-10.

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
  `results.celltype`/`results.patient` behind `results(type = )`). Both classes expose a `$`
  accessor (`slot(x, name)`) and a `show` method; validity is enforced via `setValidity()`, one
  `SpiDEFit` per bandwidth.
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
re-checking every bandwidth still converges. `.toySPE(composition = k)` plants a **between-sample
composition confound** (each sample's A cells shifted toward or away from the B-rich region, G2
shifting with that in Responders only, zero within-sample slope) — but at the default 80 cells the
niche covariate's per-sample mean is placement noise (between-sample sd 0.07 vs within 0.25), and
no estimator can see a between-sample effect the covariate does not carry. Use `n_per >= 200` and
bandwidth 50 (ratio 0.7): then `compositionTest()` sees it at t ≈ 4.6 and the nested intercept
takes the GLM's leakage to 0.00. This cost three fixture redesigns to learn; do not "fix" it by
lowering thresholds.

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
(see the *What was tried and rejected* report, `research/reports/benchmarks/spiDE-rejected.Rmd`); `re.prop=1` restores the
reproducible, all-cell path. No seed is set internally (set one externally).
`.blockedInference()` (`R/inference.R`) then uses
the **full** penalised covariance `(X'WX + Λ)⁻¹`, a per-gene **dispersion scale** (the
quasi-likelihood dispersion by default from 0.99.18, `dispersion = "ql"`; the working Pearson
dispersion as `"pearson"`; never the NB `psi`), and a reference df from `SpiDEFit@df` (see
`df.method` below) — the three together are what restore calibration (see
`tests/testthat/test-mixedEffects.R`). New `SpiDEFit` slots:
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
stratum)`. Do not rebuild it without reading the *What was tried and rejected* report
(`research/reports/benchmarks/spiDE-rejected.Rmd`), which records the measurements and the two
intermediate findings that *were* correct.

### The two-stage estimator is archived (0.99.19)

`twoStageSpiDE()` and its stages (`R/twostage*.R`) left the package in 0.99.19 and live as the
standalone research package `spiDEtwostage` (`research/twostage/`, imports spiDE for the niches and
checkers, returns a `twoStageResults` list with a `results()` method on spiDE's generic). It was a
*different estimator*, not a `random` mode: patients as units, a niche slope per (sample, index)
in stage 1, precision-pooled and contrasted with `limma` in stage 2, plain BH over triplets. It is
kept for reproducing the comparison on the research site (*The two-stage estimator* report and its
rows in the canonical tables, `method == "twostage"`); the mixed-effects model is the recommended
approach because it fits all cells jointly, so thin cell types borrow strength, where the two-stage
estimator was calibrated only in the populous types. Its calibration findings (cells per subset at
r = −0.84, dropout confounded with condition in B cell / DC / Monocyte, `stage1 = "ols"` best of
its three stage-1 paths) are recorded in that report.

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

### The null tail is per-gene, and the lever is a gene filter

`research/fdr-ordering/FINDINGS.md` (addendum, 2026-08-31) localises the FDR failure. The heavy null
p-value tail that breaks BH is **per-gene scale heterogeneity**, not per-column: standardising each
gene's `t` by its own null sd takes the `p < 1e-6` excess from **59.9x to 2.6x**, where a per-column
rescale reaches only 37.3x. Per-gene `sd(null t)` spans 0.837–2.292 against 0.908–1.245 per column,
and **90.8% of `|t| > 4.89` exceedances come from the top 5% of genes**.

The affected genes are the **highest-expressed** ones (`cor(sd(t), log median SE) = -0.60`; worst are
HLA-DPA1 2.29, CDV3, AEBP1, CST3, LAPTM5). This survives `free` shuffles, so it is not spurious
spatial regression, and it survives the per-gene **Pearson working dispersion already applied** at
`R/inference.R:539` — so "add a per-gene dispersion" is not by itself the fix. Twelve candidate
cures have been measured and refuted (ledger in `research/fdr-ordering/REPORT.md` §4): the
`calculateMu()` winsorisation clamp, a patient-clustered sandwich, an edgeR-v4 QL dispersion (which
*steepened* the gradient at the time — that refutation was scored on the composition bias below and
is withdrawn; the QL scale is the 0.99.18 default, see "The dispersion rule and the SE scale"),
cascade ordering, within-gene variance misspecification (it makes the
coefficient *conservative*), spatial autocorrelation, shared-weight IRLS inefficiency, the extreme
fits themselves, a cell-type-level dispersion and a cell-level HC0 sandwich at the converged fit.

**The production fit is not at its optimum for bright genes (2026-09-03), and that is a separate,
real defect.** `SpaNorm::fitNB` fits all genes with one gene-averaged cell-weight vector, an
aggregate convergence criterion and a cross-gene coefficient clamp; for the top 5% of genes by
expression the coefficients sit **1–4 production SEs** from the gene's own penalised-NB optimum
(log-likelihood gaps of 10^4–10^6, minimum fitted `log mu` below −20 for the brightest 35), and the
edgeR `psi` estimated at that point is **1.6× too large**. Polishing *from* the production point
diverges; damped Newton from a sane start (cell-type log means, `loglib` slope 1) converges every
gene in 5–29 iterations at ~3 s per gene (`research/fdr-ordering/R/converged_null.R`). **Converging
the fit does not fix the null**: on the same shuffle grids it calibrates ordinary genes exactly
(`sd(null t)` 1.002 vs 0.929) and *lifts* the brightest band from 0.22–0.88 to **1.25–1.42**, because
the production fit's inflated Pearson `phi` had been hiding inflation behind deflation. The residual
inflation is confined to bright genes, tracks the **number of cells in the index type** (Tumor 1.36,
B cell 1.24, Fibroblast 1.15, every type under ~3,000 cells at 1.00–1.03), is present under `free`
shuffles and larger under `block`, and is not any marginal-variance quantity: per-cell-type Pearson
scale, μ-tercile Pearson ratio and the cell-level HC0 sandwich all predict 0.92–1.0. Any future
per-gene inference should first converge the fit; the honest null variance of bright genes is then
*higher*, not lower. See `research/fdr-ordering/FINDINGS.md` (2026-09-03 entry) and `REPORT.md` §5c.

**The cause (2026-09-03, `REPORT.md` §5d): the design has no (sample × cell type) intercept.**
Holding a gene's converged fit fixed and permuting the niche rows within (sample, cell type) gives a
`t` whose spread is the HC0 prediction but whose per-column **mean** runs to ±3 and equals the
real-data `t` column by column (r = 0.66). The free shuffle preserves each group's mean niche density,
and with only a shared per-sample intercept a between-sample association between a cell type's mean
niche density and its mean expression in that type loads onto the `CellType:niche` and
`CellType:condition:niche` slopes — a patient-level composition effect (S = 55) reported with a
cell-level SE. The excess is a *bias of the estimand*, not a variance, which is why every
variance-side candidate failed; it scales with cells per index type and is confined to bright genes.
Centring the niche-dependent columns within (sample, cell type) removes it (bias 1.5 → 0.09, RMS of
the null `t` 1.83 → 0.99 = HC0, in every index type) and collapses CDV3's real Tumor `t` from 1.69 to
0.82 RMS. **Validated on all twelve raw grids**: the centred, converged null is flat at 0.96–0.99 in every
expression band under `free` shuffles (gradient 1.000, per-gene spread 0.92–1.05, zero `|t| > 4.89`
exceedances against 95 for production); `block` shuffles keep 1.10–1.19 in the ~100 brightest genes,
the spatially-smooth-covariate component, so calibrate against `block` on the fixed design.
**Confirmed on the packaged fix** (job 27965091, 0.99.17): `free` is flat at 0.96-0.98 in every band
with 0 exceedances of 101,508, `block` keeps 0.988 -> 1.206 with 10, so the nested intercept does not
remove the spatial component. **REPLICATED (5 block grids, jobs 27965091/27965201): the real cohort
now exceeds its null, reversing the "no detectable niche-dependent DE" conclusion.** Real gives 97
exceedances of `|t| > 4.89` against a null range of 5-11; the 100 random control genes sit INSIDE the
null range while every expressed band sits outside it, so it is not a residual scale artefact.
Per-**gene** calibrated (the granularity this investigation established — per-index pools over genes
and under-calibrates the bright ones), **104 calls at empirical FDP <= 0.05** on the reviewed code
(84 before the mean-function fix), concentrated in B cell and Tumor. **Do not quote that count as a
set of findings**: three artefact arms (2026-09-04) show the estimator is calibrated everywhere
tested (null `sd(t)` 0.99–1.03 across bandwidths, compartment definitions and covariate sets) but
the *identity* of the calls is not stable. Splitting Plasma out of the B cell compartment keeps the
count (70) and loses every immunoglobulin call, so the earlier plasma-activity reading was an
artefact of the merge. Five imaging covariates keep 84 of 104 and drop the epithelial-in-immune
calls (KRT7, KRT17), though collagens in B cells survive and per-cell morphology cannot capture
neighbour transcript bleed. Bandwidths 10/50/70 give 75/63/61 calls of which only **4** are shared.
The one triplet robust to every perturbation is **NDRG1 in Tumor against a Fibroblast niche**. The fix is *why*: the old confound was present in the real data
and in every shuffle, so both inflated equally and matched; removing it from both leaves a
difference. Segmentation spillover — confounded with the niche covariate BY CONSTRUCTION, since both
scale with neighbour density — was tested and is not supported (called-UP genes are *depleted* in the
neighbouring type, median log2 −0.84, where spillover predicts enrichment). Before any publication:
rerun on the full transcriptome (this panel is enriched for the pathological tail), validate the
biology independently, and check imaging artefacts beyond simple spillover. Scripts:
`research/fdr-ordering/R/score_pkgfixed.R`, `R/score_real_calibrated.R`.

**Fixed in 0.99.17.** `fitSpiDE(re.celltype = TRUE)` (the default) adds the ridge-penalised
(sample × cell type) intercept block in `.buildRandomEffects()`, tagged `Random` with
`re_group = "SampleCellTypeInt"` and carrying its own `tau2` (the Schall loop and
`.satterthwaiteDF()` pick it up unchanged; the nested df is checked against `lmerTest` in
`tests/testthat/test-satterthwaite.R`). The per-gene convergence was `fitSpiDE(converge = TRUE)`
in 0.99.17–0.99.18 and is the `polishSpiDE()` stage from 0.99.19 (see "Pipeline"), recording
diagnostics in `SpiDEFit@polish`. `random = "slope"` does **not** substitute for the nested
block. Specs: `design/specs/2026-09-03-sample-celltype-intercept.md` (the defect and its
validation) and `design/specs/2026-09-04-nested-intercept-and-convergence-design.md` (the
implementation), plan in `design/plans/`.

Converging also **sharpens real signal**, which was not why it was built: on the toy fixture the
unconverged 0.99.16 fit put the planted effect at `t` 1.68 in condition mode, and in niche mode a
spurious competing niche *outranked* the true one (|t| 7.82 vs 5.63). Under the 0.99.18 defaults
(measured 2026-09-08, `data(toySpiDE)`, default bandwidths) the planted G1/A/B effect reaches `t`
6.05 at its best bandwidth (3.72 / 5.72 / 5.45 / 6.05 at 10 / 30 / 50 / 70), and in niche mode the
true niche B is the strongest at 5.57 while the former competitor C sits at 2.15. The 0.99.17
profile-`psi` figures (10.19 and 14.27) were larger because the profile dispersion (0.28 against the
moderated 3.09 for G1) shrank the Pearson-scaled SE; the QL scale sees the residual scatter
directly, so the two `psi` rules now give the same statistic (5.7-5.8 at bandwidth 30 under every
combination). Two unit tests had encoded the 0.99.16 artefacts -- one asserting the argmax of a raw
coefficient, which a near-empty gene can win, the other asserting that the spurious call survives
FDR -- and assert the statistic instead.

Three things to know about the implementation. The polish stage re-estimates the dispersion per
gene by profile ML (`psi = "profile"`, the default) and can keep edgeR's cross-gene moderated one
(`psi = "moderated"`, the cheaper rule); the two are indistinguishable on the null and in the
cohort's calls, but the moderated value is whatever the shared fit left (fifteen times the
converged value on `.toyClustered()`), and the variance-component step needs the converged one. The nested indicator block
is absorbed by a Schur complement inside `.newtonSolver()`, so the per-gene Newton cost is one
dense-column gram regardless of how many groups exist — but `.blockedInference()` still forms
a **dense** per-gene gram over the full design, so with ~660 extra columns real-cohort
inference is ~8× slower; absorbing it there is deferred to its own spec. And two corollaries
of the finding itself: the shuffle null is a complete null only for within-group slopes, so under the shipped
design it carried the same confound as the real data (why real and null were indistinguishable); and
the between-sample association is real and should be tested at the patient level, not reported as
niche-dependent DE. The two-stage estimator is within-sample by construction, which is why its Tumor
index was calibrated where the GLM's was not.

Two levers were measured on the complete null (flat BH, false calls at alpha .05):

| restriction | tests | BH .01 | BH .05 |
|---|---|---|---|
| all 12 index x 13,348 genes | 1,815,328 | 116.5 | 223.5 |
| **Tumor + Fibroblast only** | 293,656 | 127.5 | **294.0** |
| random 500-gene panel, all index | 68,000 | 1.5 | **2.5** |

**Restricting index cell types makes the null worse**, because BH's `alpha*R/m` threshold rises as
`m` falls while the pathological genes stay: the Tumor+Fibroblast restriction is a
*calibration-of-the-estimator* argument, not a multiplicity one. Per-gene calibration alone is
necessary but **not sufficient** (23% fewer false calls, 8–18 recall points lost, `P(>=1)` still 1)
because within-gene correlation remains.

**The gene filter is worth applying, but it is a 1.5x effect, not 6x.** Dropping the genes whose
shuffle `sd(t)` exceeds 1.3 cuts complete-null false calls **1.5x at alpha .05 and 1.7x at .01**
when the hot-gene list is built from shuffle grids *other* than the one scored. An earlier 6x/18x
figure was **circular** — the filter was defined on the grid it was scored on — and is withdrawn;
it reproduces exactly (6.7x) when the circularity is reinstated. Two corollaries: about two thirds
of any single grid's hot list is that grid's own noise (126 genes on their own grid vs 38 on a
five-grid average), so most of the extreme tail is **not** a stable gene property; and the filter
costs **no power** (FDP falls 0.05–0.07, TPR moves in the third decimal, worst loss 0.011). At 1.5x
it does not remove the need to fix the variance estimator.

### There is no library-size term, and that is a real gap

`fitSpiDE()` has **no offset argument** and `R/design.R` has **no size-factor handling** — per-cell
sequencing depth is modelled only through the `CellType` cell-means intercepts, which are constant
within a cell type. The YTMA v10 write-up gives the reason ("library size is constant within a
patient, and a per-patient random intercept absorbs between-patient variation"), which is right for
patient-level covariates and **wrong for library size**: `nCount_RNA` is per *cell*, median 1732,
IQR 1108–2698, max 27,553.

Measured on the shuffle null (2026-09-01, `research/fdr-ordering/`), on **raw counts** (see below):
`base` 1.207 → `+library size` **1.149**, and median `sd(null t)` 1.055 → **1.008**. That is the
**best-calibrated configuration measured anywhere in this investigation** — every other arm is off
in one direction (adjusted counts run 0.89–0.92, conservative; raw without a depth term 1.055,
liberal). The gradient is not eliminated, so depth is a major contributor and not the only one.

The slope the data want on raw counts is **0.979** (IQR 0.941–1.021, positive for all 13,348 genes),
i.e. a **conventional offset**, now empirically justified. On the SpaNorm-adjusted assay it is 0.679
— attenuated, because that adjustment has already absorbed part of the depth effect.

**Cell-type-specific size factors are NOT appropriate.** On raw counts every cell type lands at
0.939–0.983 (spread 0.044) and modelling them separately moves the gradient 1.149 → 1.144. The
0.458–0.718 spread seen on adjusted counts, which looked biologically ordered, is an artefact of how
SpaNorm's adjustment interacts with cell type: the two orderings correlate at Spearman **0.100**.

Pass it through `covariates` for now: a numeric `colData` column gives one global slope, and
`loglib * 1[celltype == k]` passed as k columns *is* the `CellType:loglib` interaction under
cell-means coding.

**The counts assay is not raw counts.** `batch_nichede_v10.R:224` sets
`counts(spe) <- 2^logcounts(spe) - 1`, a back-transform of SpaNorm's logPAC, and it produces
infeasible values — adjusted totals reach 287,929 against a true maximum library size of 27,553,
with 460 cells (0.59%) exceeding their own library size by more than 2x. It also destroys sparsity:
7 GB dense against 362 MB sparse for the same matrix. **The raw counts are not lost** —
`ytma_4_nichede.rds` still carries them and line 224 overwrites them;
`research/fdr-ordering/R/make_raw_spe.R` rebuilds a raw twin in 3.5 minutes. On raw counts **no**
cell exceeds its own library size. Prefer the raw object for any new fit.

### Which method to use

Measured, not assumed. Simulation numbers are from the structured-LS sweep
(`research/plasmode/summary/twostage_*.csv`), re-aggregated 2026-08-31 with the arms relabelled
`twostage/{nb,ols,spanorm}[lsstruct]`. The null table is now complete; some power tables still carry
NA cells, so those numbers may still move slightly.

**Default to `fitSpiDE(random = "intercept")` when samples are plentiful (S ≥ 16).** Raw power 0.660
at S = 30 against 0.451 (`ols`), 0.293 (`nb`), 0.274 (`spanorm`); FDP 0.041 at a nominal 0.05 by
S = 16; and the best real-data shuffle calibration above.

**It fails badly at small S**: FDP **0.35** against a nominal 0.05 at S = 4 (0.33 against a nominal
0.01), recovering to 0.04 by S = 16. Note this does *not* show up in null type-I (0.049 at S = 4,
respectable) — the two measure different things, and reading only the null table would have declared
it fine. Any sub-analysis that thins the patient count re-enters this regime.

**The 0.99.17 defaults cost null type-I on the synthetic benchmark (2026-09-07), and the benchmark
cannot show their benefit.** Arms `nested-converged`, `nested-only` and `converged-only` in the
canonical tables (40 reps, intercept mode, per-cell SE ~0.0016): over S ≥ 10 the shipped
celltype-response arm is at 0.058, converged-only 0.066, nested-only 0.068, both 0.071. Two
mechanisms: at S = 4 the nested block is harmless (0.037–0.043) and convergence carries the whole
excess (0.086–0.095); at S ≥ 16 convergence fades toward the shipped arm while the nested block holds
a constant +0.010. The pair gains recall (0.402 vs 0.361 TPR at S = 30) with slightly better FDP at
S ≥ 10 and FDP 0.59 vs 0.35 at S = 4. The simulator places niche cells by the same potential in every
sample, so it plants **no** between-sample composition effect: on it the nested intercept can only
cost, and the real-cohort shuffle grids are where it earned its place. Both facts hold; do not read
the synthetic null as a reason to revert. One candidate cure is measured and adopted (next
section); the nested block's Satterthwaite df is still unmeasured.

### The dispersion rule and the SE scale (0.99.18)

Four arms on the benchmark (40 reps; `research/notes/dispersion_arms.R`, `fdr-ordering/FINDINGS.md`
2026-09-08) crossed the convergence stage's `psi` rule (profile-ML vs `fitNB`'s moderated value
kept) with the standard-error scale (Pearson vs quasi-likelihood). Null type-I at S = 4/10/16/30:
profile + Pearson 0.091/0.076/0.071/0.068, moderated + Pearson 0.091/0.074/0.071/0.068, moderated +
QL 0.050/0.055/0.055/0.055, profile + QL 0.055/0.055/0.054/0.053. **The `psi` rule is irrelevant
and the scale is the whole difference**: the QL scale is the only configuration that holds the
nominal level at every S including 4. Recall at FDR .05 returns to the no-switch design's (0.360 vs
0.402 for Pearson at S = 30, whose excess is its inflated null; raw power 0.666/0.694/0.660) and
realised FDP sits far below nominal (0.014 at S = 30, 0 at S = 4 vs 0.59) — the cascade is
conservative on QL p-values, a separate lever. On the cohort the scale changes nothing (110 of ~115
calls shared, 81–83 of the base 84 kept; raw null `sd(t)` 1.058 vs 1.014, absorbed by the per-gene
calibration). Cost: moderated `psi` runs the cohort's real grid in 178 min vs 266 (profile), and
0.6× CPU at S ≥ 16; the QL pre-pass is free at cohort scale (3× at S = 4 only, where a task is 8
CPU-min). **Defaults from 0.99.18: `dispersion = "ql"`; the moderated `psi` was the 0.99.18 default and
reverted to `"profile"` in 0.99.19** (see "Pipeline", stage 2b, for why). The QL
machinery lives in SpaNorm (>= 1.7.10): `nbUnitDeviance()`, `nbDevianceMoments()`,
`qlDispersion()` with CPU and torch backends (the moments are one shared (log mu, log phi) table
per gene block; the rest is elementwise), oracle-tested there against `edgeR::glmQLFit()`; spiDE
only wires the pre-pass over gene blocks and `limma::squeezeVar()`. A fixed-effects unpolished fit
(`random = "none"`, no `polishSpiDE()`) keeps its legacy `psi` scale with a message.

**The two rules re-measured from one fit under the 0.99.20 polish (2026-09-12, arm
`polish-rules`, `research/fdr-ordering/FINDINGS.md`).** Equivalent: null 0.058 / 0.060
(profile / moderated) at S >= 10 and 0.055 / 0.052 at S = 4, TPR 0.379 / 0.366 with FDP
0.014 / 0.019 at S = 30; on the cohort's full transcriptome the per-gene calibrated
statistics agree at r = 0.999 with the same six exceedances at the threshold (a "6 vs 1"
call count was one extra null exceedance in five grids, not power). The moderated rule saves
13-19% of the polish on the benchmark (100.9 vs 87.3 min at 24,000 cells, where the fit is
17.7 and inference 4.2) and half on the cohort (651 vs 332 min on 64 workers), so the
default-rule choice is a cost choice; the fixture's "moderated psi 15x off" is a small-fixture
artefact that neither dataset reproduces (tau2 within 3% for the per-sample component and 7% for the nested one). **Classic BH over all triplets
does not reclaim the cascade's conservatism**: on identical fits the cascade is slightly less
conservative and slightly more powerful at every nominal level (S = 30: 0.379 / 0.014 vs
0.371 / 0.011 at 0.05; 0.476 / 0.101 vs 0.461 / 0.062 at 0.20). The realised FDP sits far
below nominal under both because the QL-scaled p-values are conservative; any recall lever
is in the p-values, not the multiplicity step. Tables `timing` (per design, with
`polish_seconds`) and `fdr_procedure` in the canonical set.

**The legacy niche-only design's higher simulation recall is a different estimand, not a better
test** (`research/notes/design_power_decomposition.R`: one power dataset under four libraries,
including a build from `5027da1` that differs from the niche-only build only by the
`CellType:condition` term). The build makes no difference (slopes correlate at 1.000); the term is
the whole gap. Without it, the mean of the planted effect — β times the mean niche potential, a flat
A-in-Responders shift — loads onto the uncentred slope: the niche-only A:B slope exceeds the other
design's by 0.25× the `CellType_A:condition` coefficient (r 0.88 on planted genes, 0.92 on null
genes), a ~10% larger estimate, yet its `t` is 1.5–1.6× because its SE is 0.67–0.72×: with no
intercept column the slope is a regression through the origin. Its recall (0.36 vs 0.17, Brown
pinned; concentrated in marker genes) is bought by assuming the estimand, it saturates at ~0.5 by
β = 2, and the term-restored design recovered 41% more triplets on the real cohort. Sidedness is not
the reason: two-sided Cauchy gives *more* discoveries than one-sided Brown in all four arms.

**Keep `pool.psi = TRUE`.** A paired ablation on identical datasets
(`research/plasmode/poolpsi_ablation.R`, 30 pairs) gives type-I 0.0673 pooled vs 0.0729 unpooled,
consistently lower at every S, paired p < 0.0001. Small but unambiguous.

### Supplying `psi` changes the objective (winsorisation), it does not just save time

`fitNB()`'s `winsor` caps counts against the **current fitted mu**, so supplying `psi` changes the
IRLS weights, hence mu, hence *which cells are capped* — the estimating and supplied-psi paths
maximise **different** objectives. On the archived two-stage estimator's near-singular stage-1
subsets the gap reaches 64
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

- `vignettes/spiDE-model.Rmd` — the model stated once, top-down, organised around Frisch–Waugh–Lovell,
  with the polish stage as its own chapter; `vignettes/spiDE-calibration.Rmd` reads lambda and is
  the **only** place in the package that quotes benchmark numbers — keep it that way.
- `research/twostage/` — the archived two-stage estimator as the package `spiDEtwostage`.
- `research/reports/benchmarks/` — the six reports in reading order (simulation study, the real
  cohort, two-stage estimator, combiner, spiGSEA calibration, what was tried and rejected), rendered to a static site at
  `research/docs/` by `build_site.R` (https://bhuvalab.github.io/spiDE-research/). They moved out of
  `vignettes/` because the eight built vignettes alone exceeded Bioconductor's 10 MB tarball cap;
  the package keeps the quickstart, model, two-stage and calibration vignettes.
- `research/reports/benchmarks/tables/*.rds` — the canonical benchmark tables the reports read
  (they render without the HPC runs; refreshed by `research/R/install_results.R` and
  `research/plasmode/install_twostage.R`). **One canonical table per scenario**: a new method arm is
  extra *rows*, not a parallel file that would carry a stale copy of the others.
- `research/` — a git submodule (`BhuvaLab/spiDE-research`) holding the benchmark harness and the
  written-up negative results, collected in the *What was tried and rejected* report.
- `research/fdr-ordering/` — all eight FDR procedures scored on the real-cohort shuffle null, on
  injected signal and end-to-end; the per-gene tail addendum above. `R/recover.R` recovers exact
  p-values from any stored `fdr = 1` table, so orderings are comparable without refitting.
- `research/fdr-triplet/` — the earlier study: coefficient-level p-values are conservative, the
  cascade is clean on exact uniforms, and control is nevertheless lost in the combination between
  them. Its five refuted hypotheses are listed so they are not re-run.
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
snapshot (`.claude/skills/run-benchmark-arm/scripts/freeze_snapshot.sh`). The same applies to **every
`Rscript` driver**, not only the package: R parses a script file incrementally through an 8 KB stdio
buffer, so a task that is inside its multi-hour fit reads the *next* expression from whatever the
live file holds by then, at the old byte offset. Two losses on 2026-09-07: 255 benchmark tasks wrote
their result and then died on parse garbage (the outputs were intact; the array was marked FAILED and
its `afterok` aggregation left at `DependencyNeverSatisfied`), and at 22:26 an in-place edit of
`research/fdr-ordering/R/package_fixed_design.R` killed the eight cohort tasks whose fit ended after
it, six to seven hours each. Every sbatch therefore now copies its driver to a task-private temp
file at task start and runs the copy, and `freeze_snapshot.sh` puts frozen copies of both drivers
under `$SPIDE_PKG/drivers/`, which the sbatch scripts prefer. Two facts that made the rescue
possible: a `git checkout` or `sed -i` writes a **new inode** and cannot touch a running reader,
while the Edit tool and a shell `>` rewrite **in place** and can; and a reader that has not yet
refilled its buffer (position 8192 in `/proc/<pid>/fdinfo`, visible through `srun --overlap
--cpu-bind=none` on the task's raw job id) is saved by truncating and rewriting the file in place
with the bytes it expects. Edit the harness between arrays, never during one.

### Checkers

`R/checkers.R` holds shared input validation (`checkSPE`, `checkCondition`, `checkCovariates`,
`checkNiche`, `checkCounts`) called at the top of every exported entry point — extend these rather than
duplicating validation logic in a new function.
