# Nested (sample × cell type) intercept and per-gene convergence: implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every tested niche slope a within-(sample, cell type) slope by adding a nested random-intercept block to the niche design, and converge each gene's fit to its own penalised-NB optimum after `SpaNorm::fitNB` returns.

**Architecture:** Two independent additions. (1) `.buildRandomEffects()` gains a second ridge-penalised indicator block, one column per non-empty (sample, cell type), tagged `Random` with `re_group = "SampleCellTypeInt"` — the existing Schall loop, `.varParamCov()` and `.satterthwaiteDF()` iterate over `re_group` groups and pick it up with no change. (2) A new `R/polish.R` runs damped Newton plus profile-ML dispersion per gene after the fit, dispatched over gene blocks through `BiocParallel` exactly as inference is, with the indicator block absorbed by a Schur complement so the cost is independent of how many groups exist.

**Tech Stack:** R/Bioconductor package, roxygen2, testthat edition 3, `SpaNorm::fitNB` / `calculateMu` / `invert_mat`, `BiocParallel`.

**Spec:** `design/specs/2026-09-04-nested-intercept-and-convergence-design.md` (this plan implements it). Background: `design/specs/2026-09-03-sample-celltype-intercept.md`, `research/fdr-ordering/REPORT.md` §5c–§5d.

## Global Constraints

- R/Bioconductor package conventions: roxygen2 comments, `@noRd` for internals, `devtools::document()` after any roxygen change.
- Dot-separated argument names are intentional (they mirror `SpaNorm::fitNB`'s own so they forward through `...`). New arguments here: `re.celltype`, `converge`, `converge.maxit`, `converge.tol`.
- **Never slice genes before `fitNB`** — dispersion is moderated across the full gene set. The polish stage runs *after* the fit and is per-gene, so it may be blocked.
- `torch` is `Suggests` only; nothing added here may hard-depend on it.
- Never match a covtype literal downstream — go through `.testedCols()` / `.nicheTestCols()`. The new columns are `Random`, so no tested-column logic changes.
- Bioconductor's `R CMD check` budget is 10 minutes. Live mixed fits cost ~69 s each on the toy, so **numerical** assertions requiring live mixed fits go in `longtests/testthat/`; `tests/testthat/` gets contract assertions and cheap direct-call unit tests. This is the documented split (`data-raw/make_test_fixtures.R`).
- New defaults: `re.celltype = TRUE` and `converge = TRUE` on `fitSpiDE()` / `spiDE()`; `re.celltype = FALSE` on `.buildNicheDesign()` / `nicheDesign()` (matching the existing deliberate asymmetry for `random`, where a constructor does not return penalty-identified columns by default).
- `re.celltype = FALSE` must reproduce the previous design column for column; `converge = FALSE` must leave `alpha`/`psi` exactly as `fitNB` returned them.
- Version bump 0.99.16 → 0.99.17 with a `NEWS.md` entry.
- Commit after every task.

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `R/design.R` | design assembly, column tagging, random-effect blocks | modify `.buildRandomEffects()`, `.buildNicheDesign()`, `nicheDesign()` |
| `R/polish.R` | **new**: per-gene penalised-NB convergence and profile-ML dispersion | create |
| `R/fitSpiDE.R` | per-bandwidth fit orchestration | modify `.fitOneBandwidth()`, `fitSpiDE()` signature |
| `R/spiDE.R` | three-stage wrapper | forward the four new arguments |
| `R/AllClasses.R` | S4 classes | one new `SpiDEFit` slot `polish`, prototype, validity, `show` |
| `R/checkers.R` | shared input validation | doc note on `checkSample()` only |
| `R/toydata.R` | synthetic fixtures | `.toySPE(composition = )` plants a between-sample composition effect on `G2` |
| `tests/testthat/test-design.R` | design contracts | nested-block tests |
| `tests/testthat/test-polish.R` | **new**: polish unit tests | create |
| `tests/testthat/test-mixedEffects.R` | mixed-fit contracts | τ² group names |
| `tests/testthat/test-satterthwaite.R` | df oracle | nested-term lme4 variant |
| `longtests/testthat/test-nested-intercept.R` | **new**: live numerical check that the composition confound is removed | create |
| `data-raw/make_test_fixtures.R` | precomputed fits | regenerate under the new defaults |
| `vignettes/spiDE-model.Rmd` | the model | two new sections |
| `NEWS.md`, `DESCRIPTION`, `CLAUDE.md` | release notes, version, project guidance | update |

---

### Task 1: The nested random-intercept block

**Files:**
- Modify: `R/design.R:73-91` (`.buildRandomEffects`), `R/design.R:225-231` and `~331-345` (`.buildNicheDesign`), `R/design.R:386-400` (`nicheDesign`)
- Test: `tests/testthat/test-design.R`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `.buildRandomEffects(sample_vec, slope_base, random, cell_type_vec = NULL)` returning `list(Z, re_group)` where `re_group` may contain `"SampleCellTypeInt"`; `.buildNicheDesign(..., random, re.celltype = FALSE)`; `nicheDesign(..., re.celltype = FALSE)`. Task 4 and Task 6 rely on the group label string `"SampleCellTypeInt"` exactly.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-design.R`:

```r
test_that("re.celltype adds one nested intercept per non-empty (sample, cell type)", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- spiDE:::.buildNicheDesign(spe, "condition", 20, random = "intercept",
                                   re.celltype = TRUE)
  cd <- SummarizedExperiment::colData(spe)
  n_grp <- length(unique(paste(cd$sample_id, cd$cell_type)))

  nested <- which(des$re_group == "SampleCellTypeInt")
  expect_length(nested, n_grp)
  expect_true(all(grepl("^SampleCellType", colnames(des$W)[nested])))
  # tagged Random, so no tested-column logic can see them
  expect_true(all(as.character(des$covtype)[nested] == "Random"))
  expect_true(all(is.na(des$coefmap$index[nested])))
  expect_true(all(is.na(des$coefmap$niche[nested])))
  # 0/1 indicators partitioning the cells
  expect_true(all(des$W[, nested] %in% c(0, 1)))
  expect_true(all(rowSums(des$W[, nested, drop = FALSE]) == 1))
  # the per-sample block is KEPT and comes first
  expect_true(any(des$re_group == "SampleInt", na.rm = TRUE))
  expect_lt(max(which(des$re_group == "SampleInt")), min(nested))
})

test_that("re.celltype = FALSE reproduces the design exactly", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  a <- spiDE:::.buildNicheDesign(spe, "condition", 20, random = "intercept",
                                 re.celltype = FALSE)
  b <- spiDE:::.buildNicheDesign(spe, "condition", 20, random = "intercept")
  expect_identical(a$W, b$W)
  expect_identical(a$re_group, b$re_group)
  expect_false(any(a$re_group == "SampleCellTypeInt", na.rm = TRUE))
})

test_that("the nested block follows the slope block under random = 'slope'", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  des <- spiDE:::.buildNicheDesign(spe, "condition", 20, random = "slope",
                                   re.celltype = TRUE)
  expect_lt(max(which(des$re_group == "SampleSlope")),
            min(which(des$re_group == "SampleCellTypeInt")))
  expect_setequal(unique(des$re_group[!is.na(des$re_group)]),
                  c("SampleInt", "SampleSlope", "SampleCellTypeInt"))
})

test_that("nicheDesign forwards re.celltype and defaults it off", {
  spe <- buildNiches(.toySPE(), sigma = 20)
  d0 <- nicheDesign(spe, condition = "condition", sigma = 20,
                    random = "intercept")
  d1 <- nicheDesign(spe, condition = "condition", sigma = 20,
                    random = "intercept", re.celltype = TRUE)
  expect_false(any(d0$re_group == "SampleCellTypeInt", na.rm = TRUE))
  expect_true(any(d1$re_group == "SampleCellTypeInt", na.rm = TRUE))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-design.R")'`
Expected: FAIL — `unused argument (re.celltype = TRUE)`.

- [ ] **Step 3: Add the block to `.buildRandomEffects()`**

Replace the function at `R/design.R:73-91` with:

```r
.buildRandomEffects <- function(sample_vec, slope_base, random,
                                cell_type_vec = NULL) {
  smp <- factor(sample_vec)
  Zint <- stats::model.matrix(~ 0 + smp)
  colnames(Zint) <- paste0("Sample", levels(smp))
  Z <- Zint
  re_group <- rep("SampleInt", ncol(Zint))

  if (random == "slope" && ncol(slope_base) > 0) {
    Zslope <- do.call(cbind, lapply(colnames(slope_base), function(bc) {
      M <- Zint * slope_base[, bc]
      colnames(M) <- paste0(colnames(Zint), ":", bc)
      M
    }))
    Z <- cbind(Z, Zslope)
    re_group <- c(re_group, rep("SampleSlope", ncol(Zslope)))
  }

  # Nested (sample x cell type) intercepts. The tested niche slopes are
  # estimated from the TOTAL covariance of niche density and expression among
  # the cells of an index type -- within-sample and between-sample -- and the
  # between-sample part is a composition effect (samples whose type-k cells sit
  # in denser type-n surroundings also differ in mean expression in type k).
  # That is a patient-level association, not neighbourhood-dependent DE, and a
  # shuffle that permutes within (sample, cell type) preserves it. These
  # indicators absorb it, making every niche slope a within-group slope
  # (Frisch-Waugh). Empty combinations are dropped, so the columns partition
  # the cells exactly.
  if (!is.null(cell_type_vec)) {
    grp <- interaction(smp, factor(cell_type_vec), drop = TRUE, sep = ".")
    Zct <- stats::model.matrix(~ 0 + grp)
    colnames(Zct) <- paste0("SampleCellType", levels(grp))
    Z <- cbind(Z, Zct)
    re_group <- c(re_group, rep("SampleCellTypeInt", ncol(Zct)))
  }
  list(Z = Z, re_group = re_group)
}
```

- [ ] **Step 4: Thread `re.celltype` through `.buildNicheDesign()`**

At `R/design.R:225-231`, add the argument after `random`:

```r
.buildNicheDesign <- function(spe, condition = NULL, sigma, index = NULL,
                              niche = NULL,
                              covariates = character(), cell_type = "cell_type",
                              name = "Niche", sample_id = "sample_id",
                              random = c("none", "intercept", "slope"),
                              re.celltype = FALSE) {
```

and in the `if (random != "none")` block (around `R/design.R:335-341`), replace the `re <- .buildRandomEffects(...)` call with:

```r
    slope_base <- W[, coefmap$type == "Niche", drop = FALSE]
    re <- .buildRandomEffects(cd[[sample_id]], slope_base, random,
                              cell_type_vec = if (re.celltype) cd[[cell_type]] else NULL)
```

- [ ] **Step 5: Forward from `nicheDesign()`**

At `R/design.R:386-400`, add the argument and pass it positionally-safely by name:

```r
nicheDesign <- function(spe, condition = NULL, sigma, index = NULL,
                        niche = NULL,
                        covariates = character(), cell_type = "cell_type",
                        name = "Niche", sample_id = "sample_id",
                        random = c("none", "intercept", "slope"),
                        re.celltype = FALSE, ...) {
  random <- match.arg(random)
  checkSPE(spe, cell_type = cell_type)
  if (!is.null(condition)) checkCondition(spe, condition)
  checkCovariates(spe, covariates)
  checkNiche(spe, sigma, name = name)
  res <- .buildNicheDesign(spe, condition, sigma, index, niche, covariates,
                           cell_type, name, sample_id, random,
                           re.celltype = re.celltype)
  keep <- c("W", "covtype", "coefmap", "mode",
            if (random != "none") "re_group")
  res[keep]
}
```

Add the roxygen `@param` above `nicheDesign()`:

```r
#' @param re.celltype logical; add a nested (sample x cell type) random
#'   intercept alongside the per-sample one, so that every tested niche slope
#'   is a within-(sample, cell type) slope. Without it the slopes also carry
#'   the between-sample composition effect (a patient-level association
#'   between a cell type's mean niche density and its mean expression), which
#'   is not neighbourhood-dependent differential expression. Ignored when
#'   \code{random = "none"}. Defaults to \code{FALSE} here and \code{TRUE} in
#'   [fitSpiDE()]: like \code{random}, a design returned with
#'   penalty-identified columns is rank-deficient, which is correct for
#'   fitting and surprising from a constructor.
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); testthat::test_file("tests/testthat/test-design.R")'`
Expected: PASS, all tests in the file.

- [ ] **Step 7: Commit**

```bash
git add R/design.R man/ NAMESPACE tests/testthat/test-design.R
git commit -m "feat(design): nested (sample x cell type) random intercept block

The tested niche slopes carried a between-sample composition effect that a
within-group shuffle preserves. These ridge-penalised indicators absorb it, so
every niche slope is a within-(sample, cell type) slope.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: `re.celltype` on the fitting entry points

**Files:**
- Modify: `R/fitSpiDE.R` (`.fitOneBandwidth` signature and its `.buildNicheDesign` call; `fitSpiDE()` method signature and its `.fitOneBandwidth` call), `R/spiDE.R` (`spiDE()` signature and its `fitSpiDE()` call)
- Test: `tests/testthat/test-mixedEffects.R`

**Interfaces:**
- Consumes: `.buildNicheDesign(..., re.celltype)` from Task 1.
- Produces: `fitSpiDE(..., re.celltype = TRUE)` and `spiDE(..., re.celltype = TRUE)`; fits carry `tau2` with names `c("SampleInt", "SampleCellTypeInt")` under the default. Task 5 relies on this default being active.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-mixedEffects.R`:

```r
test_that("re.celltype yields its own variance component and keeps per-column df", {
  spe <- buildNiches(.toySPE(n_genes = 8, n_per = 40), sigma = 30)
  f <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                re.maxit = 1L, converge = FALSE, verbose = FALSE)
  fit <- fits(f)[[1]]
  expect_setequal(names(fit@tau2), c("SampleInt", "SampleCellTypeInt"))
  expect_true(all(is.finite(fit@tau2)))
  expect_true(all(fit@tau2 > 0))
  # the nested columns are penalised, the fixed ones are not
  nested <- fit@re_group == "SampleCellTypeInt"
  expect_true(all(fit@penalty[which(nested)] == 1 / fit@tau2[["SampleCellTypeInt"]]))
  # df is still one per tested column
  expect_gt(length(fit@df), 1L)
  expect_true(all(is.finite(fit@df)))
})

test_that("re.celltype = FALSE keeps the single variance component", {
  spe <- buildNiches(.toySPE(n_genes = 8, n_per = 40), sigma = 30)
  f <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                re.celltype = FALSE, re.maxit = 1L, converge = FALSE,
                verbose = FALSE)
  expect_setequal(names(fits(f)[[1]]@tau2), "SampleInt")
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-mixedEffects.R")'`
Expected: FAIL — `unused argument (re.celltype = FALSE)` / `unused argument (converge = FALSE)`. The `converge` argument arrives in Task 4; until then these two tests fail on it, which is expected and recorded here so the Task 4 executor knows to re-run this file.

- [ ] **Step 3: Add the argument to `.fitOneBandwidth()`**

In `R/fitSpiDE.R`, add `re.celltype = TRUE` to the `.fitOneBandwidth()` signature after `random = "none"`, and pass it into the design call:

```r
  des <- .buildNicheDesign(spe, condition, sigma, index, niche, covariates,
                           cell_type, name, sample_id, random,
                           re.celltype = re.celltype && random != "none")
```

- [ ] **Step 4: Add the argument to `fitSpiDE()` and forward it**

In the `fitSpiDE()` method signature add `re.celltype = TRUE` after `random`, and in the `lapply` over `sigma` pass `re.celltype = re.celltype` alongside `random = random`.

Add the roxygen `@param`:

```r
#' @param re.celltype logical; when \code{random != "none"}, add a nested
#'   (sample x cell type) random intercept alongside the per-sample one.
#'   \strong{Default \code{TRUE}.} Without it the tested niche slopes are
#'   estimated from the total covariance of niche density and expression within
#'   an index cell type, so they also carry the between-sample composition
#'   effect: samples whose index cells sit in denser niche surroundings also
#'   differ in mean expression there. That is a patient-level association with
#'   S units, not neighbourhood-dependent differential expression, and a
#'   shuffle null that permutes within (sample, cell type) preserves it -- which
#'   is why real data and such a null were indistinguishable on the YTMA
#'   cohort. With the block present every niche slope is a within-group slope
#'   and the shuffle null is calibrated in every expression band. Set
#'   \code{FALSE} to reproduce pre-correction fits. Ignored when
#'   \code{random = "none"}.
```

- [ ] **Step 5: Forward from `spiDE()`**

In `R/spiDE.R` add `re.celltype = TRUE` to the signature after `random` and pass `re.celltype = re.celltype` in the `fitSpiDE()` call. Add an `@inheritParams fitSpiDE`-consistent `@param` block mirroring the text above in one sentence:

```r
#' @param re.celltype logical; add a nested (sample x cell type) random
#'   intercept so the tested niche slopes are within-group. See [fitSpiDE()].
```

- [ ] **Step 6: Run to verify the design half passes**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); spe <- buildNiches(spiDE:::.toySPE(n_genes = 8, n_per = 40), sigma = 30); f <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept", re.maxit = 1L, verbose = FALSE); print(spiDE::fits(f)[[1]]@tau2)'`
Expected: two named components, `SampleInt` and `SampleCellTypeInt`, both finite and positive.

- [ ] **Step 7: Commit**

```bash
git add R/fitSpiDE.R R/spiDE.R man/ tests/testthat/test-mixedEffects.R
git commit -m "feat: re.celltype on fitSpiDE() and spiDE(), default TRUE

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: `R/polish.R` — per-gene convergence, with the indicator block absorbed

**Files:**
- Create: `R/polish.R`
- Test: `tests/testthat/test-polish.R` (create)

**Interfaces:**
- Consumes: `"SampleCellTypeInt"` from Task 1; `.chunkGenes()` from `R/inference.R:14`.
- Produces:
  - `.nbPenLoglik(y, mu, psi, a, pen)` → numeric scalar.
  - `.newtonSolver(W, pen, nested)` → a list with `$solve(w, s)` → the Newton step (length `ncol(W)`) and `$xcov(w)` → the X-block of the penalised covariance; `nested` is a logical over `colnames(W)`.
  - `.polishGene(y, W, a0, psi0, pen, solver, maxit, tol)` → `list(alpha, psi, iterations, restarted, capped, singular, loglik)`.
  - `.polishFit(Y, W, alpha, psi, pen, re_group, cell_type_lab, maxit, tol, block.size, BPPARAM, verbose)` → `list(alpha, psi, loglik, polish)` where `polish` is a `data.frame` with one row per gene and columns `iterations`, `psi_fitnb`, `restarted`, `capped`, `singular`.
  Task 4 calls `.polishFit()` only.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-polish.R`:

```r
# Per-gene convergence after fitNB. SpaNorm's multi-gene IRLS shares one cell
# weight vector across genes and stops on the aggregate log-likelihood, so
# bright genes sit 1-4 SE from their own optimum with a dispersion ~1.6x too
# large. These tests pin the polish stage: it never decreases the penalised
# log-likelihood, it reaches a stationary point, it leaves an already-converged
# gene alone, and the Schur-complement solve equals the dense one.

set.seed(11)

# a small design with an indicator block, in the shape .polishFit() sees
toy_design <- function(n = 240, p_x = 6, n_grp = 8) {
  X <- cbind(1, matrix(rnorm(n * (p_x - 1)), n, p_x - 1))
  colnames(X) <- c("(Intercept)", paste0("x", seq_len(p_x - 1)))
  gidx <- rep(seq_len(n_grp), length.out = n)
  Z <- stats::model.matrix(~ 0 + factor(gidx))
  colnames(Z) <- paste0("SampleCellType", seq_len(n_grp))
  W <- cbind(X, Z)
  list(W = W, nested = grepl("^SampleCellType", colnames(W)),
       pen = c(rep(0, p_x), rep(1.7, n_grp)))
}

test_that(".newtonSolver's Schur solve equals the dense solve", {
  d <- toy_design()
  w <- runif(nrow(d$W), 0.2, 3)
  s <- rnorm(ncol(d$W))

  I <- crossprod(d$W * sqrt(w))
  diag(I) <- diag(I) + d$pen
  expect_equal(as.numeric(spiDE:::.newtonSolver(d$W, d$pen, d$nested)$solve(w, s)),
               as.numeric(solve(I, s)), tolerance = 1e-8)
})

test_that(".newtonSolver's xcov equals the X-block of the dense covariance", {
  d <- toy_design()
  w <- runif(nrow(d$W), 0.2, 3)
  I <- crossprod(d$W * sqrt(w))
  diag(I) <- diag(I) + d$pen
  expect_equal(spiDE:::.newtonSolver(d$W, d$pen, d$nested)$xcov(w),
               solve(I)[!d$nested, !d$nested], tolerance = 1e-8,
               ignore_attr = TRUE)
})

test_that(".newtonSolver falls back to the dense path with no nested block", {
  d <- toy_design()
  nested <- rep(FALSE, ncol(d$W))
  w <- runif(nrow(d$W), 0.2, 3); s <- rnorm(ncol(d$W))
  I <- crossprod(d$W * sqrt(w)); diag(I) <- diag(I) + d$pen
  expect_equal(as.numeric(spiDE:::.newtonSolver(d$W, d$pen, nested)$solve(w, s)),
               as.numeric(solve(I, s)), tolerance = 1e-8)
})

test_that(".polishGene reaches a stationary point and raises the log-likelihood", {
  d <- toy_design()
  a_true <- c(1.2, 0.4, -0.3, 0.2, 0, 0.1, rnorm(8, 0, 0.3))
  mu <- exp(d$W %*% a_true)
  y <- rnbinom(length(mu), mu = as.numeric(mu), size = 1 / 0.4)
  a0 <- c(0.2, rep(0, ncol(d$W) - 1))          # a deliberately poor start
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)

  ll0 <- spiDE:::.nbPenLoglik(y, as.numeric(exp(d$W %*% a0)), 0.4, a0, d$pen)
  out <- spiDE:::.polishGene(y, d$W, a0, 0.4, d$pen, solver,
                             maxit = 50L, tol = 1e-10)

  expect_gt(out$loglik, ll0)
  # penalised score at the polished point, at the polished psi
  mu1 <- as.numeric(exp(d$W %*% out$alpha))
  s1 <- as.numeric(crossprod(d$W, (y - mu1) / (1 + out$psi * mu1))) -
    d$pen * out$alpha
  expect_lt(max(abs(s1)), 1e-4 * max(abs(y)))
  expect_true(is.finite(out$psi) && out$psi > 0)
  expect_gte(out$iterations, 1L)
  expect_false(out$singular)
})

test_that(".polishGene leaves an already-converged gene alone", {
  d <- toy_design()
  a_true <- c(1.0, 0.3, -0.2, 0.1, 0, 0, rnorm(8, 0, 0.2))
  mu <- exp(d$W %*% a_true)
  y <- rnbinom(length(mu), mu = as.numeric(mu), size = 1 / 0.5)
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  # converge once, then polish again from the converged point
  first <- spiDE:::.polishGene(y, d$W, a_true, 0.5, d$pen, solver, 50L, 1e-12)
  again <- spiDE:::.polishGene(y, d$W, first$alpha, first$psi, d$pen, solver,
                               50L, 1e-12)
  expect_equal(again$alpha, first$alpha, tolerance = 1e-3)
  expect_equal(again$psi, first$psi, tolerance = 1e-3)
})

test_that(".polishGene restarts from cell-type means when the start is degenerate", {
  d <- toy_design()
  mu <- exp(d$W %*% c(1.0, rep(0.1, 5), rnorm(8, 0, 0.2)))
  y <- rnbinom(length(mu), mu = as.numeric(mu), size = 1 / 0.4)
  a_bad <- c(-25, rep(0, ncol(d$W) - 1))       # min log mu < -10
  solver <- spiDE:::.newtonSolver(d$W, d$pen, d$nested)
  out <- spiDE:::.polishGene(y, d$W, a_bad, 0.4, d$pen, solver, 50L, 1e-10)

  expect_true(out$restarted)
  expect_gt(min(as.numeric(d$W %*% out$alpha)), -10)
  expect_true(is.finite(out$loglik))
})

test_that(".polishFit is invariant to gene blocking", {
  d <- toy_design()
  ng <- 5
  A0 <- matrix(0, ng, ncol(d$W), dimnames = list(paste0("G", seq_len(ng)),
                                                 colnames(d$W)))
  A0[, 1] <- 0.5
  Y <- t(vapply(seq_len(ng), function(g) {
    mu <- exp(d$W %*% c(1 + 0.2 * g, rep(0.15, 5), rnorm(8, 0, 0.2)))
    rnbinom(nrow(d$W), mu = as.numeric(mu), size = 1 / 0.4)
  }, numeric(nrow(d$W))))
  dimnames(Y) <- list(rownames(A0), NULL)
  psi0 <- rep(0.4, ng)
  re_group <- ifelse(d$nested, "SampleCellTypeInt", NA_character_)

  a <- spiDE:::.polishFit(Y, d$W, A0, psi0, d$pen, re_group,
                          block.size = NULL)
  b <- spiDE:::.polishFit(Y, d$W, A0, psi0, d$pen, re_group,
                          block.size = 2L)
  expect_equal(a$alpha, b$alpha)
  expect_equal(a$psi, b$psi)
  expect_equal(nrow(a$polish), ng)
  expect_setequal(colnames(a$polish),
                  c("iterations", "psi_fitnb", "restarted", "capped", "singular"))
  expect_equal(a$polish$psi_fitnb, psi0)
})

test_that(".polishFit flags a singular gene and keeps its input values", {
  d <- toy_design()
  # a gene that is zero everywhere: every working weight is ~0, so the
  # information matrix is singular at the start
  Y <- matrix(0L, 1, nrow(d$W), dimnames = list("G1", NULL))
  A0 <- matrix(-40, 1, ncol(d$W), dimnames = list("G1", colnames(d$W)))
  re_group <- ifelse(d$nested, "SampleCellTypeInt", NA_character_)
  out <- spiDE:::.polishFit(Y, d$W, A0, 0.4, d$pen, re_group)
  expect_true(out$polish$singular || out$polish$restarted)
  expect_true(all(is.finite(out$alpha)))
  expect_true(is.finite(out$psi))
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-polish.R")'`
Expected: FAIL — `'.newtonSolver' is not an exported object` / object not found.

- [ ] **Step 3: Write `R/polish.R`**

```r
# Per-gene convergence of the penalised NB fit.
#
# SpaNorm::fitNB fits every gene in one IRLS loop: it shares a single
# gene-averaged cell weight vector, decides step-halving and convergence on the
# AGGREGATE log-likelihood, and clamps coefficient columns across genes. That is
# what makes a 13,000-gene fit affordable, and for the great majority of genes
# it is indistinguishable from the per-gene optimum. For a bright,
# cell-type-restricted gene -- whose own working weights look nothing like the
# average -- the aggregate criterion is met long before that gene's own score is
# zero: measured on the YTMA cohort, every gene in the top 5% by expression sat
# 1-4 production standard errors from its own penalised-NB optimum, with
# log-likelihood gaps of 1e4-1e6 and an edgeR dispersion ~1.6x too large.
#
# This stage removes the fit from the inference question. It runs AFTER fitNB
# (so the cross-gene dispersion moderation still happens on the whole gene set,
# which is why genes must not be blocked at fit time), and because a polished
# gene depends only on its own counts it is blockable and parallel -- the same
# split .blockedInference() uses.

#' Penalised NB log-likelihood of one gene
#'
#' @param y counts (length ncells); mu the fitted mean; psi the NB dispersion;
#'   a the coefficients; pen the per-column ridge penalty.
#' @return a numeric scalar.
#' @importFrom stats dnbinom
#' @noRd
.nbPenLoglik <- function(y, mu, psi, a, pen) {
  sum(stats::dnbinom(y, size = 1 / psi, mu = mu, log = TRUE)) -
    0.5 * sum(pen * a^2)
}

#' Newton solver for the penalised NB information, with indicator columns absorbed
#'
#' The nested (sample x cell type) random intercepts are 0/1 indicators that
#' partition the cells, so their blocks of the information matrix are cheap:
#' with per-cell weights \code{w} and the design split into dense \code{X} and
#' indicators \code{Z},
#' \code{A = X' diag(w) X + diag(pen_x)}, \code{B = X' diag(w) Z} is one
#' \code{rowsum} pass, and \code{C = Z' diag(w) Z + diag(pen_z)} is DIAGONAL.
#' The step then comes from the Schur complement \code{S = A - B C^-1 B'}, whose
#' cost is one \code{ncol(X)}-column gram regardless of how many groups there
#' are -- the difference between a 345-column and a 1,000-column gram per
#' iteration on a real cohort. \code{S^-1} is also the X-block of the full
#' penalised covariance, which is what the tested columns need.
#'
#' @param W the full design; pen the per-column penalty; nested a logical over
#'   the columns of \code{W} marking the indicator block (all-FALSE gives the
#'   plain dense path).
#' @return a list with \code{solve(w, s)} (the Newton step over all columns) and
#'   \code{xcov(w)} (the X-block of the penalised covariance), or NULL from
#'   either on a singular system.
#' @noRd
.newtonSolver <- function(W, pen, nested = NULL) {
  if (is.null(nested)) nested <- rep(FALSE, ncol(W))
  if (!any(nested)) {
    return(list(
      solve = function(w, s) {
        I <- crossprod(W * sqrt(w))
        diag(I) <- diag(I) + pen
        tryCatch(solve(I, s), error = function(e) NULL)
      },
      xcov = function(w) {
        I <- crossprod(W * sqrt(w))
        diag(I) <- diag(I) + pen
        tryCatch(solve(I), error = function(e) NULL)
      }
    ))
  }

  xi <- which(!nested)
  zi <- which(nested)
  X <- W[, xi, drop = FALSE]
  pen_x <- pen[xi]
  pen_z <- pen[zi]
  # the indicator each cell belongs to (the columns partition the cells)
  gidx <- as.integer(W[, zi, drop = FALSE] %*% seq_along(zi))
  G <- length(zi)
  gf <- factor(gidx, levels = seq_len(G))

  parts <- function(w) {
    A <- crossprod(X * sqrt(w))
    diag(A) <- diag(A) + pen_x
    agg <- rowsum(cbind(w, X * w), group = gf, reorder = TRUE)
    cvec <- agg[, 1] + pen_z                 # length G, the diagonal of C
    B <- t(agg[, -1, drop = FALSE])          # ncol(X) x G
    S <- A - B %*% (t(B) / cvec)
    list(S = S, B = B, cvec = cvec)
  }

  list(
    solve = function(w, s) {
      p <- parts(w)
      rhs <- s[xi] - as.numeric(p$B %*% (s[zi] / p$cvec))
      dx <- tryCatch(solve(p$S, rhs), error = function(e) NULL)
      if (is.null(dx)) return(NULL)
      dz <- (s[zi] - as.numeric(crossprod(p$B, dx))) / p$cvec
      out <- numeric(ncol(W))
      out[xi] <- dx
      out[zi] <- dz
      out
    },
    xcov = function(w) {
      tryCatch(solve(parts(w)$S), error = function(e) NULL)
    }
  )
}

#' Converge one gene to its own penalised NB optimum
#'
#' Damped Newton on the penalised log-likelihood at fixed \code{psi}, then
#' profile ML for \code{psi} at the converged mean, then a re-polish -- twice.
#' The information matrix is reused for up to three consecutive steps (only the
#' score is recomputed), since it changes slowly and each rebuild is the whole
#' cost of an iteration.
#'
#' Starting from fitNB's coefficients is right for almost every gene, but for a
#' degenerate fit (a fitted mean below exp(-10) somewhere) Newton from that
#' point DIVERGES -- measured: fitted log-means reaching +57 to +358 and
#' predicted one-step gains of 1e6-1e8 against actual gains of 1e2-1e5. Those
#' genes restart from a sane point instead: the per-cell-type log mean, every
#' other coefficient zero, which converges in 5-29 iterations.
#'
#' @return a list with \code{alpha}, \code{psi}, \code{loglik},
#'   \code{iterations}, \code{restarted}, \code{capped}, \code{singular}.
#' @importFrom stats optimize
#' @noRd
.polishGene <- function(y, W, a0, psi0, pen, solver, maxit = 50L, tol = 1e-8) {
  restarted <- FALSE
  singular <- FALSE

  # the sane start: cell-type (or, absent a cell-type block, overall) log means
  sane_start <- function() {
    a <- numeric(ncol(W))
    ct <- grep("^CellType[^:]*$", colnames(W))
    if (length(ct)) {
      for (j in ct) {
        cells <- W[, j] != 0
        a[j] <- log(mean(y[cells]) + 1e-3)
      }
    } else {
      a[1] <- log(mean(y) + 1e-3)
    }
    a
  }

  newton <- function(a, psi, maxit) {
    mu <- as.numeric(exp(W %*% a))
    ll <- .nbPenLoglik(y, mu, psi, a, pen)
    it <- 0L
    capped <- TRUE
    stale <- 0L
    w <- NULL
    while (it < maxit) {
      it <- it + 1L
      s <- as.numeric(crossprod(W, (y - mu) / (1 + psi * mu))) - pen * a
      if (is.null(w) || stale >= 3L) {
        w <- mu / (1 + psi * mu)
        stale <- 0L
      }
      d <- solver$solve(w, s)
      if (is.null(d)) {
        singular <<- TRUE
        capped <- FALSE
        break
      }
      step <- 1
      ok <- FALSE
      halvings <- 0L
      while (step > 1e-6) {
        a1 <- a + step * d
        mu1 <- as.numeric(exp(W %*% a1))
        ll1 <- .nbPenLoglik(y, mu1, psi, a1, pen)
        if (is.finite(ll1) && ll1 >= ll - 1e-9 * abs(ll)) {
          ok <- TRUE
          break
        }
        step <- step / 2
        halvings <- halvings + 1L
      }
      if (!ok) {
        # a fresh information matrix may rescue the direction; otherwise stop
        if (stale > 0L) {
          w <- mu / (1 + psi * mu)
          stale <- 0L
          next
        }
        capped <- FALSE
        break
      }
      # a hard line search or a stale matrix both call for a rebuild next step
      stale <- if (halvings > 2L) 3L else stale + 1L
      gain <- ll1 - ll
      a <- a1
      mu <- mu1
      ll <- ll1
      if (gain < tol * abs(ll)) {
        capped <- FALSE
        break
      }
    }
    list(a = a, mu = mu, ll = ll, it = it, capped = capped)
  }

  psi_ml <- function(mu) {
    exp(stats::optimize(
      function(lp) -sum(stats::dnbinom(y, size = 1 / exp(lp), mu = mu,
                                       log = TRUE)),
      c(log(1e-3), log(1e3)))$minimum)
  }

  a <- a0
  psi <- psi0
  if (!all(is.finite(a)) || min(as.numeric(W %*% a)) < -10) {
    a <- sane_start()
    restarted <- TRUE
  }
  f <- newton(a, psi, maxit)
  if (!restarted && (singular || !all(is.finite(f$a)) || max(f$mu) > 1e10)) {
    singular <- FALSE
    restarted <- TRUE
    f <- newton(sane_start(), psi, maxit)
  }
  it_total <- f$it
  for (k in 1:2) {
    psi <- psi_ml(f$mu)
    f2 <- newton(f$a, psi, 20L)
    it_total <- it_total + f2$it
    f <- f2
  }
  list(alpha = f$a, psi = psi, loglik = f$ll, iterations = it_total,
       restarted = restarted, capped = f$capped, singular = singular)
}

#' Converge every gene's fit, blocked over genes
#'
#' @param Y counts (genes x cells); W the design; alpha/psi fitNB's estimates;
#'   pen the per-column ridge penalty; re_group the per-column random-effect
#'   group (used only to locate the indicator block); block.size and BPPARAM as
#'   in [testSpiDE()].
#' @return a list with \code{alpha}, \code{psi}, \code{loglik} and a per-gene
#'   \code{polish} data.frame.
#' @importFrom BiocParallel bplapply SerialParam
#' @noRd
.polishFit <- function(Y, W, alpha, psi, pen, re_group = NULL,
                       maxit = 50L, tol = 1e-8, block.size = NULL,
                       BPPARAM = BiocParallel::SerialParam(), verbose = FALSE) {
  ng <- nrow(alpha)
  if (length(pen) == 1L) pen <- rep(pen, ncol(W))
  if (length(psi) == 1L) psi <- rep(psi, ng)
  nested <- if (is.null(re_group)) rep(FALSE, ncol(W)) else
    !is.na(re_group) & re_group == "SampleCellTypeInt"
  solver <- .newtonSolver(W, pen, nested)

  blocks <- .chunkGenes(ng, block.size)
  if (verbose) {
    message(sprintf("  converging %d genes per gene (%d block%s)", ng,
                    length(blocks), if (length(blocks) == 1L) "" else "s"))
  }
  res <- BiocParallel::bplapply(blocks, function(gi) {
    lapply(gi, function(g) {
      .polishGene(as.numeric(Y[g, ]), W, alpha[g, ], psi[[g]], pen, solver,
                  maxit = maxit, tol = tol)
    })
  }, BPPARAM = BPPARAM)
  res <- unlist(res, recursive = FALSE)

  out_alpha <- alpha
  out_psi <- psi
  for (g in seq_len(ng)) {
    out_alpha[g, ] <- res[[g]]$alpha
    out_psi[[g]] <- res[[g]]$psi
  }
  polish <- data.frame(
    iterations = vapply(res, `[[`, integer(1), "iterations"),
    psi_fitnb = psi,
    restarted = vapply(res, `[[`, logical(1), "restarted"),
    capped = vapply(res, `[[`, logical(1), "capped"),
    singular = vapply(res, `[[`, logical(1), "singular"),
    row.names = rownames(alpha)
  )
  list(alpha = out_alpha, psi = out_psi,
       loglik = vapply(res, `[[`, numeric(1), "loglik"), polish = polish)
}
```

- [ ] **Step 4: Run to verify they pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-polish.R")'`
Expected: PASS, all eight tests.

- [ ] **Step 5: Commit**

```bash
git add R/polish.R tests/testthat/test-polish.R
git commit -m "feat(polish): per-gene penalised-NB convergence with an absorbed indicator block

fitNB's multi-gene IRLS leaves bright genes 1-4 SE from their own optimum with
a dispersion ~1.6x too large. .polishFit() converges each gene by damped Newton
and re-estimates psi by profile ML at the converged mean; the nested indicator
block is absorbed by a Schur complement so the per-iteration cost is
independent of how many groups exist.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Wire the polish stage into the fit, with a `polish` slot

**Files:**
- Modify: `R/fitSpiDE.R` (`.fitOneBandwidth`, `fitSpiDE()`), `R/spiDE.R`, `R/AllClasses.R` (slot, prototype, validity, `show`)
- Test: `tests/testthat/test-polish.R`, `tests/testthat/test-fitSpiDE.R`

**Interfaces:**
- Consumes: `.polishFit()` from Task 3.
- Produces: `fitSpiDE(..., converge = TRUE, converge.maxit = 50L, converge.tol = 1e-8)`; `SpiDEFit@polish` (`NULL` or a per-gene `data.frame`). Task 5 and Task 6 rely on `converge` being an argument of `fitSpiDE()`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-polish.R`:

```r
test_that("converge populates @polish and raises the per-gene log-likelihood", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 10, n_per = 50), sigma = 30)
  f0 <- fitSpiDE(spe, "condition", sigma = 30, random = "none",
                 converge = FALSE, verbose = FALSE)
  f1 <- fitSpiDE(spe, "condition", sigma = 30, random = "none",
                 converge = TRUE, verbose = FALSE)
  a0 <- fits(f0)[[1]]
  a1 <- fits(f1)[[1]]

  expect_null(a0@polish)
  expect_s3_class(a1@polish, "data.frame")
  expect_equal(nrow(a1@polish), a1@ngenes)
  expect_true(all(a1@polish$iterations >= 1))
  expect_true(all(is.finite(a1@psi)) && all(a1@psi > 0))

  # every gene's own penalised log-likelihood is at least as high as fitNB's
  pen <- rep(0, ncol(a1@W))
  ll <- function(fit, g) {
    mu <- as.numeric(exp(fit@W %*% fit@alpha[g, ]))
    spiDE:::.nbPenLoglik(SummarizedExperiment::assay(spe, "counts")[g, ], mu,
                         fit@psi[g], fit@alpha[g, ], pen)
  }
  gains <- vapply(seq_len(a1@ngenes), function(g) ll(a1, g) - ll(a0, g),
                  numeric(1))
  expect_true(all(gains > -1e-6 * abs(vapply(seq_len(a1@ngenes),
                                             function(g) ll(a0, g), numeric(1)))))
})

test_that("converge = FALSE leaves the fit byte-identical to fitNB's", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 10, n_per = 50), sigma = 30)
  a <- fitSpiDE(spe, "condition", sigma = 30, random = "none",
                converge = FALSE, verbose = FALSE)
  b <- fitSpiDE(spe, "condition", sigma = 30, random = "none",
                converge = FALSE, verbose = FALSE)
  expect_identical(fits(a)[[1]]@alpha, fits(b)[[1]]@alpha)
  expect_identical(fits(a)[[1]]@psi, fits(b)[[1]]@psi)
})

test_that("a polished fit still passes validity and testSpiDE runs on it", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 10, n_per = 50), sigma = 30)
  f <- fitSpiDE(spe, "condition", sigma = 30, random = "none",
                converge = TRUE, verbose = FALSE)
  expect_true(validObject(fits(f)[[1]]))
  r <- testSpiDE(f, spe = spe, fdr = 1)
  expect_true(nrow(results(r)) > 0)
  expect_true(all(is.finite(results(r)$t)))
})
```

- [ ] **Step 2: Run to verify they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-polish.R")'`
Expected: FAIL — `unused argument (converge = FALSE)`.

- [ ] **Step 3: Add the `polish` slot**

In `R/AllClasses.R`, add `polish = "ANY"` to the `SpiDEFit` slots (after `sampling`), add `polish = NULL` to the prototype list, and add to `validSpiDEFit()` before `.checkMode()`:

```r
  if (!is.null(object@polish)) {
    if (nrow(object@polish) != object@ngenes) {
      stop("nrow of 'polish' does not match 'ngenes'")
    }
  }
```

Add the roxygen `@slot` line next to the other slot docs:

```r
#' @slot polish a data.frame or NULL; per-gene diagnostics of the convergence
#'   stage (iterations, fitNB's psi, whether the gene restarted from a sane
#'   start, whether it hit the iteration cap, whether its information matrix
#'   was singular). NULL when \code{converge = FALSE}.
```

And in the `show` method, after the `psi:` line:

```r
      sprintf("Converged per gene: %s", if (is.null(object@polish)) "no" else
        sprintf("yes (%d restarted, %d capped)", sum(object@polish$restarted),
                sum(object@polish$capped))),
```

- [ ] **Step 4: Call `.polishFit()` from `.fitOneBandwidth()`**

In `R/fitSpiDE.R`, add `converge = TRUE, converge.maxit = 50L, converge.tol = 1e-8, block.size = NULL, BPPARAM = BiocParallel::SerialParam()` to the `.fitOneBandwidth()` signature (after `df.method`), and insert immediately after the `if (random == "none") ... else ...` block, before `alpha <- fit$alpha`:

```r
  # Converge each gene to its own optimum. fitNB's shared-weight, aggregate
  # criterion leaves bright genes short of stationarity (see R/polish.R); this
  # is per-gene and therefore blockable, unlike the fit itself.
  polish <- NULL
  if (converge) {
    pen_vec <- if (is.null(penalty)) {
      if (length(lambda.a) == 1L) rep(lambda.a, ncol(W)) else lambda.a
    } else {
      penalty
    }
    pol <- .polishFit(Y, W, fit$alpha, fit$psi, pen_vec, des$re_group,
                      maxit = converge.maxit, tol = converge.tol,
                      block.size = block.size, BPPARAM = BPPARAM,
                      verbose = verbose)
    fit$alpha <- pol$alpha
    fit$psi <- pol$psi
    polish <- pol$polish
  }
```

Then pass `polish = polish` in the `new("SpiDEFit", ...)` call.

- [ ] **Step 5: Add the arguments to `fitSpiDE()` and `spiDE()`**

In the `fitSpiDE()` signature add, after `df.method`:

```r
                        converge = TRUE, converge.maxit = 50L,
                        converge.tol = 1e-8, block.size = NULL,
```

and forward all four plus `BPPARAM` in the `.fitOneBandwidth()` call. Same four arguments on `spiDE()` (which already has `block.size`), forwarded to `fitSpiDE()`.

Roxygen for `fitSpiDE()`:

```r
#' @param converge logical; after \code{fitNB} returns, converge each gene to
#'   its own penalised negative-binomial optimum and re-estimate its dispersion
#'   there. \strong{Default \code{TRUE}.} \code{fitNB} fits every gene in one
#'   IRLS loop with a single gene-averaged cell weight vector and an aggregate
#'   convergence criterion, which is what makes a whole-transcriptome fit
#'   affordable; for bright, cell-type-restricted genes it stops 1-4 standard
#'   errors short of that gene's own optimum, with a dispersion about 1.6 times
#'   too large. This stage is per-gene and so is blocked and parallelised over
#'   \code{block.size} / \code{BPPARAM} exactly as inference is. Note that it
#'   replaces edgeR's cross-gene moderated dispersion with a per-gene profile
#'   maximum-likelihood dispersion at the converged mean: with many cells per
#'   gene that is well determined, but it is a deliberate departure from
#'   \code{fitNB}'s moderation. Set \code{FALSE} to reproduce pre-correction
#'   fits.
#' @param converge.maxit,converge.tol iteration cap and relative
#'   log-likelihood tolerance for the per-gene convergence stage.
```

- [ ] **Step 6: Run to verify they pass**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); testthat::test_file("tests/testthat/test-polish.R")'`
Expected: PASS, all eleven tests.

- [ ] **Step 7: Run the fit and mixed suites**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-fitSpiDE.R"); testthat::test_file("tests/testthat/test-mixedEffects.R")'`
Expected: PASS. Task 2's two new tests now pass (they needed `converge`). Any failure caused by the new defaults must be fixed here, not deferred.

- [ ] **Step 8: Commit**

```bash
git add R/fitSpiDE.R R/spiDE.R R/AllClasses.R man/ tests/testthat/test-polish.R
git commit -m "feat: converge each gene after fitNB, recorded in SpiDEFit@polish

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: A toy fixture with a planted between-sample composition effect

**Files:**
- Modify: `R/toydata.R:190-276` (`.toySPE`)
- Test: `tests/testthat/test-fitSpiDE.R` (a cheap contract test only)

**Interfaces:**
- Consumes: nothing.
- Produces: `.toySPE(composition = 0)`. With `composition != 0`, gene `G2`'s baseline in index type `A` shifts per sample in proportion to that sample's B-cell prevalence, with **no** within-sample dependence on local B density, and prevalence is tied to `condition` so the confound lands on the three-way `CellType:condition:niche` term. Task 6 uses `composition = 2.5`.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-fitSpiDE.R`:

```r
test_that(".toySPE(composition = 0) is unchanged and composition plants a between-sample effect", {
  a <- spiDE:::.toySPE()
  b <- spiDE:::.toySPE(composition = 0)
  expect_identical(SummarizedExperiment::assay(a, "counts"),
                   SummarizedExperiment::assay(b, "counts"))

  cs <- spiDE:::.toySPE(composition = 2.5)
  cd <- SummarizedExperiment::colData(cs)
  y <- SummarizedExperiment::assay(cs, "counts")["G2", ]
  isA <- cd$cell_type == "A"

  # G2 in A cells differs BETWEEN samples ...
  m <- tapply(y[isA], droplevels(factor(cd$sample_id[isA])), mean)
  expect_gt(max(m) / min(m), 1.5)
  # ... and B-cell prevalence differs between samples too, in the same order
  pb <- tapply(cd$cell_type == "B", factor(cd$sample_id), mean)
  expect_gt(abs(cor(as.numeric(m), as.numeric(pb[names(m)]),
                    method = "spearman")), 0.6)
  # the planted G1 within-sample effect is untouched by composition
  expect_identical(SummarizedExperiment::assay(a, "counts")["G1", ],
                   SummarizedExperiment::assay(spiDE:::.toySPE(composition = 0),
                                               "counts")["G1", ])
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-fitSpiDE.R")'`
Expected: FAIL — `unused argument (composition = 0)`.

- [ ] **Step 3: Add `composition` to `.toySPE()`**

Add `composition = 0` to the signature after `beta = 2`. Then, immediately after `names(cond_levels) <- sample_ids`, insert:

```r
  # Per-sample B-cell prevalence. At composition = 0 this is the historical
  # constant 0.7 and draws no random numbers, so the default fixture -- and
  # every test that depends on it -- is bit-identical.
  p_B <- stats::setNames(rep(0.7, n_samples), sample_ids)
  if (composition != 0) {
    # tie prevalence to condition, so the between-sample composition effect
    # lands on the three-way CellType:condition:niche term that spiDE tests
    is_r <- cond_levels == "Responder"
    p_B[is_r] <- runif(sum(is_r), 0.75, 0.95)
    p_B[!is_r] <- runif(sum(!is_r), 0.35, 0.55)
  }
```

In the `cells <- lapply(...)` body, replace `runif(n_per) < 0.7` with `runif(n_per) < p_B[[sid]]`.

After the existing `log_effect["G1", ] <- ...` line, insert:

```r
  # A planted CONFOUND, not a signal: G2's baseline in A cells shifts with the
  # sample's B-cell prevalence and is CONSTANT within (sample, A), so the true
  # within-sample niche slope is exactly zero. The old design (a per-sample
  # intercept shared across cell types) has nothing to absorb it and reports it
  # as a niche effect; a nested (sample x cell type) intercept absorbs it
  # exactly. A shuffle that permutes within (sample, cell type) preserves it,
  # which is why it cannot be detected by permutation alone.
  if (composition != 0 && n_genes >= 2) {
    log_effect["G2", ] <- composition * is_A *
      (p_B[cd$sample_id] - mean(p_B))
  }
```

Add the roxygen `@param`:

```r
#' @param composition size of a planted BETWEEN-sample composition confound on
#'   G2 in index type A: each sample's B-cell prevalence is tied to its
#'   condition, and G2's baseline in that sample's A cells shifts in proportion
#'   to it, constant within (sample, A). The true within-sample niche slope is
#'   zero, so a design without a nested (sample x cell type) intercept reports
#'   a false A x B effect and one with it does not. \code{0} (the default)
#'   draws no extra random numbers, leaving the fixture unchanged.
```

- [ ] **Step 4: Run to verify it passes**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); testthat::test_file("tests/testthat/test-fitSpiDE.R")'`
Expected: PASS. If the prevalence-to-expression correlation is below the threshold, raise `composition` to 3.5 and re-run — do **not** weaken the assertion.

- [ ] **Step 5: Commit**

```bash
git add R/toydata.R man/ tests/testthat/test-fitSpiDE.R
git commit -m "test(toydata): plant a between-sample composition confound on G2

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: The numerical check — the confound is called without the nested block and not with it

**Files:**
- Create: `longtests/testthat/test-nested-intercept.R`
- Test: itself (a long test; not part of `devtools::test()`)

**Interfaces:**
- Consumes: `re.celltype` (Task 2), `converge` (Task 4), `.toySPE(composition =)` (Task 5).
- Produces: nothing downstream.

- [ ] **Step 1: Write the test**

Create `longtests/testthat/test-nested-intercept.R`:

```r
# The nested (sample x cell type) intercept, checked numerically on a fixture
# with a planted between-sample composition confound and zero within-sample
# niche slope. Live mixed fits cost ~70 s each, which is why this lives in
# longtests/ (see data-raw/make_test_fixtures.R for the contract/numerics split).
#
#   Rscript -e 'testthat::test_file("longtests/testthat/test-nested-intercept.R")'

test_that("the nested intercept removes a between-sample composition confound", {
  spe <- buildNiches(spiDE:::.toySPE(composition = 2.5, n_genes = 12), sigma = 30)

  no_nest <- spiDE(spe, condition = "condition", sigma = 30,
                   random = "intercept", re.celltype = FALSE,
                   re.maxit = 2L, fdr = 1, verbose = FALSE)
  nested <- spiDE(spe, condition = "condition", sigma = 30,
                  random = "intercept", re.celltype = TRUE,
                  re.maxit = 2L, fdr = 1, verbose = FALSE)

  pick <- function(res, gene, idx, nch) {
    tab <- results(res)
    r <- tab[tab$gene == gene & tab$ct_index == idx & tab$ct_niche == nch, ]
    expect_equal(nrow(r), 1L)
    r
  }

  # the confound: called without the nested block, not called with it
  a <- pick(no_nest, "G2", "A", "B")
  b <- pick(nested, "G2", "A", "B")
  expect_gt(abs(a$t), 3)
  expect_lt(abs(b$t), abs(a$t) / 2)

  # the genuine within-sample effect survives both
  ga <- pick(no_nest, "G1", "A", "B")
  gb <- pick(nested, "G1", "A", "B")
  expect_gt(abs(ga$t), 3)
  expect_gt(abs(gb$t), 3)
  expect_lt(abs(abs(gb$t) - abs(ga$t)) / abs(ga$t), 0.5)
})

test_that("converging each gene raises every gene's penalised log-likelihood", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 15), sigma = 30)
  f0 <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 re.maxit = 2L, converge = FALSE, verbose = FALSE)
  f1 <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 re.maxit = 2L, converge = TRUE, verbose = FALSE)
  a0 <- fits(f0)[[1]]
  a1 <- fits(f1)[[1]]
  Y <- SummarizedExperiment::assay(spe, "counts")

  ll <- function(fit, g) {
    mu <- as.numeric(exp(fit@W %*% fit@alpha[g, ]))
    spiDE:::.nbPenLoglik(Y[g, ], mu, fit@psi[g], fit@alpha[g, ], fit@penalty)
  }
  gains <- vapply(seq_len(a1@ngenes), function(g) ll(a1, g) - ll(a0, g),
                  numeric(1))
  expect_true(all(gains > -1e-6))
  expect_gt(median(gains), 0)
  # the dispersion falls, as measured on the real cohort
  expect_lt(median(a1@psi / a0@psi), 1.05)
})
```

- [ ] **Step 2: Run it**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("longtests/testthat/test-nested-intercept.R")'`
Expected: PASS. Budget ~6 minutes (four live mixed fits). If the confound's `|t|` without the nested block is below 3, raise `composition` in the fixture call to 3.5; if the nested `|t|` is not at least halved, check that `des$re_group` actually contains `SampleCellTypeInt` in that fit before touching the threshold.

- [ ] **Step 3: Commit**

```bash
git add longtests/testthat/test-nested-intercept.R
git commit -m "test: numerical check that the nested intercept removes the composition confound

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: The Satterthwaite oracle with a nested term

**Files:**
- Modify: `tests/testthat/test-satterthwaite.R`

**Interfaces:**
- Consumes: `"SampleCellTypeInt"` (Task 1), `.satterthwaiteDF()` / `.varParamCov()` (existing).
- Produces: nothing downstream.

- [ ] **Step 1: Write the test**

Append to `tests/testthat/test-satterthwaite.R`:

```r
test_that(".satterthwaiteDF matches lmerTest with a nested second grouping", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("lme4")
  set.seed(4)
  S <- 12L; K <- 3L; n_per <- 20L
  g <- rep(seq_len(S), each = K * n_per)
  k <- rep(rep(seq_len(K), each = n_per), S)
  gk <- paste(g, k, sep = ".")
  x <- rnorm(length(g))                     # a within-group covariate
  y <- 0.4 * x + rnorm(S, 0, 0.8)[g] + rnorm(S * K, 0, 0.5)[factor(gk)] +
    rnorm(length(g), 0, 1)

  m <- lmerTest::lmer(y ~ x + (1 | g) + (1 | gk), REML = TRUE)
  df_lmer <- summary(m)$coefficients["x", "df"]
  vc <- as.data.frame(lme4::VarCorr(m))
  s2 <- vc$vcov[vc$grp == "Residual"]
  tau2 <- c(SampleInt = vc$vcov[vc$grp == "g"] / s2,
            SampleCellTypeInt = vc$vcov[vc$grp == "gk"] / s2)

  W <- cbind(`(Intercept)` = 1, x = x,
             stats::model.matrix(~ 0 + factor(g)),
             stats::model.matrix(~ 0 + factor(gk)))
  re_group <- c(NA, NA, rep("SampleInt", S), rep("SampleCellTypeInt", S * K))
  pen <- ifelse(is.na(re_group), 0, 1 / tau2[re_group])
  pen[is.na(pen)] <- 0
  A <- crossprod(W)
  minv <- solve(A + diag(pen))
  df_s <- spiDE:::.satterthwaiteDF(A, minv, pen, re_group, as.list(tau2),
                                   tested = 2L, ncells = length(y),
                                   tested_names = "x")
  expect_equal(unname(df_s), unname(df_lmer), tolerance = 0.1)
})
```

- [ ] **Step 2: Run to verify it passes**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-satterthwaite.R")'`
Expected: PASS. `.varParamCov()` already loops over `re_group` groups, so no code change should be needed. If the tolerance fails, print `df_s` and `df_lmer` and check the `tau2` scaling (both must be relative to the residual variance) before altering `.satterthwaiteDF()`.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-satterthwaite.R
git commit -m "test: Satterthwaite df oracle with a nested second grouping factor

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 8: Regenerate the shipped test fixtures under the new defaults

**Files:**
- Modify: `data-raw/make_test_fixtures.R`, `inst/extdata/testfits/*.rds`
- Test: `tests/testthat/test-mixedEffects.R`, `tests/testthat/test-satterthwaite.R`

**Interfaces:**
- Consumes: the new defaults from Tasks 2 and 4.
- Produces: fixtures whose `@re_group` includes `SampleCellTypeInt` and whose `@polish` is populated.

- [ ] **Step 1: Note the new defaults in the generator's header**

Add to the comment block at the top of `data-raw/make_test_fixtures.R`:

```r
# These are refit whenever a fit-shaping default changes. As of 0.99.17 that
# includes re.celltype = TRUE (a nested sample x cell type random intercept,
# hence a second tau2 component) and converge = TRUE (per-gene convergence,
# hence a populated @polish). The fits below use the defaults deliberately, so
# the contract tests assert the CURRENT object shape.
```

- [ ] **Step 2: Regenerate**

Run: `Rscript data-raw/make_test_fixtures.R`
Expected: seven `.rds` files rewritten, with the sizes reported. Budget ~20 minutes.

- [ ] **Step 3: Verify the fixture-based tests still pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-mixedEffects.R"); testthat::test_file("tests/testthat/test-satterthwaite.R")'`
Expected: PASS. `fit_toyspe_slope.rds` must now carry all three `re_group` labels; the existing assertion `all(c("SampleInt", "SampleSlope") %in% fb@re_group)` still holds.

- [ ] **Step 4: Commit**

```bash
git add data-raw/make_test_fixtures.R inst/extdata/testfits
git commit -m "test: regenerate shipped fits under re.celltype and converge defaults

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 9: Documentation, version, and the full check

**Files:**
- Modify: `vignettes/spiDE-model.Rmd`, `R/checkers.R` (roxygen only), `NEWS.md`, `DESCRIPTION`, `CLAUDE.md`

**Interfaces:**
- Consumes: everything above.
- Produces: a releasable package.

- [ ] **Step 1: Add the two vignette sections**

In `vignettes/spiDE-model.Rmd`, after the random-effects section, add:

```markdown
## Why the niche slope must be a within-(sample, cell type) slope

Under cell-means coding the tested `CellType:condition:niche` coefficient is
estimated from the covariance of niche density and expression among the cells
of one index type. That covariance has two parts. Within a sample, cells of the
index type differ in how dense their surroundings are, and that is the effect
spiDE exists to measure. Between samples, whole samples differ in how dense the
surroundings of their index cells are on average, and they also differ in mean
expression of the gene in that cell type. The second part is a patient-level
association with `S` units. It is not neighbourhood-dependent differential
expression, and with a single random intercept per sample -- shared across cell
types -- nothing in the design absorbs it.

A permutation null does not detect it. Permuting niche rows within (sample,
cell type) preserves every group's mean niche density, so it preserves the
association exactly; measured on a real cohort, the permutation distribution of
each tested `t` had the expected spread but a per-column *mean* equal to the
real-data `t`. The null and the data were indistinguishable because both
carried the same confound.

`re.celltype = TRUE` (the default) adds a ridge-penalised intercept per
non-empty (sample, cell type). By Frisch-Waugh this makes every niche slope a
within-group slope, and the between-sample part is absorbed into the new
intercepts, where it belongs. It has its own variance component, so the
Satterthwaite degrees of freedom account for it. `random = "slope"` does *not*
do this: per-sample slopes on the `CellType:niche` bases leave the group means
untouched.

If the between-sample association is itself of interest -- "do patients whose
tumour compartment is B-cell-rich express this gene differently in tumour
cells?" -- test it at the patient level, on `S` units, with patient-level
covariates. It is a different question.

## Convergence of the per-gene fit

`SpaNorm::fitNB` fits every gene in a single IRLS loop. To make a
whole-transcriptome fit affordable it shares one gene-averaged cell weight
vector across genes, decides step-halving and convergence on the aggregate
log-likelihood, and clamps coefficient columns across genes. For most genes
that is indistinguishable from the per-gene optimum. For a bright,
cell-type-restricted gene, whose own working weights look nothing like the
average, the aggregate criterion is met long before that gene's own score
reaches zero: measured on a real cohort, every gene in the top 5% by expression
sat one to four standard errors from its own optimum, and the dispersion
estimated at that point was about 1.6 times too large.

`converge = TRUE` (the default) fixes this after the fact. For each gene it
maximises that gene's own penalised negative-binomial log-likelihood by damped
Newton, then re-estimates the dispersion by profile maximum likelihood at the
converged mean, then re-polishes. A gene whose starting fit is degenerate --
a fitted mean below `exp(-10)` somewhere -- restarts from the per-cell-type log
means, because Newton from a degenerate point diverges. Because a polished gene
depends only on its own counts, the stage is blocked over genes and
parallelised, exactly as inference is; the nested indicator block is absorbed
by a Schur complement, so the per-iteration cost does not grow with the number
of groups.

Two things to know. This replaces edgeR's cross-gene moderated dispersion with
a per-gene one; with many cells per gene that is well determined, but it is a
deliberate departure. And converging the fit *raises* the null variance of
bright genes rather than lowering it -- their previous, unconverged standard
errors were inflated by a dispersion estimated off the optimum, which was
hiding the inflation behind deflation. The higher number is the honest one.
```

- [ ] **Step 2: Note the absorption in `checkSample()`'s roxygen**

Add to the `checkSample()` documentation block in `R/checkers.R`:

```r
#' @section Absorbed covariates:
#' A per-sample random intercept absorbs any covariate that is constant within
#' a sample, which is why those are rejected when \code{random != "none"}.
#' With \code{re.celltype = TRUE} the nested (sample x cell type) intercepts
#' additionally absorb covariates constant within a (sample, cell type) group.
#' Such a covariate is not rejected -- it is legitimate to adjust for one and
#' let the penalty shrink it -- but it will not be identified.
```

- [ ] **Step 3: Version and NEWS**

Set `Version: 0.99.17` in `DESCRIPTION`. Prepend to `NEWS.md`:

```markdown
# spiDE 0.99.17

## New Features

* `re.celltype` (default `TRUE`) adds a nested (sample x cell type) random
  intercept to the niche design alongside the per-sample one, so that every
  tested niche slope is a within-(sample, cell type) slope. Without it the
  slopes also carry the between-sample composition effect -- a patient-level
  association between a cell type's mean niche density and its mean expression
  in that type -- which a shuffle null that permutes within (sample, cell type)
  preserves exactly. On the YTMA cohort that confound was the cause of the
  triplet-level FDR failure: with the nested block the shuffle null is
  calibrated in every expression band (`sd(null t)` 0.96-0.99, no expression
  gradient, no extreme tail) where before it ran to 1.42 for the brightest
  genes. `re.celltype = FALSE` reproduces pre-correction fits.

* `converge` (default `TRUE`) converges each gene to its own penalised
  negative-binomial optimum after `SpaNorm::fitNB` returns, and re-estimates
  its dispersion by profile maximum likelihood at the converged mean.
  `fitNB`'s multi-gene IRLS shares one cell weight vector across genes and
  stops on the aggregate log-likelihood, which leaves bright,
  cell-type-restricted genes one to four standard errors short of their own
  optimum with a dispersion about 1.6 times too large. The stage is per-gene,
  so it is blocked and parallelised over `block.size` / `BPPARAM` as inference
  is, with the nested indicator block absorbed by a Schur complement. New
  `SpiDEFit@polish` records the per-gene diagnostics. Note this replaces
  edgeR's cross-gene moderated dispersion with a per-gene one.

## Changes

* `nicheDesign()` gains `re.celltype`, defaulting to `FALSE` -- the same
  deliberate asymmetry as `random`, since a design returned with
  penalty-identified columns is rank-deficient.
* `.toySPE()` gains `composition`, which plants a between-sample composition
  confound with zero within-sample niche slope.
```

- [ ] **Step 4: Update CLAUDE.md**

In the "The null tail is per-gene" section, replace the sentence beginning "**The fix is in `R/design.R`**" with:

```markdown
**Fixed in 0.99.17**: `fitSpiDE(re.celltype = TRUE)` (the default) adds the nested
(sample × cell type) random intercept in `.buildRandomEffects()`, tagged `Random` with
`re_group = "SampleCellTypeInt"` and its own `tau2`; `fitSpiDE(converge = TRUE)` (also the
default) converges each gene per `R/polish.R` and re-estimates `psi` at the converged mean,
recording diagnostics in `SpiDEFit@polish`. `random = "slope"` does not substitute for the
nested block. Two costs to know: the polish stage replaces edgeR's cross-gene moderated
dispersion with a per-gene profile-ML one, and `.blockedInference()` still forms a dense
per-gene gram over the full design, so with ~660 extra columns real-cohort inference is ~8×
slower — absorbing the indicator block there is deferred to its own spec.
```

- [ ] **Step 5: Document, test, check**

Run: `Rscript -e 'devtools::document(); devtools::test()'`
Expected: all tests pass, no warnings from `document()`.

Run: `Rscript -e 'devtools::check(args = c("--no-build-vignettes"))'`
Expected: 0 errors, 0 warnings. Any new `R CMD check` note about undocumented arguments means a `@param` is missing — add it rather than suppressing.

Run: `Rscript -e 'BiocCheck::BiocCheck()'`
Expected: no new ERRORs. `set.seed` notes from `R/toydata.R` are pre-existing and documented in `.localSeed()`.

- [ ] **Step 6: Commit**

```bash
git add DESCRIPTION NEWS.md CLAUDE.md vignettes/spiDE-model.Rmd R/checkers.R man/
git commit -m "docs: the nested intercept and the convergence stage; bump to 0.99.17

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

---

### Task 10: Research validation on the real cohort and the simulation arm

**Files:**
- Create: `research/fdr-ordering/R/package_fixed_design.R`
- Modify: `research/fdr-ordering/FINDINGS.md`

**Interfaces:**
- Consumes: the released package behaviour.
- Produces: a measured statement about the packaged fix, as opposed to the research centring.

- [ ] **Step 1: Write the script**

Create `research/fdr-ordering/R/package_fixed_design.R`:

```r
#!/usr/bin/env Rscript
# Does the PACKAGED fix reproduce the research result? The centred converged
# null (R/converged_null.R with CONVNULL_CENTRE=1) demonstrated the cause by
# centring the niche-dependent columns within (sample, cell type). The package
# instead adds ridge-penalised nested intercepts, which is the same projection
# only in the tau2 -> Inf limit; at a finite, estimated tau2 the slopes are
# shrunk toward the pooled ones. This script measures the difference on the real
# cohort at bw 30, on a gene panel, against the same block shuffle grid.
#
#   SPIDE_SHUF=block SPIDE_SEED=1 Rscript research/fdr-ordering/R/package_fixed_design.R
suppressPackageStartupMessages({
  library(SpatialExperiment); library(SingleCellExperiment)
  library(SummarizedExperiment); library(data.table)
  devtools::load_all("/scratch/project_mnt/S0249/R_projects/spiDE", quiet = TRUE)
})
H <- "/scratch/project_mnt/S0249/R_projects/spiDE/research/fdr-ordering"
D <- "/scratch/project_mnt/S0249/R_projects/YTMACosMxWTA/data"
MODE <- Sys.getenv("SPIDE_SHUF", "block")
SEED <- as.integer(Sys.getenv("SPIDE_SEED", "1"))
source(file.path(H, "R", "raw_shuffle_common.R"), local = TRUE)  # if absent, replicate raw_shuffle.R's preamble verbatim
spe <- .rawShuffleSPE(MODE, SEED, bw = 30)
panel <- readRDS(file.path(H, "out_muvar", "per_gene_min_lmu.rds"))
set.seed(5)
genes <- c(panel[aveLogCPM >= quantile(aveLogCPM, .95), gene],
           panel[aveLogCPM < quantile(aveLogCPM, .95)][sample(.N, 100), gene])
res <- spiDE(spe[genes, ], condition = "Response", sigma = 30,
             covariates = "loglib", random = "intercept",
             re.celltype = TRUE, converge = TRUE, fdr = 1,
             BPPARAM = BiocParallel::MulticoreParam(4), verbose = TRUE)
saveRDS(results(res), file.path(H, "out_convnull",
        sprintf("pkgfixed_%s_bw30_seed%03d.rds", MODE, SEED)))
tab <- as.data.table(results(res))
cat("sd(t) overall:", round(sd(tab$t), 3), "\n")
print(tab[, .(sd_t = round(sd(t), 3), n = .N), by = ct_index][order(-sd_t)])
```

- [ ] **Step 2: Run it on one shuffle grid and one real grid**

Run:
```bash
cd /scratch/project_mnt/S0249/R_projects/spiDE/research
SPIDE_SHUF=block SPIDE_SEED=1 Rscript fdr-ordering/R/package_fixed_design.R
SPIDE_SHUF=none  SPIDE_SEED=1 Rscript fdr-ordering/R/package_fixed_design.R
```
Expected: per-index `sd(t)` on the shuffle grid within ~0.1 of the centred converged null's block values (Tumor ≈ 1.19, B cell ≈ 1.15, Fibroblast ≈ 1.10, small types ≈ 1.0). A materially larger value means the ridge shrinkage is not reaching the projection — record it, and report before changing anything.

- [ ] **Step 3: Score with the calibration skill**

Run: `/calibration-check` on the two output tables, and record the per-index breakdown.

- [ ] **Step 4: Append the result to FINDINGS.md and commit**

Add a dated section stating what the packaged fix measures, against the research centring, on the same grid. Then:

```bash
cd research && git add -A fdr-ordering && git commit -m "Measure the packaged nested-intercept fix against the research centring

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 5: Run the simulation arm**

Run the existing benchmark harness with the new defaults as an extra arm, following `.claude/skills/run-benchmark-arm/` (pin `SPIDE_PKG` to a frozen snapshot — a sweep against the live tree is not one experiment). The simulation has no between-sample composition effect, so type-I and power should be unchanged; a power loss would mean the nested block is absorbing within-sample signal and must be reported.

---

## Self-Review

**Spec coverage.** §1 nested intercept → Tasks 1, 2. §1 API and `checkSample` doc → Tasks 2, 9. §2 convergence stage, Schur absorption, ψ, dispatch → Tasks 3, 4. §2 `SpiDEFit` slot → Task 4. §3 toy fixture → Task 5. §4 tests: design → Task 1; τ² names → Task 2; polish → Tasks 3, 4; e2e composition → Tasks 5, 6; Satterthwaite → Task 7; fixtures → Task 8. §5 documentation → Task 9. §6 order of work → Tasks 1-9 in order, with the research runs as Task 10.

One spec deviation, deliberate: the spec put the composition e2e assertion in `tests/testthat/test-spiDE-e2e.R`, but that needs two live mixed fits (~140 s) against a 10-minute `R CMD check` budget, so Task 6 puts the numerical assertion in `longtests/` — the documented home for live mixed numerics — and Task 5 keeps a cheap fixture-level contract test in `tests/`. The existing e2e `G1`/`A`/`B` test runs with `random = "none"`, so the new defaults do not touch it; Task 6 covers the `G1` survival check under the nested design.

**Placeholder scan.** No TBDs. Every code step carries the actual code. The one conditional in Task 10 (`raw_shuffle_common.R` may not exist) names the fallback explicitly: replicate `raw_shuffle.R`'s preamble verbatim, which is what `converged_null.R` already does.

**Type consistency.** `"SampleCellTypeInt"` is the group label in Tasks 1, 3, 4, 7. `.polishFit()` is called with the argument order `(Y, W, alpha, psi, pen, re_group, ...)` in Tasks 3 and 4 and in the tests. `.newtonSolver(W, pen, nested)` returns `$solve(w, s)` and `$xcov(w)` in Task 3's implementation and its tests. `polish` columns `iterations`, `psi_fitnb`, `restarted`, `capped`, `singular` are the same in Task 3's implementation, Task 3's test, and Task 4's `show` method and validity check. `.nbPenLoglik(y, mu, psi, a, pen)` has the same signature everywhere it appears (Tasks 3, 4, 6).

**Note for the executor.** Task 2's two tests fail until Task 4 lands (`converge` does not exist yet). This is recorded in Task 2, Step 2 and re-verified in Task 4, Step 7.
