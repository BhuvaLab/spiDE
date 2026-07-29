# Per-coefficient Satterthwaite df — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the mixed-effects fit a per-coefficient Satterthwaite reference df, replacing the single global `df = S − 2`, so the ResponseNiche tests stop being over-conservative while the between-sample Response effect stays ≈ 0.99× calibrated.

**Architecture:** The Satterthwaite df is a per-tested-coefficient vector computed **once per bandwidth in the fit stage**, from the representative penalised inverse `M⁻¹ = (C'W̄C + Λ)⁻¹`. It is invariant to the per-gene dispersion `φ` (proved in the spec), so one shared length-`k` vector suffices — nothing is added to the per-gene inference hot loop. `SpiDEFit@df` generalises scalar → named length-`k` vector; both inference and FDR consumers read it through one broadcast helper. Gated by a new `df.method` argument, default `"between"` (current behaviour).

**Tech Stack:** R / Bioconductor S4, `SpaNorm::fitNB` / `invert_mat` / `calculateMu`, testthat edition 3, `lmerTest` (new Suggests, test-only oracle).

## Global Constraints

- **Never slice genes before fitting** — the whole gene set goes through one `SpaNorm::fitNB` call (dispersion moderated across genes). The df computation is post-fit and per-bandwidth; it must not re-block the fit.
- **Fit whole, infer blocked** — the per-gene inference loop in `.blockedInference()` and all GPU/tensor code must stay unchanged: `@df` is precomputed and only indexed.
- **`random = "none"` is untouched** — fixed-effects path stays bit-identical (`@df` stays `NULL`, normal reference).
- **`df.method` default is `"between"`** — the behaviour change is opt-in on this branch. `df.method` is ignored when `random = "none"`.
- **Dotted argument names** (`df.method`, `lambda.a`, `re.prop`, …) are intentional — mirror `SpaNorm::fitNB` / existing spiDE style.
- **Counts `Y` are never eagerly densified** — the representative-weight pass is gene-blocked via `.chunkGenes()`, like `.blockLoglik()`.
- Spec: [specs/2026-07-29-satterthwaite-df-design.md](2026-07-29-satterthwaite-df-design.md).

---

## File structure

- `R/inference.R` — add `.ptByCol()` (df-shape-agnostic tail helper); route the `ptail` closures in `.waldBrownGene()` and `.waldCauchyBlock()` through it. **(Task 1)**
- `R/fdr.R` — route `.nicheRecords()`'s tail computation through `.ptByCol()`, subsetting per-column df by the ResponseNiche mask. **(Task 1)**
- `R/fitSpiDE.R` — add `.repWeights()`, `.varParamCov()`, `.satterthwaiteDF()`; extend `.fitNBmixed()` / `.fitOneBandwidth()` / `fitSpiDE()` and `R/spiDE.R` with `df.method`. **(Tasks 2–4)**
- `R/AllClasses.R` — `@df` slot roxygen (scalar → per-coefficient vector). **(Task 4)**
- `tests/testthat/test-satterthwaite.R` — new: `.ptByCol` shapes, `lmerTest` oracle (intercept + slope), Response-df≈S−2, between-path identical. **(Tasks 1–4)**
- `tests/testthat/test-mixedEffects.R` — extend: `@df` shape/back-compat. **(Task 4)**
- `DESCRIPTION` — add `lmerTest`, `lme4` to Suggests. **(Task 2)**
- `inst/extdata/benchmark/null_calibration.csv` + `research`-side generator note; local calibration script. **(Task 5)**
- `vignettes/spiDE-model.Rmd`, `vignettes/spiDE-mixed-benchmark.Rmd`, `NEWS.md`. **(Task 6)**

---

## Task 1: `.ptByCol()` tail helper + route all consumers through it (behaviour-identical)

Pure refactor. `df` is still a scalar everywhere, so every result is byte-identical; this lays the seam per-column df will flow through.

**Files:**
- Modify: `R/inference.R` (add helper near top; edit `.waldBrownGene` ~L85-95, `.waldCauchyBlock` ~L342-346)
- Modify: `R/fdr.R` (`.nicheRecords` ~L84-89)
- Test: `tests/testthat/test-satterthwaite.R` (new)

**Interfaces:**
- Produces: `.ptByCol(tmat, df, lower.tail = TRUE)` — `df` may be `NULL` (→ `pnorm`), scalar (→ `pt` recycled), a length-`ncol(tmat)` vector (broadcast down columns), or a matrix matching `tmat`. `tmat` may be a matrix (genes × k) or a length-k vector (single gene).

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-satterthwaite.R
test_that(".ptByCol handles NULL / scalar / vector / matrix df", {
  tm <- matrix(c(0, 1, -1, 2), nrow = 2)         # 2 genes x 2 cols
  # NULL -> normal
  expect_equal(spiDE:::.ptByCol(tm, NULL), stats::pnorm(tm))
  # scalar -> same df everywhere
  expect_equal(spiDE:::.ptByCol(tm, 5), stats::pt(tm, df = 5))
  # per-column vector -> col 1 uses df=5, col 2 uses df=50
  ref <- cbind(stats::pt(tm[, 1], df = 5), stats::pt(tm[, 2], df = 50))
  expect_equal(spiDE:::.ptByCol(tm, c(5, 50)), ref)
  # single-gene vector input with per-column df
  expect_equal(spiDE:::.ptByCol(c(0, 2), c(5, 50)),
               c(stats::pt(0, 5), stats::pt(2, 50)))
  # lower.tail forwarded
  expect_equal(spiDE:::.ptByCol(tm, 5, lower.tail = FALSE),
               stats::pt(tm, df = 5, lower.tail = FALSE))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-satterthwaite.R")'`
Expected: FAIL — `.ptByCol` not found.

- [ ] **Step 3: Add the helper to `R/inference.R`** (after the `.chunkGenes` block, ~L20)

```r
#' Apply a t (or normal) tail to a statistic matrix, with per-column df
#'
#' \code{df} may be NULL (normal reference), a scalar (one df for all columns),
#' a length-\code{ncol(tmat)} vector (a df per tested column, broadcast down the
#' rows), or a matrix matching \code{tmat}. \code{tmat} may be a genes x k matrix
#' or a length-k vector (a single gene). Centralises the reference-distribution
#' logic so both the inference and FDR stages consume \code{@df} identically.
#' @noRd
.ptByCol <- function(tmat, df, lower.tail = TRUE) {
  if (is.null(df)) {
    return(stats::pnorm(tmat, lower.tail = lower.tail))
  }
  if (length(df) == 1L) {
    return(stats::pt(tmat, df = df, lower.tail = lower.tail))
  }
  dfx <- if (is.matrix(tmat) && !is.matrix(df)) {
    matrix(df, nrow(tmat), ncol(tmat), byrow = TRUE)
  } else {
    df               # vector tmat with length-k df, or df already a matrix
  }
  stats::pt(tmat, df = dfx, lower.tail = lower.tail)
}
```

- [ ] **Step 4: Route `.waldBrownGene` through it.** In `R/inference.R`, replace the two `ptail` definitions (the `if (is.null(W_full))` branch) so both call `.ptByCol`:

```r
  if (is.null(W_full)) {
    varcov <- SpaNorm::invert_mat(crossprod(Wsub * wt_g, Wsub))
    ptail <- function(t, lower.tail) .ptByCol(t, NULL, lower.tail)
  } else {
    info <- crossprod(W_full * wt_g, W_full) + diag(penalty)
    varcov <- SpaNorm::invert_mat(info)[sel, sel, drop = FALSE]
    ptail <- function(t, lower.tail) .ptByCol(t, df, lower.tail)
  }
```

- [ ] **Step 5: Route `.waldCauchyBlock` through it.** Replace its `ptail` (the `if (is.null(W_full))` block ~L342):

```r
  ptail <- function(t, lower.tail) {
    .ptByCol(t, if (is.null(W_full)) NULL else df, lower.tail)
  }
```

- [ ] **Step 6: Route `.nicheRecords` through it.** In `R/fdr.R` replace the `lower <- if (is.null(f@df)) ... else ...` block (~L84-88):

```r
    # per-column df aligned to the ResponseNiche subset (rn) of the tested cols;
    # a scalar/NULL @df is used as-is (recycled / normal reference).
    dfn <- if (is.null(f@df) || length(f@df) == 1L) f@df else f@df[rn]
    lower <- .ptByCol(tmat, dfn)
    pmat <- pmin(lower, 1 - lower)
```

- [ ] **Step 7: Run new test + full suite to verify pass and no regressions**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-satterthwaite.R"); devtools::test()'`
Expected: PASS; all existing tests green (byte-identical behaviour, df still scalar).

- [ ] **Step 8: Commit**

```bash
git add R/inference.R R/fdr.R tests/testthat/test-satterthwaite.R
git commit -m "Add .ptByCol tail helper; route inference + FDR through it (identical behaviour)"
```

---

## Task 2: `.satterthwaiteDF()` for the random-intercept case, validated against `lmerTest`

The statistical core. Pure functions on `(A, minv, pen, re_group, tau2, tested, ncells)`; validated by an `lmerTest` oracle under the standardised working-model mapping (spec §"the one genuinely new quantity"): `wbar = 1/σ̂²`, `pen = 1/τ̂²`, `φ = 1` ⇒ `minv_ββ` equals `lmer`'s `Cov(β̂)` and the df matches.

**Files:**
- Modify: `R/fitSpiDE.R` (add `.varParamCov`, `.satterthwaiteDF` after `.fitNBmixed`)
- Modify: `DESCRIPTION` (Suggests: `lmerTest`, `lme4`)
- Test: `tests/testthat/test-satterthwaite.R`

**Interfaces:**
- Produces:
  - `.varParamCov(A, minv, pen, re_group, ncells)` → list(`cov` = (M+1)×(M+1) covariance of `θ̂ = (φ, τ²_groups)` at φ=1, `groups` = character, `gcols` = list of full-column indices per group).
  - `.satterthwaiteDF(A, minv, pen, re_group, tau2, tested, ncells, tested_names)` → named numeric length `length(tested)`; entry `j` is `2·v_jj² / (d_j' cov d_j)`, `d_j = (v_jj, g_j·)`, `g_jm = (pen_m/τ²_m)·Σ_{k∈m}minv[j,k]²`.
- Consumes: `A = crossprod(W·√w̄)` (p×p), `minv = (A + diag(pen))⁻¹`, `pen` length p (0 fixed / `1/τ²_m` random), `re_group` length p (`NA` fixed), `tau2` named vector, `tested` = full-column indices of the Response/ResponseNiche columns.

- [ ] **Step 1: Write the failing oracle test** (guarded by `lmerTest`)

```r
test_that(".satterthwaiteDF matches lmerTest df on a Gaussian random-intercept LMM", {
  skip_if_not_installed("lmerTest")
  skip_if_not_installed("lme4")
  set.seed(1)
  n_g <- 24L; n_per <- 8L; n <- n_g * n_per
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- rnorm(n)                                  # varies WITHIN group
  u <- rnorm(n_g, 0, 0.8)
  y <- 1 + 0.5 * x + u[as.integer(g)] + rnorm(n, 0, 1.0)
  m  <- lmerTest::lmer(y ~ x + (1 | g), REML = TRUE)
  df_lmer <- summary(m)$coefficients["x", "df"]
  vc   <- as.data.frame(lme4::VarCorr(m))
  tau2 <- vc$vcov[vc$grp == "g"]
  sig2 <- vc$vcov[vc$grp == "Residual"]
  # spiDE standardised working-model inputs: weights = precision, phi = 1
  Z <- stats::model.matrix(~ 0 + g)
  W <- cbind(`(Intercept)` = 1, x = x, Z)
  reg <- c(NA, NA, rep("g", ncol(Z)))
  wbar <- rep(1 / sig2, n)
  pen  <- c(0, 0, rep(1 / tau2, ncol(Z)))
  A <- crossprod(W * sqrt(wbar))
  minv <- SpaNorm::invert_mat(A + diag(pen))
  df_s <- spiDE:::.satterthwaiteDF(A, minv, pen, reg,
                                   stats::setNames(tau2, "g"),
                                   tested = 2L, ncells = n, tested_names = "x")
  expect_equal(unname(df_s), unname(df_lmer), tolerance = 0.05)
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-satterthwaite.R", filter = "lmerTest")'`
Expected: FAIL — `.satterthwaiteDF` not found.

- [ ] **Step 3: Add `lmerTest`, `lme4` to `DESCRIPTION` Suggests**

Add both to the `Suggests:` field (comma-separated, alongside the existing entries). Then:
Run: `Rscript -e 'if(!requireNamespace("lmerTest",quietly=TRUE)) install.packages("lmerTest")'`

- [ ] **Step 4: Implement `.varParamCov` + `.satterthwaiteDF`** in `R/fitSpiDE.R` (after `.fitNBmixed`)

```r
#' Covariance of the working-model variance parameters (phi, tau2_groups)
#'
#' Reduced-form REML expected (Fisher) information at phi = 1 for the working
#' linear mixed model V = W-bar^{-1} + Z G Z', built entirely from p x p M^{-1}
#' sub-blocks (never an n x n matrix). Parameter order is (phi, groups...). The
#' phi-phi entry reduces to 0.5*(ncells - p + tr((M^{-1} Lambda)^2)); the
#' group entries use B_{a,b} = A[ca,cb] - A[ca,] M^{-1} A[,cb] = Z_a' P Z_b.
#' @noRd
.varParamCov <- function(A, minv, pen, re_group, ncells) {
  p <- ncol(A)
  groups <- unique(re_group[!is.na(re_group)])
  Mn <- length(groups)
  gcols <- lapply(groups, function(gr) which(re_group == gr))
  I <- matrix(0, Mn + 1L, Mn + 1L)
  ML <- minv * rep(pen, each = p)          # M^{-1} Lambda (scale columns by pen)
  I[1, 1] <- 0.5 * (ncells - p + sum(ML * t(ML)))   # tr((M^{-1}Lambda)^2)
  minvLminv <- minv %*% ML                 # M^{-1} Lambda M^{-1}
  for (a in seq_len(Mn)) {
    ca <- gcols[[a]]
    Aa <- A[, ca, drop = FALSE]            # U_a = A[, group a]
    MinvAa <- minv %*% Aa
    for (b in seq_len(a)) {
      cb <- gcols[[b]]
      Bab <- A[ca, cb, drop = FALSE] - crossprod(Aa, minv %*% A[, cb, drop = FALSE])
      I[a + 1L, b + 1L] <- I[b + 1L, a + 1L] <- 0.5 * sum(Bab * Bab)  # 0.5||B||_F^2
    }
    tr_Baa <- sum(diag(A)[ca]) - sum(MinvAa * Aa)
    tr_UMU <- sum((minvLminv %*% Aa) * Aa)
    I[1, a + 1L] <- I[a + 1L, 1] <- 0.5 * (tr_Baa - tr_UMU)
  }
  cov <- tryCatch(solve(I),
                  error = function(e) SpaNorm::invert_mat(I))
  list(cov = cov, groups = groups, gcols = gcols)
}

#' Per-coefficient Satterthwaite degrees of freedom (shared across genes)
#'
#' df_j = 2 v_jj^2 / (d_j' Cov(theta-hat) d_j), with v_jj = [M^{-1}]_jj, gradient
#' d_j = (v_jj, g_j1, ..., g_jM) and g_jm = (pen_m / tau2_m) * sum_{k in m} M^{-1}_{jk}^2.
#' Invariant to the per-gene dispersion phi (spec), hence one shared vector.
#' @noRd
.satterthwaiteDF <- function(A, minv, pen, re_group, tau2, tested, ncells,
                             tested_names) {
  vp <- .varParamCov(A, minv, pen, re_group, ncells)
  Vt  <- minv[tested, , drop = FALSE]      # k x p
  vjj <- diag(minv)[tested]                # k
  grad <- matrix(0, length(tested), length(vp$groups) + 1L)
  grad[, 1] <- vjj
  for (a in seq_along(vp$groups)) {
    gc <- vp$gcols[[a]]
    pen_m <- pen[gc[1]]                     # constant within a group
    Sm <- rowSums(Vt[, gc, drop = FALSE]^2)
    grad[, a + 1L] <- (pen_m / tau2[[vp$groups[a]]]) * Sm
  }
  varse2 <- rowSums((grad %*% vp$cov) * grad)
  df <- 2 * vjj^2 / varse2
  df <- pmin(pmax(df, 1), ncells)
  stats::setNames(df, tested_names)
}
```

- [ ] **Step 5: Run oracle test to verify pass**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-satterthwaite.R", filter = "lmerTest")'`
Expected: PASS (df within 5% of `lmerTest`).

- [ ] **Step 6: Add a within-vs-between limiting-behaviour test** (no lmerTest needed)

```r
test_that(".satterthwaiteDF gives large df for a within-group contrast, small for between", {
  set.seed(2)
  n_g <- 20L; n_per <- 10L; n <- n_g * n_per
  g <- factor(rep(seq_len(n_g), each = n_per))
  xin  <- rnorm(n)                                   # within-group
  xbtw <- rnorm(n_g)[as.integer(g)]                  # constant within group
  Z <- stats::model.matrix(~ 0 + g)
  W <- cbind(`(Intercept)` = 1, xin = xin, xbtw = xbtw, Z)
  reg <- c(NA, NA, NA, rep("g", ncol(Z)))
  tau2 <- 0.5; pen <- c(0, 0, 0, rep(1 / tau2, ncol(Z)))
  A <- crossprod(W); minv <- SpaNorm::invert_mat(A + diag(pen))
  df <- spiDE:::.satterthwaiteDF(A, minv, pen, reg, c(g = tau2),
                                 tested = c(2L, 3L), ncells = n,
                                 tested_names = c("xin", "xbtw"))
  expect_gt(df[["xin"]], 5 * df[["xbtw"]])           # within >> between
  expect_lt(df[["xbtw"]], n_g)                        # between ~ O(n_g), not O(n)
})
```

Run the file; expect PASS.

- [ ] **Step 7: Commit**

```bash
git add R/fitSpiDE.R DESCRIPTION tests/testthat/test-satterthwaite.R
git commit -m "Add Satterthwaite df core (.varParamCov/.satterthwaiteDF); lmerTest oracle"
```

---

## Task 3: Wire `df.method` into the fit; produce the `@df` vector for `random = "intercept"`

**Files:**
- Modify: `R/fitSpiDE.R` (`.repWeights`; `.fitNBmixed` signature + final-fit block; `.fitOneBandwidth`; `fitSpiDE` method)
- Modify: `R/spiDE.R` (`spiDE` method signature + forward)
- Test: `tests/testthat/test-satterthwaite.R`

**Interfaces:**
- Produces: `fitSpiDE(..., df.method = c("between","satterthwaite"))`; `@df` is a scalar under `"between"`, a named length-`k` vector (aligned to `@t_stat`/`@se` columns) under `"satterthwaite"`.
- Consumes: `.satterthwaiteDF` (Task 2); `.repWeights(Y, alpha, W, psi)` → length-`ncells` gene-averaged working weights.

- [ ] **Step 1: Write the failing wiring test**

```r
test_that("df.method='satterthwaite' yields a per-column @df; Response df ~ S-2", {
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 16, sd_patient = 0.7),
                     sigma = 30)
  fb <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 df.method = "between", verbose = FALSE)
  fs <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 df.method = "satterthwaite", verbose = FALSE)
  db <- fits(fb)[[1]]@df
  ds <- fits(fs)[[1]]@df
  # between: scalar S-2 = 14; satterthwaite: one df per tested column
  expect_length(db, 1L)
  ct <- as.character(fits(fs)[[1]]@covtype)
  n_tested <- sum(grepl("Response", ct))
  expect_length(ds, n_tested)
  # Response column df ~ S-2 (the 0.99x regression anchor)
  cm <- fits(fs)[[1]]@coefmap
  resp_name <- cm$covariate[ct == "Response"]
  expect_equal(unname(ds[resp_name]), 14, tolerance = 0.20)   # 16 samples -> ~14
  # niche-interaction dfs are larger than the Response df (within-sample info)
  rn_names <- cm$covariate[ct == "ResponseNiche"]
  expect_gt(stats::median(ds[rn_names]), ds[resp_name])
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-satterthwaite.R", filter = "per-column")'`
Expected: FAIL — `unused argument (df.method = ...)`.

- [ ] **Step 3: Add `.repWeights` to `R/fitSpiDE.R`** (before `.fitNBmixed`)

```r
#' Gene-averaged working weights over all cells, computed gene-block-wise
#'
#' The representative w-bar_c = mean_g mu_{c,g}/(1 + psi_g mu_{c,g}) used to form
#' the representative penalised inverse for the Satterthwaite df. Blocked over
#' genes so the genes x cells mean matrix is never fully materialised.
#' @noRd
.repWeights <- function(Y, alpha, W, psi, block.size = 2000L) {
  ng <- nrow(alpha)
  acc <- numeric(ncol(Y))
  for (gi in .chunkGenes(ng, block.size)) {
    mub <- SpaNorm::calculateMu(rep(0, length(gi)), alpha[gi, , drop = FALSE], W)
    acc <- acc + colSums(mub / (1 + psi[gi] * mub))
  }
  acc / ng
}
```

- [ ] **Step 4: Extend `.fitNBmixed` to compute the df vector.** Add `df.method = "between"` and `cols_tested = NULL` to its signature. Replace the between-patient df block (the `n_samples <- ...; df <- max(n_samples - 2, 1)` lines) with:

```r
  n_samples <- sum(re_group == "SampleInt", na.rm = TRUE)
  df_between <- max(n_samples - 2, 1)
  if (df.method == "satterthwaite" && !is.null(cols_tested)) {
    wbar <- .repWeights(Y, fit$alpha, W, fit$psi)
    A <- crossprod(W * sqrt(wbar))
    minv <- SpaNorm::invert_mat(A + diag(pen))
    tested <- which(cols_tested)
    df <- .satterthwaiteDF(A, minv, pen, re_group, tau2, tested, ncol(Y),
                           colnames(W)[tested])
  } else {
    df <- df_between
  }
  c(fit, list(penalty = pen, tau2 = tau2, df = df))
```

- [ ] **Step 5: Pass `cols_tested` / `df.method` from `.fitOneBandwidth`.** Add `df.method = "between"` to its signature; in the `random != "none"` branch, pass `cols_tested = grepl("Response", as.character(des$covtype))` and `df.method = df.method` into `.fitNBmixed(...)`.

```r
    fit <- .fitNBmixed(Y, W, des$re_group, lambda.a, winsor, backend, verbose,
                       re.maxit = re.maxit, re.tol = re.tol,
                       tau2.init = tau2.init, idx = idx,
                       re.maxit.psi = re.maxit.psi,
                       df.method = df.method,
                       cols_tested = grepl("Response", as.character(des$covtype)),
                       ...)
```

- [ ] **Step 6: Add `df.method` to `fitSpiDE` and forward it.** In the `fitSpiDE` method signature add `df.method = c("between", "satterthwaite")`; after `random <- match.arg(random)` add `df.method <- match.arg(df.method)`; pass `df.method = df.method` into the `.fitOneBandwidth(...)` call.

- [ ] **Step 7: Add `df.method` to `spiDE` and forward it.** In `R/spiDE.R` add `df.method = c("between", "satterthwaite")` to the signature, `df.method <- match.arg(df.method)` after the other `match.arg`s, and `df.method = df.method` in the `fitSpiDE(...)` call.

- [ ] **Step 8: Run wiring test + between-identical regression**

```r
test_that("df.method='between' is byte-identical to the pre-existing scalar df", {
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 12, sd_patient = 0.7),
                     sigma = 30)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  fb <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                 df.method = "between", verbose = FALSE)
  expect_equal(fits(fb)[[1]]@df, 10)                        # 12 samples -> S-2
  tb <- spiDE:::.blockedInference(fits(fb)[[1]], Y)
  expect_false(any(is.na(tb@p.combined.pos)))
})
```

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-satterthwaite.R")'`
Expected: PASS.

- [ ] **Step 9: Regenerate docs and run the full suite**

Run: `Rscript -e 'devtools::document(); devtools::test()'`
Expected: NAMESPACE/man updated for the new arg; all tests green.

- [ ] **Step 10: Commit**

```bash
git add R/fitSpiDE.R R/spiDE.R man NAMESPACE tests/testthat/test-satterthwaite.R
git commit -m "Wire df.method='satterthwaite' into the intercept fit; per-column @df vector"
```

---

## Task 4: Extend to `random = "slope"`; `@df` slot docs + back-compat tests

`.satterthwaiteDF` / `.varParamCov` already handle multiple groups. This task validates the two-component (`SampleInt` + `SampleSlope`) case against `lmerTest`, hardens the `solve(I)` fallback, and finalises slot documentation.

**Files:**
- Modify: `R/AllClasses.R` (`@df` roxygen)
- Test: `tests/testthat/test-satterthwaite.R`, `tests/testthat/test-mixedEffects.R`

**Interfaces:**
- Consumes: `.satterthwaiteDF` with `re_group` containing two groups.

- [ ] **Step 1: Write the slope oracle test** (random intercept + random slope)

```r
test_that(".satterthwaiteDF matches lmerTest with a random slope", {
  skip_if_not_installed("lmerTest"); skip_if_not_installed("lme4")
  set.seed(3)
  n_g <- 30L; n_per <- 12L; n <- n_g * n_per
  g <- factor(rep(seq_len(n_g), each = n_per))
  x <- rnorm(n)
  u0 <- rnorm(n_g, 0, 0.7); u1 <- rnorm(n_g, 0, 0.5)
  y <- 1 + 0.4 * x + u0[as.integer(g)] + u1[as.integer(g)] * x + rnorm(n, 0, 1)
  m <- lmerTest::lmer(y ~ x + (1 + x || g), REML = TRUE)   # independent int+slope
  df_lmer <- summary(m)$coefficients["x", "df"]
  vc <- as.data.frame(lme4::VarCorr(m))
  t_int <- vc$vcov[vc$grp == "g" & is.na(vc$var2) & vc$var1 == "(Intercept)"]
  t_slp <- vc$vcov[grepl("g", vc$grp) & vc$var1 == "x"]
  sig2  <- vc$vcov[vc$grp == "Residual"]
  Zi <- stats::model.matrix(~ 0 + g); Zs <- Zi * x
  W  <- cbind(`(Intercept)` = 1, x = x, Zi, Zs)
  reg <- c(NA, NA, rep("SampleInt", ncol(Zi)), rep("SampleSlope", ncol(Zs)))
  wbar <- rep(1 / sig2, n)
  pen  <- c(0, 0, rep(1 / t_int, ncol(Zi)), rep(1 / t_slp, ncol(Zs)))
  A <- crossprod(W * sqrt(wbar)); minv <- SpaNorm::invert_mat(A + diag(pen))
  df_s <- spiDE:::.satterthwaiteDF(A, minv, pen, reg,
            c(SampleInt = t_int, SampleSlope = t_slp),
            tested = 2L, ncells = n, tested_names = "x")
  expect_equal(unname(df_s), unname(df_lmer), tolerance = 0.08)
})
```

Run it; if `.varParamCov`'s `solve(I)` is singular for the slope block, the `tryCatch`→`invert_mat` fallback (already in Task 2's code) engages. Expected: PASS.

- [ ] **Step 2: End-to-end slope smoke test**

```r
test_that("df.method='satterthwaite' works for random='slope'", {
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = 16, sd_patient = 0.7),
                     sigma = 30)
  fs <- fitSpiDE(spe, "condition", sigma = 30, random = "slope",
                 df.method = "satterthwaite", verbose = FALSE)
  ds <- fits(fs)[[1]]@df
  ct <- as.character(fits(fs)[[1]]@covtype)
  expect_length(ds, sum(grepl("Response", ct)))
  expect_true(all(is.finite(ds)) && all(ds >= 1))
})
```

- [ ] **Step 3: Update `@df` slot roxygen** in `R/AllClasses.R` (replace the `@slot df` line):

```r
#' @slot df a numeric (or NULL), the Wald reference degrees of freedom. NULL for
#'   a fixed-effects fit (normal reference); a scalar between-patient `S - 2`
#'   under `df.method = "between"`; or a named per-tested-coefficient vector
#'   (aligned to the columns of `t_stat`/`se`) under `df.method = "satterthwaite"`
#'   — the Response main effect stays ~ `S - 2` while niche interactions, informed
#'   by within-sample variation, get a larger df.
```

- [ ] **Step 4: Extend `test-mixedEffects.R` back-compat block.** In the first test ("random='none' is unchanged"), after the existing slope assertions add:

```r
  # a satterthwaite fit populates @df as a per-tested-column vector
  bs <- fitSpiDE(spe, "condition", sigma = 20, random = "slope",
                 df.method = "satterthwaite", verbose = FALSE)
  fbs <- fits(bs)[[1]]
  expect_gt(length(fbs@df), 1L)
  expect_equal(length(fbs@df), sum(grepl("Response", as.character(fbs@covtype))))
```

- [ ] **Step 5: Run full suite + document**

Run: `Rscript -e 'devtools::document(); devtools::test()'`
Expected: all green.

- [ ] **Step 6: Commit**

```bash
git add R/AllClasses.R man tests/testthat/test-satterthwaite.R tests/testthat/test-mixedEffects.R
git commit -m "Satterthwaite df for random='slope'; @df slot docs + back-compat tests"
```

---

## Task 5: Calibration assessment + shipped `null_calibration.csv`

Produces the branch's evidence and the CSV the benchmark vignette reads. Small in-repo grid (`.toyClustered` at a few `S`), plus an `expect`-style assertion test.

**Files:**
- Create: `inst/scripts/make_null_calibration.R` (generator; run manually, not at build/check)
- Create: `inst/extdata/benchmark/null_calibration.csv` (committed artifact)
- Test: `tests/testthat/test-satterthwaite.R`

**Interfaces:**
- Produces: `null_calibration.csv` with columns `S, method, df.method, effect, type1` (`effect ∈ {"response","responseniche"}`).

- [ ] **Step 1: Write the calibration assertion test** (fast: two S values, few genes)

```r
test_that("satterthwaite relieves ResponseNiche over-conservatism without overshoot", {
  set.seed(11)
  type1 <- function(S, df.method) {
    spe <- buildNiches(spiDE:::.toyClustered(n_samples = S, n_genes = 60,
                                             sd_patient = 0.7, seed = S), sigma = 30)
    Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
    f <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
                  df.method = df.method, verbose = FALSE)
    ff <- spiDE:::.blockedInference(fits(f)[[1]], Y)
    ct <- as.character(ff@covtype)
    resp <- ff@coefmap$covariate[ct == "Response"]
    rn   <- ff@coefmap$covariate[ct == "ResponseNiche"]
    dfv  <- ff@df
    p_of <- function(cols) {
      tt <- ff@t_stat[, cols, drop = FALSE]
      dd <- if (length(dfv) == 1L) dfv else dfv[cols]
      2 * spiDE:::.ptByCol(-abs(tt), dd)
    }
    c(response = mean(p_of(resp) < 0.05),
      responseniche = mean(p_of(rn) < 0.05))
  }
  S <- 20L
  tb <- type1(S, "between"); ts <- type1(S, "satterthwaite")
  # Response effect stays ~nominal under both (the 0.99x anchor, loose here)
  expect_lt(tb[["response"]], 0.15)
  expect_lt(ts[["response"]], 0.15)
  # niche tests: satterthwaite is >= between (less conservative) and not blown up
  expect_gte(ts[["responseniche"]] + 1e-9, tb[["responseniche"]])
  expect_lt(ts[["responseniche"]], 0.15)          # no overshoot into anti-conservatism
})
```

Run it; expect PASS. (If `responseniche` type-I under satterthwaite exceeds ~0.05–0.10 materially, that is the **KR decision gate** tripping — record it and consult the spec before proceeding.)

- [ ] **Step 2: Write the generator** `inst/scripts/make_null_calibration.R`

```r
# Regenerate inst/extdata/benchmark/null_calibration.csv (run manually).
# Local small-grid null calibration: type-I error for the Response and
# ResponseNiche effects, fixed/intercept x between/satterthwaite, across S.
devtools::load_all()
suppressPackageStartupMessages({ library(SummarizedExperiment) })
set.seed(1)
S_grid <- c(6, 10, 16, 24)
reps   <- 10L
rows <- list()
for (S in S_grid) for (r in seq_len(reps)) {
  spe <- buildNiches(spiDE:::.toyClustered(n_samples = S, n_genes = 200,
                                           sd_patient = 0.7, seed = S * 100 + r),
                     sigma = 30)
  Y <- as.matrix(assay(spe, "counts"))
  for (m in c("none", "intercept")) for (dm in c("between", "satterthwaite")) {
    if (m == "none" && dm == "satterthwaite") next   # df.method irrelevant
    f  <- fitSpiDE(spe, "condition", sigma = 30, random = m,
                   df.method = dm, verbose = FALSE)
    ff <- spiDE:::.blockedInference(fits(f)[[1]], Y)
    ct <- as.character(ff@covtype); dfv <- ff@df
    p_of <- function(sel) {
      tt <- ff@t_stat[, ff@coefmap$covariate[ct == sel], drop = FALSE]
      dd <- if (is.null(dfv) || length(dfv) == 1L) dfv else
        dfv[ff@coefmap$covariate[ct == sel]]
      2 * spiDE:::.ptByCol(-abs(tt), dd)
    }
    for (eff in c("Response", "ResponseNiche")) {
      rows[[length(rows) + 1L]] <- data.frame(
        S = S, rep = r,
        method = c(none = "fixed", intercept = "intercept")[[m]],
        df.method = dm,
        effect = tolower(eff),
        type1 = mean(p_of(eff) < 0.05), stringsAsFactors = FALSE)
    }
  }
}
out <- do.call(rbind, rows)
agg <- aggregate(type1 ~ S + method + df.method + effect, out, mean)
dir.create("inst/extdata/benchmark", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(agg, "inst/extdata/benchmark/null_calibration.csv",
                 row.names = FALSE)
```

- [ ] **Step 3: Run the generator to produce the CSV**

Run: `Rscript inst/scripts/make_null_calibration.R`
Expected: `inst/extdata/benchmark/null_calibration.csv` written; inspect that `responseniche` type-I under `satterthwaite` is ≥ `between` and ≤ ~0.05–0.08 across `S`.

- [ ] **Step 4: Run the assertion test**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-satterthwaite.R", filter = "over-conservatism")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add inst/scripts/make_null_calibration.R inst/extdata/benchmark/null_calibration.csv tests/testthat/test-satterthwaite.R
git commit -m "Null-calibration assessment: assertion test + shipped null_calibration.csv"
```

---

## Task 6: Vignette reports + NEWS

**Files:**
- Modify: `vignettes/spiDE-model.Rmd`, `vignettes/spiDE-mixed-benchmark.Rmd`, `NEWS.md`

- [ ] **Step 1: `spiDE-model.Rmd` — rewrite the reference-distribution point.** In §"Corrected inference", replace point 3 so it says: `S − 2` is the correct df **for the between-sample Response contrast**, but the ResponseNiche contrasts carry within-sample niche information and warrant a **larger, per-coefficient** df. Under `df.method = "satterthwaite"` each tested coefficient gets `2·v_jj²/(d_j'Cov(θ̂)d_j)`; the Response column returns ≈ `S − 2` (so the between-patient df is the *special case for the main effect*, not a global constant), while niche terms move toward the residual df.

- [ ] **Step 2: `spiDE-model.Rmd` — add subsection "Per-coefficient degrees of freedom (Satterthwaite)".** Include: the df formula (copy from the spec's "Statistical core"); the **φ-invariance** argument (why one shared per-coefficient vector suffices, same logic as the shared `τ²`); the `Â = Cov(τ̂²)` REML-information quantity; and a runnable chunk fitting the toy clustered data both ways:

````markdown
```{r satt, message = FALSE}
fs <- fitSpiDE(spe, "condition", sigma = 30, random = "intercept",
               df.method = "satterthwaite", verbose = FALSE)
ff <- fits(fs)[[1]]; ct <- as.character(ff@covtype)
data.frame(effect = c("Response", "ResponseNiche (median)"),
           df = c(ff@df[ff@coefmap$covariate[ct == "Response"]],
                  median(ff@df[ff@coefmap$covariate[ct == "ResponseNiche"]])))
```
````

- [ ] **Step 3: `spiDE-model.Rmd` — add subsection "The testing framework".** Document the null-calibration simulation used to judge the fix: the clustered null (real between-patient variance, no response effect), the type-I metric at α = 0.05, the Response 0.99× regression anchor, and the requirement that ResponseNiche type-I rises toward 0.05 from below without overshooting.

- [ ] **Step 4: `spiDE-model.Rmd` — update the `@df` interpretation row** in §"Interpreting a mixed `SpiDEFit`" to describe the per-coefficient vector (Response ≈ `S − 2`, niche terms larger) under `"satterthwaite"` vs the scalar under `"between"`; and add to the "Power needs samples and effect size" caveat that per-coefficient df restores legitimate power for within-sample-informed niche terms without touching the main effect's honest between-sample df.

- [ ] **Step 5: `spiDE-mixed-benchmark.Rmd` — add a "Null calibration" section** reading the shipped CSV via the existing `bench()` helper:

````markdown
```{r calibration, fig.height = 4, fig.cap = "Null type-I error at alpha = 0.05. Satterthwaite lifts the over-conservative ResponseNiche tests toward nominal while the between-sample Response effect stays ~ nominal."}
cal <- bench("null_calibration.csv")
if (!is.null(cal)) {
  ggplot(cal, aes(S, type1, colour = interaction(method, df.method))) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = spide_pal$ref) +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~ effect) +
    labs(x = "samples (S)", y = "type-I error", colour = NULL,
         title = "Null calibration: method x df.method")
}
```
````

- [ ] **Step 6: `NEWS.md` — add a bullet** under a new top version stanza describing the `df.method = "satterthwaite"` per-coefficient reference df for mixed fits (default `"between"`, back-compatible).

- [ ] **Step 7: Build the vignettes to confirm they knit**

Run: `Rscript -e 'devtools::build_vignettes()'`
Expected: both vignettes render (the calibration figure reads the shipped CSV).

- [ ] **Step 8: Commit**

```bash
git add vignettes/spiDE-model.Rmd vignettes/spiDE-mixed-benchmark.Rmd NEWS.md
git commit -m "Document per-coefficient Satterthwaite df + testing framework in the vignettes"
```

---

## Final verification

- [ ] `Rscript -e 'devtools::document(); devtools::test()'` — all green.
- [ ] `Rscript -e 'devtools::check()'` — no new ERRORs/WARNINGs (`lmerTest`/`lme4` in Suggests are `skip_if_not_installed`-guarded).
- [ ] Confirm `git grep -n "@df" R/` shows only the vector-aware consumers (`.ptByCol`, `.nicheRecords`) and the fit-stage producer.

---

## Self-review

**Spec coverage:**
- Satterthwaite df formula + φ-invariance → Task 2 (`.satterthwaiteDF`), documented in Task 6.
- `Â = Cov(τ̂²)` via reduced-form REML info → Task 2 (`.varParamCov`); lmerTest oracle Tasks 2 & 4; slope fallback Task 4 (`tryCatch`→`invert_mat`).
- Response ≈ S−2 regression → Task 3 Step 1.
- Shared-across-genes, hot loop untouched → `@df` produced at fit stage (Task 3), consumers only index (Task 1). No per-gene df (out of scope) — honoured.
- `@df` scalar→vector representation + `.ptByCol` broadcast → Task 1; both consumers (inference + fdr) → Task 1.
- `df.method` arg, default "between", ignored for random="none" → Tasks 3 (fitSpiDE/spiDE) ; between-identical → Task 3 Step 8.
- lmerTest Suggests → Task 2 Step 3.
- Assessment (null sim, local grid + assertion), KR decision gate → Task 5.
- `null_calibration.csv` artifact → Task 5; benchmark figure → Task 6 Step 5.
- Vignette updates (Corrected-inference point 3, Satterthwaite subsection, testing framework, @df row, power caveat) → Task 6.
- Tests: lmerTest cross-check, Response df≈S−2, between identical, per-column threading, slot back-compat → Tasks 1–5.

**Placeholder scan:** none — every code step carries real code.

**Type consistency:** `.ptByCol(tmat, df, lower.tail)`, `.varParamCov(A, minv, pen, re_group, ncells)`, `.satterthwaiteDF(A, minv, pen, re_group, tau2, tested, ncells, tested_names)`, `.repWeights(Y, alpha, W, psi)` used consistently across tasks; `@df` shape contract (NULL/scalar/vector) consistent in producer (Task 3) and both consumers (Task 1).
