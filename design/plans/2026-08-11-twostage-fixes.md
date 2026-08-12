# Two-Stage Pipeline Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild `twoStageSpiDE()` per the approved spec: SpaNorm-anchored one-step stage-1 slopes (log-mean scale, joint over niches), patient pooling, and a weighted moderated limma stage 2 — with diagnostics and no deprecation shims.

**Architecture:** Stage 1 extracts the stored SpaNorm fit (biology/LS split via `wtype`), forms working residuals at the full fitted mean with NB working weights, adds the biology linear predictor back, and runs one joint WLS per (sample, index) subset. Slopes are precision-pooled per patient, a DerSimonian–Laird τ² is estimated per (index, niche), and limma tests the condition contrast with 1/(v + τ²) observation weights.

**Tech Stack:** R/Bioconductor package (roxygen2, testthat edition 3), SpaNorm (>= 1.7.7), limma (new Import), S4Vectors, Matrix.

## Global Constraints

- All work on branch `twostage-fixes` (already created from `origin/twostage-estimator`; spec committed as fa468b9).
- `"nbresid"`, `"analytic"`, `"auto"` stage-1 options are **removed outright** — no deprecation messages (unpublished code, user-directed).
- Results table schema is unchanged: `gene, ct_index, ct_niche, coef, t, p.niche, fdr.niche, DirectionNiche, bandwidth.max`.
- Dot-separated argument names (`min.cells`, `patient.covariates`, `lambda.a`, `maxit.psi`) mirror existing conventions — do not rename.
- Every new function gets roxygen with `@noRd` (internal) and is placed as specified; run `devtools::document()` only in Task 11 (no exported-API changes before then except `twoStageSpiDE` itself).
- Commit messages: sentence-style summary line (repo convention, no `feat:` prefixes), body optional, ending with `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- `docs/` is gitignored (pkgdown); commit plan/spec files with `git add -f`.
- Run tests from the package root: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/<file>")'`.
- The toy fixture is `spiDE:::.toySPE()` (8 samples, cell types A–D, planted G1/A/B effect, `condition` Responder/NonResponder, patient covariate `Age`).

## File Structure

- `R/AllClasses.R` — modify: add `diagnostics` slot to `SpiDEResults`.
- `R/twostage.R` — rewrite: exported `twoStageSpiDE()` orchestration only; helpers move out. `.nbPearsonResiduals()`, `.pearsonResiduals()` deleted; `.looksLikeCounts()` kept (input guard only).
- `R/twostage-stage1.R` — create: `.spanormComponents()`, `.stage1Epsilon()`, `.jointSlopes()`, `.nicheBasisR2()`, `.sampleSlopes()`.
- `R/twostage-stage2.R` — create: `.poolPatientSlopes()`, `.tau2DL()`, `.limmaStage2()`, `.inclusionDiagnostics()`.
- `tests/testthat/helper-twostage.R` — create: `toy_spanorm_spe()`.
- `tests/testthat/test-twostage-stage1.R`, `test-twostage-stage2.R` — create.
- `tests/testthat/test-twostage.R` — modify: existing tests updated to the new option surface + new E2E and diagnostics tests.
- `DESCRIPTION` — modify: add `limma` to Imports.
- `/Users/uqdbhuva/R_packages/spiDE-twostage-bench/scenarios.R`, `runner-adv.R` — create (outside the package; adversarial validation).

---

### Task 1: `diagnostics` slot on `SpiDEResults`

**Files:**
- Modify: `R/AllClasses.R` (the `SpiDEResults` setClass at ~line 230)
- Test: `tests/testthat/test-twostage-stage2.R` (create the file with this one test)

**Interfaces:**
- Produces: `SpiDEResults@diagnostics` — a `list`, prototype `list()`. Later tasks store `list(r2 = <data.frame>, inclusion = <data.frame>, tau2 = <data.frame>)`.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-twostage-stage2.R
test_that("SpiDEResults carries a diagnostics list slot", {
  r <- new("SpiDEResults", fits = list(), sigma = 30, condition = "condition",
           mode = "condition", index = "A", niche = "B",
           covariates = character(), coldata = S4Vectors::DataFrame(),
           gene.weights = NULL, p.cauchy.pos = NULL, p.cauchy.neg = NULL,
           results = data.frame(), fdr = 0.05, call = NULL,
           diagnostics = list(r2 = data.frame()))
  expect_identical(names(r@diagnostics), "r2")
  # default construction still works and defaults to an empty list
  r0 <- new("SpiDEResults", fits = list(), sigma = 30, condition = "c",
            mode = "condition", index = "A", niche = "B",
            covariates = character(), coldata = S4Vectors::DataFrame(),
            gene.weights = NULL, p.cauchy.pos = NULL, p.cauchy.neg = NULL,
            results = data.frame(), fdr = 0.05, call = NULL)
  expect_identical(r0@diagnostics, list())
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-twostage-stage2.R")'`
Expected: FAIL — `invalid name for slot of class "SpiDEResults": diagnostics`

- [ ] **Step 3: Implement**

In `R/AllClasses.R`, inside the `SpiDEResults` `slots = c(...)`, after `results.patient = "data.frame",` add:

```r
    # Two-stage estimator diagnostics: R2 of niche columns on the biology
    # basis, the per-index patient inclusion table, and the tau2 table.
    # Empty for the GLM path.
    diagnostics = "list",
```

and in the `prototype = list(...)` add `diagnostics = list(),`. Also add to the class roxygen block:

```r
#' @slot diagnostics a list of two-stage diagnostic tables (`r2`, `inclusion`,
#'   `tau2`); empty for the GLM path.
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: PASS. Also run the migration/e2e-adjacent files to catch prototype fallout:
`Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-twostage.R")'` — existing tests must still pass (old constructor calls omit the slot; the prototype fills it).

- [ ] **Step 5: Commit**

```bash
git add R/AllClasses.R tests/testthat/test-twostage-stage2.R
git commit -m "Add a diagnostics slot to SpiDEResults for the two-stage path

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: SpaNorm fit extraction

**Files:**
- Create: `R/twostage-stage1.R`
- Create: `tests/testthat/helper-twostage.R`
- Create: `tests/testthat/test-twostage-stage1.R`

**Interfaces:**
- Consumes: `S4Vectors::metadata(spe)$SpaNorm` — a `SpaNormFit` with `$W` (cells × p), `$alpha` (genes × p), `$gmean` (genes), `$psi` (genes), `$wtype` (factor over `"biology"`/`"ls"`/`"batch"`), `$ngenes`, `$ncells`.
- Produces: `.spanormComponents(spe)` → `list(alpha, gmean, psi, W, bio)` where `bio` is a logical over `ncol(W)` marking biology columns. Task 3 consumes this list verbatim.

- [ ] **Step 1: Write the test helper**

```r
# tests/testthat/helper-twostage.R
# A toy SPE with a hand-built SpaNormFit in metadata: design [logLS, bio.x,
# bio.y] fit with SpaNorm::fitNB, wtype marking the split. Small and fast;
# exercises the exact contract .spanormComponents() reads.
toy_spanorm_spe <- function(n_genes = 30, sigma = 30) {
  spe <- spiDE:::.toySPE(n_genes = n_genes)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  xy <- SpatialExperiment::spatialCoords(spe)
  logLS <- log(pmax(colSums(Y), 1))
  W <- cbind(logLS = logLS - mean(logLS),
             bio.x = as.numeric(scale(xy[, 1])),
             bio.y = as.numeric(scale(xy[, 2])))
  f <- SpaNorm::fitNB(Y, W, lambda.a = 0, verbose = FALSE)
  fit <- methods::new("SpaNormFit",
    ngenes = nrow(Y), ncells = ncol(Y), gene.model = "nb",
    df.tps = c(1L, 1L, 1L, 1L), sample.p = 1, lambda.a = c(0, 0),
    batch = NULL, W = W, alpha = f$alpha, gmean = f$gmean, psi = f$psi,
    wtype = factor(c("ls", "biology", "biology")),
    loglik = f$loglik,
    sampling = if (!is.null(f$sampling)) f$sampling
               else factor(rep("fit", ncol(Y))))
  S4Vectors::metadata(spe)$SpaNorm <- fit
  buildNiches(spe, sigma = sigma, verbose = FALSE)
}
```

- [ ] **Step 2: Write the failing tests**

```r
# tests/testthat/test-twostage-stage1.R
test_that(".spanormComponents extracts and validates the stored fit", {
  spe <- toy_spanorm_spe()
  comp <- spiDE:::.spanormComponents(spe)
  expect_named(comp, c("alpha", "gmean", "psi", "W", "bio"))
  expect_identical(dim(comp$W), c(ncol(spe), 3L))
  expect_identical(comp$bio, c(FALSE, TRUE, TRUE))
  expect_length(comp$psi, nrow(spe))
})

test_that(".spanormComponents errors usefully without a fit or on mismatch", {
  spe <- buildNiches(spiDE:::.toySPE(n_genes = 20), sigma = 30, verbose = FALSE)
  expect_error(spiDE:::.spanormComponents(spe), "SpaNorm::SpaNorm")
  spe2 <- toy_spanorm_spe()
  expect_error(spiDE:::.spanormComponents(spe2[, 1:10]), "does not match")
})
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-twostage-stage1.R")'`
Expected: FAIL — `'.spanormComponents' not found`

- [ ] **Step 4: Implement**

```r
# R/twostage-stage1.R
# Stage-1 machinery for the two-stage estimator: SpaNorm fit extraction, the
# one-step working response, the joint weighted slope fit, and the basis-R2
# diagnostic. See design/specs/2026-08-11-twostage-fixes-design.md.

#' Extract and validate the stored SpaNorm fit
#'
#' The fit supplies the biology/LS design split (`wtype`), the coefficients
#' and the dispersion that stage 1 linearises around. Kept as small pieces --
#' full genes x cells matrices are only materialised per (sample, index)
#' subset by .stage1Epsilon().
#' @param spe a SpatialExperiment previously normalised with SpaNorm().
#' @return list(alpha, gmean, psi, W, bio) with `bio` a logical over columns
#'   of `W` marking the biology block.
#' @noRd
.spanormComponents <- function(spe) {
  fit <- S4Vectors::metadata(spe)$SpaNorm
  if (is.null(fit)) {
    stop("no SpaNorm fit found in metadata(spe)$SpaNorm; run ",
         "SpaNorm::SpaNorm(spe) first (stage1 = \"spanorm\" needs it), or ",
         "use stage1 = \"ols\"")
  }
  if (fit$ncells != ncol(spe) || fit$ngenes != nrow(spe)) {
    stop("the stored SpaNorm fit (", fit$ngenes, " x ", fit$ncells,
         ") does not match spe (", nrow(spe), " x ", ncol(spe), "); ",
         "re-run SpaNorm::SpaNorm() on this object")
  }
  wt <- as.character(fit$wtype)
  list(alpha = fit$alpha, gmean = fit$gmean, psi = fit$psi, W = fit$W,
       bio = wt == "biology")
}
```

- [ ] **Step 5: Run tests to verify they pass**

Same command. Expected: PASS (3 tests).

- [ ] **Step 6: Commit**

```bash
git add R/twostage-stage1.R tests/testthat/helper-twostage.R tests/testthat/test-twostage-stage1.R
git commit -m "Extract and validate the stored SpaNorm fit for two-stage stage 1

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: The one-step working response

**Files:**
- Modify: `R/twostage-stage1.R`
- Test: `tests/testthat/test-twostage-stage1.R`

**Interfaces:**
- Consumes: `.spanormComponents()` output.
- Produces: `.stage1Epsilon(Y, comp, cells, epsilon)` → `list(eps, w)`, both genes × length(cells) matrices. `eps = eta_bio + z` (`"addback"`) or `z` (`"residual"`); `w = mu/(1 + psi*mu)`. Task 4 consumes `eps`/`w`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-twostage-stage1.R`:

```r
test_that(".stage1Epsilon linearises log counts minus the LS effect", {
  spe <- toy_spanorm_spe()
  comp <- spiDE:::.spanormComponents(spe)
  cells <- 1:50
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  se <- spiDE:::.stage1Epsilon(Y, comp, cells, epsilon = "addback")
  expect_identical(dim(se$eps), c(nrow(Y), 50L))
  expect_identical(dim(se$w), dim(se$eps))
  expect_true(all(se$w > 0 & is.finite(se$eps)))
  # at y = 0 the addback response is eta_bio - 1 exactly (z = -1)
  Wc <- comp$W[cells, , drop = FALSE]
  eta_bio <- comp$gmean +
    tcrossprod(comp$alpha[, comp$bio, drop = FALSE],
               Wc[, comp$bio, drop = FALSE])
  zero <- Y[, cells] == 0
  expect_equal(se$eps[zero], (eta_bio - 1)[zero])
  # "residual" drops the biology term
  sr <- spiDE:::.stage1Epsilon(Y, comp, cells, epsilon = "residual")
  expect_equal(se$eps - sr$eps, eta_bio, tolerance = 1e-12)
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-twostage-stage1.R")'`
Expected: FAIL — `'.stage1Epsilon' not found`

- [ ] **Step 3: Implement**

Append to `R/twostage-stage1.R`:

```r
#' One-step working response and weights for a cell subset
#'
#' Linearises log counts at the FULL fitted mean (the correct expansion
#' point: handles zeros, no baseline leakage) and, under "addback", returns
#' the biology component so that, net, only the LS (and batch) effect is
#' removed. Materialises genes x subset matrices only.
#' @param Y counts matrix (genes x all cells; dense).
#' @param comp the .spanormComponents() list.
#' @param cells integer/logical index of the subset's cells.
#' @param epsilon "addback" (default; eta_bio + z) or "residual" (z).
#' @return list(eps, w): the response and NB working weights, genes x cells.
#' @noRd
.stage1Epsilon <- function(Y, comp, cells, epsilon = c("addback", "residual")) {
  epsilon <- match.arg(epsilon)
  Wc <- comp$W[cells, , drop = FALSE]
  eta <- comp$gmean + tcrossprod(comp$alpha, Wc)
  mu <- exp(eta)
  z <- (Y[, cells, drop = FALSE] - mu) / mu
  w <- mu / (1 + comp$psi * mu)
  eps <- if (epsilon == "addback") {
    comp$gmean + tcrossprod(comp$alpha[, comp$bio, drop = FALSE],
                            Wc[, comp$bio, drop = FALSE]) + z
  } else {
    z
  }
  list(eps = eps, w = w)
}
```

- [ ] **Step 4: Run to verify it passes**

Same command. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/twostage-stage1.R tests/testthat/test-twostage-stage1.R
git commit -m "Add the one-step working response and NB weights for stage 1

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: Joint weighted slopes with per-gene weights

**Files:**
- Modify: `R/twostage-stage1.R`
- Test: `tests/testthat/test-twostage-stage1.R`

**Interfaces:**
- Consumes: `eps`, `w` from Task 3; a niche matrix `X` (cells × k, centred columns).
- Produces: `.jointSlopes(eps, w, X)` → `list(beta, var)`, both genes × k matrices (NA columns for niches dropped by the rank guard). Design is `[1, X]`; the intercept row is not returned. Task 6 consumes this.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-twostage-stage1.R`:

```r
test_that(".jointSlopes recovers planted slopes and matches lm() weights", {
  set.seed(7)
  n <- 200; k <- 2; G <- 5
  X <- cbind(a = rnorm(n), b = rnorm(n))
  X <- sweep(X, 2, colMeans(X))
  beta_true <- matrix(c(0.5, -0.2), nrow = G, ncol = k, byrow = TRUE)
  eps <- beta_true %*% t(X) + matrix(rnorm(G * n, sd = 0.3), G, n)
  w <- matrix(runif(G * n, 0.5, 2), G, n)
  js <- spiDE:::.jointSlopes(eps, w, X)
  expect_identical(dim(js$beta), c(G, k))
  # gene 1 must agree with a weighted lm to numerical precision
  fit <- lm(eps[1, ] ~ X, weights = w[1, ])
  expect_equal(unname(js$beta[1, ]), unname(coef(fit)[-1]), tolerance = 1e-8)
  sm <- summary(fit)
  expect_equal(unname(js$var[1, ]), unname(sm$coefficients[-1, 2]^2),
               tolerance = 1e-6)
})

test_that(".jointSlopes drops degenerate niche columns to NA", {
  set.seed(8)
  n <- 100
  X <- cbind(a = rnorm(n), b = 0)          # b has no variance
  X <- sweep(X, 2, colMeans(X))
  eps <- matrix(rnorm(3 * n), 3, n)
  w <- matrix(1, 3, n)
  js <- spiDE:::.jointSlopes(eps, w, X)
  expect_true(all(is.na(js$beta[, "b"])) && all(is.na(js$var[, "b"])))
  expect_true(all(is.finite(js$beta[, "a"])))
})
```

- [ ] **Step 2: Run to verify it fails**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-twostage-stage1.R")'`
Expected: FAIL — `'.jointSlopes' not found`

- [ ] **Step 3: Implement**

Append to `R/twostage-stage1.R`:

```r
#' Joint WLS of the working response on all niche columns, per gene
#'
#' Design [1, X] (the intercept absorbs each gene's level, so no per-gene
#' weighted centring is needed). Per-gene weights preclude one shared Gram
#' matrix, so cross-products are batched as matrix products over the
#' upper-triangle column pairs, then each gene's small (k+1) system is
#' solved in a loop. Degenerate columns (weighted variance ~ 0) are dropped
#' up front and return NA. Variance is the usual WLS sandwich-free form
#' phi_g * (X'W_gX)^{-1} with phi_g the weighted RSS over n - rank.
#' @param eps,w genes x cells response and weights (Task 3).
#' @param X cells x k centred niche matrix.
#' @return list(beta, var): genes x k, NA for dropped columns.
#' @noRd
.jointSlopes <- function(eps, w, X) {
  G <- nrow(eps); n <- ncol(eps); k <- ncol(X)
  out_b <- out_v <- matrix(NA_real_, G, k,
                           dimnames = list(rownames(eps), colnames(X)))
  keep <- apply(X, 2, function(x) stats::var(x) > 1e-10)
  if (!any(keep)) return(list(beta = out_b, var = out_v))
  D <- cbind(`(Intercept)` = 1, X[, keep, drop = FALSE])
  p <- ncol(D)
  # per-gene Gram entries: (X'W_gX)_{ij} = sum_c w_gc D_ci D_cj
  ut <- which(upper.tri(diag(p), diag = TRUE), arr.ind = TRUE)
  P <- D[, ut[, 1], drop = FALSE] * D[, ut[, 2], drop = FALSE]  # n x npairs
  Gm <- w %*% P                                                  # G x npairs
  Rhs <- (w * eps) %*% D                                         # G x p
  wrss_tot <- rowSums(w * eps^2)
  for (g in seq_len(G)) {
    A <- matrix(0, p, p)
    A[cbind(ut[, 1], ut[, 2])] <- Gm[g, ]
    A[cbind(ut[, 2], ut[, 1])] <- Gm[g, ]
    Ai <- tryCatch(solve(A), error = function(e) NULL)
    if (is.null(Ai)) next
    bg <- Ai %*% Rhs[g, ]
    df <- max(n - p, 1)
    phi <- max((wrss_tot[g] - sum(bg * Rhs[g, ])) / df, 1e-12)
    out_b[g, keep] <- bg[-1]
    out_v[g, keep] <- phi * diag(Ai)[-1]
  }
  list(beta = out_b, var = out_v)
}
```

- [ ] **Step 4: Run to verify it passes**

Same command. Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/twostage-stage1.R tests/testthat/test-twostage-stage1.R
git commit -m "Add the per-gene joint weighted niche slope fit

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: Basis-overlap R² diagnostic

**Files:**
- Modify: `R/twostage-stage1.R`
- Test: `tests/testthat/test-twostage-stage1.R`

**Interfaces:**
- Consumes: subset niche matrix `X` (cells × k), subset biology basis `B` (cells × b).
- Produces: `.nicheBasisR2(X, B)` → named numeric length k (R² of each niche column on `[1, B]`). Task 8's orchestration collects these into `diagnostics$r2`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-twostage-stage1.R`:

```r
test_that(".nicheBasisR2 is 1 in-span and ~0 orthogonal", {
  set.seed(9)
  n <- 120
  B <- cbind(rnorm(n), rnorm(n))
  in_span <- 2 * B[, 1] - B[, 2] + 3
  ortho <- residuals(lm(rnorm(n) ~ B))
  r2 <- spiDE:::.nicheBasisR2(cbind(s = in_span, o = ortho), B)
  expect_equal(unname(r2["s"]), 1, tolerance = 1e-8)
  expect_equal(unname(r2["o"]), 0, tolerance = 1e-8)
})
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL — `'.nicheBasisR2' not found`

- [ ] **Step 3: Implement**

Append to `R/twostage-stage1.R`:

```r
#' R2 of each niche column on the biology basis, within a subset
#'
#' Measures the attenuation the "residual" response would suffer and the
#' smooth-trend overlap the "addback" response is exposed to. Reported per
#' (sample, index) subset in the diagnostics.
#' @noRd
.nicheBasisR2 <- function(X, B) {
  f <- stats::lm.fit(cbind(1, B), X)
  res <- as.matrix(f$residuals)
  tss <- colSums(sweep(X, 2, colMeans(X))^2)
  r2 <- 1 - colSums(res^2) / pmax(tss, 1e-12)
  stats::setNames(pmin(pmax(r2, 0), 1), colnames(X))
}
```

- [ ] **Step 4: Run to verify it passes**  — same command, PASS.

- [ ] **Step 5: Commit**

```bash
git add R/twostage-stage1.R tests/testthat/test-twostage-stage1.R
git commit -m "Add the niche-on-biology-basis R2 diagnostic

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: `.sampleSlopes()` — the new stage-1 driver for all paths

**Files:**
- Modify: `R/twostage-stage1.R` (add `.sampleSlopes()`)
- Modify: `R/twostage.R` (delete `.patientSlopes()`, `.nbPearsonResiduals()`, `.pearsonResiduals()`)
- Test: `tests/testthat/test-twostage-stage1.R`

**Interfaces:**
- Consumes: Tasks 2–5.
- Produces: `.sampleSlopes(Y, E, comp, nm, ct, smp, idx_types, min.cells, stage1, epsilon, winsor, lambda.a, maxit.psi, backend, verbose)` → `list(beta = <list per index of genes x niches x samples arrays>, var = <same shape>, r2 = <data.frame sample, index, niche, r2>, ncells = <data.frame sample, index, n>)`. `E` is log-CPM (only for `"ols"`, else NULL); `comp` is the components list (only for `"spanorm"`, else NULL). Self-niche columns are fitted but set to NA afterwards by the caller (schema keeps all niches).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-twostage-stage1.R`:

```r
test_that(".sampleSlopes returns aligned beta/var arrays for all paths", {
  spe <- toy_spanorm_spe()
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  nm <- as.matrix(SingleCellExperiment::reducedDim(spe, "Niche30"))
  ct <- as.character(spe$cell_type); smp <- as.character(spe$sample_id)
  comp <- spiDE:::.spanormComponents(spe)
  s1 <- spiDE:::.sampleSlopes(Y, NULL, comp, nm, ct, smp,
                              idx_types = c("A", "B"), min.cells = 10,
                              stage1 = "spanorm", epsilon = "addback")
  expect_named(s1$beta, c("A", "B"))
  expect_identical(dim(s1$beta$A), dim(s1$var$A))
  expect_identical(dimnames(s1$beta$A)[[2]], colnames(nm))
  expect_true(all(s1$var$A >= 0, na.rm = TRUE))
  expect_true(all(c("sample", "index", "niche", "r2") %in% names(s1$r2)))
  # ols path: same shapes, no comp needed
  lib <- colSums(Y)
  E <- log1p(sweep(Y, 2, mean(lib) / pmax(lib, 1), "*"))
  s2 <- spiDE:::.sampleSlopes(Y, E, NULL, nm, ct, smp,
                              idx_types = "A", min.cells = 10,
                              stage1 = "ols", epsilon = "addback")
  expect_identical(dim(s2$beta$A), dim(s1$beta$A))
})

test_that("subsets below min.cells contribute NA, and are counted", {
  spe <- toy_spanorm_spe()
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  nm <- as.matrix(SingleCellExperiment::reducedDim(spe, "Niche30"))
  ct <- as.character(spe$cell_type); smp <- as.character(spe$sample_id)
  comp <- spiDE:::.spanormComponents(spe)
  s <- spiDE:::.sampleSlopes(Y, NULL, comp, nm, ct, smp, idx_types = "A",
                             min.cells = 10000L, stage1 = "spanorm",
                             epsilon = "addback")
  expect_true(all(is.na(s$beta$A)))
  expect_true(all(s$ncells$n < 10000L))
})
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL — `'.sampleSlopes' not found`

- [ ] **Step 3: Implement**

Append to `R/twostage-stage1.R`:

```r
#' Per-(sample, index) joint niche slopes for every stage-1 path
#'
#' Replaces .patientSlopes(). All paths run the same joint fit and return a
#' variance for every slope (stage 2 needs it): "spanorm" uses the one-step
#' response/weights, "ols" uses log-CPM with unit weights, "nb" fits the
#' per-subset NB GLM on [1, loglib, niches] and reads Fisher variances.
#' @return list(beta, var, r2, ncells) -- see Interfaces in the plan.
#' @noRd
.sampleSlopes <- function(Y, E, comp, nm, ct, smp, idx_types, min.cells,
                          stage1 = c("spanorm", "ols", "nb"),
                          epsilon = c("addback", "residual"),
                          winsor = 4, lambda.a = 0, maxit.psi = 2,
                          backend = "cpu", verbose = FALSE) {
  stage1 <- match.arg(stage1); epsilon <- match.arg(epsilon)
  smps <- sort(unique(smp)); niches <- colnames(nm)
  gn <- rownames(Y)
  loglib <- log(pmax(colSums(Y), 1))
  beta <- var <- stats::setNames(vector("list", length(idx_types)), idx_types)
  r2 <- list(); ncl <- list()
  for (ix in idx_types) {
    A <- V <- array(NA_real_, c(nrow(Y), length(niches), length(smps)),
                    dimnames = list(gn, niches, smps))
    for (ss in smps) {
      cells <- which(ct == ix & smp == ss)
      ncl[[length(ncl) + 1L]] <- data.frame(sample = ss, index = ix,
                                            n = length(cells))
      if (length(cells) < min.cells) next
      X <- sweep(nm[cells, , drop = FALSE], 2,
                 colMeans(nm[cells, , drop = FALSE]))
      if (stage1 == "spanorm") {
        se <- .stage1Epsilon(Y, comp, cells, epsilon)
        js <- .jointSlopes(se$eps, se$w, X)
        Bbio <- comp$W[cells, comp$bio, drop = FALSE]
        r2[[length(r2) + 1L]] <- data.frame(
          sample = ss, index = ix, niche = niches,
          r2 = as.numeric(.nicheBasisR2(X, Bbio)))
      } else if (stage1 == "ols") {
        Ec <- E[, cells, drop = FALSE]
        js <- .jointSlopes(Ec, matrix(1, nrow(Ec), ncol(Ec)), X)
      } else {                                        # "nb" reference path
        keep <- apply(X, 2, function(x) stats::var(x) > 1e-10)
        Wn <- cbind(`(Intercept)` = 1,
                    loglib = loglib[cells] - mean(loglib[cells]),
                    X[, keep, drop = FALSE])
        if (qr(Wn)$rank < ncol(Wn)) next
        f <- try(SpaNorm::fitNB(Y[, cells, drop = FALSE], Wn,
                                lambda.a = lambda.a, winsor = winsor,
                                maxit.psi = maxit.psi, backend = backend,
                                verbose = FALSE), silent = TRUE)
        if (inherits(f, "try-error")) next
        mu <- SpaNorm::calculateMu(f$gmean, f$alpha, Wn)
        js <- list(beta = matrix(NA_real_, nrow(Y), length(niches),
                                 dimnames = list(gn, niches)),
                   var = matrix(NA_real_, nrow(Y), length(niches),
                                dimnames = list(gn, niches)))
        wnb <- mu / (1 + f$psi * mu)
        for (g in seq_len(nrow(Y))) {
          info <- crossprod(Wn * wnb[g, ], Wn)
          vc <- diag(SpaNorm::invert_mat(info))
          js$beta[g, keep] <- f$alpha[g, -(1:2)]
          js$var[g, keep] <- vc[-(1:2)]
        }
      }
      A[, colnames(js$beta), ss] <- js$beta
      V[, colnames(js$var), ss] <- js$var
    }
    beta[[ix]] <- A; var[[ix]] <- V
    if (verbose) {
      message(sprintf("  %-22s %d samples >= %d cells", ix,
                      sum(apply(!is.na(A[1, , , drop = FALSE]), 3, any)),
                      min.cells))
    }
  }
  list(beta = beta, var = var,
       r2 = if (length(r2)) do.call(rbind, r2) else
         data.frame(sample = character(), index = character(),
                    niche = character(), r2 = numeric()),
       ncells = do.call(rbind, ncl))
}
```

In `R/twostage.R`, delete the functions `.patientSlopes()`, `.nbPearsonResiduals()`, `.pearsonResiduals()` entirely (keep `.looksLikeCounts()`). `twoStageSpiDE()` is now broken — that is expected until Task 11; do not run the e2e test file in this task.

- [ ] **Step 4: Run to verify it passes**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-twostage-stage1.R")'`
Expected: PASS (all stage-1 tests).

- [ ] **Step 5: Commit**

```bash
git add R/twostage-stage1.R R/twostage.R tests/testthat/test-twostage-stage1.R
git commit -m "Drive all stage-1 paths through the joint slope fit with variances

Removes .patientSlopes, .nbPearsonResiduals and .pearsonResiduals; the
nbresid/analytic/auto options are gone per the spec (unpublished code).

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: Patient pooling

**Files:**
- Create: `R/twostage-stage2.R`
- Test: `tests/testthat/test-twostage-stage2.R`

**Interfaces:**
- Consumes: one index's `beta`/`var` arrays (genes × niches × samples) from Task 6; `sample2patient` — named character, `names()` = samples, values = patients.
- Produces: `.poolPatientSlopes(B, V, sample2patient)` → `list(beta, var)` arrays genes × niches × patients (precision-weighted; single-sample patients pass through; all-NA stays NA).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-twostage-stage2.R`:

```r
test_that(".poolPatientSlopes precision-pools cores within a patient", {
  B <- array(NA_real_, c(1, 1, 3),
             dimnames = list("g", "n", c("s1", "s2", "s3")))
  V <- B
  B[1, 1, ] <- c(1, 3, 10); V[1, 1, ] <- c(1, 0.5, 2)
  s2p <- c(s1 = "p1", s2 = "p1", s3 = "p2")
  pooled <- spiDE:::.poolPatientSlopes(B, V, s2p)
  expect_identical(dimnames(pooled$beta)[[3]], c("p1", "p2"))
  expect_equal(pooled$beta["g", "n", "p1"], (1 / 1 + 3 / 0.5) / (1 + 2))
  expect_equal(pooled$var["g", "n", "p1"], 1 / 3)
  expect_equal(pooled$beta["g", "n", "p2"], 10)   # single core passes through
  # an NA core is ignored, not contagious
  B[1, 1, 2] <- NA
  pooled2 <- spiDE:::.poolPatientSlopes(B, V, s2p)
  expect_equal(pooled2$beta["g", "n", "p1"], 1)
})
```

- [ ] **Step 2: Run to verify it fails** — expected: FAIL, `'.poolPatientSlopes' not found`

- [ ] **Step 3: Implement**

```r
# R/twostage-stage2.R
# Stage-2 machinery: patient pooling of per-core slopes, the DL tau2, the
# weighted moderated limma contrast, and the dropout diagnostics.

#' Precision-weighted pooling of per-core slopes within patients
#'
#' Cores share their patient's slope, so fixed-effect (1/v) pooling is
#' correct here (the between-patient component enters later, in tau2).
#' @noRd
.poolPatientSlopes <- function(B, V, sample2patient) {
  pats <- unique(unname(sample2patient[dimnames(B)[[3]]]))
  db <- dimnames(B)
  Bp <- Vp <- array(NA_real_, c(dim(B)[1:2], length(pats)),
                    dimnames = list(db[[1]], db[[2]], pats))
  for (p in pats) {
    ss <- dimnames(B)[[3]][sample2patient[dimnames(B)[[3]]] == p]
    b <- B[, , ss, drop = FALSE]; v <- V[, , ss, drop = FALSE]
    wt <- 1 / v
    wt[is.na(b) | !is.finite(wt)] <- NA
    swt <- apply(wt, 1:2, sum, na.rm = TRUE)
    num <- apply(b * wt, 1:2, sum, na.rm = TRUE)
    est <- num / swt
    est[swt == 0] <- NA_real_
    Bp[, , p] <- est
    vp <- 1 / swt
    vp[swt == 0] <- NA_real_
    Vp[, , p] <- vp
  }
  list(beta = Bp, var = Vp)
}
```

- [ ] **Step 4: Run to verify it passes** — same command, PASS.

- [ ] **Step 5: Commit**

```bash
git add R/twostage-stage2.R tests/testthat/test-twostage-stage2.R
git commit -m "Pool per-core slopes within patients by precision

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: The DerSimonian–Laird τ²

**Files:**
- Modify: `R/twostage-stage2.R`
- Test: `tests/testthat/test-twostage-stage2.R`

**Interfaces:**
- Consumes: per-(index, niche) slope matrix `B` (genes × patients), variance matrix `V`, design matrix `X` (patients × p, from `model.matrix(~ condition + covariates)`).
- Produces: `.tau2DL(B, V, X)` → a single non-negative numeric (per-gene DL estimates, median-pooled, floored at 0). Task 9 consumes it.

- [ ] **Step 1: Write the failing tests**

```r
test_that(".tau2DL recovers heterogeneity and floors at zero", {
  set.seed(11)
  S <- 40; G <- 200
  X <- cbind(1, rep(0:1, each = S / 2))
  v <- matrix(runif(G * S, 0.05, 0.1), G, S)
  # homogeneous slopes: tau2 ~ 0
  B0 <- matrix(rnorm(G * S, sd = sqrt(v)), G, S)
  expect_lt(spiDE:::.tau2DL(B0, v, X), 0.02)
  # true tau2 = 0.5 added between patients
  u <- matrix(rnorm(G * S, sd = sqrt(0.5)), G, S)
  B1 <- B0 + u
  t2 <- spiDE:::.tau2DL(B1, v, X)
  expect_gt(t2, 0.25); expect_lt(t2, 1)
})
```

- [ ] **Step 2: Run to verify it fails** — expected: FAIL, `'.tau2DL' not found`

- [ ] **Step 3: Implement**

Append to `R/twostage-stage2.R`:

```r
#' Between-patient variance component, DerSimonian-Laird, pooled over genes
#'
#' Per gene: weighted (1/v) regression on the stage-2 design, Cochran's
#' Q from the weighted residuals, DL moment estimate
#' (Q - (S - p)) / (sum(w) - tr((X'WX)^{-1} X'W^2X)); then the median over
#' genes, floored at 0. One tau2 per (index, niche) -- matching how the
#' package shares strength across genes elsewhere.
#' @noRd
.tau2DL <- function(B, V, X) {
  p <- ncol(X)
  t2 <- apply(cbind(B, V), 1, function(row) {
    S <- length(row) / 2
    b <- row[seq_len(S)]; v <- row[S + seq_len(S)]
    ok <- is.finite(b) & is.finite(v) & v > 0
    if (sum(ok) < p + 2) return(NA_real_)
    w <- 1 / v[ok]; Xo <- X[ok, , drop = FALSE]; bo <- b[ok]
    XtW <- t(Xo * w)
    M <- tryCatch(solve(XtW %*% Xo), error = function(e) NULL)
    if (is.null(M)) return(NA_real_)
    res <- bo - Xo %*% (M %*% (XtW %*% bo))
    Q <- sum(w * res^2)
    denom <- sum(w) - sum(diag(M %*% (t(Xo * w^2) %*% Xo)))
    if (denom <= 0) return(NA_real_)
    max((Q - (sum(ok) - p)) / denom, 0)
  })
  med <- stats::median(t2, na.rm = TRUE)
  if (!is.finite(med)) 0 else med
}
```

- [ ] **Step 4: Run to verify it passes** — same command, PASS.

- [ ] **Step 5: Commit**

```bash
git add R/twostage-stage2.R tests/testthat/test-twostage-stage2.R
git commit -m "Estimate the between-patient slope variance by pooled DerSimonian-Laird

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 9: Weighted moderated limma contrast

**Files:**
- Modify: `R/twostage-stage2.R`
- Modify: `DESCRIPTION` (add `limma` to Imports, alphabetical position after `BiocParallel`)
- Test: `tests/testthat/test-twostage-stage2.R`

**Interfaces:**
- Consumes: `B`, `V` (genes × patients), τ² scalar (Task 8), design `X` with the condition contrast in column 2.
- Produces: `.limmaStage2(B, V, tau2, X)` → data.frame(gene, coef, t, p) — the moderated condition-coefficient statistics per gene; NA rows for genes with too few finite slopes.

- [ ] **Step 1: Write the failing tests**

```r
test_that(".limmaStage2 matches a weighted lm directionally and handles NA", {
  set.seed(12)
  S <- 30; G <- 50
  grp <- rep(0:1, each = S / 2)
  X <- cbind(`(Intercept)` = 1, g = grp)
  V <- matrix(0.1, G, S)
  B <- matrix(rnorm(G * S, sd = 0.4), G, S,
              dimnames = list(paste0("g", 1:G), paste0("p", 1:S)))
  B[1, grp == 1] <- B[1, grp == 1] + 2          # strong true effect in g1
  B[2, 1:3] <- NA                               # missing patients tolerated
  st <- spiDE:::.limmaStage2(B, V, tau2 = 0.05, X)
  expect_named(st, c("gene", "coef", "t", "p"))
  expect_equal(nrow(st), G)
  expect_lt(st$p[1], 1e-4)
  expect_gt(st$coef[1], 1)
  expect_true(is.finite(st$p[2]))               # NA cells downweighted, not fatal
  # moderation shares variance: nulls should not all be extreme
  expect_gt(median(st$p[-1]), 0.05)
})
```

- [ ] **Step 2: Run to verify it fails** — expected: FAIL, `'.limmaStage2' not found`

- [ ] **Step 3: Implement**

Append to `R/twostage-stage2.R`:

```r
#' Moderated condition contrast on the patient slope matrix
#'
#' Observation weights 1/(v + tau2): random-effects weighting that degrades
#' toward equal weights when between-patient variation dominates (the
#' pseudobulk pitfall guard) and does real work when subset precision varies.
#' NA slopes enter with weight 0 and value 0 (the standard limma idiom).
#' @importFrom limma lmFit eBayes
#' @noRd
.limmaStage2 <- function(B, V, tau2, X) {
  Wt <- 1 / (V + tau2)
  bad <- !is.finite(B) | !is.finite(Wt)
  Bf <- B; Bf[bad] <- 0; Wt[bad] <- 0
  fit <- limma::lmFit(Bf, design = X, weights = Wt)
  fit <- limma::eBayes(fit, robust = TRUE)
  cn <- colnames(X)[2]
  data.frame(gene = rownames(B),
             coef = fit$coefficients[, cn],
             t = fit$t[, cn],
             p = fit$p.value[, cn],
             row.names = NULL, stringsAsFactors = FALSE)
}
```

In `DESCRIPTION`, add `    limma,` to `Imports:` (after `BiocParallel,`).

- [ ] **Step 4: Run to verify it passes** — same command, PASS.

- [ ] **Step 5: Commit**

```bash
git add R/twostage-stage2.R DESCRIPTION tests/testthat/test-twostage-stage2.R
git commit -m "Test the patient-level contrast with weighted moderated limma

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 10: Dropout diagnostics

**Files:**
- Modify: `R/twostage-stage2.R`
- Test: `tests/testthat/test-twostage-stage2.R`

**Interfaces:**
- Consumes: `ncells` data.frame from Task 6 (`sample`, `index`, `n`), `sample2patient`, patient condition factor `pgrp` (names = patients), `min.cells`.
- Produces: `.inclusionDiagnostics(ncells, sample2patient, pgrp, min.cells)` → data.frame(index, patient, condition, included); warns when inclusion is condition-associated (Fisher p < 0.05) for any index.

- [ ] **Step 1: Write the failing tests**

```r
test_that(".inclusionDiagnostics warns on condition-associated dropout", {
  ncells <- data.frame(sample = paste0("s", 1:12), index = "A",
                       n = c(rep(100, 6), rep(5, 6)))
  s2p <- stats::setNames(paste0("s", 1:12), paste0("s", 1:12))
  pgrp <- stats::setNames(factor(rep(c("R", "NR"), each = 6)), paste0("s", 1:12))
  # all R included, all NR dropped -> maximally associated
  expect_warning(
    d <- spiDE:::.inclusionDiagnostics(ncells, s2p, pgrp, min.cells = 30),
    "condition")
  expect_true(all(c("index", "patient", "condition", "included") %in% names(d)))
  # balanced inclusion -> silent
  ncells$n <- rep(c(100, 5), 6)
  expect_silent(spiDE:::.inclusionDiagnostics(ncells, s2p, pgrp, min.cells = 30))
})
```

- [ ] **Step 2: Run to verify it fails** — expected: FAIL, `'.inclusionDiagnostics' not found`

- [ ] **Step 3: Implement**

Append to `R/twostage-stage2.R`:

```r
#' Per-index patient inclusion table, with a condition-association warning
#'
#' min.cells dropout is informative when cell-type abundance correlates with
#' condition; this makes it visible instead of silent.
#' @noRd
.inclusionDiagnostics <- function(ncells, sample2patient, pgrp, min.cells) {
  ncells$patient <- unname(sample2patient[ncells$sample])
  agg <- stats::aggregate(n ~ index + patient, ncells, max)
  agg$included <- agg$n >= min.cells
  agg$condition <- as.character(pgrp[agg$patient])
  out <- agg[, c("index", "patient", "condition", "included")]
  for (ix in unique(out$index)) {
    d <- out[out$index == ix, ]
    if (length(unique(d$condition)) < 2 || all(d$included) || !any(d$included)) next
    p <- stats::fisher.test(table(d$condition, d$included))$p.value
    if (p < 0.05) {
      warning("min.cells dropout for index '", ix, "' is associated with ",
              "the condition (Fisher p = ", signif(p, 2), "); the tested ",
              "patient subset is condition-biased", call. = FALSE)
    }
  }
  out
}
```

- [ ] **Step 4: Run to verify it passes** — same command, PASS.

- [ ] **Step 5: Commit**

```bash
git add R/twostage-stage2.R tests/testthat/test-twostage-stage2.R
git commit -m "Report and warn on condition-associated min.cells dropout

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 11: Rewire `twoStageSpiDE()` end-to-end

**Files:**
- Modify: `R/twostage.R` (full rewrite of the exported function; keep `.looksLikeCounts()`)
- Modify: `tests/testthat/test-twostage.R` (update to the new option surface; add E2E + diagnostics tests)
- Modify: `man/` via `devtools::document()`

**Interfaces:**
- Consumes: everything from Tasks 1–10.
- Produces: the exported signature

```r
twoStageSpiDE(spe, condition, sigma, index = NULL, niche = NULL,
              patient = NULL, patient.covariates = character(),
              assay = "counts", cell_type = "cell_type",
              sample_id = "sample_id", name = "Niche",
              min.cells = 30L, fdr = 0.05,
              stage1 = c("spanorm", "ols", "nb"),
              epsilon = c("addback", "residual"),
              winsor = 4, lambda.a = 0, maxit.psi = 2,
              backend = "cpu", verbose = TRUE)
```

returning a `SpiDEResults` with the unchanged results-table schema and populated `@diagnostics`.

- [ ] **Step 1: Update the existing tests and add the new ones**

In `tests/testthat/test-twostage.R`: replace every `spiDE:::.toySPE(...)` + `buildNiches(...)` fixture with `toy_spanorm_spe(n_genes = <same>)` (the default path now needs the stored fit), keep all five existing behavioural assertions unchanged, and append:

```r
test_that("the planted G1/A/B effect is niche-specific on the new stage 1", {
  spe <- toy_spanorm_spe(n_genes = 40)
  r <- twoStageSpiDE(spe, condition = "condition", sigma = 30, index = "A",
                     min.cells = 10, fdr = 1, verbose = FALSE)
  tb <- results(r)
  g1 <- tb[tb$gene == "G1", ]
  # B is the strongest niche association for G1 in index A (repo convention:
  # niche-specificity, not genome-wide top rank)
  expect_identical(g1$ct_niche[which.min(g1$p.niche)], "B")
})

test_that("diagnostics are populated and the patient argument nests cores", {
  spe <- toy_spanorm_spe(n_genes = 20)
  # two cores per patient: collapse the 8 samples into 4 patients
  spe$patient <- paste0("pat", (as.integer(factor(spe$sample_id)) + 1) %/% 2)
  r <- twoStageSpiDE(spe, condition = "condition", sigma = 30,
                     patient = "patient", min.cells = 10, fdr = 1,
                     verbose = FALSE)
  expect_named(r@diagnostics, c("r2", "inclusion", "tau2"))
  expect_true(all(r@diagnostics$inclusion$patient %in% unique(spe$patient)))
})

test_that("removed stage1 options are rejected by match.arg", {
  spe <- toy_spanorm_spe(n_genes = 20)
  expect_error(twoStageSpiDE(spe, condition = "condition", sigma = 30,
                             stage1 = "nbresid", verbose = FALSE))
})
```

Note for the patient test: the toy's `condition` is assigned per sample; pairing samples into patients must keep condition constant within patient — pair samples within condition, i.e. build the pairing from `order(spe$condition, spe$sample_id)` if the naive pairing errors. Implementers: check `.toySPE()`'s sample-to-condition layout and pair accordingly; the assertion that matters is that a 2-core patient runs end-to-end.

- [ ] **Step 2: Run to verify the new tests fail**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-twostage.R")'`
Expected: FAIL — `twoStageSpiDE` still has the old signature/body (or was broken by Task 6).

- [ ] **Step 3: Rewrite the orchestration**

Replace the body of `twoStageSpiDE()` in `R/twostage.R` with:

```r
twoStageSpiDE <- function(spe, condition, sigma, index = NULL, niche = NULL,
                          patient = NULL, patient.covariates = character(),
                          assay = "counts", cell_type = "cell_type",
                          sample_id = "sample_id", name = "Niche",
                          min.cells = 30L, fdr = 0.05,
                          stage1 = c("spanorm", "ols", "nb"),
                          epsilon = c("addback", "residual"),
                          winsor = 4, lambda.a = 0, maxit.psi = 2,
                          backend = "cpu", verbose = TRUE) {
  stage1 <- match.arg(stage1); epsilon <- match.arg(epsilon)
  checkSPE(spe, assay = assay, cell_type = cell_type, sample_id = sample_id)
  checkCondition(spe, condition)
  checkNiche(spe, sigma, name = name)
  cd <- SummarizedExperiment::colData(spe)
  smp <- as.character(cd[[sample_id]])
  pat <- if (is.null(patient)) smp else as.character(cd[[patient]])
  ct <- as.character(cd[[cell_type]])
  grp_cell <- as.character(cd[[condition]])
  chk <- tapply(grp_cell, pat, function(z) length(unique(z[!is.na(z)])))
  if (any(chk > 1)) {
    stop("condition '", condition, "' varies within patient; two-stage needs ",
         "a patient-level condition")
  }
  s2p <- stats::setNames(vapply(unique(smp),
                                function(s) pat[match(s, smp)], character(1)),
                         unique(smp))
  nm <- as.matrix(SingleCellExperiment::reducedDim(spe, paste0(name, sigma)))
  if (!is.null(niche)) {
    keep <- intersect(niche, colnames(nm))
    if (!length(keep)) stop("no requested niche cell types found")
    nm <- nm[, keep, drop = FALSE]
  }
  idx_types <- sort(unique(ct))
  if (!is.null(index)) idx_types <- intersect(index, idx_types)
  if (!length(idx_types)) stop("no requested index cell types found")

  Y <- as.matrix(SummarizedExperiment::assay(spe, assay))
  if (stage1 %in% c("spanorm", "nb") && !.looksLikeCounts(Y)) {
    stop("stage1 = '", stage1, "' needs counts, but assay '", assay,
         "' does not look like counts")
  }
  comp <- if (stage1 == "spanorm") .spanormComponents(spe) else NULL
  E <- NULL
  if (stage1 == "ols") {
    lib <- colSums(Y)
    E <- log1p(sweep(Y, 2, mean(pmax(lib, 1)) / pmax(lib, 1), "*"))
  }
  if (verbose) {
    message(sprintf("stage 1 (%s): %d genes x %d index types x %d niches x %d samples",
                    stage1, nrow(Y), length(idx_types), ncol(nm),
                    length(unique(smp))))
  }
  sl <- .sampleSlopes(Y, E, comp, nm, ct, smp, idx_types, min.cells,
                      stage1 = stage1, epsilon = epsilon, winsor = winsor,
                      lambda.a = lambda.a, maxit.psi = maxit.psi,
                      backend = backend, verbose = verbose)

  pats <- unique(unname(s2p))
  pgrp <- factor(stats::setNames(
    vapply(pats, function(p) grp_cell[match(p, pat)], character(1)), pats))
  Xdes <- if (length(patient.covariates)) {
    cv <- as.data.frame(lapply(patient.covariates, function(v) {
      x <- tapply(cd[[v]], pat, function(z) z[1]); x[pats]
    }), col.names = patient.covariates, stringsAsFactors = TRUE)
    stats::model.matrix(~ g + ., data = cbind(data.frame(g = pgrp[pats]), cv))
  } else {
    stats::model.matrix(~ g, data = data.frame(g = pgrp[pats]))
  }
  diag_incl <- .inclusionDiagnostics(sl$ncells, s2p, pgrp, min.cells)

  recs <- list(); tau_tab <- list()
  for (ix in idx_types) {
    pooled <- .poolPatientSlopes(sl$beta[[ix]], sl$var[[ix]], s2p)
    for (nn in colnames(nm)) {
      if (identical(nn, ix)) next            # self-interaction, as the GLM drops
      B <- pooled$beta[, nn, pats, drop = TRUE]
      V <- pooled$var[, nn, pats, drop = TRUE]
      if (is.null(dim(B))) {
        B <- matrix(B, nrow = nrow(Y), dimnames = list(rownames(Y), pats))
        V <- matrix(V, nrow = nrow(Y), dimnames = list(rownames(Y), pats))
      }
      if (all(!is.finite(B))) next
      t2 <- .tau2DL(B, V, Xdes)
      st <- .limmaStage2(B, V, t2, Xdes)
      st <- st[is.finite(st$p), , drop = FALSE]
      if (!nrow(st)) next
      tau_tab[[length(tau_tab) + 1L]] <-
        data.frame(ct_index = ix, ct_niche = nn, tau2 = t2)
      recs[[length(recs) + 1L]] <- data.frame(
        gene = st$gene, ct_index = ix, ct_niche = nn,
        coef = st$coef, t = st$t, p.niche = st$p,
        stringsAsFactors = FALSE)
    }
  }
  empty <- data.frame(gene = character(), ct_index = character(),
                      ct_niche = character(), coef = numeric(), t = numeric(),
                      p.niche = numeric(), fdr.niche = numeric(),
                      DirectionNiche = character(), bandwidth.max = numeric(),
                      stringsAsFactors = FALSE)
  diagnostics <- list(
    r2 = sl$r2, inclusion = diag_incl,
    tau2 = if (length(tau_tab)) do.call(rbind, tau_tab) else
      data.frame(ct_index = character(), ct_niche = character(),
                 tau2 = numeric()))
  if (!length(recs)) {
    if (verbose) message("stage 2: no estimable triplets")
    return(new("SpiDEResults", fits = list(), sigma = sigma,
               condition = condition, mode = "condition", index = idx_types,
               niche = colnames(nm), covariates = patient.covariates,
               coldata = cd, gene.weights = NULL, p.cauchy.pos = NULL,
               p.cauchy.neg = NULL, results = empty, fdr = fdr,
               call = match.call(), diagnostics = diagnostics))
  }
  out <- do.call(rbind, recs)
  out$fdr.niche <- stats::p.adjust(out$p.niche, "BH")
  out$DirectionNiche <- ifelse(out$t > 0, "Up", "Down")
  out$bandwidth.max <- sigma
  out <- out[order(out$p.niche), , drop = FALSE]
  rownames(out) <- NULL
  if (verbose) {
    message(sprintf("stage 2: %d triplets tested, %d at FDR %.2g",
                    nrow(out), sum(out$fdr.niche <= fdr), fdr))
  }
  new("SpiDEResults", fits = list(), sigma = sigma, condition = condition,
      mode = "condition", index = idx_types, niche = colnames(nm),
      covariates = patient.covariates, coldata = cd, gene.weights = NULL,
      p.cauchy.pos = NULL, p.cauchy.neg = NULL,
      results = out[out$fdr.niche <= fdr, , drop = FALSE],
      fdr = fdr, call = match.call(), diagnostics = diagnostics)
}
```

Update the function's roxygen: new params (`patient`, `stage1` values, `epsilon`), delete the `"nbresid"`/`"analytic"` documentation paragraphs, document `@slot`-free `diagnostics` access via `r@diagnostics` in `@return`, keep the multiplicity warning text. Then run `devtools::document()`.

- [ ] **Step 4: Run the full suite**

Run: `Rscript -e 'devtools::test()'`
Expected: all files PASS (stage1, stage2, twostage, and the untouched remainder of the suite). Fix anything the rewrite broke before committing.

- [ ] **Step 5: Commit**

```bash
git add R/twostage.R tests/testthat/test-twostage.R man NAMESPACE DESCRIPTION
git commit -m "Rewire twoStageSpiDE around the SpaNorm-anchored joint estimator

Stage 1 defaults to the one-step log-scale response with the biology
added back; slopes pool per patient; limma tests the contrast with
1/(v + tau2) weights; diagnostics land on the results object.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 12: Adversarial validation scenarios (bench directory)

**Files:**
- Create: `/Users/uqdbhuva/R_packages/spiDE-twostage-bench/scenarios.R`
- Create: `/Users/uqdbhuva/R_packages/spiDE-twostage-bench/runner-adv.R`
- Modify: `/Users/uqdbhuva/R_packages/spiDE-twostage-bench/spiDE` (worktree: check out `twostage-fixes`)

These live outside the package (not committed to the package repo). The full sweep is user-triggered; this task delivers the scenario code plus a smoke run proving mechanics.

**Interfaces:**
- Consumes: `simulate_dataset()` from `research/R/sim_model.R`; the new `twoStageSpiDE()`.
- Produces: `simulate_adversarial(scenario, S, n_per, n_genes, seed)` for `scenario` in `"maineffect"` (condition baseline shift + shared nonzero niche slope, no true interaction), `"depth"` (condition-dependent within-patient depth–density coupling, no expression effect), `"attribution"` (two correlated niches, effect through one).

- [ ] **Step 1: Point the bench worktree at the branch**

```bash
git -C /Users/uqdbhuva/R_packages/spiDE-twostage-bench/spiDE checkout twostage-fixes
```

- [ ] **Step 2: Write `scenarios.R`**

```r
# Adversarial scenarios the published benchmark cannot see. Each wraps
# simulate_dataset() from research/R/sim_model.R and post-processes counts or
# labels to plant exactly one failure mode. All return list(spe, truth).
# SpaNorm is fit on each simulated SPE (df.tps = 2, LS from colSums) so the
# default stage1 = "spanorm" path runs.
source("/Users/uqdbhuva/R_packages/spiDE/research/R/sim_model.R")

.fit_spanorm <- function(spe) {
  SpaNorm::SpaNorm(spe, sample.p = 1, df.tps = 2,
                   batch = spe$sample_id, verbose = FALSE)
}

simulate_adversarial <- function(scenario = c("maineffect", "depth", "attribution"),
                                 S, n_per, n_genes, seed) {
  scenario <- match.arg(scenario)
  sim <- simulate_dataset(S = S, n_per = n_per, n_genes = n_genes,
                          n_marker = 50, n_hk = 150, layout = "gradient",
                          beta = 0, seed = seed)
  spe <- sim$spe
  set.seed(seed + 1L)
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  nm_target <- 50L                    # genes carrying the planted structure
  gsel <- seq_len(nm_target)
  resp <- spe$condition == "Responder"
  if (scenario == "maineffect") {
    # shared niche slope in EVERY patient + 2x baseline in responders;
    # no interaction -> every discovery among gsel is a false positive
    dens <- .niche_density_of(spe, type = "B", sigma = 50)
    sc <- exp(0.4 * scale(dens)[, 1]) * ifelse(resp, 2, 1)
    Y[gsel, ] <- round(sweep(Y[gsel, , drop = FALSE], 2, sc, "*"))
  } else if (scenario == "depth") {
    # depth scales with local density, 2x more strongly in responders;
    # expression itself is untouched
    dens <- .niche_density_of(spe, type = "B", sigma = 50)
    sc <- exp(0.3 * scale(dens)[, 1] * ifelse(resp, 2, 1))
    Y <- round(sweep(Y, 2, sc, "*"))
  } else {
    # two correlated niches: relabel half of type C cells near B as B-like
    # decoy "C"; true effect (added in responders only) follows B density
    dens <- .niche_density_of(spe, type = "B", sigma = 50)
    eff <- exp(0.5 * scale(dens)[, 1] * resp)
    Y[gsel, ] <- round(sweep(Y[gsel, , drop = FALSE], 2, eff, "*"))
  }
  SummarizedExperiment::assay(spe, "counts") <- Y
  spe <- .fit_spanorm(spe)
  list(spe = spe, truth = rownames(Y)[gsel], scenario = scenario)
}

.niche_density_of <- function(spe, type, sigma) {
  spe2 <- spiDE::buildNiches(spe, sigma = sigma, verbose = FALSE)
  as.matrix(SingleCellExperiment::reducedDim(spe2, paste0("Niche", sigma)))[, type]
}
```

Implementers: `simulate_dataset()`'s return fields and cell-type labels must be checked against `research/R/sim_model.R` (lines ~121–312) before finalising — the `"attribution"` scenario in particular should induce niche correlation the way the simulator's geometry allows (e.g. co-locating two types), not via relabelling if the layout already provides correlated types. The acceptance criteria below are what matter; adjust the generator mechanics to the simulator's actual API.

- [ ] **Step 3: Write `runner-adv.R`**

```r
# One (scenario, S, block) -> one part file, mirroring runner.R's shape.
# Measured quantities: type-I at 0.05 on the planted-structure genes
# (maineffect/depth: any discovery among them is false) or TPR/decoy-rate
# (attribution: discoveries on the true niche vs the correlated decoy).
BENCH <- "/Users/uqdbhuva/R_packages/spiDE-twostage-bench"
a <- commandArgs(trailingOnly = TRUE)
scenario <- a[1]; S <- as.integer(a[2]); block <- as.integer(a[3])
reps <- ((block - 1L) * 10L + 1L):(block * 10L)
part <- sprintf("%s/parts-adv/part_%s_S%02d_b%d.rds", BENCH, scenario, S, block)
if (file.exists(part)) quit(save = "no")
suppressPackageStartupMessages({
  library(SpatialExperiment); library(SingleCellExperiment); library(SummarizedExperiment)
  devtools::load_all(file.path(BENCH, "spiDE"), quiet = TRUE)
})
source(file.path(BENCH, "scenarios.R"))
out <- list()
for (r in reps) {
  sim <- simulate_adversarial(scenario, S = S, n_per = 300L, n_genes = 2000L,
                              seed = 200000L + 977L * r + 31L * S)
  spe <- spiDE::buildNiches(sim$spe, sigma = 50, verbose = FALSE)
  res <- try(spiDE::twoStageSpiDE(spe, condition = "condition", sigma = 50,
                                  index = "A", niche = "B", fdr = 1,
                                  min.cells = 30L, verbose = FALSE), silent = TRUE)
  if (inherits(res, "try-error")) next
  tb <- spiDE::results(res)
  p <- stats::setNames(tb$p.niche, tb$gene)
  hit <- p[sim$truth]
  out[[length(out) + 1L]] <- data.frame(
    scenario = scenario, S = S, rep = r,
    type1 = mean(hit < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE)
}
if (length(out)) saveRDS(do.call(rbind, out), part)
```

- [ ] **Step 4: Smoke-run one replicate per scenario**

```bash
cd /Users/uqdbhuva/R_packages/spiDE-twostage-bench && mkdir -p parts-adv
Rscript --vanilla runner-adv.R maineffect 10 1   # first rep only is enough:
```

Expected: a part file appears and, for `maineffect` and `depth`, per-rep `type1` is near 0.05 (the whole point of the fixes); the old estimator would show gross inflation. If type1 is grossly inflated (> 0.2), STOP — the fix or the scenario has a bug; investigate before proceeding. Do NOT launch the full sweep (user-triggered; ~hours).

- [ ] **Step 5: Report**

No package commit (bench files live outside the repo). Summarise smoke numbers for the user and hand back the decision on running the full sweep.

---

## Self-Review Notes

- **Spec coverage:** stage-1 estimator (Tasks 2–4), R² diagnostic (5), all-paths joint slopes + option removal (6), patient pooling (7), τ² (8), weighted limma + Welch/lm replacement (9), dropout diagnostics (10), diagnostics slot (1), orchestration + `patient` arg + E2E (11), adversarial validation (12). Out-of-scope items (multiplicity, bandwidths, one-stage path) untouched — matches spec.
- **Known judgement calls for implementers:** `.jointSlopes()` loops genes over small solves — acceptable at panel scale per subset; `SpaNorm::invert_mat_batched` is available if profiling demands it. The `"nb"` path's per-gene Fisher loop is slow by design (reference path only). Task 12's scenario generators must be reconciled with `sim_model.R`'s actual API; the acceptance criteria, not the sketch mechanics, are normative.
- **Type consistency check:** `.sampleSlopes` consumes `comp` (Task 2 list) and produces `beta`/`var`/`r2`/`ncells` exactly as Tasks 7/10/11 consume them; `.limmaStage2` takes the τ² scalar from `.tau2DL`; design matrix `Xdes` column 2 is the condition contrast in both Tasks 8 and 9 (factor `g`, treatment coding).
