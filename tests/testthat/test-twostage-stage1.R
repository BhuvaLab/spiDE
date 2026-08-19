# stage1 = "spanorm" now needs SpaNorm::fitNB(offset=) (in SpaNorm > 1.7.7,
# added on the fitnb-offset branch). Until that lands in the loaded SpaNorm,
# the spanorm-path tests skip rather than fail, and the pure-R pieces
# (.stage1Offset, .jointSlopes, .nicheBasisR2) still run.
has_offset <- "offset" %in% names(formals(SpaNorm::fitNB))

test_that(".spanormComponents extracts and validates the stored fit", {
  spe <- toy_spanorm_spe()
  comp <- spiDE:::.spanormComponents(spe)
  expect_named(comp, c("alpha", "gmean", "psi", "W", "bio", "off"))
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

test_that(".spanormComponents flags the ls+batch offset block", {
  spe <- toy_spanorm_spe()
  comp <- spiDE:::.spanormComponents(spe)
  expect_true(is.logical(comp$off) && any(comp$off))
  expect_false(any(comp$off & comp$bio))          # blocks are disjoint
})

test_that(".stage1Offset is the ls+batch linear predictor exactly", {
  spe <- toy_spanorm_spe()
  comp <- spiDE:::.spanormComponents(spe)
  cells <- 1:50
  O <- spiDE:::.stage1Offset(comp, cells)
  expect_identical(dim(O), c(length(comp$gmean), 50L))
  expect_true(all(is.finite(O)))
  # identity: full eta minus the biology part equals gmean + offset part
  Wc <- comp$W[cells, , drop = FALSE]
  eta_full <- comp$gmean + tcrossprod(comp$alpha, Wc)
  eta_bio <- tcrossprod(comp$alpha[, comp$bio, drop = FALSE],
                        Wc[, comp$bio, drop = FALSE])
  other <- !(comp$bio | comp$off)
  eta_rest <- if (any(other)) {
    tcrossprod(comp$alpha[, other, drop = FALSE], Wc[, other, drop = FALSE])
  } else 0
  expect_equal(O, eta_full - eta_bio - eta_rest, tolerance = 1e-12)
})

test_that(".jointSlopes recovers planted slopes and matches lm() weights", {
  set.seed(7)
  n <- 200; k <- 2; G <- 5
  X <- cbind(a = rnorm(n), b = rnorm(n))
  X <- sweep(X, 2, colMeans(X))
  beta_true <- matrix(c(0.5, -0.2), nrow = G, ncol = k, byrow = TRUE)
  eps <- beta_true %*% t(X) + matrix(rnorm(G * n, sd = 0.3), G, n)
  w <- matrix(runif(G * n, 0.5, 2), G, n)
  js <- spiDE:::.jointSlopes(eps, w, X)
  expect_equal(dim(js$beta), c(G, k))
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

test_that(".nicheBasisR2 returns NA, not ~1, for a degenerate/constant column", {
  # A constant column has tss ~ 0; the old code floored tss at 1e-12 and
  # reported r2 = 1 - 0/floor = 1 (a perfect-fit false positive). The slope
  # fit (.jointSlopes()) already reports NA for such a column, so the R2
  # diagnostic should agree instead of contradicting it.
  set.seed(9)
  n <- 120
  B <- cbind(rnorm(n), rnorm(n))
  const <- rep(5, n)
  r2 <- spiDE:::.nicheBasisR2(cbind(s = 2 * B[, 1] - B[, 2] + 3, c = const), B)
  expect_true(is.na(unname(r2["c"])))
  expect_equal(unname(r2["s"]), 1, tolerance = 1e-8)   # unaffected columns unchanged
})

test_that(".sampleSlopes returns aligned beta/var arrays for all paths", {
  spe <- toy_spanorm_spe()
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  nm <- as.matrix(SingleCellExperiment::reducedDim(spe, "Niche30"))
  ct <- as.character(spe$cell_type); smp <- as.character(spe$sample_id)
  comp <- spiDE:::.spanormComponents(spe)
  skip_if_not(has_offset, "SpaNorm::fitNB() lacks offset support")
  s1 <- spiDE:::.sampleSlopes(Y, NULL, comp, nm, ct, smp,
                              idx_types = c("A", "B"), min.cells = 10,
                              stage1 = "spanorm")
  expect_named(s1$beta, c("A", "B"))
  expect_identical(dim(s1$beta$A), dim(s1$var$A))
  expect_identical(dimnames(s1$beta$A)[[2]], colnames(nm))
  expect_true(all(s1$var$A >= 0, na.rm = TRUE))
  expect_true(all(c("sample", "index", "niche", "basis", "r2") %in%
                    names(s1$r2)))
  expect_setequal(unique(s1$r2$basis), c("biology", "ls"))
  # ols path: same shapes, no comp needed
  lib <- colSums(Y)
  E <- log1p(sweep(Y, 2, mean(lib) / pmax(lib, 1), "*"))
  s2 <- spiDE:::.sampleSlopes(Y, E, NULL, nm, ct, smp,
                              idx_types = "A", min.cells = 10,
                              stage1 = "ols")
  expect_identical(dim(s2$beta$A), dim(s1$beta$A))
})

test_that(".sampleSlopes excludes the index cell type's own niche column", {
  # The one-stage GLM design drops symmetric self-interactions (an index cell
  # type against its own niche density); stage 1 must do the same, or the
  # joint fit estimates a coefficient that's untestable/uninterpretable by
  # spec. The array keeps the full niches dimension for reporting -- the self
  # column is NA, not dropped from dimnames.
  spe <- toy_spanorm_spe()
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  nm <- as.matrix(SingleCellExperiment::reducedDim(spe, "Niche30"))
  ct <- as.character(spe$cell_type); smp <- as.character(spe$sample_id)
  comp <- spiDE:::.spanormComponents(spe)
  skip_if_not(has_offset, "SpaNorm::fitNB() lacks offset support")
  s <- spiDE:::.sampleSlopes(Y, NULL, comp, nm, ct, smp, idx_types = "A",
                             min.cells = 10, stage1 = "spanorm")
  expect_identical(dimnames(s$beta$A)[[2]], colnames(nm))  # full niches dim
  expect_true(all(is.na(s$beta$A["G1", "A", ])))
  expect_true(any(is.finite(s$beta$A["G1", "B", ])))
})

test_that("subsets below min.cells contribute NA, and are counted", {
  spe <- toy_spanorm_spe()
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  nm <- as.matrix(SingleCellExperiment::reducedDim(spe, "Niche30"))
  ct <- as.character(spe$cell_type); smp <- as.character(spe$sample_id)
  comp <- spiDE:::.spanormComponents(spe)
  s <- spiDE:::.sampleSlopes(Y, NULL, comp, nm, ct, smp, idx_types = "A",
                             min.cells = 10000L, stage1 = "spanorm")
  expect_true(all(is.na(s$beta$A)))
  expect_true(all(s$ncells$n < 10000L))
})

test_that("pooled psi is a sane dispersion estimate and pooling keeps stage 1 well-formed", {
  skip_if_not(has_offset, "SpaNorm::fitNB() lacks offset support")
  skip_if_not("psi" %in% names(formals(SpaNorm::fitNB)),
              "SpaNorm::fitNB() lacks psi passthrough")
  spe <- toy_spanorm_spe()
  Y <- as.matrix(SummarizedExperiment::assay(spe, "counts"))
  nm <- as.matrix(SingleCellExperiment::reducedDim(spe, "Niche30"))
  ct <- as.character(spe$cell_type); smp <- as.character(spe$sample_id)
  comp <- spiDE:::.spanormComponents(spe)

  # WHAT THIS DOES NOT TEST, and why. An earlier version asserted that the
  # pooled and unpooled slopes correlate > 0.98, and that the toy's planted
  # G1/A/B effect survives pooling. BOTH assertions were uninformative: on this
  # fixture the two-stage estimator has no power (every |t| < 1 -- measured
  # B = 0.34, C = 0.04 unpooled; C = 0.74, B = 0.35 pooled), so "which niche
  # ranks first" is noise, and a correlation between two near-null t vectors
  # carries no signal either. A test that cannot fail for the right reason is
  # worse than no test, because it invites exactly the misreading it got:
  # its failure was taken as evidence of a bug and drove a long, wrong
  # investigation. Power-based validation belongs in the simulator sweep,
  # where effects are detectable by construction.
  #
  # WHAT THIS DOES TEST: that .pooledPsi() returns a dispersion estimate in the
  # right ballpark per sample (it is the whole point of pooling), and that
  # pooling leaves stage 1 structurally sound. Both can fail.
  pp <- spiDE:::.pooledPsi(Y, comp, ct, smp, winsor = 4, lambda.a = 0,
                           maxit.psi = 2, backend = "cpu")
  expect_equal(sort(names(pp)), sort(unique(smp)))
  expect_true(all(vapply(pp, function(z) all(is.finite(z)) && all(z > 0), TRUE)))
  expect_true(all(vapply(pp, length, 1L) == nrow(Y)))

  # pooled vs per-subset dispersion: pooling borrows across cell types within a
  # sample, so it should be CLOSE to, not equal to, a subset-specific estimate.
  # Measured ratio ~1.07; the bounds below would catch a pooling fit that
  # collapsed, exploded, or silently returned another sample's values.
  ss <- sort(unique(smp))[1]
  cells <- which(ct == "A" & smp == ss)
  O <- spiDE:::.stage1Offset(comp, cells)
  X <- sweep(nm[cells, , drop = FALSE], 2, colMeans(nm[cells, , drop = FALSE]))
  X <- X[, setdiff(colnames(X), "A"), drop = FALSE]
  f <- SpaNorm::fitNB(Y[, cells], cbind(`(Intercept)` = 1, X), offset = O,
                      maxit.psi = 2, winsor = 4, verbose = FALSE)
  ratio <- stats::median(pp[[ss]] / f$psi)
  expect_gt(ratio, 0.5)
  expect_lt(ratio, 2.0)

  a <- spiDE:::.sampleSlopes(Y, NULL, comp, nm, ct, smp, idx_types = "A",
                             min.cells = 10, stage1 = "spanorm",
                             pool.psi = FALSE)
  b <- spiDE:::.sampleSlopes(Y, NULL, comp, nm, ct, smp, idx_types = "A",
                             min.cells = 10, stage1 = "spanorm",
                             pool.psi = TRUE)
  # structural: same shape, same triplets, finite where the unpooled arm is,
  # and strictly positive variances. A broken pooling path fails these.
  expect_identical(dim(a$beta$A), dim(b$beta$A))
  expect_identical(dimnames(a$beta$A), dimnames(b$beta$A))
  fin <- is.finite(a$beta$A) & is.finite(b$beta$A)
  expect_gt(sum(fin), 0)
  expect_true(all(b$var$A[fin] > 0))
})
