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
