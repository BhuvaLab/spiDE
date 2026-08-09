test_that(".cauchyCombine of identical p-values returns that p-value", {
  # each row constant across columns (columns = bandwidths)
  p <- matrix(c(0.02, 0.3, 0.8), nrow = 3, ncol = 3)
  expect_equal(spiDE:::.cauchyCombine(p), c(0.02, 0.3, 0.8), tolerance = 1e-8)
})

test_that(".cauchyCombine ignores zero-weight columns", {
  p <- cbind(c(0.01, 0.5), c(0.9, 0.9))
  w <- cbind(c(1, 1), c(0, 0))
  # with the second column zero-weighted, result equals the first column
  expect_equal(spiDE:::.cauchyCombine(p, w), c(0.01, 0.5), tolerance = 1e-8)
})

test_that(".cauchyCombine of a single column is (near) identity", {
  p <- matrix(c(0.05, 0.4, 0.7), ncol = 1)
  expect_equal(spiDE:::.cauchyCombine(p), c(0.05, 0.4, 0.7), tolerance = 1e-8)
})

test_that(".cauchyCombine handles NA entries via na.rm", {
  p <- cbind(c(0.02, 0.6), c(NA, 0.6))
  # first gene: only the non-NA 0.02 contributes
  expect_equal(spiDE:::.cauchyCombine(p)[1], 0.02, tolerance = 1e-8)
})

test_that(".geneWeights are relative likelihoods thresholded to zero", {
  spe <- buildNiches(.toySPE(), sigma = c(10, 30))
  res <- fitSpiDE(spe, condition = "condition", random = "none", verbose = FALSE)
  w <- spiDE:::.geneWeights(fits(res), thresh = 0.1)
  expect_equal(dim(w), c(20L, 2L))
  expect_true(all(w >= 0))
  expect_true(all(w <= 1))
  # each gene's best bandwidth has weight 1 (exp(0))
  expect_true(all(matrixStats::rowMaxs(w) == 1))
  # thresholded weights are exactly zero
  expect_true(all(w == 0 | w >= 0.1))
})
